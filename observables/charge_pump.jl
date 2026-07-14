# ═══════════════════════════════════════════════════════════════════════════
# charge_pump — one-dimensional flux insertion and many-body polarization
# ═══════════════════════════════════════════════════════════════════════════

"""
    many_body_position_phases(lattice, direction; include_sublattice=true)

Return the single-site phase factors for Resta's periodic position operator in direction-α `U_α = exp(2π i / L_α * Σ_j x_{j,α} n̂_j)`.

The coordinate x_{j,α} is the α-component of the crystal coordinate of site j, i.e. `cell_int[α] / L_α` plus an optional sub-lattice offset. The returned phase for site j is `exp(2π i * x_{j,α` (used as `site_phases[i]` in `_position_phase_for_mask`).
"""
function many_body_position_phases(
    lattice::TightBinding.Real_Space_Lattice,
    direction::Int;
    include_sublattice::Bool=true,
)::Vector{ComplexF64}
    1 <= direction <= lattice.dim || error("direction must be in 1:$(lattice.dim).")
    L = Float64(lattice.sample_size[direction])
    phases = Vector{ComplexF64}(undef, lattice.n_site)
    @inbounds for (i, (cell_int, isub)) in enumerate(lattice.site_list)
        x = Float64(cell_int[direction])
        if include_sublattice
            x += Float64(lattice.sub_crys_list[isub][direction])
        end
        phases[i] = cis(2π * x / L)
    end
    return phases
end

# ═══════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════

@inline function _position_phase_for_mask(m::Mask, site_phases::Vector{ComplexF64})::ComplexF64
    z = COMPLEX_ONE
    tmp = m
    @inbounds while tmp != 0
        lsb = tmp & -tmp
        idx = trailing_zeros(lsb) + 1
        z *= site_phases[idx]
        tmp ⊻= lsb
    end
    return z
end

"Build the projected position-operator matrix `⟨χ_to|U_d|χ_from⟩` between two symmetry sectors."
function _position_operator_matrix(
    basis_to::Symmetry_Sector_Basis,
    basis_from::Symmetry_Sector_Basis,
    site_phases::Vector{ComplexF64},
    particle_statistics::Particle_Statistics,
)::SparseMatrixCSC{ComplexF64,Int}
    G = basis_from.symmetry_group
    @assert basis_to.symmetry_group.name == G.name
    @assert group_order(basis_to.symmetry_group) == group_order(G)
    @assert basis_to.symmetry_group.n_site == G.n_site

    n_to = length(basis_to.representative_mask_list)
    n_from = length(basis_from.representative_mask_list)
    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    sizehint!(Is, min(n_to, n_from))
    sizehint!(Js, min(n_to, n_from))
    sizehint!(Vs, min(n_to, n_from))

    χ_to = basis_to.irrep.values
    χ_from = basis_from.irrep.values
    inv_nG = 1.0 / group_order(G)

    @inbounds for (col, repr) in enumerate(basis_from.representative_mask_list)
        row = _basis_index(basis_to, repr)
        row == 0 && continue

        elem = zero(ComplexF64)
        for (gidx, op) in enumerate(G.operations)
            shifted, _ = apply_operation_to_mask(repr, op, particle_statistics)
            elem += χ_to[gidx] * conj(χ_from[gidx]) *
                    _position_phase_for_mask(shifted, site_phases)
        end
        elem *= inv_nG

        if abs(elem) > 1e-13
            push!(Is, row)
            push!(Js, col)
            push!(Vs, elem)
        end
    end

    return sparse(Is, Js, Vs, n_to, n_from)
end

function _default_polarization_direction(dim::Int, flux_direction::Int)::Int
    dim == 1 && return flux_direction
    return mod1(flux_direction + 1, dim)
end

function _irrep_for_label(ed_data::Symmetry_Resolved_ED_Data, label)
    idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
    idx === nothing && error("No irrep with label $(repr(label)).")
    return ed_data.irrep_list[idx], idx
end

"Unwrap phase branches along the flux path using optimal permutation matching at each step."
function _tracked_phase_branches(raw_phases::Matrix{Float64})::Matrix{Float64}
    nθ, nbranch = size(raw_phases)
    tracked = similar(raw_phases)
    order0 = sortperm(raw_phases[1, :])
    tracked[1, :] .= raw_phases[1, order0]
    prev = copy(tracked[1, :])

    for iθ in 2:nθ
        best_perm = _best_phase_permutation(view(raw_phases, iθ, :), prev)
        for b in 1:nbranch
            φ = raw_phases[iθ, best_perm[b]]
            tracked[iθ, b] = φ + round(prev[b] - φ)
        end
        prev .= tracked[iθ, :]
    end
    return tracked
end

"Brute-force optimal assignment of `nbranch` phases to minimize discontinuity."
function _best_phase_permutation(raw, prev::Vector{Float64})::Vector{Int}
    nbranch = length(prev)
    used = falses(nbranch)
    current = zeros(Int, nbranch)
    best = collect(1:nbranch)
    best_cost = Ref(Inf)

    function visit(depth::Int, cost::Float64)
        if depth > nbranch
            if cost < best_cost[]
                best_cost[] = cost
                best .= current
            end
            return nothing
        end
        for j in 1:nbranch
            used[j] && continue
            δ = raw[j] - prev[depth]
            δ -= round(δ)
            next_cost = cost + δ^2
            next_cost >= best_cost[] && continue
            used[j] = true
            current[depth] = j
            visit(depth + 1, next_cost)
            used[j] = false
        end
        return nothing
    end

    visit(1, 0.0)
    return best
end

"Normalize a low-energy manifold to explicit `(sector, level)` state specifications."
function _normalize_manifold_state_specs(sector_labels, nev_per_sector::Int, manifold_states)
    nev_per_sector > 0 || error("nev_per_sector must be positive.")
    is_identity = sector_labels == :identity
    base_labels = is_identity ? [:identity] :
                  sector_labels isa Tuple ? [Tuple(Int.(sector_labels))] :
                  [Tuple(Int.(label)) for label in sector_labels]

    specs = if manifold_states === nothing
        [(sector=label, level=level) for label in base_labels for level in 1:nev_per_sector]
    else
        is_identity && error("Explicit manifold_states are only supported for symmetry-sector ED.")
        normalized = NamedTuple[]
        for state in manifold_states
            sector_raw = hasproperty(state, :sector) ? getproperty(state, :sector) : state[1]
            level_raw = hasproperty(state, :level) ? getproperty(state, :level) : state[2]
            sector = Tuple(Int.(sector_raw))
            level = Int(level_raw)
            level > 0 || error("Manifold levels must be positive; got $level in sector $sector.")
            sector in base_labels || error(
                "Manifold state ($sector, level $level) is absent from sector_labels=$base_labels.")
            push!(normalized, (sector=sector, level=level))
        end
        isempty(normalized) && error("manifold_states must not be empty.")
        length(unique((state.sector, state.level) for state in normalized)) == length(normalized) ||
            error("manifold_states contains duplicate (sector, level) entries.")
        normalized
    end

    labels = unique(state.sector for state in specs)
    required_levels = Dict{Any,Int}()
    for state in specs
        required_levels[state.sector] = max(get(required_levels, state.sector, 0), state.level)
    end
    return labels, specs, required_levels, is_identity
end

# ═══════════════════════════════════════════════════════════════════════════
# flux_charge_pump — scan many-body polarization P(θ)
# ═══════════════════════════════════════════════════════════════════════════
"""
Core Method to Compute the Many-Body Charge Pump Under a Twisted Boundary Condition
---
Reads the full symmetry-resolved ED data from canonical per-θ checkpoints
(produced by [`ed_scan!`](@ref) in flux-scan mode), projects Resta's periodic
position operator `exp(2π i X/L)` into the low-energy manifold for the requested
sectors, and unwraps the phase branches.

- Args:
    - `model::Real_Space_Second_Quantized_Model`
    - `sector_labels`: `:identity` or a `Vector` of sector tuples
- Named Args:
    - `filling_fraction::Rational{Int}`: particles per flattened vertex
    - `flux_direction::Int=1`
    - `polarization_direction::Int`: default transverse to flux
    - `twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`
    - `nev_per_sector::Int=1`: low-lying states per sector for the legacy
      one-size-fits-all manifold specification
    - `manifold_states=nothing`: optional explicit `(sector, level)` states;
      required when multiple manifold states occupy the same momentum sector
    - `include_sublattice::Bool=true`
    - `checkpoint_dir::String=\"checkpoints\"`: directory for per-θ checkpoints
    - `fig_path::Union{Nothing,String}=nothing`
    - `overwrite::Bool=false`
- Returns:
    - `res::NamedTuple`: `pumped_charges`, `pumped_charge_trajectories`, `polarizations`, ...
"""
function flux_charge_pump(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    flux_direction::Int=1,
    polarization_direction::Int=_default_polarization_direction(model.lattice.dim, flux_direction),
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    nev_per_sector::Int=1,
    manifold_states=nothing,
    include_sublattice::Bool=true,
    checkpoint_dir::String="checkpoints",
    fig_path::Union{Nothing,String}=nothing,
    overwrite::Bool=false,
)
    # ── Resolve every state in the low-energy manifold explicitly.  This is
    # essential when, e.g., three FCI states are levels 1:3 of one momentum
    # sector rather than level 1 of three distinct sectors.
    labels, state_specs, required_levels, is_identity =
        _normalize_manifold_state_specs(sector_labels, nev_per_sector, manifold_states)

    dim = model.lattice.dim
    1 <= flux_direction <= dim || error("flux_direction must be in 1:$dim.")
    1 <= polarization_direction <= dim || error("polarization_direction must be in 1:$dim.")

    nstates = length(state_specs)
    max_level = maximum(values(required_levels))
    energies = fill(NaN, length(twisted_phases_over_2π_list), length(labels), max_level)
    polarization_eigenvalues = Matrix{ComplexF64}(undef, length(twisted_phases_over_2π_list), nstates)

    # ── Ensure all per-θ ED checkpoints exist (shared with spectrum flow) ──
    flux0 = zeros(Float64, model.lattice.dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_ed_data = build_ed_data(model; filling_fraction=filling_fraction,
        symmetry_group=build_translation_group(model.lattice, flux0))

    checkpoint_paths = ed_scan!(init_ed_data;
        nev=max(max_level, 2),
        mode=:matrix,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
        flux_direction=flux_direction,
        checkpoint_dir=checkpoint_dir,
        overwrite=overwrite,
        scanned_sectors=(is_identity ? nothing : labels),
    )

    for (iθ, ckpt_path) in enumerate(checkpoint_paths)
        ed_data = load_checkpoint(ckpt_path)

        # ── Build sector bases and extract eigenvectors ──
        bases = Dict{Any,Symmetry_Sector_Basis}()
        eigvecs = Dict{Any,Matrix{ComplexF64}}()
        for (ilabel, label) in enumerate(labels)
            irrep_idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
            if irrep_idx === nothing || !haskey(ed_data.ed_scan_res, irrep_idx)
                @warn "Sector $label not found in checkpoint $ckpt_path"
                continue
            end
            vecs = ed_data.ed_scan_res[irrep_idx][2]
            nv = required_levels[label]
            size(vecs, 2) >= nv || error(
                "Sector $label in $ckpt_path has $(size(vecs, 2)) eigenvectors; need $nv.")
            eigvecs[label] = vecs[:, 1:nv]
            vals = ed_data.ed_scan_res[irrep_idx][1]
            length(vals) >= nv || error(
                "Sector $label in $ckpt_path has $(length(vals)) eigenvalues; need $nv.")
            energies[iθ, ilabel, 1:nv] .= vals[1:nv]

            irrep = ed_data.irrep_list[irrep_idx]
            bases[label] = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
        end

        # ── Build position operator matrix projected into the low-energy manifold ──
        site_phases = many_body_position_phases(
            model.lattice, polarization_direction;
            include_sublattice=include_sublattice,
        )
        P = zeros(ComplexF64, nstates, nstates)

        position_blocks = Dict{Tuple{Any,Any},Any}()
        for (row, state_to) in enumerate(state_specs), (col, state_from) in enumerate(state_specs)
            label_to, lev_to = state_to.sector, state_to.level
            label_from, lev_from = state_from.sector, state_from.level
            block_key = (label_to, label_from)
            Ublock = get!(position_blocks, block_key) do
                _position_operator_matrix(
                    bases[label_to], bases[label_from], site_phases,
                    model.particle_statistics)
            end
            P[row, col] = dot(eigvecs[label_to][:, lev_to],
                Ublock * eigvecs[label_from][:, lev_from])
        end

        ev = eigen(P).values
        order = sortperm(angle.(ev))
        polarization_eigenvalues[iθ, :] .= ev[order]
    end

    # ── Unwrap phase branches and compute pumped charge ──
    raw_phases = angle.(polarization_eigenvalues) ./ (2π)
    polarizations = _tracked_phase_branches(raw_phases)
    pumped_charge_trajectories = polarizations .- reshape(polarizations[1, :], 1, :)
    pumped_charges = pumped_charge_trajectories[end, :]

    # ── Plot ──
    if !isnothing(fig_path)
        mkpath(dirname(fig_path))
        fig = Figure(size=(800, 600))
        sectors_str = is_identity ? "Full Hilbert Space" : join(repr.(labels), ", ")
        ax = Axis(fig[1, 1];
            xlabel="Inserted Flux [2π] along Direction-$(flux_direction)",
            ylabel="Pumped Charge ΔQ",
            title="$(model.lattice.sample_size)-sample Charge Pump — sectors: $sectors_str"
        )
        for b in 1:nstates
            state = state_specs[b]
            lbl = is_identity ? "branch $b" :
                  "sector $(repr(state.sector)), level $(state.level)"
            lines!(ax, twisted_phases_over_2π_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), linewidth=2, label=lbl)
            scatter!(ax, twisted_phases_over_2π_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), markersize=6)
        end
        axislegend(ax; position=:lt)
        save(fig_path, fig)
        @info "charge pump plot saved to `$fig_path`"
    end

    @info "Pumped charges for manifold states" state_specs pumped_charges

    res = (; twisted_phases_over_2π_list, energies, sector_labels=labels, flux_direction,
        polarization_direction, nev_per_sector=max_level, manifold_states=state_specs,
        polarization_eigenvalues,
        polarizations, pumped_charge_trajectories, pumped_charges,
        include_sublattice, is_identity, fig_path, checkpoint_paths)

    return res
end
