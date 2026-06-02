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
    statistics::Statistics,
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
            shifted, _ = apply_operation_to_mask(repr, op, statistics)
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
    - `nev_per_sector::Int=1`: low-lying states per sector
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
    include_sublattice::Bool=true,
    checkpoint_dir::String="checkpoints",
    fig_path::Union{Nothing,String}=nothing,
    overwrite::Bool=false,
)
    # ── Resolve sector labels ──
    is_identity = sector_labels == :identity
    labels = is_identity ? [:identity] :
             sector_labels isa Tuple ? [Tuple(Int.(sector_labels))] :
             [Tuple(Int.(l)) for l in sector_labels]

    dim = model.lattice.dim
    1 <= flux_direction <= dim || error("flux_direction must be in 1:$dim.")
    1 <= polarization_direction <= dim || error("polarization_direction must be in 1:$dim.")

    nstates = length(labels) * nev_per_sector
    energies = fill(NaN, length(twisted_phases_over_2π_list), length(labels), nev_per_sector)
    polarization_eigenvalues = Matrix{ComplexF64}(undef, length(twisted_phases_over_2π_list), nstates)

    # ── Ensure all per-θ ED checkpoints exist (shared with spectrum flow) ──
    flux0 = zeros(Float64, model.lattice.dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_ed_data = build_ed_data(model; filling_fraction=filling_fraction,
        symmetry_group=build_translation_group(model.lattice, flux0))

    checkpoint_paths = ed_scan!(init_ed_data;
        nev=max(nev_per_sector, 2),  # need at least 1 eigenvector
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
            nv = min(size(vecs, 2), nev_per_sector)
            eigvecs[label] = vecs[:, 1:nv]
            vals = ed_data.ed_scan_res[irrep_idx][1]
            energies[iθ, ilabel, 1:min(length(vals), nev_per_sector)] .= vals[1:min(length(vals), nev_per_sector)]

            irrep = ed_data.irrep_list[irrep_idx]
            bases[label] = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
        end

        # ── Build position operator matrix projected into the low-energy manifold ──
        site_phases = many_body_position_phases(
            model.lattice, polarization_direction;
            include_sublattice=include_sublattice,
        )
        P = zeros(ComplexF64, nstates, nstates)

        for (ito_label, label_to) in enumerate(labels)
            basis_to = bases[label_to]
            vecs_to = eigvecs[label_to]
            for (ifrom_label, label_from) in enumerate(labels)
                basis_from = bases[label_from]
                vecs_from = eigvecs[label_from]
                Ublock = _position_operator_matrix(
                    basis_to, basis_from, site_phases, model.statistics,
                )
                for lev_to in 1:size(vecs_to, 2), lev_from in 1:size(vecs_from, 2)
                    row = (ito_label - 1) * nev_per_sector + lev_to
                    col = (ifrom_label - 1) * nev_per_sector + lev_from
                    P[row, col] = dot(vecs_to[:, lev_to], Ublock * vecs_from[:, lev_from])
                end
            end
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
            lbl = is_identity ? "branch $b" : "sector $(repr(labels[mod1(b, length(labels))])) branch $(div(b-1, length(labels))+1)"
            lines!(ax, twisted_phases_over_2π_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), linewidth=2, label=lbl)
            scatter!(ax, twisted_phases_over_2π_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), markersize=6)
        end
        axislegend(ax; position=:lt)
        save(fig_path, fig)
        @info "charge pump plot saved to `$fig_path`"
    end

    @info "Pumped Charges for sectors $sector_labels: $pumped_charges"

    res = (; twisted_phases_over_2π_list, energies, sector_labels=labels, flux_direction,
        polarization_direction, nev_per_sector, polarization_eigenvalues,
        polarizations, pumped_charge_trajectories, pumped_charges,
        include_sublattice, is_identity, fig_path, checkpoint_paths)

    return res
end
