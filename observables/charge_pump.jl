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
Insert a flux `θ` along `flux_direction`, diagonalise the low-energy many-body states, project Resta's periodic position operator `exp(2π i X/L)` into that manifold, and unwrap the phase branches to obtain the pumped charge.

This is a one-dimensional flux-cylinder diagnostic (not a 2D many-body Chern number). At fractional filling a single-sector expectation value can vanish by translation symmetry; passing the full topological multiplet (e.g. `[(0,0), (1,0)]` for the `[2,3]` Haldane FCI) allows the position operator to connect the sectors.

- Args:
    - `model::Real_Space_Second_Quantized_Model`: the second quantized model
    - `sector_labels`: `:identity` (full Hilbert space) or a `Vector` of sector labels, e.g. `[(0,0), (1,0)]`
- Named Args:
    - `filling_fraction::Rational{Int}`: particles per flattened vertex
    - `flux_direction::Int=1`: direction along which flux is threaded
    - `polarization_direction::Int`: direction of the polarisation measurement (default: transverse to flux)
    - `twisted_phases_list::Vector{Float64}=collect(range(0.0,1.0;length=9))`: list of flux values
    - `nev_per_sector::Int=1`: number of low-lying states per symmetry sector
    - `include_sublattice::Bool=true`: include sublattice offset in position coordinate
    - `fig_path::Union{Nothing,String}=nothing`: path to save SVG/PNG figure
    - `checkpoint_path::Union{Nothing,String}=nothing`: path to JLD2 checkpoint
- Returns:
    - `res::NamedTuple`: `twisted_phases_list`, `energies`, `pumped_charge_trajectories`, `pumped_charges`, ...
"""
function flux_charge_pump(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    flux_direction::Int=1,
    polarization_direction::Int=_default_polarization_direction(model.lattice.dim, flux_direction),
    twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    nev_per_sector::Int=1,
    include_sublattice::Bool=true,
    fig_path::Union{Nothing,String}=nothing,
    checkpoint_path::Union{Nothing,String}=nothing,
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
    energies = fill(NaN, length(twisted_phases_list), length(labels), nev_per_sector)
    polarization_eigenvalues = Matrix{ComplexF64}(undef, length(twisted_phases_list), nstates)
    polarization_matrices = Vector{Matrix{ComplexF64}}(undef, length(twisted_phases_list))

    # ── Initialise model at θ=0, build ED data once ──
    flux0 = zeros(Float64, dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_group = is_identity ? build_identity_group(model.lattice.n_site) :
                 build_translation_group(model.lattice, flux0)
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=init_group)

    for (iθ, θ_val) in enumerate(twisted_phases_list)
        flux = zeros(Float64, dim)
        flux[flux_direction] = θ_val

        # ── In-place update: bilinear terms + symmetry group + catalog ──
        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux)

        if is_identity
            ed_data.symmetry_group = build_identity_group(model.lattice.n_site)
        else
            active_group = build_translation_group(model.lattice, flux)
            ed_data.symmetry_group = active_group
            update_orbit_stabilizer_phases!(ed_data.orbit_catalog, active_group, model.statistics)
        end

        # ── Diagonalise & build sector bases ──
        bases = Dict{Any,Symmetry_Sector_Basis}()
        eigvecs = Dict{Any,Matrix{ComplexF64}}()
        for (ilabel, label) in enumerate(labels)
            vals, vecs = ed_scan_at_irrep_matrix!(label, ed_data; nev=nev_per_sector)
            energies[iθ, ilabel, 1:length(vals)] .= vals
            irrep, _ = _irrep_for_label(ed_data, label)
            bases[label] = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
            eigvecs[label] = vecs
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

        polarization_matrices[iθ] = P
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
        ax = Axis(fig[1, 1];
            xlabel="Inserted Flux [2π] along Direction-$(flux_direction)",
            ylabel="Pumped Charge ΔQ",
            title="$(model.lattice.sample_size)-sample Charge Pump ($(is_identity ? "Full" : "Sector-Resolved"))"
        )
        for b in 1:nstates
            lines!(ax, twisted_phases_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), linewidth=2, label="branch $b")
            scatter!(ax, twisted_phases_list, pumped_charge_trajectories[:, b];
                color=Makie.Cycled(b), markersize=6)
        end
        axislegend(ax; position=:lt)
        save(fig_path, fig)
        @info "charge pump plot saved to `$fig_path`"
    end

    # ── Checkpoint ──
    if !isnothing(checkpoint_path)
        mkpath(dirname(checkpoint_path))
        result = (; twisted_phases_list, energies, sector_labels=labels, flux_direction,
            polarization_direction, nev_per_sector, polarization_eigenvalues,
            polarization_matrices, polarizations, pumped_charge_trajectories, pumped_charges,
            include_sublattice, is_identity, fig_path)
        @save checkpoint_path result
    end

    res = (; twisted_phases_list, energies, sector_labels=labels, flux_direction,
        polarization_direction, nev_per_sector, polarization_eigenvalues,
        polarization_matrices, polarizations, pumped_charge_trajectories, pumped_charges,
        include_sublattice, is_identity, fig_path)

    return res
end
