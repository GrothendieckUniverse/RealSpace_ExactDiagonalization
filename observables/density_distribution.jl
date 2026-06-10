# ═══════════════════════════════════════════════════════════════════════════
# density_distribution — real-space occupation of graph vertices
#
# After extensive analysis (see `doc/observables.ipynb`), we have simplified
# this module: the symmetry-resolved single-sector density operator is
# mathematically guaranteed to give uniform ⟨n_i⟩ for any translation-invariant
# Hamiltonian with PBC, regardless of the underlying phase.
#
# Therefore this module retains ONLY the full real-space ED path (identity
# group).  It requires **open boundary conditions** — the only setting where
# a single ground-state eigenvector can display genuine real-space density
# modulation.
#
# For diagnosing charge order / superfluidity in PBC systems, use:
#   - `observables/static_structure_factor.jl`   → connected S(q)
#   - `observables/off_diagonal_long_range_order.jl` → ρ(k) = ⟨a†_i a_j⟩
# ═══════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════
# Fock-diagonal density kernels (shared with static_structure_factor)
# ═══════════════════════════════════════════════════════════════════════════

@inline function _density_value_for_mask(m::Mask, weights::AbstractVector{<:Number})::ComplexF64
    value = zero(ComplexF64)
    tmp = m
    @inbounds while tmp != 0
        lsb = tmp & -tmp
        idx = trailing_zeros(lsb) + 1
        value += weights[idx]
        tmp ⊻= lsb
    end
    return value
end

function _one_vertex_density_weights(n_site::Int, vertex::Int)::Vector{ComplexF64}
    weights = zeros(ComplexF64, n_site)
    weights[vertex] = COMPLEX_ONE
    return weights
end

function _flavor_density_weights(n_site::Int, density_flavor::Function)::Vector{ComplexF64}
    weights = zeros(ComplexF64, n_site)
    @inbounds for i in 1:n_site
        weights[i] = density_flavor(i) ? COMPLEX_ONE : zero(ComplexF64)
    end
    return weights
end

# ═══════════════════════════════════════════════════════════════════════════
# Full real-space ED occupation — OBC only
# ═══════════════════════════════════════════════════════════════════════════

"""
    vertices_occupation_distribution_full_ed(
        model;
        target_eigval_idx = 1,
        filling_fraction,
        flavor_filter = i -> true,
        ed_mode = :matrix,
        ed_data = nothing,
    ) -> Vector{Float64}

Compute `⟨n_i⟩` directly from a **full real-space ED eigenvector**, using the
identity symmetry group (each symmetry orbit = one raw Fock mask).

## ⚠️ Requires open boundary conditions

For translation-invariant Hamiltonians with periodic boundary conditions, a
single symmetry-sector ground state gives identically uniform density — a
mathematical consequence of `[H, T] = 0` and the non-degeneracy of the sector.
This function asserts `!any(pbc_indicator)` to prevent misuse.

To diagnose charge order or superfluidity in PBC systems, use:
- `static_structure_factor` (connected `S(q)`) in `observables/static_structure_factor.jl`
- `off_diagonal_long_range_order` (`ρ(k) = ⟨a_i^† a_j⟩`) in `observables/off_diagonal_long_range_order.jl`
"""
function vertices_occupation_distribution_full_ed(
    model::Real_Space_Second_Quantized_Model;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    flavor_filter::Function=i -> true,
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::Vector{Float64}
    @assert !any(model.lattice.pbc_indicator) "vertices_occupation_distribution_full_ed " *
        "requires open boundary conditions.  For PBC systems, use `static_structure_factor` " *
        "or `off_diagonal_long_range_order` instead (see `doc/observables.ipynb`)."

    n_site = model.lattice.n_site

    if ed_data === nothing
        flux0 = zeros(Float64, model.lattice.dim)
        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
        G = build_identity_group(n_site)
        ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=G)
        ed_scan!(ed_data; nev=max(target_eigval_idx, 2), mode=ed_mode)
    end

    irrep_idx = findfirst(irrep -> irrep.label == :identity, ed_data.irrep_list)
    irrep_idx === nothing &&
        error("Full real-space ED sector :identity was not found in ed_data.")

    haskey(ed_data.ed_scan_res, irrep_idx) ||
        error("Full real-space ED sector :identity was not scanned.")

    _, vecs = ed_data.ed_scan_res[irrep_idx]
    target_eigval_idx <= size(vecs, 2) ||
        error("target_eigval_idx=$target_eigval_idx exceeds the $(size(vecs, 2)) " *
              "computed full-ED eigenvectors.")

    basis = build_symmetry_sector_basis(ed_data.orbit_catalog, ed_data.irrep_list[irrep_idx])
    c = vecs[:, target_eigval_idx]
    length(c) == length(basis.representative_mask_list) ||
        error("Full-ED eigenvector length does not match identity-basis dimension.")

    density = zeros(Float64, n_site)
    @inbounds for (col, mask) in enumerate(basis.representative_mask_list)
        w = abs2(c[col])
        w == 0.0 && continue

        tmp = mask
        while tmp != 0
            lsb = tmp & -tmp
            vertex = trailing_zeros(lsb) + 1
            if flavor_filter(vertex)
                density[vertex] += w
            end
            tmp ⊻= lsb
        end
    end

    return density
end
