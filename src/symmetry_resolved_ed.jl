# ============================================================================
# Symmetry-Resolved Exact Diagonalization — XDiag-inspired High-Performance Engine
#
# Supports hard-core bosons and fermions on arbitrary real-space lattices.
#
# Two modes (following XDiag):
#   • "matrix" mode  — precompute sparse H, then diagonalize (fast, memory-heavy)
#   • "matrixfree" mode — on-the-fly H|ψ⟩ via multithreaded Lanczos (~1.5× slower,
#                          near-zero memory overhead beyond basis vectors)
#
# Architecture:
#   1. Symmetry group G → Symmetry_Operation (perm + U(1) phases)
#   2. Gosper's hack → enumerate all bitmasks at fixed particle number
#   3. Orbit-stabilizer decomposition → Symmetry_Orbit_Catalog
#   4. 1D irreps → filter orbits → Symmetry_Sector_Basis
#   5a. [matrix]  build sparse CSC block → Arpack
#   5b. [matrixfree] CanonicalMap for O(1) lookup + Threads.@threads apply_H! → KrylovKit
#
# Key references:
#   - XDiag: https://github.com/awietek/xdiag  /  arXiv:2505.02901
#   - D.N. Sheng et al., Phys. Rev. Lett. 107, 146803 (2011)
# ============================================================================

using Base.Threads
using Distributed
using LinearAlgebra, SparseArrays, Arpack, KrylovKit
using Printf, JLD2
using MLStyle

# ═══════════════════════════════════════════════════════════════════════════
# Parallelism strategy (HPC-first):
#   - Launch: julia -p N  (distributed workers)
#   - Default: BLAS threads = 1 (don't compete with distributed workers)
#   - Hamiltonian construction: Distributed.pmap across workers
#   - Diagonalization: temporarily BLAS.set_num_threads(nprocs()), restore to 1
#   - Matrix-free H|ψ⟩: Threads.@threads (shared-memory) + optional pmap columns
# ═══════════════════════════════════════════════════════════════════════════

BLAS.set_num_threads(1)  # HPC default: one BLAS thread per process

# Import bitwise operations from the submodule
using .BitWise_Operations: Mask, COMPLEX_ONE, bitmask_of_site,
    occupy_site_for_mask, empty_site_for_mask,
    is_site_occupied, is_site_empty, n_occupied_for_mask,
    filled_site_iter_for_mask, empty_site_iter_for_mask,
    encode_configuration_to_bit_mask, decode_bit_mask_to_configuration, decode_bit_mask_to_configuration!

# ═══════════════════════════════════════════════════════════════════════════
# 1. Symmetry Operation — a single group element
# ═══════════════════════════════════════════════════════════════════════════

"""
Struct `Symmetry_Operation{Group_Label}` for One Single Group Action on the Occupation Basis
---
following (for both boson and fermions) `U_g |n₁,…,n_N⟩ = (∏_{i∈occ} η_g(i)) * sgn_g(occ) * |n_{π(1)},…,n_{π(N)}⟩`

In the second quantized form of the occupation basis, it acts as `U_g b_i^† U_g^{-1} = η_g(i) b_{perm[i]}^†`,
- Fields:
    - `label::Group_Label`: the generic `Group_Label`-typed label of the group element, depending on the symmetry transformation we are considering
    - `perm::Vector{Int}`: the permutation of the whole vertices, which is a vector of length `N` with each element being an integer from `1` to `N`
    - `perm_phases::Vector{ComplexF64}`: the associated phase factor for each vertex, which is a vector of length `N` with each element being a complex number. This is used to capture the nontrivial U(1) phase factor in the symmetry action
"""
struct Symmetry_Operation{Group_Label}
    label::Group_Label
    perm::Vector{Int}               # site permutation (1-based)
    perm_phases::Vector{ComplexF64} # U(1) phase η_g(i) per site
end

function Symmetry_Operation(label::Group_Label, perm::Vector{Int}; perm_phases=nothing) where {Group_Label}
    n_site = length(perm)
    phases = perm_phases === nothing ? fill(COMPLEX_ONE, n_site) : ComplexF64.(perm_phases)
    @assert length(phases) == n_site
    return Symmetry_Operation{Group_Label}(label, perm, phases)
end

# ═══════════════════════════════════════════════════════════════════════════
# 2. Finite Symmetry Group
# ═══════════════════════════════════════════════════════════════════════════

abstract type Abstract_Symmetry_Group end
"""
Struct `Finite_Symmetry_Group`
---
as a collection of `Symmetry_Operation`s together with some metadata.
- Fields:
    - `name::String`: group name
    - `n_site::Int`: number of vertices in the graph
    - `operations::Vector{<:Symmetry_Operation}`: vector of symmetry operations, each of which is a group element acting on the occupation basis
    - `identity_idx::Int`: index of the identity operation in the `operations` vector
"""
struct Finite_Symmetry_Group <: Abstract_Symmetry_Group
    name::String
    n_site::Int
    operations::Vector{<:Symmetry_Operation}
    identity_idx::Int
end

"constructor for `Finite_Symmetry_Group` with optional `identity_idx` argument indicating the linear-index of the identity group operation"
function Finite_Symmetry_Group(name::String, ops::Vector{<:Symmetry_Operation}; identity_idx::Int=1)
    @assert !isempty(ops)
    n_site = length(ops[1].perm)
    @assert all(length(op.perm) == n_site && length(op.perm_phases) == n_site for op in ops)
    @assert ops[identity_idx].perm == collect(1:n_site) "Identity must have trivial permutation"
    return Finite_Symmetry_Group(name, n_site, ops, identity_idx)
end

@inline group_order(G::Finite_Symmetry_Group)::Int = length(G.operations)

function ensure_distributed_workers_loaded!()
    nprocs() == 1 && return nothing
    for w in workers()
        remotecall_fetch(Core.eval, w, Main, :(using RealSpace_ExactDiagonalization))
        remotecall_fetch(Core.eval, w, Main, :(using LinearAlgebra; LinearAlgebra.BLAS.set_num_threads(1)))
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
# 3. Group action on a bitmask — O(k) with bit tricks
# ═══════════════════════════════════════════════════════════════════════════

"""
    apply_operation_to_mask(m, op, statistics) -> (new_mask, phase)

Apply g ∈ G to a Fock state.

Bosons:  α = ∏ η_g(i)
Fermions: α = ∏ η_g(i) × (-1)^{#inversions}, tracked via count_ones(new_mask >> π(i))
"""
@inline function apply_operation_to_mask(m::Mask, op::Symmetry_Operation, ::Bosonic)::Tuple{Mask,ComplexF64}
    tmp = m
    new_mask = zero(Mask)
    phase = COMPLEX_ONE
    @inbounds @fastmath while tmp != 0
        lsb = tmp & -tmp
        idx = trailing_zeros(lsb) + 1
        new_mask |= Mask(1) << (op.perm[idx] - 1)
        phase *= op.perm_phases[idx]
        tmp ⊻= lsb
    end
    return new_mask, phase
end

@inline function apply_operation_to_mask(m::Mask, op::Symmetry_Operation, ::Fermionic)::Tuple{Mask,ComplexF64}
    tmp = m
    new_mask = zero(Mask)
    phase = COMPLEX_ONE
    parity = Int8(1)
    @inbounds @fastmath while tmp != 0
        lsb = tmp & -tmp
        idx = trailing_zeros(lsb) + 1
        p = op.perm[idx]
        if isodd(count_ones(new_mask >> p))
            parity = -parity
        end
        new_mask |= Mask(1) << (p - 1)
        phase *= op.perm_phases[idx]
        tmp ⊻= lsb
    end
    return new_mask, phase * parity
end

"Convenience 2-argument wrapper matching the legacy `apply_operation_to_mask` API (defaults to Bosonic statistics)."
@inline function apply_operation_to_mask(m::Mask, op::Symmetry_Operation)::Tuple{Mask,ComplexF64}
    return apply_operation_to_mask(m, op, Bosonic())
end

# ═══════════════════════════════════════════════════════════════════════════
# 4. Canonical representative (O(|G|) — used only during precomputation)
# ═══════════════════════════════════════════════════════════════════════════

function get_canonical_representative(m::Mask, G::Finite_Symmetry_Group, stats::Particle_Statistics)::Tuple{Mask,Int,ComplexF64}
    repr = m
    best_g = G.identity_idx
    best_amp = COMPLEX_ONE
    @inbounds for (g_idx, g) in enumerate(G.operations)
        shifted, α = apply_operation_to_mask(m, g, stats)
        if shifted < repr
            repr = shifted
            best_g = g_idx
            best_amp = α
        end
    end
    return repr, best_g, best_amp
end

# ═══════════════════════════════════════════════════════════════════════════
# 5. Gosper's hack
# ═══════════════════════════════════════════════════════════════════════════

@inline function _gosper_next(x::Mask)::Mask
    c = x & -x
    r = x + c
    return (((r ⊻ x) >> 2) ÷ c) | r
end

@inline function _first_combination_mask(n_filled::Int)::Mask
    n_filled == 0 && return zero(Mask)
    return (one(Mask) << n_filled) - one(Mask)
end

# ═══════════════════════════════════════════════════════════════════════════
# 6. Symmetry Orbit Catalog
# ═══════════════════════════════════════════════════════════════════════════

mutable struct Symmetry_Orbit_Catalog
    symmetry_group::Finite_Symmetry_Group
    representative_mask_list::Vector{Mask}
    stabilizer_order_list::Vector{Int}
    stabilizer_g_indices_list::Vector{Vector{Int}}
    stabilizer_phases_list::Vector{Vector{ComplexF64}}
end

"""
    update_orbit_stabilizer_phases!(catalog, new_group, statistics)

In-place update of the stabilizer phases when the symmetry group is replaced by
a gauge-covariant (e.g. flux-aware) copy.  The orbit *partition* (representative
masks and stabilizer group-element indices) depends only on the permutation part
of the operations, which is unchanged, so only the accumulated U(1) phases need
recomputation.  This avoids re-running Gosper's hack at every flux point.
"""
function update_orbit_stabilizer_phases!(catalog::Symmetry_Orbit_Catalog,
    new_group::Finite_Symmetry_Group, statistics::Particle_Statistics)
    catalog.symmetry_group == new_group && return catalog  # no-op
    catalog.symmetry_group = new_group
    @inbounds for orbit_idx in eachindex(catalog.representative_mask_list)
        gidxs = catalog.stabilizer_g_indices_list[orbit_idx]
        phases = catalog.stabilizer_phases_list[orbit_idx]
        for (j, gidx) in enumerate(gidxs)
            _, α = apply_operation_to_mask(
                catalog.representative_mask_list[orbit_idx],
                new_group.operations[gidx], statistics)
            phases[j] = α
        end
    end
    return catalog
end

function build_symmetry_orbit_catalog(;
    second_quantized_model::Second_Quantized_Model,
    n_filled::Int,
    symmetry_group::Finite_Symmetry_Group,
    statistics::Particle_Statistics,
)::Symmetry_Orbit_Catalog
    n_site = symmetry_group.n_site
    @assert n_site == second_quantized_model.lattice.n_site

    n_total = binomial(n_site, n_filled)
    n_orbits_est = cld(n_total, group_order(symmetry_group))

    repr_list = Mask[]
    stab_order_list = Int[]
    stab_gidx_list = Vector{Int}[]
    stab_phase_list = Vector{ComplexF64}[]
    sizehint!(repr_list, n_orbits_est)
    sizehint!(stab_order_list, n_orbits_est)
    sizehint!(stab_gidx_list, n_orbits_est)
    sizehint!(stab_phase_list, n_orbits_est)

    seen = Set{Mask}()
    sizehint!(seen, n_total)

    nG = group_order(symmetry_group)
    orbit_masks = Vector{Mask}(undef, nG)
    orbit_amps = Vector{ComplexF64}(undef, nG)

    print("\tBuilding symmetry-orbit catalog (n_filled=$n_filled, |G|=$nG) ... ")

    @assert 0 <= n_filled <= n_site

    res = @timed begin
        if n_filled == 0 || n_filled == n_site
            m = _first_combination_mask(n_filled)
            gidx_stab = Int[]
            phase_stab = ComplexF64[]
            @inbounds for gidx in 1:nG
                shifted, α = apply_operation_to_mask(m, symmetry_group.operations[gidx], statistics)
                shifted == m || continue
                push!(gidx_stab, gidx)
                push!(phase_stab, α)
            end
            push!(repr_list, m)
            push!(stab_order_list, length(gidx_stab))
            push!(stab_gidx_list, gidx_stab)
            push!(stab_phase_list, phase_stab)
        else
            x = _first_combination_mask(n_filled)
            upper = one(Mask) << n_site
            while x < upper
                m = x
                if m in seen
                    x = _gosper_next(x)
                    continue
                end

                min_mask = typemax(Mask)
                @inbounds for gidx in 1:nG
                    shifted, α = apply_operation_to_mask(m, symmetry_group.operations[gidx], statistics)
                    orbit_masks[gidx] = shifted
                    orbit_amps[gidx] = α
                    shifted < min_mask && (min_mask = shifted)
                end
                @assert min_mask == m "Gosper ordering violation"

                @inbounds for shifted in orbit_masks
                    push!(seen, shifted)
                end

                gidx_stab = Int[]
                phase_stab = ComplexF64[]
                @inbounds for gidx in 1:nG
                    if orbit_masks[gidx] == m
                        push!(gidx_stab, gidx)
                        push!(phase_stab, orbit_amps[gidx])
                    end
                end

                push!(repr_list, m)
                push!(stab_order_list, length(gidx_stab))
                push!(stab_gidx_list, gidx_stab)
                push!(stab_phase_list, phase_stab)

                x = _gosper_next(x)
            end
        end
    end

    n_orbits = length(repr_list)
    printstyled("Done. $(n_orbits) orbits (reduction $(round(n_orbits/n_total*100, digits=1))%).  t=$(round(res.time, digits=3))s\n", bold=true)
    return Symmetry_Orbit_Catalog(symmetry_group, repr_list, stab_order_list, stab_gidx_list, stab_phase_list)
end

# ═══════════════════════════════════════════════════════════════════════════
# 7. 1D Irreps
# ═══════════════════════════════════════════════════════════════════════════

"""
Struct `OneDim_Irrep{Irrep_Label}` for _One Single_ One-dimensional Irreducible Representation (Irrep)
---
labeled by a generic typed `label::Irrep_Label`.
- Fields:
    - `label::Irrep_Label`: the `Irrep_Label`-typed label of the irrep, can be integer, tuples etc., depending on the symmetry transformation we are considering
    - `values::Vector{ComplexF64}`: the irrep values for each group element, which is a vector of length equal to the order of the group, with each element being a complex number representing the irrep value χ(g) for the g-th group operation
"""
struct OneDim_Irrep{Irrep_Label}
    label::Irrep_Label
    values::Vector{ComplexF64}   # χ(g) for each g ∈ G
end
"constructor for `OneDim_Irrep`"
OneDim_Irrep(label::Irrep_Label, values::Vector{T}) where {Irrep_Label,T<:Number} = OneDim_Irrep{Irrep_Label}(label, ComplexF64.(collect(values)))

# ═══════════════════════════════════════════════════════════════════════════
# 8. Symmetry Sector Basis
# ═══════════════════════════════════════════════════════════════════════════

"""
Struct `Symmetry_Sector_Basis` to Store the Symmetry Sector Basis for a Given 1D Irrep
---
- Fields:
    - `irrep::OneDim_Irrep`: the 1D irrep for which this symmetry sector basis is constructed
    - `symmetry_group::Finite_Symmetry_Group`: the symmetry group under consideration
    - `representative_mask_list::Vector{Mask}`: the list of representative masks for each orbit that belongs to this irrep/sector, inherited from the `Symmetry_Orbit_Catalog`
    - `stabilizer_group_order_list::Vector{Int}`: the list of stabilizer group orders for each representative mask in this sector, inherited from the `Symmetry_Orbit_Catalog`
    - `representative_mask_to_mask_idx_map::Dict{Mask,Int}`: `Hashmap<representative_mask, idx in representative_mask_list>`
"""
struct Symmetry_Sector_Basis
    irrep::OneDim_Irrep
    symmetry_group::Finite_Symmetry_Group
    representative_mask_list::Vector{Mask}
    stabilizer_order_list::Vector{Int}
    representative_mask_to_mask_idx_map::Dict{Mask,Int}
end

@inline function _basis_index(basis::Symmetry_Sector_Basis, repr::Mask)::Int
    return get(basis.representative_mask_to_mask_idx_map, repr, 0)
end

@inline function _is_orbit_compatible(catalog::Symmetry_Orbit_Catalog, orbit_idx::Int,
    irrep::OneDim_Irrep; atol::Float64=1e-12)::Bool
    gidxs = catalog.stabilizer_g_indices_list[orbit_idx]
    phases = catalog.stabilizer_phases_list[orbit_idx]
    @inbounds for j in eachindex(gidxs)
        !isapprox(irrep.values[gidxs[j]], phases[j]; atol=atol) && return false
    end
    return true
end

function build_symmetry_sector_basis(catalog::Symmetry_Orbit_Catalog, irrep::OneDim_Irrep;
    atol::Float64=1e-12)::Symmetry_Sector_Basis
    repr_list = Mask[]
    stab_order_list = Int[]
    for i in eachindex(catalog.representative_mask_list)
        if _is_orbit_compatible(catalog, i, irrep; atol=atol)
            push!(repr_list, catalog.representative_mask_list[i])
            push!(stab_order_list, catalog.stabilizer_order_list[i])
        end
    end
    representative_mask_to_mask_idx_map = Dict{Mask,Int}(m => idx for (idx, m) in enumerate(repr_list))
    return Symmetry_Sector_Basis(irrep, catalog.symmetry_group, repr_list, stab_order_list, representative_mask_to_mask_idx_map)
end

# ═══════════════════════════════════════════════════════════════════════════
# 9. CanonicalMap — O(1) canonical-representative lookup (all modes)
# ═══════════════════════════════════════════════════════════════════════════

"""
    CanonicalMap

Lazily-populated cache: scattered_mask → (repr_mask, g_idx, α_g).

Used uniformly across ALL ED modes (matrix, distributed-matrix, and matrix-free)
to avoid repeated O(|G|) canonicalization of the same scattered masks.

- **Matrix mode**: rows are built single-threaded; the cache warms naturally
  via `get!`, providing 3–10× speedup over uncached canonicalization.
- **Distributed matrix mode**: each worker maintains its own cache copy;
  columns are processed single-threaded per worker — no synchronisation needed.
- **Matrix-free mode**: MUST be pre-populated via `populate_canonical_map!`
  before the multi-threaded Lanczos loop; after population all `get_canonical`
  calls are read-only Dict lookups, which are thread-safe without locks.
"""
struct CanonicalMap
    symmetry_group::Finite_Symmetry_Group
    statistics::Particle_Statistics
    cache::Dict{Mask,Tuple{Mask,Int,ComplexF64}}   # scattered mask → canonical data
end

"O(1) canonical-representative lookup (falls back to O(|G|) on cache miss)"
@inline function get_canonical(cmap::CanonicalMap, m::Mask)::Tuple{Mask,Int,ComplexF64}
    return get!(cmap.cache, m) do
        get_canonical_representative(m, cmap.symmetry_group, cmap.statistics)
    end
end

"""
    populate_canonical_map!(cmap, basis, bilinear_terms)

Single-threaded cache warmup: iterates every representative mask × every valid
hopping move, calling `get_canonical` to populate the Dict.

**Required before any multi-threaded use** (e.g., `apply_hamiltonian!` with
`Threads.@threads`) because `get!` on a shared Dict is NOT thread-safe during
concurrent writes.  After this call, all subsequent `get_canonical` lookups
are pure Dict reads, which are safe across threads.

For single-threaded matrix construction, this call is optional (the cache
warms naturally during the column loop) but harmless.
"""
function populate_canonical_map!(cmap::CanonicalMap, basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}})
    @inbounds for repr_mask in basis.representative_mask_list
        for (i_from, i_to, _) in bilinear_terms
            if is_site_occupied(repr_mask, i_from) && is_site_empty(repr_mask, i_to)
                new_mask = empty_site_for_mask(repr_mask, i_from)
                new_mask = occupy_site_for_mask(new_mask, i_to)
                get_canonical(cmap, new_mask)  # populates cache if not present
            end
        end
    end
    return cmap
end

# ═══════════════════════════════════════════════════════════════════════════
# 10. project_to_sector — unified irrep projection via CanonicalMap
# ═══════════════════════════════════════════════════════════════════════════

"""
    project_to_sector(m, basis, cmap) -> Union{Nothing, Tuple{Int, ComplexF64}}

Project a scattered Fock mask `m` onto the symmetry-sector basis.

1. Obtain canonical representative `(repr, g_idx, α)` via `CanonicalMap`
   (O(1) cached, O(|G|) on first encounter).
2. Look up the representative in the basis Dict.
3. Return `(row_index, α · χ(g_idx)ˣ)` or `nothing` if the orbit does not
   belong to this irrep sector.

Used uniformly by matrix, distributed-matrix, and matrix-free modes.
"""
@inline function project_to_sector(m::Mask, basis::Symmetry_Sector_Basis, cmap::CanonicalMap)
    repr, g_idx, α = get_canonical(cmap, m)
    idx = _basis_index(basis, repr)
    idx == 0 && return nothing
    return (idx, α * conj(basis.irrep.values[g_idx]))
end

"Legacy wrapper: project a raw mask (not necessarily a canonical representative) to the unnormalized symmetry sector basis. Returns `(repr_idx, coeff)` or `nothing`."
@inline function project_to_unnormalized_sector(m::Mask, basis::Symmetry_Sector_Basis, stats=Bosonic())
    repr, g_idx, α = get_canonical_representative(m, basis.symmetry_group, stats)
    idx = _basis_index(basis, repr)
    idx == 0 && return nothing
    return (idx, α * conj(basis.irrep.values[g_idx]))
end

# ═══════════════════════════════════════════════════════════════════════════
# 11. Matrix-free Hamiltonian application — the core Lanczos operation
# ═══════════════════════════════════════════════════════════════════════════

@inline hopping_phase_for_stats(::Bosonic, ::Mask, ::Int, ::Int)::ComplexF64 = COMPLEX_ONE

@inline function hopping_phase_for_stats(::Fermionic, m::Mask, i_from::Int, i_to::Int)::ComplexF64
    i_from == i_to && return COMPLEX_ONE
    lo = min(i_from, i_to)
    hi = max(i_from, i_to)
    between = hi - lo <= 1 ? zero(Mask) : (((one(Mask) << (hi - lo - 1)) - one(Mask)) << lo)
    return isodd(count_ones(m & between)) ? -COMPLEX_ONE : COMPLEX_ONE
end

function _matrixfree_buffers(n::Int)
    return [zeros(ComplexF64, n) for _ in 1:Threads.maxthreadid()]
end

struct MatrixFreeHamiltonian
    basis::Symmetry_Sector_Basis
    bilinear_terms::Vector{Tuple{Int,Int,ComplexF64}}
    density_terms::Vector{Tuple{Int,Int,ComplexF64}}
    cmap::CanonicalMap
    y_threads::Vector{Vector{ComplexF64}}
end

function MatrixFreeHamiltonian(
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap,
)
    n = length(basis.representative_mask_list)
    bilin = [(i, j, ComplexF64(t)) for (i, j, t) in bilinear_terms]
    density = [(i, j, ComplexF64(v)) for (i, j, v) in density_terms]
    return MatrixFreeHamiltonian(basis, bilin, density, cmap, _matrixfree_buffers(n))
end

"""
    apply_hamiltonian!(y, x, basis, bilinear_terms, density_terms, cmap)

Compute y = H·x on-the-fly.  No sparse matrix is ever stored.

Uses `Threads.@threads` with per-thread accumulation buffers for race-free
shared-memory parallelism (equivalent to XDiag's OpenMP backend).

The `cmap::CanonicalMap` must be pre-populated via `populate_canonical_map!`
before the first call.
"""
function apply_hamiltonian!(y::Vector{ComplexF64}, x::Vector{ComplexF64},
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap)
    y_threads = _matrixfree_buffers(length(x))
    return apply_hamiltonian!(y, x, basis, bilinear_terms, density_terms, cmap, y_threads)
end

function apply_hamiltonian!(y::Vector{ComplexF64}, x::Vector{ComplexF64},
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap,
    y_threads::Vector{Vector{ComplexF64}})
    n = length(x)
    @assert length(y) == n

    max_tid = Threads.maxthreadid()
    length(y_threads) < max_tid && error("not enough thread-local buffers")
    @inbounds for tid in 1:max_tid
        fill!(y_threads[tid], 0)
    end

    Threads.@threads :static for col in 1:n
        tid = Threads.threadid()
        yt = y_threads[tid]
        x_col = x[col]
        x_col == 0 && continue

        repr_mask = basis.representative_mask_list[col]
        stab_col = basis.stabilizer_order_list[col]
        inv_sqrt_stab_col = 1.0 / sqrt(stab_col)

        # ── Diagonal (density-density interactions) ──
        H_diag = zero(ComplexF64)
        for (i, j, V) in density_terms
            if is_site_occupied(repr_mask, i) && is_site_occupied(repr_mask, j)
                H_diag += V
            end
        end
        yt[col] += H_diag * x_col

        # ── Off-diagonal (hopping) ──
        for (i_from, i_to, t) in bilinear_terms
            if is_site_occupied(repr_mask, i_from) && is_site_empty(repr_mask, i_to)
                new_mask = empty_site_for_mask(repr_mask, i_from)
                new_mask = occupy_site_for_mask(new_mask, i_to)

                proj = project_to_sector(new_mask, basis, cmap)
                proj === nothing && continue

                row, coeff = proj
                stab_row = basis.stabilizer_order_list[row]
                hop_phase = hopping_phase_for_stats(cmap.statistics, repr_mask, i_from, i_to)
                # H_{row,col} = t × coeff × √(|Stab(row)| / |Stab(col)|)
                H_elem = t * hop_phase * coeff * sqrt(stab_row) * inv_sqrt_stab_col
                yt[row] += H_elem * x_col
            end
        end
    end

    # ── Reduce thread-local buffers → y ──
    fill!(y, 0)
    @inbounds for tid in 1:max_tid
        yt = y_threads[tid]
        for i in 1:n
            y[i] += yt[i]
        end
    end
    return y
end

function (H::MatrixFreeHamiltonian)(x::Vector{ComplexF64})
    y = similar(x)
    apply_hamiltonian!(y, x, H.basis, H.bilinear_terms, H.density_terms, H.cmap, H.y_threads)
    return y
end

"""
    hamiltonian_linear_operator(basis, bilinear_terms, density_terms, cmap)

Return a closure `H_op(x) -> y` suitable for iterative eigensolvers (KrylovKit, Arpack).
"""
function hamiltonian_linear_operator(basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap)
    n = length(basis.representative_mask_list)
    return MatrixFreeHamiltonian(basis, bilinear_terms, density_terms, cmap), n
end

# ═══════════════════════════════════════════════════════════════════════════
# 12. Sparse-matrix construction (matrix mode)
# ═══════════════════════════════════════════════════════════════════════════

function build_ed_Hamiltonian_symmetry_block(
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap,
)::SparseMatrixCSC{ComplexF64,Int}
    sector_dim = length(basis.representative_mask_list)
    print("\tBuilding H block (matrix mode) @ irrep $(basis.irrep.label) (dim=$sector_dim) ... ")

    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    est_nnz = 1 + 4 * basis.symmetry_group.n_site
    sizehint!(Is, sector_dim * est_nnz)
    sizehint!(Js, sector_dim * est_nnz)
    sizehint!(Vs, sector_dim * est_nnz)

    res = @timed begin
        @inbounds for (col, repr_mask) in enumerate(basis.representative_mask_list)
            stab_col = basis.stabilizer_order_list[col]

            # Diagonal
            H_diag = zero(ComplexF64)
            for (i, j, V) in density_terms
                if is_site_occupied(repr_mask, i) && is_site_occupied(repr_mask, j)
                    H_diag += V
                end
            end
            push!(Is, col)
            push!(Js, col)
            push!(Vs, H_diag)

            # Off-diagonal
            for (i_from, i_to, t) in bilinear_terms
                if is_site_occupied(repr_mask, i_from) && is_site_empty(repr_mask, i_to)
                    new_mask = empty_site_for_mask(repr_mask, i_from)
                    new_mask = occupy_site_for_mask(new_mask, i_to)
                    proj = project_to_sector(new_mask, basis, cmap)
                    proj === nothing && continue
                    row, coeff = proj
                    stab_row = basis.stabilizer_order_list[row]
                    hop_phase = hopping_phase_for_stats(cmap.statistics, repr_mask, i_from, i_to)
                    H_elem = t * hop_phase * coeff * sqrt(stab_row / stab_col)
                    push!(Is, row)
                    push!(Js, col)
                    push!(Vs, H_elem)
                end
            end
        end
    end

    H = sparse(Is, Js, Vs, sector_dim, sector_dim)
    dropzeros!(H)
    nnz_H = nnz(H)
    sparsity = sector_dim == 0 ? 0.0 : nnz_H / (sector_dim^2)
    printstyled("Done. nnz=$nnz_H, sparsity=$(round(sparsity,digits=6)). t=$(round(res.time,digits=3))s\n", bold=true)
    return H
end

# ═══════════════════════════════════════════════════════════════════════════
# 13. Diagonalization helpers
# ═══════════════════════════════════════════════════════════════════════════

function diagonalize_block_dense(H::SparseMatrixCSC{ComplexF64,Int}; nev::Int=5)
    n = size(H, 1)
    n == 0 && return Float64[], Matrix{ComplexF64}(undef, 0, 0)
    H_dense = Hermitian(Matrix(H))
    vals, vecs = eigen(H_dense)
    k = min(nev, length(vals))
    return real.(vals[1:k]), vecs[:, 1:k]
end

function diagonalize_block_arpack(H::SparseMatrixCSC{ComplexF64,Int}; nev::Int=5)
    n = size(H, 1)
    n == 0 && return Float64[], Matrix{ComplexF64}(undef, 0, 0)
    # Temporarily set BLAS threads = nprocs() for diagonalization (HPC: one per worker)
    blas_threads = max(1, nprocs())
    BLAS.set_num_threads(blas_threads)
    try
        if n <= 500
            return diagonalize_block_dense(H; nev=nev)
        else
            vals, vecs, _ = Arpack.eigs(H; nev=nev, which=:SR)
            return real.(vals), Matrix(vecs)
        end
    finally
        BLAS.set_num_threads(1)  # restore HPC default
    end
end

"Diagonalize using matrix-free Lanczos (KrylovKit)"
function diagonalize_block_matrixfree(H_op!, n::Int; nev::Int=5)
    n == 0 && return Float64[], Matrix{ComplexF64}(undef, 0, 0)
    # Temporarily set BLAS threads = nprocs() for diagonalization
    blas_threads = max(1, nprocs())
    BLAS.set_num_threads(blas_threads)
    try
        if n <= 500
            # Build dense matrix for small sectors
            H_dense = Matrix{ComplexF64}(undef, n, n)
            x = zeros(ComplexF64, n)
            for j in 1:n
                x[j] = 1.0
                y = H_op!(x)
                H_dense[:, j] .= y
                x[j] = 0.0
            end
            H_herm = Hermitian(H_dense)
            vals, vecs = eigen(H_herm)
            k = min(nev, length(vals))
            return real.(vals[1:k]), vecs[:, 1:k]
        else
            # Use KrylovKit for large sectors
            x0 = randn(ComplexF64, n)
            x0 ./= norm(x0)
            vals, vecs, info = KrylovKit.eigsolve(H_op!, x0, nev, :SR; tol=1e-10, maxiter=300)
            return real.(vals), hcat(vecs...)
        end
    finally
        BLAS.set_num_threads(1)  # restore HPC default
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# 14. High-level ED data structure
# ═══════════════════════════════════════════════════════════════════════════
abstract type ED_Data end


"""
Struct `Symmetry_Resolved_ED_Data` to Store All Data for Symmetry-resolved ED
---
- Fields:
    - `second_quantized_model::Real_Space_Second_Quantized_Model`: the second quantized model containing the lattice and the Hamiltonian parameters
    - `n_filled::Int`: the number of filled particles (set bits in the occupation basis)
    - `filling_fraction::Rational{Int}`: the filling fraction, defined as `n_filled / n_site`
    - `symmetry_group::Finite_Symmetry_Group`: the symmetry group used for symmetry resolution
    - `irrep_list::Vector{<:OneDim_Irrep}`: the list of 1D irreps for which we want to resolve the ED
    - `orbit_catalog::Symmetry_Orbit_Catalog`: the catalog of symmetry orbits for the given model and symmetry
    - `sector_dims::Vector{Int}`: the dimensions of each symmetry sector corresponding to each irrep, which is equal to the number of representative masks in the sector basis for that irrep
    - `ed_scan_res::Dict{Int, Tuple{Vector{Float64},Matrix{ComplexF64}}}`: a hashmap `Dict<irrep_idx, (eigvals, eigvecs)>` to store the ED results for each symmetry sector, where `irrep_idx` is the index of the irrep in `irrep_list`, and `(eigvals, eigvecs)` is the tuple of eigenvalues and eigenvectors obtained from the ED scan for that sector
"""
mutable struct Symmetry_Resolved_ED_Data <: ED_Data
    second_quantized_model::Real_Space_Second_Quantized_Model
    n_filled::Int
    filling_fraction::Rational{Int}
    symmetry_group::Finite_Symmetry_Group
    irrep_list::Vector{OneDim_Irrep}
    orbit_catalog::Symmetry_Orbit_Catalog
    sector_dims::Vector{Int}
    ed_scan_res::Dict{Int,Tuple{Vector{Float64},Matrix{ComplexF64}}}
end

# ═══════════════════════════════════════════════════════════════════════════
# 15. ED scan — matrix mode (stores sparse H)
# ═══════════════════════════════════════════════════════════════════════════

function ed_scan_at_irrep_matrix!(irrep_label, ed_data::Symmetry_Resolved_ED_Data; nev::Int=5)
    irrep_idx = findfirst(irrep -> irrep.label == irrep_label, ed_data.irrep_list)
    @assert irrep_idx !== nothing

    irrep = ed_data.irrep_list[irrep_idx]
    statistics = ed_data.second_quantized_model.statistics
    basis = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
    ed_data.sector_dims[irrep_idx] = length(basis.representative_mask_list)

    # ── Create CanonicalMap (shared across all matrix-construction paths) ──
    cmap = CanonicalMap(ed_data.symmetry_group, statistics, Dict{Mask,Tuple{Mask,Int,ComplexF64}}())

    # Use distributed construction when workers are available (HPC)
    if nprocs() > 1 && ed_data.sector_dims[irrep_idx] > 500
        H = build_ed_Hamiltonian_symmetry_block_distributed(basis, ed_data.second_quantized_model.bilinear_terms,
            ed_data.second_quantized_model.density_density_terms, cmap)
    else
        H = build_ed_Hamiltonian_symmetry_block(basis, ed_data.second_quantized_model.bilinear_terms,
            ed_data.second_quantized_model.density_density_terms, cmap)
    end
    vals, vecs = diagonalize_block_arpack(H; nev=nev)
    ed_data.ed_scan_res[irrep_idx] = (vals, vecs)
    H = nothing
    GC.gc(true)  # force GC after each sector to avoid memory pressure accumulation
    return vals, vecs
end

# ═══════════════════════════════════════════════════════════════════════════
# 16. ED scan — matrix-free mode (on-the-fly Lanczos, XDiag-style)
# ═══════════════════════════════════════════════════════════════════════════

function ed_scan_at_irrep_matrixfree!(irrep_label, ed_data::Symmetry_Resolved_ED_Data;
    nev::Int=5, use_distributed::Bool=false)
    irrep_idx = findfirst(irrep -> irrep.label == irrep_label, ed_data.irrep_list)
    @assert irrep_idx !== nothing

    irrep = ed_data.irrep_list[irrep_idx]
    statistics = ed_data.second_quantized_model.statistics
    basis = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
    sector_dim = length(basis.representative_mask_list)
    ed_data.sector_dims[irrep_idx] = sector_dim

    if sector_dim == 0
        ed_data.ed_scan_res[irrep_idx] = (Float64[], Matrix{ComplexF64}(undef, 0, 0))
        return Float64[], Matrix{ComplexF64}(undef, 0, 0)
    end

    print("\tMatrix-free mode @ irrep $(irrep.label) (dim=$sector_dim, threads=$(Threads.nthreads())) ... ")

    bilinear = ed_data.second_quantized_model.bilinear_terms
    density = ed_data.second_quantized_model.density_density_terms

    res = @timed begin
        # Build and populate CanonicalMap
        cmap = CanonicalMap(ed_data.symmetry_group, statistics, Dict{Mask,Tuple{Mask,Int,ComplexF64}}())
        populate_canonical_map!(cmap, basis, bilinear)

        # Create the linear operator. The distributed variant avoids @everywhere;
        # workers are loaded through ensure_distributed_workers_loaded!().
        H_op, n = if use_distributed && nprocs() > 1
            hamiltonian_linear_operator_distributed(basis, bilinear, density, cmap)
        else
            hamiltonian_linear_operator(basis, bilinear, density, cmap)
        end
        vals, vecs = diagonalize_block_matrixfree(H_op, n; nev=nev)
    end
    printstyled("Done. t=$(round(res.time,digits=3))s\n", bold=true)

    ed_data.ed_scan_res[irrep_idx] = (vals, vecs)
    cmap = nothing
    GC.gc(true)  # force GC after each sector (CanonicalMap can be large)
    return vals, vecs
end

# ═══════════════════════════════════════════════════════════════════════════
# 17. Distributed matrix construction (HPC: pmap across workers)
# ═══════════════════════════════════════════════════════════════════════════

"""
    build_ed_Hamiltonian_symmetry_block_distributed(basis, bilinear, density, stats)
        -> SparseMatrixCSC

Build the Hamiltonian matrix using `Distributed.pmap`: each worker processes a
chunk of the representative masks in parallel.  Used on HPC clusters where
the matrix is too large to build on a single node.
"""
function build_ed_Hamiltonian_symmetry_block_distributed(
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap,
)::SparseMatrixCSC{ComplexF64,Int}
    sector_dim = length(basis.representative_mask_list)
    if sector_dim == 0
        return spzeros(ComplexF64, Int, 0, 0)
    end
    if nprocs() == 1
        return build_ed_Hamiltonian_symmetry_block(basis, bilinear_terms, density_terms, cmap)
    end
    ensure_distributed_workers_loaded!()

    nw = nworkers()
    chunk_size = cld(sector_dim, nw)
    col_ranges = [r[1]:r[end] for r in Iterators.partition(1:sector_dim, chunk_size)]
    pool = CachingPool(workers())

    print("\tBuilding H block (distributed, $nw workers) @ irrep $(basis.irrep.label) (dim=$sector_dim) ... ")

    res = @timed begin
        results = pmap(pool, col_ranges) do rng
            Is = Int[]
            Js = Int[]
            Vs = ComplexF64[]
            est_nnz = 1 + 4 * basis.symmetry_group.n_site
            sizehint!(Is, length(rng) * est_nnz)
            sizehint!(Js, length(rng) * est_nnz)
            sizehint!(Vs, length(rng) * est_nnz)
            @inbounds for col in rng
                repr_mask = basis.representative_mask_list[col]
                stab_col = basis.stabilizer_order_list[col]
                H_diag = zero(ComplexF64)
                for (i, j, V) in density_terms
                    if is_site_occupied(repr_mask, i) && is_site_occupied(repr_mask, j)
                        H_diag += V
                    end
                end
                push!(Is, col)
                push!(Js, col)
                push!(Vs, H_diag)
                for (i_from, i_to, t) in bilinear_terms
                    if is_site_occupied(repr_mask, i_from) && is_site_empty(repr_mask, i_to)
                        new_mask = empty_site_for_mask(repr_mask, i_from)
                        new_mask = occupy_site_for_mask(new_mask, i_to)
                        proj = project_to_sector(new_mask, basis, cmap)
                        proj === nothing && continue
                        row, coeff = proj
                        stab_row = basis.stabilizer_order_list[row]
                        hop_phase = hopping_phase_for_stats(cmap.statistics, repr_mask, i_from, i_to)
                        H_elem = t * hop_phase * coeff * sqrt(stab_row / stab_col)
                        push!(Is, row)
                        push!(Js, col)
                        push!(Vs, H_elem)
                    end
                end
            end
            return (Is, Js, Vs)
        end
    end

    Is = reduce(vcat, getindex.(results, 1))
    Js = reduce(vcat, getindex.(results, 2))
    Vs = reduce(vcat, getindex.(results, 3))
    H = sparse(Is, Js, Vs, sector_dim, sector_dim)
    dropzeros!(H)
    GC.gc(true)

    nnz_H = nnz(H)
    sparsity = sector_dim == 0 ? 0.0 : nnz_H / (sector_dim^2)
    printstyled("Done. nnz=$nnz_H, sparsity=$(round(sparsity,digits=6)). t=$(round(res.time,digits=3))s\n", bold=true)
    return H
end

# ═══════════════════════════════════════════════════════════════════════════
# 17c. Distributed matrix-free H|ψ⟩ (HPC: pmap column chunks)
# ═══════════════════════════════════════════════════════════════════════════

"""
    apply_hamiltonian_distributed!(y, x, basis, bilinear, density, cmap)

Distributed variant of `apply_hamiltonian!`: each worker computes its chunk
of columns, partial results are reduced on the master.

The `cmap` is broadcast to all workers so they can look up canonical
representatives without network round-trips.
"""
function apply_hamiltonian_distributed!(y::Vector{ComplexF64}, x::Vector{ComplexF64},
    basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap)
    n = length(x)
    if nprocs() == 1
        return apply_hamiltonian!(y, x, basis, bilinear_terms, density_terms, cmap)
    end
    nw = nworkers()
    ensure_distributed_workers_loaded!()

    chunk_size = cld(n, nw)
    col_ranges = [r[1]:r[end] for r in Iterators.partition(1:n, chunk_size)]
    pool = CachingPool(workers())

    results = pmap(pool, col_ranges) do rng
        y_part = zeros(ComplexF64, n)
        @inbounds for col in rng
            x_col = x[col]
            x_col == 0 && continue
            repr_mask = basis.representative_mask_list[col]
            stab_col = basis.stabilizer_order_list[col]
            inv_sqrt_stab_col = 1.0 / sqrt(stab_col)

            H_diag = zero(ComplexF64)
            for (i, j, V) in density_terms
                if is_site_occupied(repr_mask, i) && is_site_occupied(repr_mask, j)
                    H_diag += V
                end
            end
            y_part[col] += H_diag * x_col

            for (i_from, i_to, t) in bilinear_terms
                if is_site_occupied(repr_mask, i_from) && is_site_empty(repr_mask, i_to)
                    new_mask = empty_site_for_mask(repr_mask, i_from)
                    new_mask = occupy_site_for_mask(new_mask, i_to)
                    proj = project_to_sector(new_mask, basis, cmap)
                    proj === nothing && continue
                    row, coeff = proj
                    stab_row = basis.stabilizer_order_list[row]
                    hop_phase = hopping_phase_for_stats(cmap.statistics, repr_mask, i_from, i_to)
                    H_elem = t * hop_phase * coeff * sqrt(stab_row) * inv_sqrt_stab_col
                    y_part[row] += H_elem * x_col
                end
            end
        end
        return y_part
    end

    fill!(y, 0)
    for yp in results
        y .+= yp
    end
    return y
end

"""
    hamiltonian_linear_operator_distributed(basis, bilinear, density, cmap)

Return a distributed linear operator `H_op(x) -> y` using `pmap` column chunks.
"""
function hamiltonian_linear_operator_distributed(basis::Symmetry_Sector_Basis,
    bilinear_terms::Vector{<:Tuple{Int,Int,<:Number}},
    density_terms::Vector{<:Tuple{Int,Int,<:Number}},
    cmap::CanonicalMap)
    n = length(basis.representative_mask_list)
    function H_op(x::Vector{ComplexF64})
        y = similar(x)
        fill!(y, 0)
        apply_hamiltonian_distributed!(y, x, basis, bilinear_terms, density_terms, cmap)
        return y
    end
    return H_op, n
end

# ═══════════════════════════════════════════════════════════════════════════
# 17d. Checkpoint support (HPC preemption resilience)
# ═══════════════════════════════════════════════════════════════════════════

"Save ED data to a checkpoint file using JLD2"
function save_checkpoint(ed_data::Symmetry_Resolved_ED_Data, path::String)::Nothing
    mkpath(dirname(path))
    @save path ed_data
    return nothing
end

"Load ED data from a checkpoint file"
function load_checkpoint(path::String)::Symmetry_Resolved_ED_Data
    @load path ed_data
    return ed_data
end

# ═══════════════════════════════════════════════════════════════════════════
# 17e. Flux-insertion checkpoint — standard filename + scan over twisted phases
# ═══════════════════════════════════════════════════════════════════════════

"""
Generate the Universal Checkpoint Filename for ED Scan
---
for both conventional scan and flux-insertion scan. The format reads `{tb_model.model_name}_{sample_size}_ν_graph={num}_{den}_twisted_phases_over_2π_{twisted_phases_over_2π}_params={params_short}.jld2`.

Here `params_short` just round ALL values of `params` to 3 digits.
"""
function ed_scan_checkpoint_filename(
    model::Real_Space_Second_Quantized_Model,
    twisted_phases_over_2π::AbstractVector{<:Real},
    filling_fraction::Rational{Int},
)::String
    params_short = Dict{keytype(model.params),valtype(model.params)}()
    for (k, v) in model.params
        params_short[k] = round(v; digits=3)
    end

    return "$(model.tb_model.model_name)_$(model.lattice.sample_size)_ν_graph=$(numerator(filling_fraction))_$(denominator(filling_fraction))_twisted_phases_over_2π=$(twisted_phases_over_2π)_params=$(params_short).jld2"
end

"_Internal_: scan sectors of `ed_data` for optionally provided `scanned_sectors`"
function _ed_scan_sectors!(ed_data::Symmetry_Resolved_ED_Data;
    nev::Int=5,
    mode::Symbol=:matrix,
    checkpoint_path::Union{String,Nothing}=nothing,
    use_distributed::Bool=true,
    scanned_sectors::Union{Nothing,Vector{<:Tuple{Int,Int}}}=nothing,
    overwrite::Bool=false,
)::Union{String,Nothing}
    # Resume from existing checkpoint if not overwriting
    if checkpoint_path !== nothing && isfile(checkpoint_path) && !overwrite
        @info "Loading existing checkpoint: $(checkpoint_path)"
        loaded = load_checkpoint(checkpoint_path)
        ed_data.ed_scan_res = loaded.ed_scan_res
        return checkpoint_path
    end
    n_total = length(ed_data.irrep_list)
    n_done = length(ed_data.ed_scan_res)
    n_done > 0 && println("[ED scan] $(n_done)/$(n_total) sectors already computed; resuming.")

    for (irrep_idx, irrep) in enumerate(ed_data.irrep_list)
        haskey(ed_data.ed_scan_res, irrep_idx) && continue
        # If sector filter is provided, skip sectors not in the list
        if scanned_sectors !== nothing && !(irrep.label in scanned_sectors)
            continue
        end
        println("[ED scan] Sector $(irrep_idx)/$(n_total) — irrep $(irrep.label)  [mode=$mode]")
        flush(stdout)
        if mode == :matrixfree
            ed_scan_at_irrep_matrixfree!(irrep.label, ed_data; nev=nev, use_distributed=use_distributed)
        elseif mode == :matrix
            ed_scan_at_irrep_matrix!(irrep.label, ed_data; nev=nev)
        else
            error("Unknown ED scan mode: $mode. Use :matrix or :matrixfree.")
        end

        if checkpoint_path !== nothing
            save_checkpoint(ed_data, checkpoint_path)
            println("    [checkpoint → $(checkpoint_path)]")
        end
    end
    return checkpoint_path
end

"""
Unified API for ED Scan 
---
with support of flux-insertion scan and checkpoint resume.
- Args:
    - `ed_data::ED_Data`: the ED data structure.
- Named Args:
    - `nev::Int=5`: number of eigenvalues/eigenvectors to compute for each sector.
    - `mode::Symbol=:matrix`: the ED mode, can be `:matrix`
    - `use_distributed::Bool=true`: enable distributed computation via `pmap`.
    - `scanned_sectors::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing`: filter for specific sector labels.
    - `checkpoint_path::Union{String,Nothing}=nothing`: path to the checkpoint file to resume from. If a checkpoint file exists, resume scanning from the last saved state.
    - `flux_direction::Int=1`: for flux-insertion scans, the direction of the twisted boundary condition
    - `twisted_phases_over_2π_list::Union{Nothing,Vector{Float64}}=nothing`: the list of twisted phase (flux) values [in units of 2π] to scan over. If `nothing`, fallback to conventional ED scan (at fixed twisted boundary conditions)
Returns:
    - `Union{Nothing, Vector{String}}`:
        - If `twisted_phases_over_2π_list === nothing`, returns the (single-element) checkpoint paths (if provided) or `nothing` (if not).
        - If `twisted_phases_over_2π_list` is provided, returns the list of checkpoint paths saved.
"""
function ed_scan!(ed_data::Symmetry_Resolved_ED_Data;
    nev::Int=5,
    mode::Symbol=:matrix,
    use_distributed::Bool=true,
    scanned_sectors::Union{Nothing,Vector{<:Tuple{Int,Int}}}=nothing, # sector-specific scan support
    # checkpoint support
    checkpoint_path::Union{String,Nothing}=nothing,
    checkpoint_dir::String="checkpoints",
    overwrite::Bool=false,
    # ── flux-scan extensions ──
    flux_direction::Int=1,
    twisted_phases_over_2π_list::Union{Nothing,Vector{Float64}}=nothing,
)::Union{Nothing,Vector{String}}
    @assert mode in [:matrix, :matrixfree]

    if twisted_phases_over_2π_list === nothing # fallback to conventional ED scan
        ckpt_path = _ed_scan_sectors!(ed_data; nev, mode, checkpoint_path, use_distributed, scanned_sectors, overwrite)
        if isnothing(ckpt_path)
            return nothing
        else
            return [ckpt_path]
        end
    end

    # ── Flux-scan mode ──
    model = ed_data.second_quantized_model
    filling_fraction = ed_data.filling_fraction
    mkpath(checkpoint_dir)
    dim = model.lattice.dim

    checkpoint_paths = Vector{String}()
    for θ_val in twisted_phases_over_2π_list
        flux = zeros(Float64, dim)
        flux[flux_direction] = θ_val

        ckpt_name = ed_scan_checkpoint_filename(model, flux, filling_fraction)
        ckpt_path = joinpath(checkpoint_dir, ckpt_name)

        if isfile(ckpt_path) && !overwrite
            @info "Skipping θ=$θ_val — checkpoint exists: $ckpt_name"
            push!(checkpoint_paths, ckpt_path)
            continue
        end

        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux)
        active_group = build_translation_group(model.lattice, flux)
        θ_ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=active_group)
        _ed_scan_sectors!(θ_ed_data; nev, mode, use_distributed, scanned_sectors, overwrite)
        save_checkpoint(θ_ed_data, ckpt_path)
        @info "Saved ED scan @ θ=$θ_val → $ckpt_name"
        push!(checkpoint_paths, ckpt_path)
    end

    return checkpoint_paths
end

# ═══════════════════════════════════════════════════════════════════════════
# 18. Pre-built symmetry groups
# ═══════════════════════════════════════════════════════════════════════════

function build_identity_group(n_site::Int)::Finite_Symmetry_Group
    id_op = Symmetry_Operation(:identity, collect(1:n_site))
    return Finite_Symmetry_Group("identity", [id_op]; identity_idx=1)
end

"""
Build the (Flux-aware) Lattice Translation Group ℤ^{L₁}×ℤ^{L₂}.
---
When `θ` is omitted or all zeros, the operations are bare permutations (no per-site phases).  This is the fast path used for ordinary symmetry-resolved ED.

When `θ` is supplied (or read from `lattice.twisted_phases_over_2π`), each site that crosses a periodic boundary under translation acquires the gauge-covariant phase `g(x)/g(Tx)` with `g(x)=exp(i 2π θ⋅x/L)`.  The resulting group commutes with the flux-inserted Hamiltonian `H(θ)` while keeping the irrep labels (i.e. the ordinary momentum labels) unchanged — the entire flux physics is captured in the per-site phases of the group operations themselves.
"""
function build_translation_group(
    lattice::TightBinding.Real_Space_Lattice,
    θ::AbstractVector{<:Real}=Float64[],
)::Finite_Symmetry_Group
    n_site = lattice.n_site
    L = lattice.sample_size

    # Resolve θ: explicit argument, then lattice built-in, then zeros
    if isempty(θ)
        θ_use = isempty(lattice.twisted_phases_over_2π) ?
                zeros(Float64, lattice.dim) : Float64.(lattice.twisted_phases_over_2π)
    else
        length(θ) == lattice.dim || error("θ must have length $(lattice.dim).")
        θ_use = Float64.(θ)
    end

    Lf = Float64.(L)

    ops = Symmetry_Operation{Tuple{Int,Int}}[]
    for dx in 0:(L[1]-1), dy in 0:(L[2]-1)
        shift = [dx, dy]
        perm = Vector{Int}(undef, n_site)
        phases = Vector{ComplexF64}(undef, n_site)
        @inbounds for (i, (cell_int, isub)) in enumerate(lattice.site_list)
            shifted_cell = mod.(cell_int .+ shift, L)
            perm[i] = lattice.site_to_index_map[(shifted_cell, isub)]
            phases[i] = cis(2π * dot(θ_use, (cell_int .- shifted_cell) ./ Lf))
        end
        push!(ops, Symmetry_Operation((dx, dy), perm; perm_phases=phases))
    end
    return Finite_Symmetry_Group("translations", ops; identity_idx=1)
end

# ═══════════════════════════════════════════════════════════════════════════
# 19. Pre-built irrep lists
# ═══════════════════════════════════════════════════════════════════════════

function build_identity_irrep_list()::Vector{OneDim_Irrep}
    return [OneDim_Irrep(:identity, [1.0 + 0.0im])]
end

function build_translation_irrep_list(G::Finite_Symmetry_Group,
    lattice::TightBinding.Real_Space_Lattice)::Vector{OneDim_Irrep{Tuple{Int,Int}}}
    @assert G.name == "translations" && group_order(G) == prod(lattice.sample_size)
    L1, L2 = lattice.sample_size
    irrep_list = OneDim_Irrep{Tuple{Int,Int}}[]
    for k1 in 0:(L1-1), k2 in 0:(L2-1)
        χ_list = ComplexF64[]
        sizehint!(χ_list, group_order(G))
        for op in G.operations
            dx, dy = op.label
            push!(χ_list, cis(2π * (k1 * dx / L1 + k2 * dy / L2)))
        end
        push!(irrep_list, OneDim_Irrep((k1, k2), χ_list))
    end
    return irrep_list
end

function build_irrep_list(G::Finite_Symmetry_Group, lattice::TightBinding.Real_Space_Lattice)::Vector{<:OneDim_Irrep}
    if G.name == "identity"
        return build_identity_irrep_list()
    elseif G.name == "translations"
        return build_translation_irrep_list(G, lattice)
    else
        error("Irrep construction for group '$(G.name)' is not yet implemented.")
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# 20. Convenience constructors
# ═══════════════════════════════════════════════════════════════════════════
"""
Core Constructor for `Symmetry_Resolved_ED_Data <: ED_Data`
---
by internally call `build_symmetry_orbit_catalog`, `build_irrep_list`, etc.

- Args:
    - `second_quantized_model::Real_Space_Second_Quantized_Model`
- Named Args:
    - `filling_fraction::Rational{Int}=1//2`: the filling per FLATTENED VERTICES!
    - `symmetry_group::Finite_Symmetry_Group=build_identity_group(second_quantized_model

**Note:** Here `filling_fraction` is the *particle number per _flatten_ed vertex*. So we ALWAYA have `n_filled = Int(filling_fraction * n_site)`. For example:
- For a spinful Hubbard model with 2×(Lx·Ly) vertices, the so-called "half-filling" means site-filling, i.e., one particle per vertex → `filling_fraction=1//2`.
- For most fractional Chern insulator communities, such as the bose Hubbard model over Haldane honeycomb lattice, or the spinless fermoinic Hubbard model over checkerboard lattice, the filling is per band!!!
"""
function build_ed_data(second_quantized_model::Real_Space_Second_Quantized_Model;
    filling_fraction::Rational{Int}=1 // 2,
    symmetry_group::Finite_Symmetry_Group,
)::Symmetry_Resolved_ED_Data
    n_site = second_quantized_model.lattice.n_site
    n_filled = Int(filling_fraction * n_site)
    @assert denominator(filling_fraction) * n_filled == numerator(filling_fraction) * n_site

    irrep_list = build_irrep_list(symmetry_group, second_quantized_model.lattice)
    catalog = build_symmetry_orbit_catalog(; second_quantized_model=second_quantized_model, n_filled=n_filled, symmetry_group=symmetry_group, statistics=second_quantized_model.statistics)
    sector_dims = zeros(Int, length(irrep_list))
    ed_scan_res = Dict{Int,Tuple{Vector{Float64},Matrix{ComplexF64}}}()

    return Symmetry_Resolved_ED_Data(second_quantized_model, n_filled, filling_fraction, symmetry_group,
        irrep_list, catalog, sector_dims, ed_scan_res)
end

function full_ed(second_quantized_model::Real_Space_Second_Quantized_Model, n_filled::Int; nev::Int=5)
    symmetry = build_identity_group(second_quantized_model.lattice.n_site)
    ed_data = build_ed_data(second_quantized_model; filling_fraction=n_filled // second_quantized_model.lattice.n_site, symmetry_group=symmetry)
    ed_scan!(ed_data; nev=nev)
    return ed_data.ed_scan_res[1]
end

# ═══════════════════════════════════════════════════════════════════════════
# 21. Utility: print & plot spectrum
# ═══════════════════════════════════════════════════════════════════════════

function print_spectrum(ed_data::Symmetry_Resolved_ED_Data; shift_to_zero::Bool=true)
    scanned = sort(collect(keys(ed_data.ed_scan_res)))
    isempty(scanned) && return
    n_irrep = length(ed_data.irrep_list)
    # Determine max nev across all scanned sectors (may differ per sector)
    nev = maximum(length(ed_data.ed_scan_res[i][1]) for i in scanned; init=0)
    spec = fill(NaN, n_irrep, nev)
    for irrep_idx in scanned
        vals = ed_data.ed_scan_res[irrep_idx][1]
        spec[irrep_idx, 1:length(vals)] = vals
    end
    if shift_to_zero
        finite = spec[.!isnan.(spec)]
        !isempty(finite) && (spec .-= minimum(finite))
    end
    println("ED Spectrum (rows = irrep, cols = eigenvalue #):")
    for irrep_idx in 1:n_irrep
        label = ed_data.irrep_list[irrep_idx].label
        vals_str = join([isnan(spec[irrep_idx, e]) ? "  NaN" : @sprintf("%.8f", spec[irrep_idx, e]) for e in 1:nev], "  ")
        println("  [$irrep_idx] $(repr(label)): $vals_str")
    end
    return spec
end

function plot_spectrum(ed_data::Symmetry_Resolved_ED_Data; shift_to_zero::Bool=true)
    scanned = sort(collect(keys(ed_data.ed_scan_res)))
    isempty(scanned) && return
    n_irrep = length(ed_data.irrep_list)
    # Determine max nev across all scanned sectors (may differ per sector)
    nev = maximum(length(ed_data.ed_scan_res[i][1]) for i in scanned; init=0)
    spec = fill(NaN, n_irrep, nev)
    for irrep_idx in scanned
        vals = ed_data.ed_scan_res[irrep_idx][1]
        spec[irrep_idx, 1:length(vals)] = vals
    end
    if shift_to_zero
        finite = spec[.!isnan.(spec)]
        !isempty(finite) && (spec .-= minimum(finite))
    end
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1]; xlabel="Irrep index", ylabel="E",
        title="ED Spectrum — $(ed_data.second_quantized_model.lattice.sample_size), ν=$(ed_data.filling_fraction)",
        xticks=0:2:(n_irrep-1), xminorticksvisible=true, yminorticksvisible=true)
    for k in 1:n_irrep, e in 1:nev
        val = spec[k, e]
        !isnan(val) && scatter!(ax, k - 1, val; color=:royalblue1, markersize=14, alpha=0.75,
            strokecolor=:blue, strokewidth=0.5)
    end
    return fig, ax
end



"""
Plot the ED Spectrum Resolved by Symmetry Sectors (1D Irreps)
---
- Args:
    - `ed_data::Symmetry_Resolved_Hard_Core_ED_Data`
- Named Args:
    - `shift_to_zero::Bool=true`: whether to shift the spectrum such that the minimum eigenvalue is zero
    - `show_spec::Bool=true`: whether to print the numerical values of the ED spectrum
    - `save_plot_path::Union{String,Nothing}`: absolute path to save the generated plot, default to `nothing`
"""
function plot_ed_scan_res(
    ed_data::Symmetry_Resolved_ED_Data;
    shift_to_zero::Bool=true,
    save_plot_path::Union{Nothing,String}=nothing,
    show_spec::Bool=true,
)::Tuple{CairoMakie.Figure,CairoMakie.Axis}
    scanned_indices = sort(collect(keys(ed_data.ed_scan_res)))
    nev = isempty(scanned_indices) ? 1 : length(ed_data.ed_scan_res[scanned_indices[1]][1])
    spec = reduce(hcat, [
        haskey(ed_data.ed_scan_res, irrep_idx) ?
        ed_data.ed_scan_res[irrep_idx][1] : fill(NaN, nev)
        for irrep_idx in eachindex(ed_data.irrep_list)
    ])' # the eigvals of one symmetry sector/irrep is `spec[<irrep_idx>,:]`; unscanned sectors filled with NaN

    show_spec && begin
        @info "Raw ED spectrum (symmetry-resolved):"
        show(stdout, "text/plain", spec)
        println()
    end

    finite_vals = spec[.!isnan.(spec)]
    if !isempty(finite_vals)
        shift_to_zero && (spec .-= minimum(finite_vals))
    end
    show_spec && begin
        @info "Shifted ED spectrum (symmetry-resolved):"
        show(stdout, "text/plain", spec)
        println()
    end

    sample_size = ed_data.second_quantized_model.lattice.sample_size
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1];
        xlabel="symmetry sector index",
        ylabel="E [unit]",
        title="ED Spectrum of $(sample_size) Sample",
        xticks=0:2:size(spec, 1),
        xminorticks=IntervalsBetween(2),
        xminorticksvisible=true,
        yminorticksvisible=true,
    )

    for k in axes(spec, 1), e in axes(spec, 2)
        val = spec[k, e]
        if !isnan(val)
            scatter!(ax, k - 1, val; color=:royalblue1, markersize=14, alpha=0.75,
                strokecolor=:blue, strokewidth=0.5)
        end
    end

    if !isnothing(save_plot_path)
        save(save_plot_path, fig)
        @info "Saved plot to $save_plot_path"
    end

    return fig, ax
end
