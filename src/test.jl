module Test

using LinearAlgebra
using SparseArrays, Arpack
using KrylovKit
using Printf
using Random
using MLStyle
using TightBinding
using RealSpace_ExactDiagonalization

# Import symmetry-resolved ED machinery
using RealSpace_ExactDiagonalization: Symmetry_Operation, Finite_Symmetry_Group,
    Symmetry_Orbit_Catalog, OneDim_Irrep, Symmetry_Sector_Basis,
    Symmetry_Resolved_ED_Data, CanonicalMap,
    build_identity_group, build_translation_group,
    build_identity_irrep_list, build_translation_irrep_list, build_irrep_list,
    build_symmetry_orbit_catalog, build_symmetry_sector_basis,
    build_ed_Hamiltonian_symmetry_block, build_ed_data,
    ed_scan!, ed_scan_at_irrep_matrix!, ed_scan_at_irrep_matrixfree!,
    diagonalize_block_dense, diagonalize_block_arpack, diagonalize_block_matrixfree,
    print_spectrum, plot_spectrum,
    apply_operation_to_mask, get_canonical_representative,
    project_to_sector, apply_hamiltonian!, populate_canonical_map!

# Import bitwise operations
using RealSpace_ExactDiagonalization.BitWise_Operations: Mask,
    is_site_occupied, is_site_empty,
    occupy_site_for_mask, empty_site_for_mask,
    n_occupied_for_mask

"parameters from D. N. Sheng's bosonic FCI [PhysRevLett.107.146803]"
const params_DNSheng = Dict(
    "t" => 1.0, # nearest-neighbor hopping
    "t′" => 0.60, # next-nearest-neighbor hopping
    "t′′" => -0.58, # next-next-nearest-neighbor hopping
    "ϕ_over_2π" => 0.2, # flux per 2π
    "V1" => 0.0, # nearest neighbor density-density interaction strength
    "V2" => 0.0, # next-nearest neighbor density-density interaction strength
)


"build `TightBinding.Real_Space_TightBinding_Model` over Haldane honeycomb lattice following D. N. Sheng's bosonic FCI [PhysRevLett.107.146803]"
function build_bose_hubbard_real_space_tb_model(; sample_size::Vector{Int}=[4, 3], params::Dict{String,<:Number}=params_DNSheng)::TightBinding.Real_Space_TightBinding_Model
    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [1 / 2, sqrt(3) / 2]],
        sub_crys_list=[[0.0, 0.0], [1 / 3, 1 / 3]],
        lattice_name="Haldane_Honeycomb",
        pbc_indicator=[true, true]
    )

    real_space_tb_model = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="test")

    t = params["t"]
    t′ = params["t′"]
    t′′ = params["t′′"]
    ϕ_over_2π = params["ϕ_over_2π"]

    # see [PhysRevLett.134.076601] for illustration of hoppings
    # add nearest neighbor hoppings
    # between sublattices
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([0, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([0, -1], 2)) => -t; is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([-1, 0], 2)) => -t; is_hermitian=true)

    # add next-nearest neighbor hoppings
    # for sublattice 1
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 1), ([-1, 1], 1)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    # for sublattice 2
    add_hopping_term!(real_space_tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([0, 0], 2), ([-1, 1], 2)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)


    # add next-next-nearest neighbor hoppings
    # between sublattices
    add_hopping_term!(real_space_tb_model, (([0, 0], 2), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([1, 0], 1), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(real_space_tb_model, (([1, 0], 2), ([0, 1], 1)) => -t′′; is_hermitian=true)

    n_sub = real_space_tb_model.lattice.n_sub
    n_site = real_space_tb_model.lattice.n_site
    @info "Initialize `TightBinding.Real_Space_TightBinding_Model` with\n\t\t`sample_size`: $(sample_size), `n_sub`: $n_sub --- in total `n_site`: $n_site."

    real_space_tb_model
end


"model-specific constructor for `ShortRange_Real_Space_Second_Quantized_Model`"
function initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
    tb_model::TightBinding.Real_Space_TightBinding_Model;
    params::Dict{String,<:Number},
    statistics::RealSpace_ExactDiagonalization.Statistics
)::RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model
    lattice = tb_model.lattice
    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((site_from, site_to), t) in tb_model.full_hopping_map
        i_from = lattice.site_to_index_map[site_from]
        i_to = lattice.site_to_index_map[site_to]
        push!(bilinear_terms, (i_from, i_to, ComplexF64(t)))
    end

    density_density_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    V1 = get(params, "V1", 0.0)
    V2 = get(params, "V2", 0.0)
    V3 = get(params, "V3", 0.0)
    if V1 != 0.0
        for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
            (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)

            cell_shift = cell_to .- cell_from
            if (mod.(cell_shift, lattice.sample_size) == [0, 0] && sub_from == 1 && sub_to == 2) ||
               (mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 2 && sub_to == 1) ||
               (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 2 && sub_to == 1)

                push!(density_density_terms, (i_from, i_to, ComplexF64(V1)))
            end
        end
    end
    if V2 != 0.0
        for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
            (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)

            cell_shift = cell_to .- cell_from
            if (mod.(cell_shift, lattice.sample_size) == [0, 0] && sub_from == 1 && sub_to == 2) ||
               (mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 2 && sub_to == 1) ||
               (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 2 && sub_to == 1)

                push!(density_density_terms, (i_from, i_to, ComplexF64(V2)))
            end
        end
    end
    if V3 != 0.0
        for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
            (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)

            cell_shift = cell_to .- cell_from
            if (mod.(cell_shift, lattice.sample_size) == [0, 0] && sub_from == 1 && sub_to == 2) ||
               (mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 2 && sub_to == 1) ||
               (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 2 && sub_to == 1)

                push!(density_density_terms, (i_from, i_to, ComplexF64(V3)))
            end
        end
    end

    second_quantized_model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        params,
        lattice,
        tb_model,
        statistics,
        bilinear_terms,
        density_density_terms,
    )

    @info "Initialize `ShortRange_Real_Space_Second_Quantized_Model <: Second_Quantized_Model` with\n\t\t`params`: $(params)"

    return second_quantized_model
end


# ============================================================================
# FCI Test: Haldane honeycomb lattice, 2×3 unit cells, 3 hard-core bosons
#
# Expected: two nearly degenerate ground states at energies ≈ -7.1633 and -7.1638
# (characteristic of ν=1/2 bosonic FCI on a torus)
# ============================================================================

"""
    test_FCI_ground_state(; sample_size=[2,3], n_filled=3, use_translation_symmetry=true)

Run ED on the Haldane honeycomb model at ν=1/2 filling of the lower band.

- Builds the tight-binding model and second-quantized model
- Optionally uses translation symmetry to block-diagonalize
- Returns the lowest eigenvalues across all sectors

Reference values (D.N. Sheng, Phys. Rev. Lett. 107, 146803):
    E₀ ≈ -7.1633, E₁ ≈ -7.1638  (two nearly degenerate FCI ground states)
"""
function test_FCI_ground_state(;
    sample_size::Vector{Int}=[2, 3],
    n_filled::Int=3,
    use_translation_symmetry::Bool=true,
    params::Dict{String,<:Number}=params_DNSheng,
)
    @info "=== FCI Ground State Test ==="
    @info "  lattice: $(sample_size) (×2 sublattices = $(2*prod(sample_size)) sites)"
    @info "  particles: $n_filled hard-core bosons"
    @info "  params: $params"

    # Build tight-binding model
    tb_model = build_bose_hubbard_real_space_tb_model(;
        sample_size=sample_size,
        params=params,
    )

    # Build second-quantized model
    model = initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
        tb_model;
        params,
        statistics=RealSpace_ExactDiagonalization.Bosonic(),
    )

    n_site = model.lattice.n_site
    filling_fraction = n_filled // n_site
    @info "  n_site=$n_site, filling_fraction=$filling_fraction"

    # Build symmetry group
    if use_translation_symmetry
        @assert all(model.lattice.pbc_indicator) "Translation symmetry requires PBC in all directions"
        symmetry = build_translation_group(model.lattice)
        @info "  Using translation symmetry: |G| = $(length(symmetry.operations))"
    else
        symmetry = build_identity_group(n_site)
        @info "  Using identity group (no symmetry resolution)"
    end

    # Build ED data
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry=symmetry)

    @info "  Orbit catalog: $(length(ed_data.orbit_catalog.representative_mask_list)) orbits"
    @info "  Irrep list: $(length(ed_data.irrep_list)) irreps"

    # Run ED scan
    ed_scan!(ed_data; nev=5)

    # Collect all eigenvalues across all sectors
    all_vals = Float64[]
    all_labels = Tuple{Int,Any,Int}[]  # (irrep_idx, irrep_label, eigenvalue_idx)
    for (irrep_idx, (vals, _)) in ed_data.ed_scan_res
        for (e_idx, v) in enumerate(vals)
            push!(all_vals, v)
            push!(all_labels, (irrep_idx, ed_data.irrep_list[irrep_idx].label, e_idx))
        end
    end

    # Sort by energy
    perm = sortperm(all_vals)
    all_vals = all_vals[perm]
    all_labels = all_labels[perm]

    # Print results
    println("\n" * "="^70)
    println("FCI Ground State Test Results")
    println("="^70)
    println("Lowest 10 eigenvalues across all symmetry sectors:")
    for i in 1:min(10, length(all_vals))
        irrep_idx, label, e_idx = all_labels[i]
        println("  E$i = $(@sprintf("%.8f", all_vals[i]))  (irrep $irrep_idx: $(repr(label)), eigenstate #$e_idx)")
    end

    # Check FCI signature: two nearly degenerate ground states
    if length(all_vals) >= 2
        gap = all_vals[2] - all_vals[1]
        println("\n  ΔE₁₂ = $(@sprintf("%.6f", gap))")
        if gap < 0.01
            println("  ✓ Two nearly degenerate ground states detected (FCI signature)")
        end
    end
    if length(all_vals) >= 3
        gap_to_excited = all_vals[3] - all_vals[2]
        println("  ΔE₂₃ = $(@sprintf("%.6f", gap_to_excited))")
    end

    println("="^70)

    return ed_data, all_vals, all_labels
end


"""
    validate_matrix_matrixfree_agreement(; sample_size=[2,3], n_filled=prod(sample_size) ÷ 2, nev=3, atol=1e-9)

Build the Haldane hard-core boson problem once and solve every translation
sector with both sparse-matrix and matrix-free modes. Throws if any sector
differs by more than `atol`.
"""
function validate_matrix_matrixfree_agreement(;
    sample_size::Vector{Int}=[2, 3],
    n_filled::Int=prod(sample_size) ÷ 2,
    params::Dict{String,<:Number}=params_DNSheng,
    nev::Int=3,
    atol::Float64=1e-9,
)
    tb_model = build_bose_hubbard_real_space_tb_model(; sample_size, params)
    model = initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
        tb_model; params, statistics=RealSpace_ExactDiagonalization.Bosonic(),
    )
    symmetry = build_translation_group(model.lattice)
    filling_fraction = n_filled // model.lattice.n_site

    ed_matrix = build_ed_data(model; filling_fraction, symmetry)
    ed_matrixfree = build_ed_data(model; filling_fraction, symmetry)
    ed_scan!(ed_matrix; nev, mode=:matrix)
    ed_scan!(ed_matrixfree; nev, mode=:matrixfree)

    maxdiff = 0.0
    for i in eachindex(ed_matrix.irrep_list)
        vals_matrix = ed_matrix.ed_scan_res[i][1]
        vals_matrixfree = ed_matrixfree.ed_scan_res[i][1]
        diff = maximum(abs.(vals_matrix .- vals_matrixfree); init=0.0)
        maxdiff = max(maxdiff, diff)
        @assert diff <= atol "matrix/matrix-free mismatch in sector $(ed_matrix.irrep_list[i].label): $diff"
    end
    println("matrix/matrix-free agreement: sample_size=$sample_size, maxdiff=$maxdiff")
    return maxdiff, ed_matrix, ed_matrixfree
end


"""
    sanity_check_bose_hubbard_2x3(; atol=1e-6)

Regression test for the hard-core boson Haldane-Hubbard model on a 2x3
periodic honeycomb sample at three particles. The two lowest states should be
the FCI pair near -7.1638 and -7.1634.
"""
function sanity_check_bose_hubbard_2x3(; atol::Float64=1e-6)
    _, vals, labels = test_FCI_ground_state(; sample_size=[2, 3], n_filled=3, use_translation_symmetry=true)
    @assert length(vals) >= 2
    @assert abs(vals[1] - (-7.16380536)) <= atol "unexpected boson FCI E0=$(vals[1])"
    @assert abs(vals[2] - (-7.16337536)) <= atol "unexpected boson FCI E1=$(vals[2])"
    return vals[1:2], labels[1:2]
end


@inline function _linear_square_site(x::Int, y::Int, Lx::Int)::Int
    return x + (y - 1) * Lx
end

function _square_lattice_undirected_bonds_obc(Lx::Int, Ly::Int)
    bonds = Tuple{Int,Int}[]
    for y in 1:Ly, x in 1:Lx
        i = _linear_square_site(x, y, Lx)
        x < Lx && push!(bonds, (i, _linear_square_site(x + 1, y, Lx)))
        y < Ly && push!(bonds, (i, _linear_square_site(x, y + 1, Lx)))
    end
    return bonds
end

function _fixed_weight_masks(n_site::Int, n_filled::Int)
    masks = Mask[]
    if n_filled == 0
        push!(masks, zero(Mask))
        return masks
    end
    x = (one(Mask) << n_filled) - one(Mask)
    upper = one(Mask) << n_site
    while x < upper
        push!(masks, x)
        c = x & -x
        r = x + c
        x = (((r ⊻ x) >> 2) ÷ c) | r
    end
    return masks
end

@inline function _fermion_hop_sign(m::Mask, from::Int, to::Int)::Float64
    from == to && return 1.0
    lo = min(from, to)
    hi = max(from, to)
    between = hi - lo <= 1 ? zero(Mask) : (((one(Mask) << (hi - lo - 1)) - one(Mask)) << lo)
    return isodd(count_ones(m & between)) ? -1.0 : 1.0
end

function _apply_spinful_hubbard_3x4!(
    y::Vector{Float64},
    x::Vector{Float64},
    masks::Vector{Mask},
    mask_to_idx::Dict{Mask,Int},
    directed_bonds::Vector{Tuple{Int,Int}},
    U::Float64,
)
    nm = length(masks)
    fill!(y, 0.0)
    @inbounds for idn in 1:nm
        mdn = masks[idn]
        col_offset = (idn - 1) * nm
        for iup in 1:nm
            col = iup + col_offset
            amp = x[col]
            amp == 0.0 && continue
            mup = masks[iup]

            y[col] += U * count_ones(mup & mdn) * amp

            for (from, to) in directed_bonds
                if is_site_occupied(mup, from) && is_site_empty(mup, to)
                    new_up = occupy_site_for_mask(empty_site_for_mask(mup, from), to)
                    row = mask_to_idx[new_up] + col_offset
                    y[row] += -_fermion_hop_sign(mup, from, to) * amp
                end
                if is_site_occupied(mdn, from) && is_site_empty(mdn, to)
                    new_dn = occupy_site_for_mask(empty_site_for_mask(mdn, from), to)
                    row = iup + (mask_to_idx[new_dn] - 1) * nm
                    y[row] += -_fermion_hop_sign(mdn, from, to) * amp
                end
            end
        end
    end
    return y
end

"""
    sanity_check_spinful_fermi_hubbard_3x4(; atol=1e-7)

Reference check against the ExactDiagonalization.jl documentation example:
3x4 open square lattice, spin-1/2 Fermi Hubbard model, t=1, U=8, half filling
with N_up=N_down=6. The documented ground-state energy is
-4.913259209075605.

Source: https://quantum-many-body.github.io/ExactDiagonalization.jl/dev/examples/HubbardModel/
"""
function sanity_check_spinful_fermi_hubbard_3x4(;
    Lx::Int=3,
    Ly::Int=4,
    U::Float64=8.0,
    atol::Float64=1e-7,
    tol::Float64=1e-10,
)
    n_spatial = Lx * Ly
    n_each_spin = n_spatial ÷ 2
    masks = _fixed_weight_masks(n_spatial, n_each_spin)
    mask_to_idx = Dict{Mask,Int}(m => i for (i, m) in enumerate(masks))
    undirected = _square_lattice_undirected_bonds_obc(Lx, Ly)
    directed = Tuple{Int,Int}[]
    sizehint!(directed, 2 * length(undirected))
    for (i, j) in undirected
        push!(directed, (i, j))
        push!(directed, (j, i))
    end

    dim = length(masks)^2
    Random.seed!(20240528)
    x0 = randn(dim)
    x0 ./= norm(x0)
    y = similar(x0)
    H = v -> begin
        _apply_spinful_hubbard_3x4!(y, v, masks, mask_to_idx, directed, U)
        copy(y)
    end

    vals, _, info = KrylovKit.eigsolve(H, x0, 1, :SR; tol=tol, maxiter=400, krylovdim=50)
    E0 = real(vals[1])
    reference = -4.913259209075605
    @assert abs(E0 - reference) <= atol "spinful Hubbard E0=$E0 differs from reference $reference; Krylov info=$info"
    println("spinful Hubbard 3x4 sanity check: E0=$(@sprintf("%.15f", E0))")
    return E0, info
end


function build_spinless_fermion_chain_model(; L::Int=4, t::Float64=1.0)
    lattice = TightBinding.initialize_real_space_lattice(;
        sample_size=[L, 1],
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0]],
        lattice_name="Open_Chain",
        pbc_indicator=[false, false],
    )
    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="spinless_chain")

    bilinear_terms = Tuple{Int,Int,ComplexF64}[]
    for x in 0:(L-2)
        i = lattice.site_to_index_map[([x, 0], 1)]
        j = lattice.site_to_index_map[([x + 1, 0], 1)]
        push!(bilinear_terms, (i, j, ComplexF64(-t)))
        push!(bilinear_terms, (j, i, ComplexF64(-t)))
    end

    return RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("t" => t),
        lattice,
        tb_model,
        RealSpace_ExactDiagonalization.Fermionic(),
        bilinear_terms,
        Tuple{Int,Int,ComplexF64}[],
    )
end


"""
    sanity_check_spinless_fermion_chain(; L=4, n_filled=2)

Exercise the generic `Fermionic()` ED path on a free open chain. The many-body
ground energy must equal the sum of the lowest `n_filled` single-particle
eigenvalues.
"""
function sanity_check_spinless_fermion_chain(; L::Int=4, n_filled::Int=2, atol::Float64=1e-10)
    model = build_spinless_fermion_chain_model(; L)
    vals, _ = full_ed(model, n_filled; nev=1)

    h1 = zeros(Float64, L, L)
    for i in 1:(L-1)
        h1[i, i + 1] = -1.0
        h1[i + 1, i] = -1.0
    end
    expected = sum(sort(eigvals(Symmetric(h1)))[1:n_filled])
    @assert abs(vals[1] - expected) <= atol "spinless fermion chain E0=$(vals[1]) differs from free-fermion value $expected"
    println("spinless fermion chain sanity check: E0=$(@sprintf("%.15f", vals[1]))")
    return vals[1], expected
end


"""
    run_haldane_size_sweep(; sample_sizes=[[2,3],[2,4],[2,5],[3,4],[4,4]], mode=:matrix, nev=3, max_sectors=nothing)

Run a translation-resolved size sweep for the hard-core-boson Haldane model.
Use `max_sectors` for quick smoke tests; leave it as `nothing` for the full
scan over every momentum sector.
"""
function run_haldane_size_sweep(;
    sample_sizes::Vector{Vector{Int}}=[[2, 3], [2, 4], [2, 5], [3, 4], [4, 4]],
    params::Dict{String,<:Number}=params_DNSheng,
    mode::Symbol=:matrix,
    nev::Int=3,
    max_sectors::Union{Nothing,Int}=nothing,
    use_distributed::Bool=false,
)
    rows = NamedTuple[]
    for sample_size in sample_sizes
        n_filled = prod(sample_size) ÷ 2
        tb_model = build_bose_hubbard_real_space_tb_model(; sample_size, params)
        model = initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
            tb_model; params, statistics=RealSpace_ExactDiagonalization.Bosonic(),
        )
        symmetry = build_translation_group(model.lattice)
        ed_data = build_ed_data(model; filling_fraction=n_filled // model.lattice.n_site, symmetry)

        n_run = max_sectors === nothing ? length(ed_data.irrep_list) : min(max_sectors, length(ed_data.irrep_list))
        t0 = time()
        for i in 1:n_run
            irrep = ed_data.irrep_list[i]
            if mode == :matrix
                ed_scan_at_irrep_matrix!(irrep.label, ed_data; nev)
            elseif mode == :matrixfree
                ed_scan_at_irrep_matrixfree!(irrep.label, ed_data; nev, use_distributed)
            else
                error("Unknown mode: $mode")
            end
        end
        elapsed = time() - t0
        dims = ed_data.sector_dims[1:n_run]
        row = (
            sample_size=sample_size,
            n_site=model.lattice.n_site,
            n_filled=n_filled,
            full_dim=binomial(model.lattice.n_site, n_filled),
            n_orbits=length(ed_data.orbit_catalog.representative_mask_list),
            n_sectors_run=n_run,
            max_sector_dim=maximum(dims; init=0),
            elapsed=elapsed,
        )
        push!(rows, row)
        println(row)
        GC.gc(true)
    end
    return rows
end


"precompilation/demo call"
function demo_run(;
    sample_size::Vector{Int}=[3, 2],
    params::Dict{String,<:Number}=params_DNSheng,
)
    tb_model = build_bose_hubbard_real_space_tb_model(;
        sample_size=sample_size,
        params=params,
    )

    second_quantized_model = initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
        tb_model;
        params,
        statistics=RealSpace_ExactDiagonalization.Bosonic(),
    )
end


# demo_run(; sample_size=[3, 2])

end
