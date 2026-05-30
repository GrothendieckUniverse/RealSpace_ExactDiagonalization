# XDiag Comparison: FCI validation + single-sector benchmark
# Usage: julia --project=. -p 8 benchmark/xdiag_compare.jl
#
# Compares XDiag.jl performance against RealSpace_ExactDiagonalization
# on the Haldane honeycomb bose-Hubbard model at ν=1/2 per band.

using XDiag
using RealSpace_ExactDiagonalization
using RealSpace_ExactDiagonalization.BitWise_Operations: Mask
using TightBinding, LinearAlgebra, Printf, CairoMakie
using Distributed


const PARAMS = Dict(
    "t" => 1.0,        # nearest-neighbour hopping
    "t′" => 0.60,       # next-nearest-neighbour hopping
    "t′′" => -0.58,       # next-next-nearest-neighbour hopping
    "ϕ_over_2π" => 0.2,  # flux per 2π (time-reversal breaking)
    "V1" => 0.0,         # NN density interaction
    "V2" => 0.0,         # NNN density interaction
)
const NEV = 3
const N_SECTORS = 3

# ═══════════════════════════════════════════════════════════════════════════
# Build Haldane honeycomb model (reusable, replaces removed test.jl helpers)
# ═══════════════════════════════════════════════════════════════════════════

"Build the Haldane honeycomb bose-Hubbard model for given sample_size."
function build_haldane_model(ss::Vector{Int})
    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=ss,
        brav_vec_list=[[1.0, 0.0], [1 / 2, sqrt(3) / 2]],
        sub_crys_list=[[0.0, 0.0], [1 / 3, 1 / 3]],
        lattice_name="Haldane_Honeycomb",
        pbc_indicator=[true, true],
    )
    lattice = r_data
    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="haldane")

    t, t′, t′′ = PARAMS["t"], PARAMS["t′"], PARAMS["t′′"]
    ϕ = PARAMS["ϕ_over_2π"]
    sϕ = 2π * ϕ

    # Nearest-neighbour (inter-sublattice, real)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, -1], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 0], 2)) => -t; is_hermitian=true)

    # Next-nearest-neighbour (intra-sublattice, complex)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′ * exp(-im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 1], 1)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′ * exp(-im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([-1, 1], 2)) => -t′ * exp(-im * sϕ); is_hermitian=true)

    # Next-next-nearest-neighbour (inter-sublattice, real)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 1), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t′′; is_hermitian=true)

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((sf, st), tamp) in tb_model.full_hopping_map
        push!(bilinear_terms, (lattice.site_to_index_map[sf], lattice.site_to_index_map[st], ComplexF64(tamp)))
    end
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()  # V1=V2=0 for benchmark

    model = RealSpace_ExactDiagonalization.Real_Space_Second_Quantized_Model(
        PARAMS, lattice, tb_model, Bosonic(), bilinear_terms, density_terms,
    )
    return model
end

# ═══════════════════════════════════════════════════════════════════════════
# Build Haldane honeycomb OpSum for XDiag's Spinhalf block
#
# Hard-core bosons → spin-1/2 via Matsubara-Matsuda:
#   b_i^† → S_i^+,  b_i → S_i^-,  n_i → S_i^z + 1/2
#
# H = Σ t_ij S_i^+ S_j^- + Σ V_ij (S_i^z+½)(S_j^z+½)
#
# We use XDiag's built-in operators: "SdotS", "SzSz", "Sz"
# Note: S_i^+ S_j^- + S_j^+ S_i^- = 2*(SdotS - SzSz)
#       → t_ij S_i^+ S_j^- + h.c. = 2*Re(t_ij)*(SdotS - SzSz)_ij
#                                    + i*Im(t_ij)*(S_i^+ S_j^- - S_j^+ S_i^-)
# For the imaginary part we use explicit Matrix operators.
# ═══════════════════════════════════════════════════════════════════════════

const Sp = [0.0 1.0; 0.0 0.0]   # S^+
const Sm = [0.0 0.0; 1.0 0.0]   # S^-

function build_xdiag_opsum(model::RealSpace_ExactDiagonalization.Real_Space_Second_Quantized_Model)
    ops = OpSum()

    # Hopping terms: directly use Matrix operator for each S_i^+ S_j^-
    # S_i^+ S_j^- = kron(Sp, Sm) as 4x4 matrix on sites [i,j]
    for (i, j, t) in model.bilinear_terms
        M = t .* kron(ComplexF64.(Sp), ComplexF64.(Sm))
        ops += Op("Matrix", [i, j], M)
    end

    # Density-density: Σ V_ij (S_i^z+½)(S_j^z+½) = V*SzSz + V/2*(Sz_i+Sz_j) + const
    for (i, j, V) in model.density_density_terms
        Vr = real(V)
        if abs(Vr) > 1e-14
            ops += Vr * Op("SzSz", [i, j])
            ops += Vr / 2 * Op("Sz", i)
            ops += Vr / 2 * Op("Sz", j)
        end
    end

    return ops
end

"Build translation permutation group for XDiag from TightBinding lattice"
function build_xdiag_translation_group(lattice::TightBinding.Real_Space_Lattice)
    n_site = lattice.n_site
    L1, L2 = lattice.sample_size
    perms = Permutation[]
    for dx in 0:(L1-1), dy in 0:(L2-1)
        perm_vec = Vector{Int}(undef, n_site)
        @inbounds for (i, (cell_int, isub)) in enumerate(lattice.site_list)
            shifted_cell = mod.(cell_int .+ [dx, dy], lattice.sample_size)
            perm_vec[i] = lattice.site_to_index_map[(shifted_cell, isub)]
        end
        push!(perms, Permutation(perm_vec))
    end
    return PermutationGroup(perms)
end

"Build XDiag irrep list for translation symmetry"
function build_xdiag_translation_irreps(group::PermutationGroup, lattice::TightBinding.Real_Space_Lattice)
    L1, L2 = lattice.sample_size
    nG = length(group)
    irreps = Representation[]
    for k1 in 0:(L1-1), k2 in 0:(L2-1)
        chi = [cis(2π * (k1 * dx / L1 + k2 * dy / L2)) for (dx, dy) in Iterators.product(0:L1-1, 0:L2-1)]
        # Build character table: each irrep is a row
        push!(irreps, Representation(group, [chi]))
    end
    return irreps
end

# ═══════════════════════════════════════════════════════════════════════════
# FCI Validation: [2,3] → E₀≈-7.1638, E₁≈-7.1634
# ═══════════════════════════════════════════════════════════════════════════

function validate_xdiag_fci()
    println("\n=== XDiag FCI Validation [2,3] ===")
    ss = [2, 3]
    model = build_haldane_model(ss)
    n_site = model.lattice.n_site
    n_filled = prod(ss) ÷ 2  # 3 particles

    ops = build_xdiag_opsum(model)
    group = build_xdiag_translation_group(model.lattice)
    L1, L2 = model.lattice.sample_size

    # Full Hilbert space (no symmetry)
    block_full = Spinhalf(n_site, n_filled)
    e0_full = eigval0(ops, block_full)
    println("  XDiag (full): E₀ = $(round(e0_full, digits=6))")

    # Symmetry-resolved
    e0_sym = Inf
    e1_sym = Inf
    for k1 in 0:(L1-1), k2 in 0:(L2-1)
        chi = ComplexF64[cis(2π * (k1 * dx / L1 + k2 * dy / L2))
                         for dx in 0:(L1-1) for dy in 0:(L2-1)]
        irrep = Representation(group, chi)
        block = Spinhalf(n_site, n_filled, irrep)
        dim(block) == 0 && continue
        res = eigs_lanczos(ops, block; neigvals=min(2, dim(block)))
        for ev in res.eigenvalues
            if ev < e0_sym
                e1_sym = e0_sym
                e0_sym = ev
            elseif ev < e1_sym
                e1_sym = ev
            end
        end
    end
    println("  XDiag (k-resolved): E₀ = $(round(e0_sym, digits=6)), E₁ = $(round(e1_sym, digits=6))")
    println("  Splitting: $(round(e1_sym - e0_sym, digits=6))")

    # Compare with our code
    symmetry_group = build_translation_group(model.lattice)
    ed_data = build_ed_data(model; filling_fraction=n_filled // n_site, symmetry_group=symmetry_group)
    ed_scan!(ed_data; nev=3, mode=:matrix)
    all_vals = sort!(reduce(vcat, [v for (_, (v, _)) in ed_data.ed_scan_res]))
    println("  Our code (k-resolved): E₀ = $(round(all_vals[1], digits=6)), E₁ = $(round(all_vals[2], digits=6))")
    println("  Splitting: $(round(all_vals[2] - all_vals[1], digits=6))")

    ok = abs(e0_sym - all_vals[1]) < 1e-6 && abs(e1_sym - all_vals[2]) < 1e-6
    println(ok ? "  ✓ XDiag FCI validation PASSED" : "  ✗ MISMATCH!")
    return ok
end

# ═══════════════════════════════════════════════════════════════════════════
# Single-sector benchmark for XDiag at given sample size
# ═══════════════════════════════════════════════════════════════════════════

function bench_xdiag_one(ss::Vector{Int})
    nu = prod(ss)
    ns = 2 * nu
    nf = nu ÷ 2
    sl = "$(ss[1])x$(ss[2])"
    println("\n", repeat("=", 60))
    println("XDiag: $sl -> $ns sites, $nf particles")
    println(repeat("=", 60))

    model = build_haldane_model(ss)

    ops = build_xdiag_opsum(model)
    group = build_xdiag_translation_group(model.lattice)
    L1, L2 = model.lattice.sample_size
    nG = L1 * L2  # Translation group order

    full_dim = binomial(ns, nf)
    println("  Full dim: $full_dim, |G|=$nG")

    # === Matrix mode: build CSR matrix + Lanczos ===
    irreps_data = [(k1, k2) for k1 in 0:(L1-1) for k2 in 0:(L2-1)]
    times_mat = Float64[]
    dims_mat = Int[]
    for idx in 1:min(N_SECTORS, length(irreps_data))
        k1, k2 = irreps_data[idx]
        chi = [cis(2π * (k1 * dx / L1 + k2 * dy / L2))
               for dx in 0:(L1-1) for dy in 0:(L2-1)]
        irrep = Representation(group, chi)
        block = Spinhalf(ns, nf, irrep)
        d = dim(block)
        d == 0 && continue

        GC.gc(true)

        t0 = time()

        # Build CSR matrix + diagonalize (XDiag's "matrix" mode)
        res = eigs_lanczos(ops, block; neigvals=NEV)

        t = time() - t0
        push!(times_mat, t)
        push!(dims_mat, d)
        println("    k=($k1,$k2) dim=$d  $(round(t,digits=3))s")
        flush(stdout)
    end
    avg_mat = isempty(times_mat) ? NaN : sum(times_mat) / length(times_mat)
    println("  Matrix avg: $(round(avg_mat,digits=3))s")

    return (sample_size=ss, n_sites=ns, n_filled=nf, full_dim=full_dim,
        nG=nG, max_dim=maximum(dims_mat; init=0),
        avg_time=avg_mat, mode=:xdiag)
end

# ═══════════════════════════════════════════════════════════════════════════
# Single-sector benchmark for our code (same conditions)
# ═══════════════════════════════════════════════════════════════════════════

function bench_ours_one(ss::Vector{Int}; mode::Symbol=:matrix)
    nu = prod(ss)
    ns = 2 * nu
    nf = nu ÷ 2
    sl = "$(ss[1])x$(ss[2])"
    println("\n", repeat("=", 60))
    println("Ours [$mode]: $sl -> $ns sites, $nf particles")
    println(repeat("=", 60))

    model = build_haldane_model(ss)
    symmetry_group = build_translation_group(model.lattice)
    nG = length(symmetry_group.operations)
    ed_data = build_ed_data(model; filling_fraction=nf // ns, symmetry_group=symmetry_group)

    full_dim = binomial(ns, nf)
    n_orbits = length(ed_data.orbit_catalog.representative_mask_list)
    println("  Full dim: $full_dim, orbits: $n_orbits ($(round(n_orbits/full_dim*100,digits=1))%), |G|=$nG")

    times = Float64[]
    dims = Int[]
    for i in 1:min(N_SECTORS, length(ed_data.irrep_list))
        haskey(ed_data.ed_scan_res, i) && continue
        irrep = ed_data.irrep_list[i]
        basis = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
        dim = length(basis.representative_mask_list)
        ed_data.sector_dims[i] = dim
        dim == 0 && continue

        GC.gc(true)

        t0 = time()

        if mode == :matrix
            cmap = CanonicalMap(symmetry_group, model.statistics, Dict{Mask,Tuple{Mask,Int,ComplexF64}}())
            H = build_ed_Hamiltonian_symmetry_block(basis, model.bilinear_terms,
                model.density_density_terms, cmap)
            vals, _ = diagonalize_block_arpack(H; nev=NEV)
            H = nothing
        else
            cmap = CanonicalMap(symmetry_group, model.statistics, Dict{Mask,Tuple{Mask,Int,ComplexF64}}())
            populate_canonical_map!(cmap, basis, model.bilinear_terms)
            H_op, n = hamiltonian_linear_operator(basis, model.bilinear_terms,
                model.density_density_terms, cmap)
            vals, _ = diagonalize_block_matrixfree(H_op, n; nev=NEV)
            cmap = nothing
        end

        t = time() - t0
        GC.gc(true)

        push!(times, t)
        push!(dims, dim)
        ed_data.ed_scan_res[i] = (vals, Matrix{ComplexF64}(undef, 0, 0))
        println("    k=$(irrep.label) dim=$dim  $(round(t,digits=3))s")
        flush(stdout)
    end
    avg = isempty(times) ? NaN : sum(times) / length(times)
    println("  Avg: $(round(avg,digits=3))s")
    return (sample_size=ss, n_sites=ns, n_filled=nf, full_dim=full_dim,
        n_orbits=n_orbits, nG=nG, max_dim=maximum(dims; init=0),
        avg_time=avg, mode=mode)
end

# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

println("XDiag vs RealSpace_ExactDiagonalization Benchmark")
println("Workers: $(nprocs()), BLAS: $(BLAS.get_num_threads())\n")

# 1. FCI validation
validate_xdiag_fci()

# 2. Benchmark suite
SS = [[2, 3], [2, 4], [2, 5], [3, 4]]

results = []
for ss in SS
    push!(results, bench_xdiag_one(ss))
    push!(results, bench_ours_one(ss; mode=:matrix))
    push!(results, bench_ours_one(ss; mode=:matrixfree))
    GC.gc(true)
end

# Summary
println("\n", repeat("=", 100))
println(rpad("Lattice", 8), rpad("Sites", 6), rpad("Np", 4), rpad("Full dim", 14),
    rpad("Mode", 14), rpad("Max dim", 14), rpad("Avg time", 12))
println(repeat("-", 80))
for r in results
    sl = "$(r.sample_size[1])x$(r.sample_size[2])"
    println(rpad(sl, 8), rpad(r.n_sites, 6), rpad(r.n_filled, 4),
        rpad(r.full_dim, 14), rpad(string(r.mode), 14),
        rpad(r.max_dim, 14), rpad("$(round(r.avg_time,digits=3))s", 12))
end

println("\nDone.")
