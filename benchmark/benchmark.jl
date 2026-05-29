#!/usr/bin/env julia
# ============================================================================
# Benchmark: Single-Sector Symmetry-Resolved ED Timings
#
# Three models benchmarked in both matrix-construct and matrix-free modes:
#   1. Spin:    S=½ Heisenberg chain  (N_site = [20,22,24,26,28,30])
#   2. Bosonic: Haldane honeycomb FCI  (sample_size = [[2,3],[2,4],[2,5],[3,4],[2,7]])
#   3. Fermionic: Hubbard on square   (sample_size = [[2,3],[2,4],[2,5],[3,4],[2,7]])
#
# Each: JIT-warmup → time one sector (sector index 1) → CSV output
#
# Usage:
#   julia --project=. -p 8 benchmark/benchmark.jl
# ============================================================================

using Distributed

# ---- ensure enough workers (fallback to single-process if none) ----
const REQUIRED_PROCS = 8
if nprocs() < REQUIRED_PROCS + 1
    try
        addprocs(REQUIRED_PROCS + 1 - nprocs(); exeflags="--project=$(Base.active_project())")
    catch
        @warn "Could not add workers; running single-process."
    end
end

using RealSpace_ExactDiagonalization
using TightBinding
using LinearAlgebra, Printf, Dates
using CairoMakie

# ═══════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════

const OUTDIR = @__DIR__
const RESULT_DIR = joinpath(OUTDIR, "benchmark_data")
const FIG_DIR = joinpath(OUTDIR, "figures")
const MODES = [:matrix, :matrixfree]
const SECTOR_IDX = 1             # only benchmark the first sector
const NEV = 1             # we only need the lowest eigenvalue

mkpath(RESULT_DIR)
mkpath(FIG_DIR)
BLAS.set_num_threads(1)

# ═══════════════════════════════════════════════════════════════════════════
# Model builders
# ═══════════════════════════════════════════════════════════════════════════

"Build Heisenberg chain ED data for a given N (half-filling, PBC)"
function build_heisenberg_ed(N::Int)
    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=[N, 1],
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0]],
        lattice_name="Heisenberg_Chain",
        pbc_indicator=[true, false],
    )
    lattice = r_data
    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="heisenberg")

    J = 1.0
    # NOTE: add_hopping_term! already expands via translation symmetry → only ONE representative bond needed.
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => J / 2; is_hermitian=true)

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((sf, st), t) in tb_model.full_hopping_map
        push!(bilinear_terms, (lattice.site_to_index_map[sf], lattice.site_to_index_map[st], ComplexF64(t)))
    end

    # Density terms: J·n_i·n_j on bonds  +  −J/2·n_i on each site
    #   S_i·S_j = n_i n_j + ½(b†b+h.c.) − ½ n_i − ½ n_j + ¼
    #   Σ_{bonds}(−½ n_i−½ n_j) = −Σ_i n_i = −N_e
    #   Σ_{bonds} ¼ = N/4
    # →  H_spin = J·[bond] + J/2·[hop] − J·N/4  (at half-filling)
    # The −J·N/4 is absorbed via on-site (i,i,−J/2):  Σ_i (−J/2)·n_i = −J·N/4  ✓
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for x in 0:(N-1)
        i = lattice.site_to_index_map[([x, 0], 1)]
        j = lattice.site_to_index_map[([mod(x + 1, N), 0], 1)]
        push!(density_terms, (i, j, ComplexF64(J)))
    end
    for x in 0:(N-1)
        i = lattice.site_to_index_map[([x, 0], 1)]
        push!(density_terms, (i, i, ComplexF64(-J / 2)))
    end

    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("J" => J), lattice, tb_model, Bosonic(), bilinear_terms, density_terms,
    )
    n_site = lattice.n_site
    n_filled = N ÷ 2
    symmetry = build_translation_group(lattice)
    return build_ed_data(model; filling_fraction=n_filled // n_site, symmetry=symmetry)
end

"Build Haldane honeycomb FCI ED data"
function build_haldane_ed(sample_size::Vector{Int})
    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [1 / 2, sqrt(3) / 2]],
        sub_crys_list=[[0.0, 0.0], [1 / 3, 1 / 3]],
        lattice_name="Haldane_Honeycomb",
        pbc_indicator=[true, true],
    )
    lattice = r_data
    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="haldane")

    t, t′, t′′, ϕ = 1.0, 0.60, -0.58, 0.2
    sϕ = 2π * ϕ
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, -1], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′ * exp(-im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 1], 1)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′ * exp(-im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′ * exp(im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([-1, 1], 2)) => -t′ * exp(-im * sϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 1), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t′′; is_hermitian=true)

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((sf, st), t) in tb_model.full_hopping_map
        push!(bilinear_terms, (lattice.site_to_index_map[sf], lattice.site_to_index_map[st], ComplexF64(t)))
    end
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()  # V1=V2=0 for benchmark

    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("t" => t, "t′" => t′, "t′′" => t′′, "ϕ" => ϕ), lattice, tb_model, Bosonic(), bilinear_terms, density_terms,
    )
    n_filled = prod(sample_size) ÷ 2
    symmetry = build_translation_group(lattice)
    return build_ed_data(model; filling_fraction=n_filled // lattice.n_site, symmetry=symmetry)
end

"Build spinful Fermi-Hubbard ED data on square lattice"
function build_hubbard_ed(sample_size::Vector{Int})
    Lx, Ly = sample_size

    # spatial lattice (one sublattice)
    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0]],
        lattice_name="Square",
        pbc_indicator=[true, true],
    )
    tb_spinless = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="hubbard_spinless")
    # NOTE: add_hopping_term! already expands via translation symmetry → only ONE representative bond per direction.
    # x-neighbour: cell (0,0) → cell (1,0)
    add_hopping_term!(tb_spinless, (([0, 0], 1), ([1, 0], 1)) => ComplexF64(-1.0); is_hermitian=true)
    # y-neighbour: cell (0,0) → cell (0,1)
    add_hopping_term!(tb_spinless, (([0, 0], 1), ([0, 1], 1)) => ComplexF64(-1.0); is_hermitian=true)

    # spinful lattice: 2 "sublattices" = ↑, ↓
    r_data_spinful = TightBinding.initialize_real_space_lattice(;
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0], [0.0, 0.0]],
        lattice_name="Square_Hubbard",
        pbc_indicator=[true, true],
    )
    lattice = r_data_spinful
    n_site = lattice.n_site   # 2*Lx*Ly

    _sv(x, y, s) = lattice.site_to_index_map[([x, y], s)]  # s=1↑, s=2↓

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((sf, st), tamp) in tb_spinless.full_hopping_map
        x1, y1 = sf[1]
        x2, y2 = st[1]
        push!(bilinear_terms, (_sv(x1, y1, 1), _sv(x2, y2, 1), ComplexF64(tamp)))
        push!(bilinear_terms, (_sv(x1, y1, 2), _sv(x2, y2, 2), ComplexF64(tamp)))
    end
    density_terms = Tuple{Int,Int,ComplexF64}[]
    for y in 0:(Ly-1), x in 0:(Lx-1)
        push!(density_terms, (_sv(x, y, 1), _sv(x, y, 2), ComplexF64(8.0)))
    end

    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="Hubbard")
    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("t" => 1.0, "U" => 8.0), lattice, tb_model, Fermionic(), bilinear_terms, density_terms,
    )
    n_spatial = prod(sample_size)
    n_filled = n_spatial   # half-filling: N_up = N_down = n_spatial/2

    # build translation symmetry
    ops = Vector{Symmetry_Operation{Tuple{Int,Int}}}()
    for (δx, δy) in Iterators.product(0:(Lx-1), 0:(Ly-1))
        perm = Vector{Int}(undef, n_site)
        for (cell_int, isub) in lattice.site_list
            shifted = (mod.(cell_int .+ [δx, δy], sample_size), isub)
            i_old = lattice.site_to_index_map[(cell_int, isub)]
            i_new = lattice.site_to_index_map[shifted]
            perm[i_old] = i_new
        end
        push!(ops, Symmetry_Operation((δx, δy), perm))
    end
    symmetry = Finite_Symmetry_Group("translations", ops; identity_idx=1)
    return build_ed_data(model; filling_fraction=n_filled // n_site, symmetry=symmetry)
end

# ═══════════════════════════════════════════════════════════════════════════
# Benchmark helpers
# ═══════════════════════════════════════════════════════════════════════════

"Run a single-sector ED and return (elapsed_s, sector_dim, lowest_energy)"
function time_single_sector!(ed_data, mode::Symbol, sector_index::Int=1)
    irrep = ed_data.irrep_list[sector_index]
    GC.gc(true)
    t0 = time()
    if mode == :matrix
        vals, _ = ed_scan_at_irrep_matrix!(irrep.label, ed_data; nev=NEV)
    elseif mode == :matrixfree
        vals, _ = ed_scan_at_irrep_matrixfree!(irrep.label, ed_data; nev=NEV, use_distributed=false)
    else
        error("unknown mode $mode")
    end
    elapsed = time() - t0
    dim = ed_data.sector_dims[sector_index]
    GC.gc(true)
    return elapsed, dim, vals[1]
end

"JIT warmup — run a small system twice in each mode to trigger compilation"
function warmup!()
    println("=== JIT Warmup ===\n")
    for mode in MODES
        for rep in 1:2
            ed = build_haldane_ed([2, 3])
            time_single_sector!(ed, mode, 1)
            println("  warmup haldane [2,3] mode=$mode rep=$rep/2 done")
        end
    end
    # also warmup the distributed branch
    ed = build_haldane_ed([2, 4])
    time_single_sector!(ed, :matrix, 1)
    GC.gc(true)
    println("  warmup distributed matrix branch done")
    println("=== Warmup complete ===\n")
end

# ═══════════════════════════════════════════════════════════════════════════
# Run benchmarks
# ═══════════════════════════════════════════════════════════════════════════

struct BenchRow
    model::String
    label::String
    n_site::Int
    n_filled::Int
    full_dim::Int
    n_orbits::Int
    n_group::Int
    sector_dim::Int
    mode::String
    elapsed_s::Float64
    energy::Float64
end

function run_benchmarks()
    warmup!()

    all_rows = BenchRow[]

    # ---- 1. Heisenberg chain ----
    println("\n" * "="^70)
    println("  Model 1: Spin-½ Heisenberg Chain (translation symmetry)")
    println("="^70)
    for N in [20, 22, 24, 26, 28]
        for mode in MODES
            ed = build_heisenberg_ed(N)
            t_elapsed, dim, e0 = time_single_sector!(ed, mode)
            push!(all_rows, BenchRow(
                "Heisenberg", "N=$N", ed.second_quantized_model.lattice.n_site, ed.n_filled,
                binomial(ed.second_quantized_model.lattice.n_site, ed.n_filled),
                length(ed.orbit_catalog.representative_mask_list),
                length(ed.symmetry.operations), dim, string(mode), t_elapsed, e0,
            ))
            println("  N=$N  mode=$(rpad(mode,11))  dim=$dim  t=$(round(t_elapsed,digits=4))s")
        end
    end

    # ---- 2. Bosonic Haldane FCI ----
    println("\n" * "="^70)
    println("  Model 2: Bosonic Haldane FCI (translation symmetry)")
    println("="^70)
    for ss in [[2, 3], [2, 4], [2, 5], [3, 4], [2, 7], [4, 4]]
        for mode in MODES
            ed = build_haldane_ed(ss)
            t_elapsed, dim, e0 = time_single_sector!(ed, mode)
            push!(all_rows, BenchRow(
                "Haldane_Boson", "$(ss[1])×$(ss[2])", ed.second_quantized_model.lattice.n_site, ed.n_filled,
                binomial(ed.second_quantized_model.lattice.n_site, ed.n_filled),
                length(ed.orbit_catalog.representative_mask_list),
                length(ed.symmetry.operations), dim, string(mode), t_elapsed, e0,
            ))
            println("  $(ss[1])×$(ss[2])  mode=$(rpad(mode,11))  dim=$dim  t=$(round(t_elapsed,digits=4))s")
        end
    end

    # ---- 3. Fermionic Hubbard ----
    println("\n" * "="^70)
    println("  Model 3: Spinful Fermi-Hubbard (translation symmetry)")
    println("="^70)
    for ss in [[2, 2], [2, 3], [2, 4], [2, 5], [3, 4], [2, 7]]
        for mode in MODES
            ed = build_hubbard_ed(ss)
            t_elapsed, dim, e0 = time_single_sector!(ed, mode)
            push!(all_rows, BenchRow(
                "Hubbard_Fermion", "$(ss[1])×$(ss[2])", ed.second_quantized_model.lattice.n_site, ed.n_filled,
                binomial(ed.second_quantized_model.lattice.n_site, ed.n_filled),
                length(ed.orbit_catalog.representative_mask_list),
                length(ed.symmetry.operations), dim, string(mode), t_elapsed, e0,
            ))
            println("  $(ss[1])×$(ss[2])  mode=$(rpad(mode,11))  dim=$dim  t=$(round(t_elapsed,digits=4))s")
        end
    end

    return all_rows
end

# ═══════════════════════════════════════════════════════════════════════════
# Write CSV
# ═══════════════════════════════════════════════════════════════════════════

function write_csv(path::String, rows::Vector{BenchRow})
    open(path, "w") do io
        println(io, "model,label,n_site,n_filled,full_dim,n_orbits,n_group,sector_dim,mode,elapsed_s,energy")
        for r in rows
            println(io, join((
                    r.model, r.label, r.n_site, r.n_filled, r.full_dim, r.n_orbits, r.n_group,
                    r.sector_dim, r.mode,
                    @sprintf("%.9f", r.elapsed_s), @sprintf("%.15f", r.energy),
                ), ","))
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════

function plot_model(rows::Vector{BenchRow}, model_name::String, outfile::String;
    x_axis::Symbol=:n_site)
    rs = [r for r in rows if r.model == model_name]
    isempty(rs) && return

    labels_mat = unique(r.label for r in rs if r.mode == "matrix")
    n = length(labels_mat)
    xs = 1:n

    colors = Dict("matrix" => :royalblue3, "matrixfree" => :darkorange2)
    mode_labels = Dict("matrix" => "matrix construction", "matrixfree" => "matrix-free")

    # Build x-axis tick labels with sector dim info
    dim_by_label = Dict(r.label => r.sector_dim for r in rs if r.mode == "matrix")
    tick_lbls = ["$lbl\nD=$(dim_by_label[lbl])" for lbl in labels_mat]

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1];
        xlabel=(x_axis == :n_site ? "System size (sector dimension below)" : "Hilbert-space dimension"),
        ylabel="Single-sector ED time (s)",
        title="$model_name — One-Sector Benchmark",
        yscale=log10,
        xticks=(xs, tick_lbls),
        xticklabelrotation=0.3,
    )

    for mode in ["matrix", "matrixfree"]
        mrs = [r for r in rs if r.mode == mode]
        sort!(mrs; by=r -> findfirst(==(r.label), labels_mat))
        ys = [r.elapsed_s for r in mrs]
        scatter!(ax, xs, ys; markersize=14, color=colors[mode], label=mode_labels[mode])
        lines!(ax, xs, ys; linewidth=2, color=colors[mode])
    end
    axislegend(ax; position=:lt)
    save(outfile, fig)
    println("  → $outfile")
end

function plot_scaling(rows::Vector{BenchRow}, model_name::String, outfile::String)
    rs = [r for r in rows if r.model == model_name]
    isempty(rs) && return

    colors = Dict("matrix" => :royalblue3, "matrixfree" => :darkorange2)
    mode_labels = Dict("matrix" => "matrix construction", "matrixfree" => "matrix-free")

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1];
        xlabel="Sector Hilbert-space dimension",
        ylabel="Single-sector ED time (s)",
        title="$model_name — Scaling with Sector Dimension",
        xscale=log10, yscale=log10,
    )

    for mode in ["matrix", "matrixfree"]
        mrs = sort([r for r in rs if r.mode == mode]; by=r -> r.sector_dim)
        xs = [r.sector_dim for r in mrs]
        ys = [r.elapsed_s for r in mrs]
        scatter!(ax, xs, ys; markersize=12, color=colors[mode], label=mode_labels[mode])
        lines!(ax, xs, ys; linewidth=2, color=colors[mode])
    end
    axislegend(ax; position=:lt)
    save(outfile, fig)
    println("  → $outfile")
end

# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

println("="^70)
println("  Single-Sector ED Benchmark")
# println("  Timestamp: $TIMESTAMP")
println("  Workers: $(nworkers())  |  Threads: $(Threads.nthreads())")
println("="^70)

rows = run_benchmarks()

# Write CSV
raw_path = joinpath(RESULT_DIR, "benchmark_raw.csv")
write_csv(raw_path, rows)
println("\nResults → $raw_path")

# Also symlink/overwrite as "latest"
latest_path = joinpath(RESULT_DIR, "benchmark_raw_latest.csv")
cp(raw_path, latest_path; force=true)

# Plot
println("\n--- Generating plots ---")
plot_model(rows, "Heisenberg", joinpath(FIG_DIR, "heisenberg_1D.svg"))
plot_model(rows, "Haldane_Boson", joinpath(FIG_DIR, "bose_hubbard_2D.svg"))
plot_model(rows, "Hubbard_Fermion", joinpath(FIG_DIR, "spinful_fermi_hubbard_2D.svg"))

plot_scaling(rows, "Heisenberg", joinpath(FIG_DIR, "heisenberg_1D_scaling.svg"))
plot_scaling(rows, "Haldane_Boson", joinpath(FIG_DIR, "bose_hubbard_2D_scaling.svg"))
plot_scaling(rows, "Hubbard_Fermion", joinpath(FIG_DIR, "spinful_fermi_hubbard_2D_scaling.svg"))

println("\nDone. All plots saved to $FIG_DIR")
