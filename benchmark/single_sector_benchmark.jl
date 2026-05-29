#!/usr/bin/env julia

using Distributed

const REQUIRED_PROCS = 8
if nprocs() < REQUIRED_PROCS + 1
    addprocs(REQUIRED_PROCS + 1 - nprocs(); exeflags="--project=$(Base.active_project())")
end

using LinearAlgebra
using Printf
using Dates
using CairoMakie
using RealSpace_ExactDiagonalization

include("../src/test.jl")
const ModelBench = Main.Test

const SAMPLE_SIZES = [[2, 3], [2, 4], [2, 5], [3, 4], [2, 7], [4, 4]]
const MODES = [:matrix, :matrixfree]
const MAX_SECTORS_PER_SAMPLE = 3
const NEV = 1
const OUTDIR = @__DIR__
const RESULT_DIR = joinpath(OUTDIR, "results")
const FIG_DIR = joinpath(OUTDIR, "figures")

mkpath(RESULT_DIR)
mkpath(FIG_DIR)
BLAS.set_num_threads(1)

function build_haldane_ed(sample_size::Vector{Int})
    n_filled = prod(sample_size) ÷ 2
    tb = ModelBench.build_bose_hubbard_real_space_tb_model(;
        sample_size,
        params=ModelBench.params_DNSheng,
    )
    model = ModelBench.initialize_second_quantized_model_for_Haldane_honeycomb_lattice(
        tb; params=ModelBench.params_DNSheng, statistics=Bosonic(),
    )
    symmetry = build_translation_group(model.lattice)
    return build_ed_data(model; filling_fraction=n_filled // model.lattice.n_site, symmetry)
end

function run_sector!(ed_data, mode::Symbol, sector_index::Int)
    irrep = ed_data.irrep_list[sector_index]
    if mode == :matrix
        return ed_scan_at_irrep_matrix!(irrep.label, ed_data; nev=NEV)
    elseif mode == :matrixfree
        return ed_scan_at_irrep_matrixfree!(irrep.label, ed_data; nev=NEV, use_distributed=false)
    else
        error("unknown mode $mode")
    end
end

function warmup!()
    println("JIT warmup")
    for mode in MODES, rep in 1:2
        ed_data = build_haldane_ed([2, 3])
        run_sector!(ed_data, mode, 1)
        GC.gc(true)
        println("  warmup mode=$mode rep=$rep/2 complete")
    end

    # Pre-specialize the distributed matrix construction branch so the [2,5]
    # timed sectors do not pay worker-side compilation.
    ed_data = build_haldane_ed([2, 5])
    run_sector!(ed_data, :matrix, 1)
    GC.gc(true)
    println("  warmup distributed matrix branch complete")
end

function benchmark_one_sector!(ed_data, sample_size::Vector{Int}, mode::Symbol, sector_index::Int)
    irrep = ed_data.irrep_list[sector_index]
    GC.gc(true)
    t0 = time()
    vals, _ = run_sector!(ed_data, mode, sector_index)
    elapsed = time() - t0
    dim = ed_data.sector_dims[sector_index]
    return (
        sample_size="$(sample_size[1])x$(sample_size[2])",
        n_unit_cell=prod(sample_size),
        n_site=ed_data.model.lattice.n_site,
        n_filled=ed_data.n_filled,
        full_dim=binomial(ed_data.model.lattice.n_site, ed_data.n_filled),
        n_orbits=length(ed_data.orbit_catalog.representative_mask_list),
        sector_index=sector_index,
        sector_label=repr(irrep.label),
        sector_dim=dim,
        mode=String(mode),
        elapsed_s=elapsed,
        energy=vals[1],
        nprocs=nprocs(),
        nworkers=nworkers(),
        nthreads=Threads.nthreads(),
    )
end

function write_raw_csv(path::String, rows)
    header = [
        "sample_size", "n_unit_cell", "n_site", "n_filled", "full_dim", "n_orbits",
        "sector_index", "sector_label", "sector_dim", "mode",
        "elapsed_s", "energy", "nprocs", "nworkers", "nthreads",
    ]
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join((
                r.sample_size, r.n_unit_cell, r.n_site, r.n_filled, r.full_dim, r.n_orbits,
                r.sector_index, r.sector_label, r.sector_dim, r.mode,
                @sprintf("%.9f", r.elapsed_s), @sprintf("%.15f", r.energy),
                r.nprocs, r.nworkers, r.nthreads,
            ), ","))
        end
    end
end

function summarize(rows)
    summaries = NamedTuple[]
    keys = unique((r.sample_size, r.mode) for r in rows)
    for (sample_size, mode) in keys
        rs = [r for r in rows if r.sample_size == sample_size && r.mode == mode]
        times = [r.elapsed_s for r in rs]
        push!(summaries, (
            sample_size=sample_size,
            n_unit_cell=rs[1].n_unit_cell,
            n_site=rs[1].n_site,
            n_filled=rs[1].n_filled,
            full_dim=rs[1].full_dim,
            n_orbits=rs[1].n_orbits,
            mean_sector_dim=sum(r.sector_dim for r in rs) / length(rs),
            max_sector_dim=maximum(r.sector_dim for r in rs),
            mode=mode,
            n_sectors=length(rs),
            mean_elapsed_s=sum(times) / length(times),
            min_elapsed_s=minimum(times),
            max_elapsed_s=maximum(times),
            energy=rs[end].energy,
        ))
    end
    sort!(summaries, by=r -> (r.n_site, r.mode))
    return summaries
end

function write_summary_csv(path::String, rows)
    open(path, "w") do io
        println(io, "sample_size,n_unit_cell,n_site,n_filled,full_dim,n_orbits,mean_sector_dim,max_sector_dim,mode,n_sectors,mean_elapsed_s,min_elapsed_s,max_elapsed_s,energy")
        for r in rows
            println(io, join((
                r.sample_size, r.n_unit_cell, r.n_site, r.n_filled, r.full_dim, r.n_orbits,
                @sprintf("%.3f", r.mean_sector_dim), r.max_sector_dim, r.mode, r.n_sectors,
                @sprintf("%.9f", r.mean_elapsed_s),
                @sprintf("%.9f", r.min_elapsed_s),
                @sprintf("%.9f", r.max_elapsed_s),
                @sprintf("%.15f", r.energy),
            ), ","))
        end
    end
end

function plot_results(summaries)
    site_order = [2 * prod(ss) for ss in SAMPLE_SIZES]
    dim_by_site = Dict(r.n_site => r.max_sector_dim for r in summaries)
    tick_labels = ["$(n)\nD=$(get(dim_by_site, n, 0))" for n in site_order]

    colors = Dict("matrix" => :royalblue3, "matrixfree" => :darkorange2)
    labels = Dict("matrix" => "matrix construction", "matrixfree" => "matrix-free")

    fig = Figure(size=(860, 500))
    ax = Axis(fig[1, 1];
        xlabel="number of sites\n(single-sector dimension shown below tick)",
        ylabel="single-sector ED time (s)",
        title="Haldane hard-core bosons, half-filling per band",
        yscale=log10,
        xticks=(site_order, tick_labels),
    )

    for mode in ["matrix", "matrixfree"]
        rs = [r for r in summaries if r.mode == mode]
        xs = [r.n_site for r in rs]
        ys = [r.mean_elapsed_s for r in rs]
        scatter!(ax, xs, ys; markersize=14, color=colors[mode], label=labels[mode])
        lines!(ax, xs, ys; linewidth=2, color=colors[mode])
    end
    axislegend(ax; position=:lt)
    save(joinpath(FIG_DIR, "single_sector_time_vs_sites.svg"), fig)

    fig2 = Figure(size=(860, 500))
    ax2 = Axis(fig2[1, 1];
        xlabel="single-sector Hilbert-space dimension",
        ylabel="single-sector ED time (s)",
        title="Haldane hard-core bosons, half-filling per band",
        xscale=log10,
        yscale=log10,
    )
    for mode in ["matrix", "matrixfree"]
        rs = [r for r in summaries if r.mode == mode]
        xs = [r.mean_sector_dim for r in rs]
        ys = [r.mean_elapsed_s for r in rs]
        scatter!(ax2, xs, ys; markersize=14, color=colors[mode], label=labels[mode])
        lines!(ax2, xs, ys; linewidth=2, color=colors[mode])
    end
    axislegend(ax2; position=:lt)
    save(joinpath(FIG_DIR, "single_sector_time_vs_sector_dim.svg"), fig2)
end

println("Single-sector ED benchmark")
println("processes=$(nprocs()), workers=$(nworkers()), threads=$(Threads.nthreads()), sectors/sample≤$MAX_SECTORS_PER_SAMPLE, nev=$NEV")
warmup!()

rows = NamedTuple[]
for sample_size in SAMPLE_SIZES
    for mode in MODES
        ed_data = build_haldane_ed(sample_size)
        n_sectors = min(MAX_SECTORS_PER_SAMPLE, length(ed_data.irrep_list))
        for sector_index in 1:n_sectors
            println("\nBenchmark sample_size=$sample_size mode=$mode sector=$sector_index/$n_sectors")
            row = benchmark_one_sector!(ed_data, sample_size, mode, sector_index)
            push!(rows, row)
            println(row)
        end
    end
end

timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
raw_path = joinpath(RESULT_DIR, "single_sector_raw_$timestamp.csv")
summary_path = joinpath(RESULT_DIR, "single_sector_summary_$timestamp.csv")
latest_raw = joinpath(RESULT_DIR, "single_sector_raw_latest.csv")
latest_summary = joinpath(RESULT_DIR, "single_sector_summary_latest.csv")
summaries = summarize(rows)
write_raw_csv(raw_path, rows)
write_raw_csv(latest_raw, rows)
write_summary_csv(summary_path, summaries)
write_summary_csv(latest_summary, summaries)
plot_results(summaries)

println("\nSummary")
for r in summaries
    println(r)
end
println("\nWrote:")
println("  $raw_path")
println("  $summary_path")
println("  $(joinpath(FIG_DIR, "single_sector_time_vs_sites.svg"))")
println("  $(joinpath(FIG_DIR, "single_sector_time_vs_sector_dim.svg"))")
