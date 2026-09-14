#!/usr/bin/env julia
# Render stored Julia timings and the imported Python comparison snapshot.
# Usage: julia --project=. benchmark/plot_benchmark.jl [julia_csv] [python_csv]
using CairoMakie
using Printf

const MODEL_TITLES = Dict(
    "Heisenberg" => "Spin-½ Heisenberg chain",
    "Haldane_Boson" => "Bosonic Haldane FCI",
    "Hubbard_Fermion" => "Spinful Fermi-Hubbard",
)
const FIGURE_STEMS = Dict(
    "Heisenberg" => "heisenberg_1D", "Haldane_Boson" => "bose_hubbard_2D",
    "Hubbard_Fermion" => "spinful_fermi_hubbard_2D",
)
const MODE_COLORS = Dict("matrix" => :royalblue3, "matrixfree" => :darkorange2)

function read_benchmarks(path)
    rows = NamedTuple[]
    lines = readlines(path)
    headers = Symbol.(split(strip(first(lines)), ','))
    for line in lines[2:end]
        isempty(strip(line)) && continue
        data = Dict(zip(headers, split(strip(line), ',')))
        push!(rows, (model=data[:model], label=data[:label], mode=data[:mode],
            n_site=parse(Int, data[:n_site]), sector_dim=parse(Int, data[:sector_dim]),
            elapsed_s=parse(Float64, data[:elapsed_s]), energy=parse(Float64, data[:energy])))
    end
    all(r.elapsed_s > 0 && isfinite(r.elapsed_s) for r in rows) || error("Invalid timing in $path")
    return rows
end

function plot_model(rows, model, outdir)
    selected = sort(filter(r -> r.model == model, rows); by=r -> (r.n_site, r.label, r.mode))
    isempty(selected) && return
    labels = unique(r.label for r in selected)
    indices = Dict(label => i for (i,label) in enumerate(labels))
    dimensions = Dict(r.label => r.sector_dim for r in selected)
    for scaling in (false, true)
        fig = Figure(size=(900, 520))
        ax = if scaling
            Axis(fig[1,1]; title=MODEL_TITLES[model] * " — Julia",
                xlabel="Sector Hilbert-space dimension", ylabel="Single-sector ED time (s)",
                xscale=log10, yscale=log10)
        else
            Axis(fig[1,1]; title=MODEL_TITLES[model] * " — Julia",
                xlabel="Sample size and sector dimension", ylabel="Single-sector ED time (s)",
                yscale=log10, xticklabelsize=11,
                xticks=(collect(1:length(labels)), ["$label\nD=$(dimensions[label])" for label in labels]))
        end
        for mode in ("matrix", "matrixfree")
            data = filter(r -> r.mode == mode, selected)
            isempty(data) && continue
            scaling && sort!(data; by=r -> r.sector_dim)
            x = scaling ? [r.sector_dim for r in data] : [indices[r.label] for r in data]
            scatterlines!(ax, x, [r.elapsed_s for r in data]; color=MODE_COLORS[mode],
                linestyle=mode == "matrix" ? :solid : :dash, marker=mode == "matrix" ? :circle : :rect,
                markersize=9, linewidth=2, label=mode)
        end
        axislegend(ax; position=:lt)
        path = joinpath(outdir, FIGURE_STEMS[model] * (scaling ? "_scaling" : "") * ".svg")
        save(path, fig)
        println(path)
    end
end

function plot_comparison(julia_rows, python_rows, outdir)
    fig = Figure(size=(1500, 510))
    for (index, model) in enumerate(("Heisenberg", "Haldane_Boson", "Hubbard_Fermion"))
        ax = Axis(fig[1,index]; title=MODEL_TITLES[model], xlabel="Sector dimension",
            ylabel=index == 1 ? "Single-sector ED time (s)" : "", xscale=log10, yscale=log10)
        for (language, rows, color) in (("Julia", julia_rows, :royalblue3), ("Python", python_rows, :darkorange2))
            for mode in ("matrix", "matrixfree")
                data = sort(filter(r -> r.model == model && r.mode == mode, rows); by=r -> r.sector_dim)
                isempty(data) && continue
                scatterlines!(ax, [r.sector_dim for r in data], [r.elapsed_s for r in data];
                    color, linestyle=mode == "matrix" ? :solid : :dash,
                    marker=mode == "matrix" ? :circle : :rect, markersize=7,
                    linewidth=2, label="$language $mode")
            end
        end
        index == 1 && axislegend(ax; position=:lt, labelsize=11)
    end
    Label(fig[0,1:3], "Stored benchmark scans — Julia and Python"; fontsize=22)
    path = joinpath(outdir, "benchmark_comparison.svg")
    save(path, fig)
    println(path)
end

function main(args=ARGS)
    julia_path = length(args) >= 1 ? args[1] : joinpath(@__DIR__, "benchmark_data", "benchmark_raw_latest.csv")
    python_path = length(args) >= 2 ? args[2] : joinpath(@__DIR__, "benchmark_data", "benchmark_raw_py_scan.csv")
    julia_rows = read_benchmarks(julia_path)
    outdir = joinpath(@__DIR__, "figures")
    mkpath(outdir)
    for model in ("Heisenberg", "Haldane_Boson", "Hubbard_Fermion")
        plot_model(julia_rows, model, outdir)
    end
    if isfile(python_path)
        plot_comparison(julia_rows, read_benchmarks(python_path), outdir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
