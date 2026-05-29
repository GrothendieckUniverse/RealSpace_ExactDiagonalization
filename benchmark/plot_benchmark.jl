#!/usr/bin/env julia
# ============================================================================
# Plot benchmark results from CSV files.
#
# Usage:
#   julia --project=. benchmark/plot_benchmark.jl [csv_file]
#
# If no csv_file is provided, reads benchmark/benchmark_data/benchmark_raw_latest.csv
# ============================================================================

using CairoMakie
using Printf

# ═══════════════════════════════════════════════════════════════════════════
# Helpers
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

function read_csv(path::String)::Vector{BenchRow}
    rows = BenchRow[]
    open(path, "r") do io
        readline(io)  # skip header
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ",")
            push!(rows, BenchRow(
                parts[1], parts[2],
                parse(Int, parts[3]), parse(Int, parts[4]),
                parse(Int, parts[5]), parse(Int, parts[6]),
                parse(Int, parts[7]), parse(Int, parts[8]),
                parts[9],
                parse(Float64, parts[10]), parse(Float64, parts[11]),
            ))
        end
    end
    return rows
end

# ═══════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════

const COLORS = Dict("matrix" => :royalblue3, "matrixfree" => :darkorange2)
const MODE_LABELS = Dict("matrix" => "matrix construction", "matrixfree" => "matrix-free")
const MODEL_NAMES = Dict(
    "Heisenberg" => "Spin-½ Heisenberg Chain",
    "Haldane_Boson" => "Bosonic Haldane FCI",
    "Hubbard_Fermion" => "Spinful Fermi-Hubbard",
)

function plot_model_times(rows::Vector{BenchRow}, model_key::String, outdir::String)
    rs = [r for r in rows if r.model == model_key]
    isempty(rs) && return

    labels_mat = unique(r.label for r in rs if r.mode == "matrix")
    n = length(labels_mat)
    xs = 1:n
    dim_by_label = Dict(r.label => r.sector_dim for r in rs if r.mode == "matrix")
    tick_lbls = ["$lbl\nD=$(dim_by_label[lbl])" for lbl in labels_mat]
    title = get(MODEL_NAMES, model_key, model_key)

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1];
        xlabel="System size (sector dimension below)",
        ylabel="Single-sector ED time (s)",
        title="$title — One-Sector Benchmark",
        yscale=log10,
        xticks=(xs, tick_lbls),
        xticklabelrotation=0.3,
    )
    for mode in ["matrix", "matrixfree"]
        mrs = sort([r for r in rs if r.mode == mode]; by=r -> findfirst(==(r.label), labels_mat))
        ys = [r.elapsed_s for r in mrs]
        scatter!(ax, xs, ys; markersize=14, color=COLORS[mode], label=MODE_LABELS[mode])
        lines!(ax, xs, ys; linewidth=2, color=COLORS[mode])
    end
    axislegend(ax; position=:lt)
    outpath = joinpath(outdir, "benchmark_$(model_key).svg")
    save(outpath, fig)
    println("  → $outpath")
end

function plot_model_scaling(rows::Vector{BenchRow}, model_key::String, outdir::String)
    rs = [r for r in rows if r.model == model_key]
    isempty(rs) && return
    title = get(MODEL_NAMES, model_key, model_key)

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1];
        xlabel="Sector Hilbert-space dimension",
        ylabel="Single-sector ED time (s)",
        title="$title — Scaling with Sector Dimension",
        xscale=log10, yscale=log10,
    )
    for mode in ["matrix", "matrixfree"]
        mrs = sort([r for r in rs if r.mode == mode]; by=r -> r.sector_dim)
        xs = [r.sector_dim for r in mrs]
        ys = [r.elapsed_s for r in mrs]
        scatter!(ax, xs, ys; markersize=12, color=COLORS[mode], label=MODE_LABELS[mode])
        lines!(ax, xs, ys; linewidth=2, color=COLORS[mode])
    end
    axislegend(ax; position=:lt)
    outpath = joinpath(outdir, "scaling_$(model_key).svg")
    save(outpath, fig)
    println("  → $outpath")
end

# ═══════════════════════════════════════════════════════════════════════════
# Combined summary plot (all three models in one figure)
# ═══════════════════════════════════════════════════════════════════════════

function plot_combined_summary(rows::Vector{BenchRow}, outdir::String)
    fig = Figure(size=(1200, 400))
    for (i, model_key) in enumerate(["Heisenberg", "Haldane_Boson", "Hubbard_Fermion"])
        rs = [r for r in rows if r.model == model_key]
        isempty(rs) && continue
        labels_mat = unique(r.label for r in rs if r.mode == "matrix")
        n = length(labels_mat)
        xs = 1:n
        dim_by_label = Dict(r.label => r.sector_dim for r in rs if r.mode == "matrix")
        tick_lbls = ["$lbl\nD=$(dim_by_label[lbl])" for lbl in labels_mat]
        title = get(MODEL_NAMES, model_key, model_key)

        ax = Axis(fig[1, i];
            xlabel="System",
            ylabel=i == 1 ? "Time (s)" : "",
            title=title,
            yscale=log10,
            xticks=(xs, tick_lbls),
            xticklabelrotation=0.3,
        )
        for mode in ["matrix", "matrixfree"]
            mrs = sort([r for r in rs if r.mode == mode]; by=r -> findfirst(==(r.label), labels_mat))
            ys = [r.elapsed_s for r in mrs]
            scatter!(ax, xs, ys; markersize=10, color=COLORS[mode], label=MODE_LABELS[mode])
            lines!(ax, xs, ys; linewidth=2, color=COLORS[mode])
        end
        i == 1 && axislegend(ax; position=:lt)
    end
    outpath = joinpath(outdir, "benchmark_combined.svg")
    save(outpath, fig)
    println("  → $outpath")
end

# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

function main()
    outdir = joinpath(@__DIR__, "figures")
    mkpath(outdir)

    csv_path = if length(ARGS) >= 1
        ARGS[1]
    else
        joinpath(@__DIR__, "benchmark_data/", "benchmark_raw_latest.csv")
    end

    if !isfile(csv_path)
        error("CSV file not found: $csv_path\n  Run benchmark/benchmark.jl first.")
    end

    println("Reading: $csv_path")
    rows = read_csv(csv_path)
    println("  $(length(rows)) rows, $(length(unique(r.model for r in rows))) models\n")

    # for model_key in ["Heisenberg", "Haldane_Boson", "Hubbard_Fermion"]
    #     plot_model_times(rows, model_key, outdir)
    #     plot_model_scaling(rows, model_key, outdir)
    # end
    plot_combined_summary(rows, outdir)
    println("\nAll plots saved to: $outdir")
end

main()
