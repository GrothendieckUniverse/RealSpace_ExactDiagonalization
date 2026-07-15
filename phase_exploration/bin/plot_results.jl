#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Render figures from generated CSV data (no diagonalization is performed).

Usage:
  julia --project=. phase_exploration/bin/plot_results.jl \\
    [--kind all|sweep|ed-spectra|structure|diagnostics|charge-gap] \\
    [--geometries '3x4;3x5;3x6'] [--max-ranks 12] [--max-flow-curves 10]
"""

opts, _ = parse_cli(ARGS)
if haskey(opts, "help") || haskey(opts, "h")
    println(HELP)
    exit()
end
kind = lowercase(get(opts, "kind", "all"))
samples = if haskey(opts, "geometries")
    [parse_geometry(text) for text in split(opts["geometries"], ';')]
else
    STUDY_GEOMETRIES
end
max_ranks = parse(Int, get(opts, "max-ranks", "12"))
max_flow_curves = parse(Int, get(opts, "max-flow-curves", "10"))

if kind == "all"
    plot_sweep_results(; samples=samples, max_ranks=max_ranks)
    plot_ed_spectrum_results(; samples=samples)
    plot_structure_factor_results(; samples=samples)
    plot_diagnostic_results(; samples=samples, max_flow_curves=max_flow_curves)
    plot_charge_gap_results()
elseif kind == "sweep"
    plot_sweep_results(; samples=samples, max_ranks=max_ranks)
elseif kind in ("ed-spectra", "ed_spectra", "spectrum")
    plot_ed_spectrum_results(; samples=samples)
elseif kind == "structure"
    plot_structure_factor_results(; samples=samples)
elseif kind == "diagnostics"
    plot_diagnostic_results(; samples=samples, max_flow_curves=max_flow_curves)
elseif kind in ("charge-gap", "charge_gap")
    plot_charge_gap_results()
else
    error("Unknown plot kind `$kind`.\n$HELP")
end
