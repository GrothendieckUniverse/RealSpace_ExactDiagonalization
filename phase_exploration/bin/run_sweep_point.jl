#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Generate one independently resumable t′′ sweep point.

Usage:
  julia --project=. phase_exploration/bin/run_sweep_point.jl \\
    --geometry 3x5 --x -1.0 [--task all] [--mode auto] [--nev 10] \\
    [--dense-resolution 101] [--overwrite false]

`--x` is the numerator. Files record both x and the physical t′′=x/(2+2√2).
Task choices are all, spectrum, and structure.
"""

opts, _ = parse_cli(ARGS)
if haskey(opts, "help") || haskey(opts, "h")
    println(HELP)
    exit()
end
haskey(opts, "geometry") || error("Missing --geometry.\n$HELP")
haskey(opts, "x") || error("Missing --x.\n$HELP")
sample = parse_geometry(opts["geometry"])
x = parse(Float64, opts["x"])
task = Symbol(lowercase(get(opts, "task", "all")))
mode = get(opts, "mode", "auto") == "auto" ? mode_for(sample, :sweep) :
       CheckerboardPhaseStudy.parse_mode(opts["mode"])
nev = parse(Int, get(opts, "nev", "10"))
dense_resolution = parse(Int, get(opts, "dense-resolution", "101"))
overwrite = parse_bool(get(opts, "overwrite", "false"))

run_sweep_point(sample, x; task=task, nev=nev, mode=mode,
    dense_resolution=dense_resolution, overwrite=overwrite)
