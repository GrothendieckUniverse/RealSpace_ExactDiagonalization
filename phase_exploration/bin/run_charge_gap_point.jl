#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Compute Δc=E0(N+1)+E0(N-1)-2E0(N) for one phase and geometry.

Usage:
  julia --project=. phase_exploration/bin/run_charge_gap_point.jl \\
    --phase FCI --geometry 3x5 (--tpp -0.15 | --x -0.724264...) [--mode auto] \\
    [--nev 2] [--refresh false] [--overwrite false]

`--tpp` is the physical hopping; `--x` is its numerator. Exactly one is
required because each phase has multiple configured characteristic points.
"""

opts, _ = parse_cli(ARGS)
if haskey(opts, "help") || haskey(opts, "h")
    println(HELP)
    exit()
end
haskey(opts, "phase") || error("Missing --phase.\n$HELP")
haskey(opts, "geometry") || error("Missing --geometry.\n$HELP")
sample = parse_geometry(opts["geometry"])
phase = opts["phase"]
haskey(opts, "x") && haskey(opts, "tpp") && error("Choose only one of --tpp and --x.\n$HELP")
x = if haskey(opts, "tpp")
    numerator_at_tpp(parse(Float64, opts["tpp"]))
elseif haskey(opts, "x")
    parse(Float64, opts["x"])
else
    nothing
end
mode = get(opts, "mode", "auto") == "auto" ? mode_for(sample, :charge_gap) :
       CheckerboardPhaseStudy.parse_mode(opts["mode"])
nev = parse(Int, get(opts, "nev", "2"))
refresh = parse_bool(get(opts, "refresh", "false"))
overwrite = parse_bool(get(opts, "overwrite", "false"))

run_charge_gap_point(phase, sample; x=x, mode=mode, nev=nev, refresh=refresh,
    overwrite=overwrite)
