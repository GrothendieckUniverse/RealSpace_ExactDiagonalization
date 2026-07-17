#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Generate phase diagnostics (structure factor, all-sector flow, charge pump,
spatial ES, and PES).

Usage:
  julia --project=. phase_exploration/bin/run_diagnostic_point.jl \\
    --phase FCI --geometry 3x5 (--tpp -0.15 | --x -0.724264...) \\
    [--mode auto] [--observables structure,flow,pump,spatial_es,pes] \\
    [--manifold-size 3] [--zero-nev 10] [--flow-nev 4] [--flow-cycles 1] \\
    [--flow-steps 21] [--pump-steps 21] \\
    [--pes-na 2] [--dense-resolution 101] \\
    [--refresh false] [--overwrite false]

`--tpp` is the physical hopping. `--x` is the numerator in
t′′=x/(2+2sqrt(2)). Exactly one is required because every phase now has
multiple characteristic points in src/config.jl.
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
manifold_size = haskey(opts, "manifold-size") ? parse(Int, opts["manifold-size"]) : nothing
mode = get(opts, "mode", "auto") == "auto" ? mode_for(sample, :diagnostics) :
       CheckerboardPhaseStudy.parse_mode(opts["mode"])
observables = Symbol.(lowercase.(strip.(split(get(opts, "observables", "structure,flow,pump,spatial_es,pes"), ','))))
zero_nev = parse(Int, get(opts, "zero-nev", "10"))
flow_nev = parse(Int, get(opts, "flow-nev", "4"))
flow_steps = parse(Int, get(opts, "flow-steps", "21"))
pump_steps = parse(Int, get(opts, "pump-steps", "21"))
flow_cycles = parse(Float64, get(opts, "flow-cycles", "1"))
pes_na = parse(Int, get(opts, "pes-na", "2"))
dense_resolution = parse(Int, get(opts, "dense-resolution", "101"))
refresh = parse_bool(get(opts, "refresh", "false"))
overwrite = parse_bool(get(opts, "overwrite", "false"))

run_phase_diagnostics(phase, sample;
    x=x,
    manifold_size=manifold_size,
    mode=mode,
    observables=observables,
    zero_nev=zero_nev,
    flow_nev=flow_nev,
    flow_flux_values=collect(range(0.0, flow_cycles; length=flow_steps)),
    pump_flux_values=collect(range(0.0, 1.0; length=pump_steps)),
    n_particles_a=pes_na,
    dense_resolution=dense_resolution,
    refresh=refresh,
    overwrite=overwrite)
