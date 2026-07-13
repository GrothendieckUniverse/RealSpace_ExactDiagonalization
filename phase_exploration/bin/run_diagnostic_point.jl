#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Generate phase diagnostics (all-sector flow, charge pump, spatial ES, PES).

Usage:
  julia --project=. phase_exploration/bin/run_diagnostic_point.jl \\
    --phase FCI --geometry 3x5 [--x -1.0] [--mode auto] \\
    [--observables flow,pump,spatial_es,pes] [--manifold-size 3] \\
    [--zero-nev 10] [--flow-nev 3] [--flow-steps 25] [--pump-steps 17] \\
    [--pes-na 2] [--overwrite false]

Without --x, the common deep-phase point in src/config.jl is used.
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
x = haskey(opts, "x") ? parse(Float64, opts["x"]) : nothing
manifold_size = haskey(opts, "manifold-size") ? parse(Int, opts["manifold-size"]) : nothing
mode = get(opts, "mode", "auto") == "auto" ? mode_for(sample, :diagnostics) :
       CheckerboardPhaseStudy.parse_mode(opts["mode"])
observables = Symbol.(lowercase.(strip.(split(get(opts, "observables", "flow,pump,spatial_es,pes"), ','))))
zero_nev = parse(Int, get(opts, "zero-nev", "10"))
flow_nev = parse(Int, get(opts, "flow-nev", "3"))
flow_steps = parse(Int, get(opts, "flow-steps", "25"))
pump_steps = parse(Int, get(opts, "pump-steps", "17"))
flow_cycles = parse(Float64, get(opts, "flow-cycles", "3"))
pes_na = parse(Int, get(opts, "pes-na", "2"))
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
    overwrite=overwrite)
