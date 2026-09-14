#!/usr/bin/env julia

include(joinpath(@__DIR__, "_bootstrap.jl"))

const HELP = """
Generate phase diagnostics (structure factor, all-sector flow, charge pump,
and FCI particle entanglement spectrum).

Usage:
  julia --project=. phase_exploration/bin/run_diagnostic_point.jl \\
    --phase FCI --geometry 3x5 (--tpp -0.15 | --x -0.724264...) \\
    [--mode auto] [--observables structure,flow,pump[,pes]] \\
    [--manifold-size 3] [--zero-nev 10] [--flow-nev 4] [--flow-cycles 1] \\
    [--flow-steps 17] [--pump-steps 17] [--pump-cycles 1] [--flux-direction 1] [--polarization-direction 2] \\
    [--pes-na 2] [--dense-resolution 101] \\
    [--refresh false] [--overwrite false]

`--tpp` is the physical hopping. `--x` is the numerator in
t′′=x/(2+2sqrt(2)). Exactly one is required because every phase now has
multiple characteristic points in src/config.jl.
Flux cycles are in units of 2π; the default is 17 points including 0 and 1.
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
observables = haskey(opts, "observables") ?
    Symbol.(lowercase.(strip.(split(opts["observables"], ',')))) :
    default_diagnostic_observables(phase)
zero_nev = parse(Int, get(opts, "zero-nev", "10"))
flow_nev = parse(Int, get(opts, "flow-nev", "4"))
flow_steps = parse(Int, get(opts, "flow-steps", "17"))
pump_steps = parse(Int, get(opts, "pump-steps", "17"))
flow_cycles = parse(Float64, get(opts, "flow-cycles", "1"))
pump_cycles = parse(Float64, get(opts, "pump-cycles", "1"))
flow_flux = CheckerboardPhaseStudy.diagnostic_flux_grid(flow_steps, flow_cycles)
pump_flux = CheckerboardPhaseStudy.diagnostic_flux_grid(pump_steps, pump_cycles)
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
    flux_direction=parse(Int,get(opts,"flux-direction","1")),
    polarization_direction=parse(Int,get(opts,"polarization-direction","2")),
    flow_flux_values=flow_flux,
    pump_flux_values=pump_flux,
    n_particles_a=pes_na,
    dense_resolution=dense_resolution,
    refresh=refresh,
    overwrite=overwrite)
