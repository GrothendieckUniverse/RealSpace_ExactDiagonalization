# End-to-end flux CSV coverage, restart, history, and rendering checks.
# Run: julia --project=. --threads=2 test/phase_flux_protocol.jl
using Test
const PHASE_TEST_ROOT = mktempdir()
cp(joinpath(@__DIR__, "..", "phase_exploration", "src"), joinpath(PHASE_TEST_ROOT, "src"))
include(joinpath(PHASE_TEST_ROOT, "src", "CheckerboardPhaseStudy.jl"))
const PE = CheckerboardPhaseStudy
PE.run_phase_diagnostics(:FCI, (3,3); x=PE.numerator_at_tpp(-0.3), observables=[:flow,:pump], zero_nev=4, flow_nev=4, mode=:matrix)
@testset "HPC flux protocol and resumed pump" begin
    grid = PE.diagnostic_flux_grid()
    @test grid == collect(0.0:1/16:1.0)
    @test_throws ArgumentError PE.diagnostic_flux_grid(1)
    @test_throws ArgumentError PE.diagnostic_flux_grid(17, Inf)
    @test_throws ArgumentError PE.diagnostic_flux_grid(17, -1)
    dir = joinpath(PE.RESULT_ROOT,"diagnostics/FCI/3x3/tpp_m0p3000")
    flow = PE.read_simple_csv(joinpath(dir,"spectrum_flow.csv"))
    pump = PE.read_simple_csv(joinpath(dir,"charge_pump.csv"))
    @test sort(unique(PE.csv_float(r.flux_over_2pi) for r in flow)) == grid
    @test sort(unique(PE.csv_float(r.flux_over_2pi) for r in pump)) == grid
    @test length(flow) == 17*9*4
    @test length(pump) == 17*3
    for theta in grid
        rows = filter(r -> PE.csv_float(r.flux_over_2pi) == theta,flow)
        @test length(unique((r.k1,r.k2) for r in rows)) == 9
    end
    @test all(PE.csv_float(r.manifold_gap)>1e-9 for r in pump)
    @test all(PE.csv_float(r.min_position_singular_value)>1e-10 for r in pump)
    @test all(length(split(r.selected_states,';'))==3 for r in pump)
    charges = [PE.csv_float(r.pumped_charge) for r in pump if PE.csv_float(r.flux_over_2pi)==1.0]
    @test sum(charges) ≈ 1 atol=1e-9
    println("PUMP_BRANCHES=",charges,"; TOTAL=",sum(charges))
    println("MIN_GAP=",minimum(PE.csv_float(r.manifold_gap) for r in pump))
    println("MIN_POSITION_SV=",minimum(PE.csv_float(r.min_position_singular_value) for r in pump))
    flow_path = joinpath(dir, "spectrum_flow.csv")
    pump_path = joinpath(dir, "charge_pump.csv")
    @test PE.flux_output_complete(flow_path, grid; flux_direction=1, sample=(3,3), levels=4)
    @test PE.flux_output_complete(pump_path, grid; flux_direction=1, polarization_direction=2, manifold_size=3)
    @test !PE.flux_output_complete(flow_path, grid; flux_direction=2, sample=(3,3), levels=4)
    @test !PE.flux_output_complete(flow_path, grid[1:end-1]; flux_direction=1, sample=(3,3), levels=4)
    @test !PE.flux_output_complete(flow_path, grid; flux_direction=1, sample=(3,3), levels=5)
    @test !PE.flux_output_complete(pump_path, grid; flux_direction=1, polarization_direction=1, manifold_size=3)
    @test !PE.flux_output_complete(pump_path, grid; flux_direction=1, polarization_direction=2, manifold_size=2)
    old_flow, old_pump = read(flow_path), read(pump_path)
    before_mtime = stat(flow_path).mtime
    PE.run_phase_diagnostics(:FCI, (3,3); x=PE.numerator_at_tpp(-0.3), observables=[:flow,:pump], zero_nev=4,flow_nev=4,mode=:matrix)
    @test stat(flow_path).mtime == before_mtime
    PE.run_phase_diagnostics(:FCI, (3,3); x=PE.numerator_at_tpp(-0.3), observables=[:flow,:pump], zero_nev=4,flow_nev=4,mode=:matrix,refresh=true)
    archive_dirs = readdir(joinpath(dir,"history"); join=true)
    @test any(isfile(joinpath(d,"spectrum_flow.csv")) && read(joinpath(d,"spectrum_flow.csv")) == old_flow && read(joinpath(d,"charge_pump.csv")) == old_pump for d in archive_dirs)
    @test read(flow_path) == old_flow
    @test read(pump_path) == old_pump
    outputs = PE.plot_diagnostic_results(phases=[:FCI],samples=[(3,3)])
    @test length(outputs)==3
    @test all(isfile,outputs)
end
