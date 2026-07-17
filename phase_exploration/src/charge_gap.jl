function ground_energy_for_particle_number(sample, x, n_particles;
    mode, nev, checkpoint, overwrite, seed_checkpoint=nothing)
    _, ed_data, _, _ = build_checkerboard_problem(sample, x; n_particles=n_particles)
    seed_paths = seed_checkpoint === nothing ? String[] : [String(seed_checkpoint)]
    ed_data = scan_with_resume!(ed_data; nev=nev, mode=mode,
        checkpoint_path=checkpoint, seed_checkpoint_paths=seed_paths,
        overwrite=overwrite)
    table = spectrum_table(ed_data)
    return table[1].energy, table[1].sector
end

function run_charge_gap_point(phase_name, sample::Tuple{Int,Int};
    x::Union{Nothing,Real}=nothing,
    mode::Symbol=mode_for(sample, :charge_gap),
    nev::Int=2,
    refresh::Bool=false,
    overwrite::Bool=false,
)
    spec = phase_spec(phase_name)
    phase = String(spec.name)
    xvalue = require_phase_numerator(phase_name, x)
    n0 = default_particle_number(sample)
    outdir = phase_result_point_dir("charge_gap", phase, sample, xvalue)
    ckptdir = joinpath(CHECKPOINT_ROOT, "charge_gap", phase, geometry_tag(sample),
        "x_$(tpp_tag(xvalue))")
    summary_path = joinpath(outdir, "charge_gap.csv")
    if isfile(summary_path) && !refresh && !overwrite
        @info "Charge-gap point already complete; skipping" phase sample summary_path
        return summary_path
    end

    energies = Dict{Int,Float64}()
    sectors = Dict{Int,Tuple{Int,Int}}()
    for np in (n0 - 1, n0, n0 + 1)
        checkpoint = joinpath(ckptdir, "N_$(np).jld2")
        seed_checkpoint = np == n0 ? joinpath(CHECKPOINT_ROOT, "sweep",
            geometry_tag(sample), "x_$(tpp_tag(xvalue))", "zero_flux.jld2") : nothing
        energy, sector = ground_energy_for_particle_number(sample, xvalue, np;
            mode=mode, nev=nev, checkpoint=checkpoint, overwrite=overwrite,
            seed_checkpoint=seed_checkpoint)
        energies[np] = energy
        sectors[np] = sector
    end
    gap = energies[n0 + 1] + energies[n0 - 1] - 2 * energies[n0]

    ensure_parent(summary_path)
    open(summary_path, "w") do io
        println(io, "phase,L1,L2,n_unit_cells,n_sites,N0,tpp_numerator,tpp_actual,solver_mode,E_Nminus,ground_k1_Nminus,ground_k2_Nminus,E_N,ground_k1_N,ground_k2_N,E_Nplus,ground_k1_Nplus,ground_k2_Nplus,charge_gap")
        @printf(io, "%s,%d,%d,%d,%d,%d,%.16g,%.16g,%s,%.16g,%d,%d,%.16g,%d,%d,%.16g,%d,%d,%.16g\n",
            phase, sample[1], sample[2], prod(sample), 2 * prod(sample), n0,
            xvalue, actual_tpp(xvalue), mode,
            energies[n0 - 1], sectors[n0 - 1][1], sectors[n0 - 1][2],
            energies[n0], sectors[n0][1], sectors[n0][2],
            energies[n0 + 1], sectors[n0 + 1][1], sectors[n0 + 1][2], gap)
    end
    @info "Completed charge-gap point" phase sample xvalue tpp_actual=actual_tpp(xvalue) gap
    return summary_path
end
