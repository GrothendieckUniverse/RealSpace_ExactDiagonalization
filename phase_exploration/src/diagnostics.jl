function write_flow_csv(path, checkpoint_paths, flux_values, sectors, nev)
    ensure_parent(path)
    rows = NamedTuple[]
    for (itheta, checkpoint) in enumerate(checkpoint_paths)
        ed_data = load_checkpoint(checkpoint)
        local_rows = NamedTuple[]
        for sector in sectors
            idx = findfirst(irrep -> irrep.label == sector, ed_data.irrep_list)
            idx === nothing && continue
            haskey(ed_data.ed_scan_res, idx) || continue
            values = ed_data.ed_scan_res[idx][1]
            for level in 1:min(nev, length(values))
                push!(local_rows, (sector=sector, level=level, energy=Float64(values[level])))
            end
        end
        e0 = minimum(row.energy for row in local_rows)
        for row in local_rows
            push!(rows, (theta=Float64(flux_values[itheta]), sector=row.sector,
                level=row.level, energy=row.energy, shifted=row.energy - e0))
        end
    end
    open(path, "w") do io
        println(io, "flux_over_2pi,k1,k2,level,energy,energy_minus_flux_ground")
        for row in rows
            @printf(io, "%.16g,%d,%d,%d,%.16g,%.16g\n",
                row.theta, row.sector[1], row.sector[2], row.level, row.energy, row.shifted)
        end
    end
    return path
end

function write_pump_csv(path, pump)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "flux_over_2pi,branch,polarization,pumped_charge")
        for itheta in eachindex(pump.twisted_phases_over_2π_list),
            branch in axes(pump.pumped_charge_trajectories, 2)
            @printf(io, "%.16g,%d,%.16g,%.16g\n",
                pump.twisted_phases_over_2π_list[itheta], branch,
                pump.polarizations[itheta, branch],
                pump.pumped_charge_trajectories[itheta, branch])
        end
    end
    return path
end

function spatial_partition_a(model)
    sample = Tuple(Int.(model.lattice.sample_size))
    direction = sample[2] >= sample[1] ? 2 : 1
    cut = fld(sample[direction], 2)
    partition = Int[]
    for (site, (cell, _)) in enumerate(model.lattice.site_list)
        cell[direction] < cut && push!(partition, site)
    end
    isempty(partition) && error("Spatial cut produced an empty subsystem for $sample.")
    return partition, direction, cut
end

function write_spatial_es_csv(path, result)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "N_A,level,probability,entanglement_energy,dim_A,dim_B")
        for row in result.levels
            @printf(io, "%d,%d,%.16g,%.16g,%d,%d\n", row.n_a, row.level,
                row.probability, row.entanglement_energy, row.dim_a, row.dim_b)
        end
    end
    return path
end

function write_pes_csv(path, result)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "k1,k2,level,probability,entanglement_energy,sector_dim")
        for row in result.levels
            @printf(io, "%d,%d,%d,%.16g,%.16g,%d\n",
                row.momentum[1], row.momentum[2], row.level,
                row.probability, row.entanglement_energy, row.sector_dim)
        end
    end
    return path
end

function pes_gap_summary(levels)
    energies = sort([row.entanglement_energy for row in levels])
    length(energies) < 2 && return (largest_gap=NaN, levels_below=length(energies))
    gaps = diff(energies)
    idx = argmax(gaps)
    return (largest_gap=gaps[idx], levels_below=idx)
end

function run_phase_diagnostics(phase_name, sample::Tuple{Int,Int};
    x::Union{Nothing,Real}=nothing,
    manifold_size::Union{Nothing,Int}=nothing,
    mode::Symbol=mode_for(sample, :diagnostics),
    observables::Vector{Symbol}=[:flow, :pump, :spatial_es, :pes],
    zero_nev::Int=10,
    flow_nev::Int=3,
    flow_flux_values::Vector{Float64}=collect(range(0.0, 3.0; length=61)),
    pump_flux_values::Vector{Float64}=collect(range(0.0, 1.0; length=17)),
    flux_direction::Int=1,
    polarization_direction::Int=2,
    n_particles_a::Int=2,
    refresh::Bool=false,
    overwrite::Bool=false,
)
    valid = Set([:flow, :pump, :spatial_es, :pes])
    all(obs in valid for obs in observables) || error("Unknown diagnostic in $observables.")
    spec = phase_spec(phase_name)
    xvalue = isnothing(x) ? spec.numerator : Float64(x)
    nmanifold = isnothing(manifold_size) ? spec.manifold_size : manifold_size
    phase = String(spec.name)
    outdir = joinpath(RESULT_ROOT, "diagnostics", phase, geometry_tag(sample))
    checkpoint_point = "x_$(tpp_tag(xvalue))"
    ckpt_zero = joinpath(CHECKPOINT_ROOT, "diagnostics", phase, geometry_tag(sample),
        checkpoint_point, "zero_flux.jld2")
    ckpt_flux = joinpath(CHECKPOINT_ROOT, "diagnostics", phase, geometry_tag(sample),
        checkpoint_point, "flux")
    output_for = Dict(
        :flow => joinpath(outdir, "spectrum_flow.csv"),
        :pump => joinpath(outdir, "charge_pump.csv"),
        :spatial_es => joinpath(outdir, "spatial_entanglement_spectrum.csv"),
        :pes => joinpath(outdir, "particle_entanglement_spectrum.csv"),
    )
    # `refresh` rebuilds derived CSVs while reusing compatible ED checkpoints.
    # `overwrite` additionally discards and recomputes those checkpoints.
    todo = [obs for obs in observables if refresh || overwrite || !isfile(output_for[obs])]
    if isempty(todo) && isfile(joinpath(outdir, "summary.csv")) &&
       isfile(joinpath(outdir, "zero_flux_spectrum.csv"))
        @info "Requested diagnostics already complete; skipping" phase sample outdir
        return outdir
    end

    model, ed_data, n_particles, filling = build_checkerboard_problem(sample, xvalue)
    ed_data = scan_with_resume!(ed_data; nev=zero_nev, mode=mode,
        checkpoint_path=ckpt_zero, overwrite=overwrite)
    table = spectrum_table(ed_data)
    manifold_states = lowest_manifold_states(table; count=nmanifold)
    manifold_sectors = unique(row.sector for row in manifold_states)
    manifold_nev = maximum(row.level for row in manifold_states)
    all_sectors = all_sector_labels(ed_data)
    write_spectrum_csv(joinpath(outdir, "zero_flux_spectrum.csv"), table, sample, xvalue)

    union_flux = sort(unique(vcat(
        (:flow in todo ? flow_flux_values : Float64[]),
        (:pump in todo ? pump_flux_values : Float64[]))))
    flux_path_map = Dict{Float64,String}()
    if !isempty(union_flux)
        paths = ensure_flux_checkpoints(model, filling, union_flux;
            flux_direction=flux_direction,
            nev=max(flow_nev, manifold_nev, 2),
            mode=mode,
            checkpoint_dir=ckpt_flux,
            sectors=all_sectors,
            overwrite=overwrite)
        for (theta, path) in zip(union_flux, paths)
            flux_path_map[theta] = path
        end
    end

    if :flow in todo
        flow_paths = [flux_path_map[theta] for theta in flow_flux_values]
        write_flow_csv(output_for[:flow],
            flow_paths, flow_flux_values, all_sectors, flow_nev)
    end

    pump = nothing
    if :pump in todo
        # The checkpoints above are solver-mode agnostic. The toolbox pump call
        # sees complete files and only performs the polarization projection.
        pump = flux_charge_pump(model, manifold_sectors;
            filling_fraction=filling,
            flux_direction=flux_direction,
            polarization_direction=polarization_direction,
            twisted_phases_over_2π_list=pump_flux_values,
            manifold_states=manifold_states,
            checkpoint_dir=ckpt_flux,
            overwrite=false)
        write_pump_csv(output_for[:pump], pump)
    end

    spatial = nothing
    partition = Int[]
    cut_direction = 0
    cut_cell = 0
    if :spatial_es in todo
        partition, cut_direction, cut_cell = spatial_partition_a(model)
        spatial = entanglement_spectrum(model, table[1].sector;
            partition_a=partition,
            filling_fraction=filling,
            ed_mode=mode,
            ed_data=ed_data)
        write_spatial_es_csv(output_for[:spatial_es], spatial)
    end

    pes = nothing
    pes_summary = (largest_gap=NaN, levels_below=0)
    if :pes in todo
        n_particles_a <= n_particles ||
            error("PES N_A=$n_particles_a exceeds total N=$n_particles.")
        pes = particle_entanglement_spectrum(model, manifold_sectors;
            n_particles_a=n_particles_a,
            filling_fraction=filling,
            manifold_states=manifold_states,
            ed_mode=mode,
            ed_data=ed_data)
        write_pes_csv(output_for[:pes], pes)
        pes_summary = pes_gap_summary(pes.levels)
    end

    # Preserve summaries when a long run is resumed with only missing
    # observables left to compute.
    prior_path = joinpath(outdir, "summary.csv")
    prior = isfile(prior_path) ? read_simple_csv(prior_path)[1] : nothing
    if pump === nothing && isfile(output_for[:pump])
        pump_rows = read_simple_csv(output_for[:pump])
        final_flux = maximum(csv_float(row.flux_over_2pi) for row in pump_rows)
        pump_charges = [csv_float(row.pumped_charge) for row in pump_rows
            if isapprox(csv_float(row.flux_over_2pi), final_flux; atol=1e-12)]
        pump_sum_existing = sum(pump_charges)
        pump_abs_sum_existing = sum(abs, pump_charges)
    else
        pump_sum_existing = pump === nothing ? NaN : sum(pump.pumped_charges)
        pump_abs_sum_existing = pump === nothing ? NaN : sum(abs, pump.pumped_charges)
    end
    if pes === nothing && isfile(output_for[:pes])
        pes_rows = read_simple_csv(output_for[:pes])
        pes_energies = sort([csv_float(row.entanglement_energy) for row in pes_rows])
        if length(pes_energies) >= 2
            gaps = diff(pes_energies)
            gap_idx = argmax(gaps)
            pes_summary = (largest_gap=gaps[gap_idx], levels_below=gap_idx)
        end
        npes_existing = length(pes_rows)
    else
        npes_existing = pes === nothing ? 0 : length(pes.levels)
    end
    if spatial === nothing && prior !== nothing
        cut_direction = csv_int(prior.spatial_cut_direction)
        cut_cell = csv_int(prior.spatial_cut_cells)
    end

    ensure_parent(joinpath(outdir, "summary.csv"))
    open(joinpath(outdir, "summary.csv"), "w") do io
        println(io, "phase,L1,L2,n_sites,n_particles,tpp_numerator,tpp_actual,solver_mode,manifold_size,manifold_sectors,pumped_charge_sum,pumped_charge_abs_sum,spatial_cut_direction,spatial_cut_cells,PES_NA,PES_levels,PES_largest_gap,PES_levels_below_largest_gap,manifold_state_levels")
        sector_text = join(["$(row.sector[1]):$(row.sector[2])" for row in manifold_states], ';')
        state_text = join(["$(row.sector[1]):$(row.sector[2]):$(row.level)" for row in manifold_states], ';')
        @printf(io, "%s,%d,%d,%d,%d,%.16g,%.16g,%s,%d,%s,%.16g,%.16g,%d,%d,%d,%d,%.16g,%d,%s\n",
            phase, sample[1], sample[2], model.lattice.n_site, n_particles,
            xvalue, actual_tpp(xvalue), mode, nmanifold, sector_text,
            pump_sum_existing, pump_abs_sum_existing, cut_direction, cut_cell, n_particles_a,
            npes_existing, pes_summary.largest_gap, pes_summary.levels_below, state_text)
    end
    @info "Completed phase diagnostics" phase sample xvalue tpp_actual=actual_tpp(xvalue) observables outdir
    return outdir
end
