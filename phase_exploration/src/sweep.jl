function run_sweep_point(sample::Tuple{Int,Int}, x::Real;
    task::Symbol=:all,
    nev::Int=10,
    mode::Symbol=mode_for(sample, :sweep),
    dense_resolution::Int=101,
    overwrite::Bool=false,
)
    task in (:all, :spectrum, :structure) ||
        error("Sweep task must be all, spectrum, or structure; got $task.")
    x = Float64(x)
    point_dir = result_point_dir("sweep", sample, x)
    checkpoint = joinpath(CHECKPOINT_ROOT, "sweep", geometry_tag(sample),
        "x_$(tpp_tag(x))", "zero_flux.jld2")
    spectrum_path = joinpath(point_dir, "spectrum.csv")
    metrics_path = joinpath(point_dir, "structure_metrics.csv")
    gsd_allowed_path = joinpath(point_dir, "structure_fci_gsd_allowed.csv")
    gsd_dense_path = joinpath(point_dir, "structure_fci_gsd_dense.csv")
    gsd_metrics_path = joinpath(point_dir, "structure_fci_gsd_metrics.csv")
    gsd_state_metrics_path = joinpath(point_dir, "structure_fci_gsd_state_metrics.csv")
    structure_outputs = [joinpath(point_dir, "structure_allowed.csv"),
        joinpath(point_dir, "structure_dense.csv"), metrics_path,
        gsd_allowed_path, gsd_dense_path,
        gsd_metrics_path, gsd_state_metrics_path]

    requested_outputs_exist =
        (task == :spectrum && isfile(spectrum_path)) ||
        (task == :structure && all(isfile, structure_outputs)) ||
        (task == :all && isfile(spectrum_path) && all(isfile, structure_outputs))
    if requested_outputs_exist && !overwrite
        @info "Sweep point already complete; skipping" sample x task point_dir
        return point_dir
    end

    model, ed_data, _, filling = build_checkerboard_problem(sample, x)
    ed_data = scan_with_resume!(ed_data; nev=nev, mode=mode,
        checkpoint_path=checkpoint, overwrite=overwrite)
    table = spectrum_table(ed_data)

    if task in (:all, :spectrum)
        write_spectrum_csv(spectrum_path, table, sample, x)
    end

    if task in (:all, :structure)
        ground_sector = table[1].sector
        qx, qy, allowed = structure_factor_allowed_momenta(
            model, ground_sector;
            filling_fraction=filling,
            ed_mode=mode,
            ed_data=ed_data,
        )
        kx, ky, dense = compute_structure_factor_map(
            model, ground_sector;
            filling_fraction=filling,
            k_resolution=dense_resolution,
            ed_data=ed_data,
        )
        metrics = sf_metrics(qx, qy, allowed)
        write_allowed_sf_csv(joinpath(point_dir, "structure_allowed.csv"), qx, qy, allowed)
        write_dense_sf_csv(joinpath(point_dir, "structure_dense.csv"), kx, ky, dense)
        write_sf_metrics_csv(metrics_path, metrics, sample, x, ground_sector)

        reference_specs = get(FCI_REFERENCE_MANIFOLD, sample,
            Tuple{Tuple{Int,Int},Int}[])
        isempty(reference_specs) && error("No FCI reference manifold is configured for $sample.")
        reference_states = NamedTuple[]
        reference_ranks = Int[]
        for spec in reference_specs
            rank = findfirst(row -> (row.sector, row.level) == spec, table)
            rank === nothing && error("FCI reference state $spec is absent at $sample, x=$x.")
            push!(reference_states, table[rank])
            push!(reference_ranks, rank)
        end

        state_metric_rows = NamedTuple[]
        for (state, rank) in zip(reference_states, reference_ranks)
            state_qx, state_qy, state_values = structure_factor_allowed_momenta(
                model, state.sector;
                target_eigval_idx=state.level,
                filling_fraction=filling,
                ed_mode=mode,
                ed_data=ed_data)
            push!(state_metric_rows, (state=state, rank=rank,
                metrics=sf_metrics(state_qx, state_qy, state_values)))
        end
        write_sf_state_metrics_csv(gsd_state_metrics_path, state_metric_rows, sample, x;
            selection="fixed_fci_reference_states")

        gsd_qx, gsd_qy, gsd_allowed = structure_factor_manifold_allowed_momenta(
            model, reference_states;
            filling_fraction=filling,
            ed_mode=mode,
            ed_data=ed_data)
        gsd_kx, gsd_ky, gsd_dense = compute_structure_factor_manifold_average_map(
            model, reference_states;
            filling_fraction=filling,
            k_resolution=dense_resolution,
            ed_mode=mode,
            ed_data=ed_data)
        gsd_metrics = sf_metrics(gsd_qx, gsd_qy, gsd_allowed)
        write_allowed_sf_csv(gsd_allowed_path, gsd_qx, gsd_qy, gsd_allowed)
        write_dense_sf_csv(gsd_dense_path, gsd_kx, gsd_ky, gsd_dense)
        write_manifold_sf_metrics_csv(gsd_metrics_path, gsd_metrics, sample, x,
            reference_states; selection="fixed_fci_reference_states")
    end

    write_key_values(joinpath(point_dir, "run_metadata.csv"), [
        "geometry" => geometry_tag(sample),
        "tpp_numerator" => x,
        "tpp_actual" => actual_tpp(x),
        "solver_mode" => mode,
        "nev_per_sector" => nev,
        "dense_k_resolution" => dense_resolution,
        "task" => task,
        "structure_ensembles" => "absolute_ground_state;fixed_fci_reference_projector",
    ])
    @info "Completed sweep point" sample x tpp_actual=actual_tpp(x) mode point_dir
    return point_dir
end
