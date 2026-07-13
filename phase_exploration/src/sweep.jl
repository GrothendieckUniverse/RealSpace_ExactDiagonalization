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

    requested_outputs_exist =
        (task == :spectrum && isfile(spectrum_path)) ||
        (task == :structure && isfile(metrics_path)) ||
        (task == :all && isfile(spectrum_path) && isfile(metrics_path))
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
    end

    write_key_values(joinpath(point_dir, "run_metadata.csv"), [
        "geometry" => geometry_tag(sample),
        "tpp_numerator" => x,
        "tpp_actual" => actual_tpp(x),
        "solver_mode" => mode,
        "nev_per_sector" => nev,
        "dense_k_resolution" => dense_resolution,
        "task" => task,
    ])
    @info "Completed sweep point" sample x tpp_actual=actual_tpp(x) mode point_dir
    return point_dir
end
