function existing_sweep_points(sample::Tuple{Int,Int})
    root = joinpath(RESULT_ROOT, "sweep", geometry_tag(sample))
    isdir(root) || return String[]
    dirs = [joinpath(root, entry) for entry in readdir(root) if isdir(joinpath(root, entry))]
    function point_x(dir)
        metrics = joinpath(dir, "structure_metrics.csv")
        spectrum = joinpath(dir, "spectrum.csv")
        if isfile(metrics)
            return csv_float(read_simple_csv(metrics)[1].tpp_actual)
        elseif isfile(spectrum)
            return csv_float(read_simple_csv(spectrum)[1].tpp_actual)
        end
        return Inf
    end
    sort!(dirs; by=point_x)
    return dirs
end

function sweep_spectrum_series(sample::Tuple{Int,Int}; max_ranks::Int=12)
    series = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for dir in existing_sweep_points(sample)
        path = joinpath(dir, "spectrum.csv")
        isfile(path) || continue
        rows = read_simple_csv(path)
        for row in rows
            rank = csv_int(row.rank)
            rank <= max_ranks || continue
            push!(get!(series, rank, Tuple{Float64,Float64}[]),
                (csv_float(row.tpp_actual), csv_float(row.energy_minus_E0)))
        end
    end
    for values in values(series)
        sort!(values; by=first)
    end
    return series
end

function sweep_metric_series(sample::Tuple{Int,Int}, field::Symbol)
    values = Tuple{Float64,Float64}[]
    for dir in existing_sweep_points(sample)
        path = joinpath(dir, "structure_metrics.csv")
        isfile(path) || continue
        row = read_simple_csv(path)[1]
        push!(values, (csv_float(row.tpp_actual), csv_float(getproperty(row, field))))
    end
    sort!(values; by=first)
    return values
end

function draw_spectrum_axis!(ax, sample; max_ranks=12, legend=false)
    series = sweep_spectrum_series(sample; max_ranks=max_ranks)
    for rank in sort(collect(keys(series)))
        points = series[rank]
        isempty(points) && continue
        x = first.(points)
        y = last.(points)
        lines!(ax, x, y; color=Makie.Cycled(rank), linewidth=1.8,
            label="rank $rank")
        scatter!(ax, x, y; color=Makie.Cycled(rank), markersize=5)
    end
    legend && !isempty(series) && axislegend(ax; position=:lt, nbanks=2)
    return ax
end

function draw_metric_axis!(ax, sample, field; color=:royalblue3)
    points = sweep_metric_series(sample, field)
    isempty(points) && return ax
    lines!(ax, first.(points), last.(points); color=color, linewidth=2)
    scatter!(ax, first.(points), last.(points); color=color, markersize=7)
    return ax
end

function plot_sweep_results(; samples=STUDY_GEOMETRIES, max_ranks::Int=12)
    outdir = joinpath(FIGURE_ROOT, "sweep")
    mkpath(outdir)
    outputs = String[]

    for sample in samples
        isempty(existing_sweep_points(sample)) && continue
        tag = geometry_tag(sample)

        fig_s = Figure(size=(850, 580))
        ax_s = Axis(fig_s[1, 1]; xlabel="physical t′′ = x/(2+2√2)", ylabel="E - E₀",
            title="Lowest all-sector spectrum ranks — $tag")
        draw_spectrum_axis!(ax_s, sample; max_ranks=max_ranks, legend=true)
        path_s = joinpath(outdir, "spectrum_ranks_$(tag).svg")
        save(path_s, fig_s)
        push!(outputs, path_s)

        fig_m = Figure(size=(760, 520))
        ax_m = Axis(fig_m[1, 1]; xlabel="physical t′′ = x/(2+2√2)", ylabel="max |S(q)|",
            title="Structure-factor maximum — $tag")
        draw_metric_axis!(ax_m, sample, :max_abs_S; color=:darkorange2)
        path_m = joinpath(outdir, "max_abs_structure_factor_$(tag).svg")
        save(path_m, fig_m)
        push!(outputs, path_m)

        fig_r = Figure(size=(760, 520))
        ax_r = Axis(fig_r[1, 1]; xlabel="physical t′′ = x/(2+2√2)",
            ylabel="max |S(q)/mean(S(q))|",
            title="Normalized structure-factor maximum — $tag")
        draw_metric_axis!(ax_r, sample, :max_abs_S_over_mean_S; color=:seagreen4)
        path_r = joinpath(outdir, "max_abs_structure_factor_over_mean_$(tag).svg")
        save(path_r, fig_r)
        push!(outputs, path_r)
    end

    active = [sample for sample in samples if !isempty(existing_sweep_points(sample))]
    if !isempty(active)
        ncol = 2
        nrow = cld(length(active), ncol)
        fig = Figure(size=(760 * ncol, 520 * nrow))
        for (idx, sample) in enumerate(active)
            row, col = fldmod1(idx, ncol)
            ax = Axis(fig[row, col]; xlabel="physical t′′", ylabel="E - E₀",
                title=geometry_tag(sample))
            draw_spectrum_axis!(ax, sample; max_ranks=max_ranks)
        end
        Label(fig[0, 1:ncol], "Lowest all-sector spectrum ranks (x-axis is the physical hopping)"; fontsize=20)
        path = joinpath(outdir, "spectrum_ranks_all_geometries.svg")
        save(path, fig)
        push!(outputs, path)

        for (field, label, filename, color) in [
            (:max_abs_S, "max |S(q)|", "max_abs_structure_factor_all_geometries.svg", :darkorange2),
            (:max_abs_S_over_mean_S, "max |S(q)/mean(S(q))|", "max_abs_structure_factor_over_mean_all_geometries.svg", :seagreen4),
        ]
            figm = Figure(size=(760 * ncol, 500 * nrow))
            for (idx, sample) in enumerate(active)
                row, col = fldmod1(idx, ncol)
                ax = Axis(figm[row, col]; xlabel="physical t′′", ylabel=label,
                    title=geometry_tag(sample))
                draw_metric_axis!(ax, sample, field; color=color)
            end
            Label(figm[0, 1:ncol], "$label versus physical t′′"; fontsize=20)
            path = joinpath(outdir, filename)
            save(path, figm)
            push!(outputs, path)
        end
    end
    @info "Wrote sweep figures" count=length(outputs) outdir
    return outputs
end

function dense_sf_arrays(path)
    rows = read_simple_csv(path)
    kx = sort(unique(csv_float(row.kx) for row in rows))
    ky = sort(unique(csv_float(row.ky) for row in rows))
    ix = Dict(x => i for (i, x) in enumerate(kx))
    iy = Dict(y => i for (i, y) in enumerate(ky))
    values = fill(NaN, length(kx), length(ky))
    for row in rows
        values[ix[csv_float(row.kx)], iy[csv_float(row.ky)]] = csv_float(row.S_q)
    end
    return kx, ky, values
end

function draw_square_bz!(ax)
    lines!(ax, [-π, π, π, -π, -π], [-π, -π, π, π, -π]; color=:black, linewidth=2)
end

function plot_structure_factor_results(; samples=STUDY_GEOMETRIES)
    outputs = String[]
    for sample in samples
        outdir = joinpath(FIGURE_ROOT, "structure_factor", geometry_tag(sample))
        for point in existing_sweep_points(sample)
            allowed_path = joinpath(point, "structure_allowed.csv")
            dense_path = joinpath(point, "structure_dense.csv")
            metrics_path = joinpath(point, "structure_metrics.csv")
            all(isfile, (allowed_path, dense_path, metrics_path)) || continue
            metrics = read_simple_csv(metrics_path)[1]
            xactual = csv_float(metrics.tpp_actual)
            xnumerator = csv_float(metrics.tpp_numerator)
            allowed = read_simple_csv(allowed_path)
            qx = [csv_float(row.qx) for row in allowed]
            qy = [csv_float(row.qy) for row in allowed]
            sq = [csv_float(row.S_q) for row in allowed]
            kx, ky, dense = dense_sf_arrays(dense_path)
            lo = min(minimum(sq), minimum(dense))
            hi = max(maximum(sq), maximum(dense))

            fig = Figure(size=(1120, 500))
            ax1 = Axis(fig[1, 1]; xlabel="qₓ", ylabel="qᵧ", title="finite-torus momenta",
                aspect=DataAspect())
            sc = scatter!(ax1, qx, qy; color=sq, colorrange=(lo, hi), colormap=:viridis,
                markersize=30)
            draw_square_bz!(ax1)
            ax2 = Axis(fig[1, 2]; xlabel="qₓ", ylabel="qᵧ", title="dense q grid",
                aspect=DataAspect())
            heatmap!(ax2, kx, ky, dense; colorrange=(lo, hi), colormap=:viridis)
            draw_square_bz!(ax2)
            Colorbar(fig[1, 3], sc; label="connected S(q)")
            Label(fig[0, 1:2], @sprintf("%s, physical t′′ = %.8f (x = %.3f)",
                geometry_tag(sample), xactual, xnumerator); fontsize=18)
            mkpath(outdir)
            path = joinpath(outdir, "x_$(tpp_tag(xnumerator))_allowed_and_dense.svg")
            save(path, fig)
            push!(outputs, path)
        end
    end
    @info "Wrote structure-factor maps" count=length(outputs)
    return outputs
end

function group_rows(rows, fields::Tuple)
    groups = Dict{Tuple,Vector{Any}}()
    for row in rows
        key = Tuple(getproperty(row, field) for field in fields)
        push!(get!(groups, key, Any[]), row)
    end
    return groups
end

function plot_diagnostic_results(; phases=[:AHC, :FCI, :CDW], samples=STUDY_GEOMETRIES)
    outputs = String[]
    for phase_symbol in phases, sample in samples
        phase = String(phase_symbol)
        datadir = joinpath(RESULT_ROOT, "diagnostics", phase, geometry_tag(sample))
        isdir(datadir) || continue
        outdir = joinpath(FIGURE_ROOT, "diagnostics", phase, geometry_tag(sample))
        mkpath(outdir)
        title_prefix = "$phase — $(geometry_tag(sample))"

        flow_path = joinpath(datadir, "spectrum_flow.csv")
        if isfile(flow_path)
            rows = read_simple_csv(flow_path)
            groups = group_rows(rows, (:k1, :k2, :level))
            fig = Figure(size=(900, 600))
            ax = Axis(fig[1, 1]; xlabel="inserted flux / 2π", ylabel="E - E₀(θ)",
                title="$title_prefix spectrum flow (all momentum sectors)")
            sectors = sort(unique((csv_int(row.k1), csv_int(row.k2)) for row in rows))
            sector_index = Dict(k => i for (i, k) in enumerate(sectors))
            for (key, values) in groups
                sort!(values; by=row -> csv_float(row.flux_over_2pi))
                sector = (csv_int(key[1]), csv_int(key[2]))
                level = csv_int(key[3])
                color = Makie.Cycled(sector_index[sector])
                alpha = 1 / sqrt(level)
                lines!(ax, [csv_float(row.flux_over_2pi) for row in values],
                    [csv_float(row.energy_minus_flux_ground) for row in values];
                    color=color, alpha=alpha, linewidth=1.5)
            end
            path = joinpath(outdir, "spectrum_flow.svg")
            save(path, fig)
            push!(outputs, path)
        end

        pump_path = joinpath(datadir, "charge_pump.csv")
        if isfile(pump_path)
            rows = read_simple_csv(pump_path)
            groups = group_rows(rows, (:branch,))
            fig = Figure(size=(820, 560))
            ax = Axis(fig[1, 1]; xlabel="inserted flux / 2π", ylabel="pumped charge ΔQ",
                title="$title_prefix charge pump")
            for (key, values) in sort(collect(groups); by=x -> csv_int(x[1][1]))
                branch = csv_int(key[1])
                sort!(values; by=row -> csv_float(row.flux_over_2pi))
                x = [csv_float(row.flux_over_2pi) for row in values]
                y = [csv_float(row.pumped_charge) for row in values]
                lines!(ax, x, y; color=Makie.Cycled(branch), linewidth=2, label="branch $branch")
                scatter!(ax, x, y; color=Makie.Cycled(branch), markersize=6)
            end
            axislegend(ax; position=:lt)
            path = joinpath(outdir, "charge_pump.svg")
            save(path, fig)
            push!(outputs, path)
        end

        spatial_path = joinpath(datadir, "spatial_entanglement_spectrum.csv")
        if isfile(spatial_path)
            rows = read_simple_csv(spatial_path)
            fig = Figure(size=(760, 540))
            ax = Axis(fig[1, 1]; xlabel="particles in spatial/orbital region A (N_A)",
                ylabel="ξ = -log(λ)", title="$title_prefix spatial-orbital ES")
            x = [csv_int(row.N_A) for row in rows]
            y = [csv_float(row.entanglement_energy) for row in rows]
            scatter!(ax, x, y; color=y, colormap=:viridis, markersize=7)
            path = joinpath(outdir, "spatial_entanglement_spectrum.svg")
            save(path, fig)
            push!(outputs, path)
        end

        pes_path = joinpath(datadir, "particle_entanglement_spectrum.csv")
        if isfile(pes_path)
            rows = read_simple_csv(pes_path)
            momenta = sort(unique((csv_int(row.k1), csv_int(row.k2)) for row in rows))
            index = Dict(k => i for (i, k) in enumerate(momenta))
            x = [index[(csv_int(row.k1), csv_int(row.k2))] for row in rows]
            y = [csv_float(row.entanglement_energy) for row in rows]
            fig = Figure(size=(1000, 560))
            ax = Axis(fig[1, 1]; xlabel="subsystem momentum sector", ylabel="ξ = -log(λ)",
                title="$title_prefix momentum-resolved PES",
                xticks=(1:length(momenta), ["($(k[1]),$(k[2]))" for k in momenta]),
                xticklabelrotation=π / 3)
            scatter!(ax, x, y; color=y, colormap=:viridis, markersize=7)
            path = joinpath(outdir, "particle_entanglement_spectrum.svg")
            save(path, fig)
            push!(outputs, path)
        end
    end
    @info "Wrote diagnostic figures" count=length(outputs)
    return outputs
end

function charge_gap_rows(phase)
    root = joinpath(RESULT_ROOT, "charge_gap", String(phase))
    isdir(root) || return NamedTuple[]
    rows = NamedTuple[]
    for entry in readdir(root)
        path = joinpath(root, entry, "charge_gap.csv")
        isfile(path) && append!(rows, read_simple_csv(path))
    end
    sort!(rows; by=row -> csv_int(row.n_sites))
    return rows
end

function plot_charge_gap_results(; phases=[:AHC, :FCI, :CDW])
    outputs = String[]
    outdir = joinpath(FIGURE_ROOT, "charge_gap")
    mkpath(outdir)
    for phase in phases
        rows = charge_gap_rows(phase)
        length(rows) >= 2 || continue
        x = [1 / csv_int(row.n_sites) for row in rows]
        y = [csv_float(row.charge_gap) for row in rows]
        design = hcat(ones(length(x)), x)
        fit = design \ y
        intercept, slope = fit
        residual = sqrt(mean((design * fit .- y).^2))
        labels = ["$(row.L1)x$(row.L2)" for row in rows]

        fig = Figure(size=(760, 540))
        ax = Axis(fig[1, 1]; xlabel="1 / N_sites", ylabel="charge gap Δc",
            title="$(String(phase)) charge-gap finite-size scaling")
        line_x = collect(range(0.0, maximum(x); length=200))
        lines!(ax, line_x, intercept .+ slope .* line_x; color=:firebrick3, linewidth=2,
            label=@sprintf("linear: Δ∞ = %.6f", intercept))
        scatter!(ax, x, y; color=:royalblue3, markersize=12)
        scatter!(ax, [0.0], [intercept]; color=:black, marker=:diamond, markersize=14)
        for i in eachindex(x)
            text!(ax, x[i], y[i]; text=labels[i], offset=(7, 7), fontsize=12)
        end
        axislegend(ax; position=:lt)
        path = joinpath(outdir, "$(lowercase(String(phase)))_charge_gap_scaling.svg")
        save(path, fig)
        push!(outputs, path)

        fit_path = joinpath(RESULT_ROOT, "charge_gap", String(phase), "finite_size_fit.csv")
        ensure_parent(fit_path)
        open(fit_path, "w") do io
            println(io, "phase,tpp_numerator,tpp_actual,n_sizes,intercept_at_inverse_size_zero,slope,rms_residual")
            @printf(io, "%s,%.16g,%.16g,%d,%.16g,%.16g,%.16g\n", String(phase),
                csv_float(rows[1].tpp_numerator), csv_float(rows[1].tpp_actual), length(rows),
                intercept, slope, residual)
        end
    end
    @info "Wrote charge-gap scaling figures" count=length(outputs)
    return outputs
end

function plot_all_results()
    return vcat(
        plot_sweep_results(),
        plot_structure_factor_results(),
        plot_diagnostic_results(),
        plot_charge_gap_results(),
    )
end
