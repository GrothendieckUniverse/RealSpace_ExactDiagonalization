format_tpp(value::Real) = @sprintf("%.2f", Float64(value))
format_tpp_ticks(values) = format_tpp.(values)
const FOCUS_STATE_COLORS = [:firebrick3, :royalblue3, :seagreen4]
const OTHER_FLOW_COLOR = :gray35
const SpectrumSweepPoint = NamedTuple{
    (:tpp, :energy, :sector, :level),Tuple{Float64,Float64,Tuple{Int,Int},Int}}
const SweepGapPoint = NamedTuple{
    (:tpp, :fci_manifold_width, :roton_gap),Tuple{Float64,Float64,Float64}}

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

function sweep_spectrum_series(sample::Tuple{Int,Int}; max_ranks::Int=20)
    series = Dict{Int,Vector{SpectrumSweepPoint}}()
    references = Set(get(FCI_REFERENCE_MANIFOLD, sample,
        Tuple{Tuple{Int,Int},Int}[]))
    for dir in existing_sweep_points(sample)
        path = joinpath(dir, "spectrum.csv")
        isfile(path) || continue
        rows = read_simple_csv(path)
        for row in rows
            rank = csv_int(row.rank)
            state = ((csv_int(row.k1), csv_int(row.k2)),
                csv_int(row.level_in_sector))
            (rank <= max_ranks || state in references) || continue
            point = (tpp=csv_float(row.tpp_actual),
                energy=csv_float(row.energy_minus_E0),
                sector=state[1], level=state[2])
            push!(get!(series, rank, SpectrumSweepPoint[]), point)
        end
    end
    for values in values(series)
        sort!(values; by=point -> point.tpp)
    end
    return series
end

"Zero-flux FCI-manifold width and roton gap from global energies E₀ ≤ E₁ ≤ ... ."
function sweep_gap_series(sample::Tuple{Int,Int})
    points = SweepGapPoint[]

    for dir in existing_sweep_points(sample)
        path = joinpath(dir, "spectrum.csv")
        isfile(path) || continue
        rows = read_simple_csv(path)
        sort!(rows; by=row -> csv_int(row.rank))
        length(rows) >= 4 || continue
        ranked_energies = [csv_float(row.energy) for row in rows]
        e0, e2, e3 = ranked_energies[1], ranked_energies[3], ranked_energies[4]
        push!(points, (
            tpp=csv_float(rows[1].tpp_actual),
            fci_manifold_width=e2 - e0,
            roton_gap=e3 - e2,
        ))
    end
    sort!(points; by=point -> point.tpp)
    return points
end

function sweep_metric_series(sample::Tuple{Int,Int}, field::Symbol; ensemble::Symbol=:ground)
    values = Tuple{Float64,Float64}[]
    filename = if ensemble in (:ground, :ground_total)
        "structure_metrics.csv"
    elseif ensemble == :ground_aa
        "structure_ground_aa_metrics.csv"
    elseif ensemble == :ground_ab
        "structure_ground_ab_metrics.csv"
    elseif ensemble in (:fci_gsd, :fci_projector)
        "structure_fci_gsd_metrics.csv"
    else
        error("Unknown structure-factor ensemble $ensemble.")
    end
    for dir in existing_sweep_points(sample)
        path = joinpath(dir, filename)
        isfile(path) || continue
        row = read_simple_csv(path)[1]
        push!(values, (csv_float(row.tpp_actual), csv_float(getproperty(row, field))))
    end
    sort!(values; by=first)
    return values
end

function metric_series_in_window(series, window=FCI_PROJECTOR_PLOT_WINDOW)
    lower, upper = window
    return [point for point in series if lower <= first(point) <= upper]
end

function draw_ground_sublattice_max_axis!(ax, sample)
    plotted = false
    for (ensemble, label, color) in [
        (:ground_aa, "max |Sᴬᴬ(q)|", :darkorange2),
        (:ground_ab, "max |Re Sᴬᴮ(q)|", :mediumpurple3),
    ]
        series = sweep_metric_series(sample, :max_abs_S; ensemble=ensemble)
        isempty(series) && continue
        lines!(ax, first.(series), last.(series); color=color, linewidth=2.0,
            label=label)
        scatter!(ax, first.(series), last.(series); color=color, markersize=6)
        plotted = true
    end
    switches = sweep_ground_state_switches(sample)
    isempty(switches) || vlines!(ax, switches; color=(:black, 0.22),
        linestyle=:dot, linewidth=1.0)
    plotted && axislegend(ax; position=:rt, labelsize=10)
    return ax
end

function draw_fci_projector_zoom_axis!(ax, sample;
    window=FCI_PROJECTOR_PLOT_WINDOW, label=geometry_tag(sample), color=:royalblue3)
    series = metric_series_in_window(
        sweep_metric_series(sample, :max_abs_S; ensemble=:fci_projector), window)
    isempty(series) && return false
    lines!(ax, first.(series), last.(series); color=color, linewidth=2.2, label=label)
    scatter!(ax, first.(series), last.(series); color=color, markersize=7)
    xlims!(ax, window...)
    return true
end

function sweep_ground_state_switches(sample::Tuple{Int,Int})
    switches = Float64[]
    previous = nothing
    for dir in existing_sweep_points(sample)
        path = joinpath(dir, "spectrum.csv")
        isfile(path) || continue
        rows = read_simple_csv(path)
        isempty(rows) && continue
        ground = rows[argmin([csv_int(row.rank) for row in rows])]
        state = (csv_int(ground.k1), csv_int(ground.k2), csv_int(ground.level_in_sector))
        tpp = csv_float(ground.tpp_actual)
        previous !== nothing && state != previous && push!(switches, tpp)
        previous = state
    end
    return switches
end

function draw_spectrum_axis!(ax, sample; max_ranks=20, series=nothing)
    series === nothing && (series = sweep_spectrum_series(sample; max_ranks=max_ranks))
    references = get(FCI_REFERENCE_MANIFOLD, sample, Tuple{Tuple{Int,Int},Int}[])
    length(references) <= length(FOCUS_STATE_COLORS) || error(
        "FCI manifold has $(length(references)) states, but only " *
        "$(length(FOCUS_STATE_COLORS)) highlight colors are configured.")
    reference_colors = Dict(reference => FOCUS_STATE_COLORS[index]
        for (index, reference) in enumerate(references))
    all_points = [point for points in values(series) for point in points]
    state_groups = Dict{Tuple{Tuple{Int,Int},Int},Vector{SpectrumSweepPoint}}()
    for point in all_points
        push!(get!(state_groups, (point.sector, point.level), SpectrumSweepPoint[]), point)
    end
    sweep_index = Dict(tpp => index
        for (index, tpp) in enumerate(sort!(unique(point.tpp for point in all_points))))

    function draw_state!(points, color; focused=false)
        sort!(points; by=point -> point.tpp)
        segments = Vector{Vector{SpectrumSweepPoint}}()
        segment = SpectrumSweepPoint[]
        previous_index = nothing
        for point in points
            point_index = sweep_index[point.tpp]
            if !isnothing(previous_index) && point_index != previous_index + 1
                push!(segments, segment)
                segment = SpectrumSweepPoint[]
            end
            push!(segment, point)
            previous_index = point_index
        end
        isempty(segment) || push!(segments, segment)
        for values in segments
            length(values) >= 2 || continue
            lines!(ax, [point.tpp for point in values], [point.energy for point in values];
                color=color, linewidth=focused ? 1.6 : 1.0,
                alpha=focused ? 1.0 : 0.82)
        end
        scatter!(ax, [point.tpp for point in points], [point.energy for point in points];
            color=color, markersize=focused ? 12 : 10,
            alpha=focused ? 1.0 : 0.82,
            strokecolor=focused ? :black : color,
            strokewidth=focused ? 0.35 : 0)
    end

    # Plot every non-FCI level first. A gray roton crossing below the colored
    # reference states is then directly visible instead of being hidden by a
    # rank color or momentum-dependent marker.
    background_states = [state for state in keys(state_groups)
                         if !haskey(reference_colors, state)]
    sort!(background_states; by=state -> (state[1]..., state[2]))
    for state in background_states
        draw_state!(state_groups[state], OTHER_FLOW_COLOR)
    end
    for reference in references
        haskey(state_groups, reference) || continue
        draw_state!(state_groups[reference], reference_colors[reference]; focused=true)
    end
    return ax
end

function spectrum_legend!(position, sample=nothing)
    references = isnothing(sample) ? nothing : get(FCI_REFERENCE_MANIFOLD, sample, nothing)
    elements = [MarkerElement(marker=:circle, color=color, strokecolor=:black,
        strokewidth=0.35, markersize=10) for color in FOCUS_STATE_COLORS]
    labels = isnothing(references) ?
        ["FCI reference state $index" for index in eachindex(FOCUS_STATE_COLORS)] :
        ["FCI sector ($(state[1][1]), $(state[1][2])), level $(state[2])"
         for state in references]
    push!(elements, MarkerElement(marker=:circle, color=OTHER_FLOW_COLOR,
        strokewidth=0, markersize=9))
    push!(labels, "other levels")
    return Legend(position, elements, labels; framevisible=true,
        labelsize=11, patchsize=(20, 12))
end

function draw_metric_axis!(ax, sample, field; color=:royalblue3,
    compare_fci_gsd::Bool=true)
    ground = sweep_metric_series(sample, field; ensemble=:ground)
    if !isempty(ground)
        lines!(ax, first.(ground), last.(ground); color=:gray35, linewidth=1.6,
            linestyle=:dash, label="absolute ground state")
        scatter!(ax, first.(ground), last.(ground); color=:gray35, markersize=6)
    end
    averaged = compare_fci_gsd ?
        sweep_metric_series(sample, field; ensemble=:fci_gsd) : Tuple{Float64,Float64}[]
    if !isempty(averaged)
        lines!(ax, first.(averaged), last.(averaged); color=color, linewidth=2.2,
            label="FCI-manifold projector")
        scatter!(ax, first.(averaged), last.(averaged); color=color, markersize=7)
    end
    switches = sweep_ground_state_switches(sample)
    isempty(switches) || vlines!(ax, switches; color=(:black, 0.22),
        linestyle=:dot, linewidth=1.0)
    (!isempty(ground) || !isempty(averaged)) && axislegend(ax; position=:rt, labelsize=9)
    return ax
end


function draw_peak_wavevector_axis!(ax, sample, component::Symbol)
    field = component == :qx ? :peak_qx : component == :qy ? :peak_qy :
            error("Peak component must be qx or qy.")
    ground = sweep_metric_series(sample, field; ensemble=:ground)
    averaged = sweep_metric_series(sample, field; ensemble=:fci_gsd)
    if !isempty(ground)
        lines!(ax, first.(ground), last.(ground); color=:gray35, linestyle=:dash,
            linewidth=1.6, label="absolute ground state")
        scatter!(ax, first.(ground), last.(ground); color=:gray35, markersize=6)
    end
    if !isempty(averaged)
        lines!(ax, first.(averaged), last.(averaged); color=:royalblue3,
            linewidth=2.0, label="FCI-manifold projector")
        scatter!(ax, first.(averaged), last.(averaged); color=:royalblue3, markersize=7)
    end
    switches = sweep_ground_state_switches(sample)
    isempty(switches) || vlines!(ax, switches; color=(:black, 0.22),
        linestyle=:dot, linewidth=1.0)
    axislegend(ax; position=:rt, labelsize=9)
    return ax
end

function draw_zero_flux_gap_axis!(ax, points)
    for (field, label, color, style) in [
        (:fci_manifold_width, "E₂−E₀: FCI manifold width", :mediumpurple3, :solid),
        (:roton_gap, "E₃−E₂: roton gap", :firebrick3, :solid),
    ]
        y = [getproperty(point, field) for point in points]
        lines!(ax, [point.tpp for point in points], y; color=color,
            linestyle=style, linewidth=2, label=label)
        scatter!(ax, [point.tpp for point in points], y; color=color, markersize=6)
    end
    axislegend(ax; position=:lt, labelsize=10)
    return ax
end

function draw_competing_point_guides!(ax)
    sweep_min, sweep_max = extrema(actual_tpp.(SWEEP_NUMERATORS))
    guide_values = sort(filter(value -> sweep_min <= value <= sweep_max,
        vcat(Base.values(CHARACTERISTIC_TPP_VALUES)...)))
    vlines!(ax, guide_values; color=(:black, 0.35), linestyle=:dashdot, linewidth=1.2)
    return ax
end

function plot_sweep_results(; samples=STUDY_GEOMETRIES, max_ranks::Int=20)
    outdir = joinpath(FIGURE_ROOT, "sweep")
    mkpath(outdir)
    outputs = String[]

    for sample in samples
        isempty(existing_sweep_points(sample)) && continue
        tag = geometry_tag(sample)

        series = sweep_spectrum_series(sample; max_ranks=max_ranks)
        fig_s = Figure(size=(980, 580))
        ax_s = Axis(fig_s[1, 1]; xlabel="t′′", ylabel="E - E₀",
            title="Lowest all-sector spectrum — $tag", xtickformat=format_tpp_ticks)
        draw_spectrum_axis!(ax_s, sample; max_ranks=max_ranks, series=series)
        spectrum_legend!(fig_s[1, 2], sample)
        colgap!(fig_s.layout, 12)
        path_s = joinpath(outdir, "spectrum_ranks_$(tag).svg")
        save(path_s, fig_s)
        push!(outputs, path_s)

        fig_m = Figure(size=(760, 520))
        ax_m = Axis(fig_m[1, 1]; xlabel="t′′", ylabel="max |S(q)|",
            title="Ground-state sublattice structure-factor maxima — $tag",
            xtickformat=format_tpp_ticks)
        draw_ground_sublattice_max_axis!(ax_m, sample)
        path_m = joinpath(outdir,
            "max_abs_ground_sublattice_structure_factors_$(tag).svg")
        save(path_m, fig_m)
        push!(outputs, path_m)

        fig_mp = Figure(size=(760, 520))
        ax_mp = Axis(fig_mp[1, 1]; xlabel="t′′", ylabel="max |Sᶠᶜⁱ_proj(q)|",
            title="FCI-manifold projector maximum — $tag",
            xtickformat=format_tpp_ticks)
        if draw_fci_projector_zoom_axis!(ax_mp, sample;
            label="FCI-manifold projector")
            axislegend(ax_mp; position=:rt, labelsize=10)
            path_mp = joinpath(outdir,
                "max_abs_fci_manifold_projector_zoom_$(tag).svg")
            save(path_mp, fig_mp)
            push!(outputs, path_mp)
        end

        fig_r = Figure(size=(760, 520))
        ax_r = Axis(fig_r[1, 1]; xlabel="t′′",
            ylabel="max |S(q)/mean(S(q))|",
            title="Normalized structure-factor maximum — $tag",
            xtickformat=format_tpp_ticks)
        draw_metric_axis!(ax_r, sample, :max_abs_S_over_mean_S; color=:seagreen4)
        path_r = joinpath(outdir, "max_abs_structure_factor_over_mean_$(tag).svg")
        save(path_r, fig_r)
        push!(outputs, path_r)

        fig_q = Figure(size=(920, 720))
        ax_qx = Axis(fig_q[1, 1]; xlabel="t′′", ylabel="peak qₓ",
            title="Wavevector selected by max |S(q)| — $tag",
            xtickformat=format_tpp_ticks)
        ax_qy = Axis(fig_q[2, 1]; xlabel="t′′", ylabel="peak qᵧ",
            xtickformat=format_tpp_ticks)
        draw_peak_wavevector_axis!(ax_qx, sample, :qx)
        draw_peak_wavevector_axis!(ax_qy, sample, :qy)
        path_q = joinpath(outdir, "structure_factor_peak_wavevector_$(tag).svg")
        save(path_q, fig_q)
        push!(outputs, path_q)

        gap_points = sweep_gap_series(sample)
        if !isempty(gap_points)
            fig_g = Figure(size=(760, 520))
            ax_g = Axis(fig_g[1, 1]; xlabel="t′′", ylabel="energy difference",
                title="FCI manifold width and roton gap — $tag",
                xtickformat=format_tpp_ticks)
            draw_zero_flux_gap_axis!(ax_g, gap_points)
            draw_competing_point_guides!(ax_g)
            path_g = joinpath(outdir, "zero_flux_gap_diagnostics_$(tag).svg")
            save(path_g, fig_g)
            push!(outputs, path_g)
        end
    end

    active = [sample for sample in samples if !isempty(existing_sweep_points(sample))]
    if !isempty(active)
        ncol = 2
        nrow = cld(length(active), ncol)
        series_by_sample = Dict(sample => sweep_spectrum_series(sample;
            max_ranks=max_ranks) for sample in active)
        fig = Figure(size=(760 * ncol + 340, 520 * nrow))
        for (idx, sample) in enumerate(active)
            row, col = fldmod1(idx, ncol)
            ax = Axis(fig[row, col]; xlabel="t′′", ylabel="E - E₀",
                title=geometry_tag(sample), xtickformat=format_tpp_ticks)
            draw_spectrum_axis!(ax, sample; max_ranks=max_ranks,
                series=series_by_sample[sample])
        end
        Label(fig[0, 1:ncol], "Lowest all-sector spectrum versus t′′"; fontsize=20)
        spectrum_legend!(fig[1:nrow, ncol+1])
        colgap!(fig.layout, 12)
        path = joinpath(outdir, "spectrum_ranks_all_geometries.svg")
        save(path, fig)
        push!(outputs, path)

        for (field, label, filename, color) in [
            (:max_abs_S_over_mean_S, "max |S(q)/mean(S(q))|", "max_abs_structure_factor_over_mean_all_geometries.svg", :seagreen4),
        ]
            figm = Figure(size=(760 * ncol, 500 * nrow))
            for (idx, sample) in enumerate(active)
                row, col = fldmod1(idx, ncol)
                ax = Axis(figm[row, col]; xlabel="t′′", ylabel=label,
                    title=geometry_tag(sample), xtickformat=format_tpp_ticks)
                draw_metric_axis!(ax, sample, field; color=color)
            end
            Label(figm[0, 1:ncol], "$label versus physical t′′"; fontsize=20)
            path = joinpath(outdir, filename)
            save(path, figm)
            push!(outputs, path)
        end

        geometry_colors = [:royalblue3, :darkorange2, :seagreen4,
            :mediumpurple3, :firebrick3]
        for (ensemble, component_label, filename) in [
            (:ground_aa, "max |Sᴬᴬ(q)|",
                "max_abs_ground_SAA_all_geometries.svg"),
            (:ground_ab, "max |Re Sᴬᴮ(q)|",
                "max_abs_ground_SAB_all_geometries.svg"),
        ]
            fig_component = Figure(size=(800, 540))
            ax_component = Axis(fig_component[1, 1]; xlabel="t′′",
                ylabel=component_label,
                title="$component_label — absolute ground state",
                xtickformat=format_tpp_ticks)
            plotted = false
            for (index, sample) in enumerate(active)
                series = sweep_metric_series(sample, :max_abs_S; ensemble=ensemble)
                isempty(series) && continue
                color = geometry_colors[mod1(index, length(geometry_colors))]
                lines!(ax_component, first.(series), last.(series); color=color,
                    linewidth=2.0, label=geometry_tag(sample))
                scatter!(ax_component, first.(series), last.(series); color=color,
                    markersize=6)
                plotted = true
            end
            if plotted
                axislegend(ax_component; position=:rt)
                path = joinpath(outdir, filename)
                save(path, fig_component)
                push!(outputs, path)
            end
        end

        fig_projector = Figure(size=(800, 540))
        ax_projector = Axis(fig_projector[1, 1]; xlabel="t′′",
            ylabel="max |Sᶠᶜⁱ_proj(q)|",
            title="FCI-manifold projector maximum — all geometries",
            xtickformat=format_tpp_ticks)
        projector_plotted = false
        for (index, sample) in enumerate(active)
            color = geometry_colors[mod1(index, length(geometry_colors))]
            projector_plotted |= draw_fci_projector_zoom_axis!(
                ax_projector, sample; label=geometry_tag(sample), color=color)
        end
        if projector_plotted
            axislegend(ax_projector; position=:rt)
            path = joinpath(outdir,
                "max_abs_fci_manifold_projector_zoom_all_geometries.svg")
            save(path, fig_projector)
            push!(outputs, path)
        end


        gap_active = [sample for sample in active if !isempty(sweep_gap_series(sample))]
        if !isempty(gap_active)
            gap_ncol = 2
            gap_nrow = cld(length(gap_active), gap_ncol)
            figg = Figure(size=(760 * gap_ncol, 430 * gap_nrow))
            for (idx, sample) in enumerate(gap_active)
                row, col = fldmod1(idx, gap_ncol)
                points = sweep_gap_series(sample)
                ax = Axis(figg[row, col]; xlabel="t′′", ylabel="energy difference",
                    title=geometry_tag(sample),
                    xtickformat=format_tpp_ticks)
                draw_zero_flux_gap_axis!(ax, points)
                draw_competing_point_guides!(ax)
            end
            Label(figg[0, 1:gap_ncol], "Zero-flux FCI manifold width and roton gap";
                fontsize=20)
            path = joinpath(outdir, "zero_flux_gap_diagnostics_all_geometries.svg")
            save(path, figg)
            push!(outputs, path)
        end
    end
    @info "Wrote sweep figures" count=length(outputs) outdir
    return outputs
end

"Render the zero-twist ED spectrum stored at every available sweep point."
function plot_ed_spectrum_results(; samples=STUDY_GEOMETRIES)
    outputs = String[]
    for sample in samples
        outdir = joinpath(FIGURE_ROOT, "ed_spectra", geometry_tag(sample))
        for point in existing_sweep_points(sample)
            spectrum_path = joinpath(point, "spectrum.csv")
            isfile(spectrum_path) || continue
            rows = read_simple_csv(spectrum_path)
            isempty(rows) && continue

            xnumerator = csv_float(rows[1].tpp_numerator)
            xactual = csv_float(rows[1].tpp_actual)
            irrep_indices = [csv_int(row.irrep_index) - 1 for row in rows]
            energies = [csv_float(row.energy_minus_E0) for row in rows]
            n_irrep = prod(sample)

            fig = Figure(size=(760, 520))
            ax = Axis(fig[1, 1];
                xlabel="Irrep index",
                ylabel="E - E₀",
                title="ED Spectrum — $(geometry_tag(sample)), ν=1/6\nt′′ = $(format_tpp(xactual))",
                xticks=0:2:(n_irrep-1),
                xminorticksvisible=true,
                yminorticksvisible=true,
            )
            scatter!(ax, irrep_indices, energies; color=:royalblue1, markersize=14,
                alpha=0.75, strokecolor=:blue, strokewidth=0.5)

            mkpath(outdir)
            path = joinpath(outdir, "x_$(tpp_tag(xnumerator)).svg")
            save(path, fig)
            push!(outputs, path)
        end
    end
    @info "Wrote zero-twist ED spectrum figures" count=length(outputs)
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
            for (ensemble, allowed_name, dense_name, metrics_name, suffix, title) in [
                (:ground, "structure_allowed.csv", "structure_dense.csv",
                    "structure_metrics.csv", "", "absolute ground state"),
                (:fci_gsd, "structure_fci_gsd_allowed.csv", "structure_fci_gsd_dense.csv",
                    "structure_fci_gsd_metrics.csv", "_fci_gsd_average",
                    "FCI-manifold projector structure factor"),
            ]
                allowed_path = joinpath(point, allowed_name)
                dense_path = joinpath(point, dense_name)
                metrics_path = joinpath(point, metrics_name)
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
                ax1 = Axis(fig[1, 1]; xlabel="qₓ", ylabel="qᵧ",
                    title="finite-torus momenta", aspect=DataAspect())
                sc = scatter!(ax1, qx, qy; color=sq, colorrange=(lo, hi),
                    colormap=:viridis, markersize=30)
                draw_square_bz!(ax1)
                ax2 = Axis(fig[1, 2]; xlabel="qₓ", ylabel="qᵧ", title="dense q grid",
                    aspect=DataAspect())
                heatmap!(ax2, kx, ky, dense; colorrange=(lo, hi), colormap=:viridis)
                draw_square_bz!(ax2)
                Colorbar(fig[1, 3], sc; label="connected S(q)")
                Label(fig[0, 1:2],
                    "$(geometry_tag(sample)), t′′ = $(format_tpp(xactual)) — $title";
                    fontsize=18)
                mkpath(outdir)
                path = joinpath(outdir,
                    "x_$(tpp_tag(xnumerator))$(suffix)_allowed_and_dense.svg")
                save(path, fig)
                push!(outputs, path)
            end
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

function manifold_state_specs(datadir)
    path = joinpath(datadir, "summary.csv")
    isfile(path) || return NamedTuple[]
    rows = read_simple_csv(path)
    isempty(rows) && return NamedTuple[]
    specs = NamedTuple[]
    for state in split(rows[1].manifold_state_levels, ';')
        fields = split(state, ':')
        length(fields) == 3 || continue
        push!(specs, (sector=(parse(Int, fields[1]), parse(Int, fields[2])),
            level=parse(Int, fields[3])))
    end
    return specs
end

function diagnostic_tpp(datadir)
    for filename in ("summary.csv", "zero_flux_spectrum.csv")
        path = joinpath(datadir, filename)
        isfile(path) || continue
        rows = read_simple_csv(path)
        isempty(rows) || return csv_float(rows[1].tpp_actual)
    end
    return nothing
end

function candidate_point_is_active(phase, tpp)
    tpp === nothing && return false
    return any(value -> isapprox(tpp, value; atol=1e-10, rtol=0.0),
        characteristic_tpp_values(phase))
end

charge_pump_point_is_active(phase, tpp) = candidate_point_is_active(phase, tpp)

function existing_diagnostic_points(phase, sample::Tuple{Int,Int})
    root = joinpath(RESULT_ROOT, "diagnostics", String(phase), geometry_tag(sample))
    isdir(root) || return String[]
    points = String[]
    candidate_point_is_active(phase, diagnostic_tpp(root)) &&
        push!(points, root) # legacy flat layout
    for entry in readdir(root)
        dir = joinpath(root, entry)
        isdir(dir) && candidate_point_is_active(phase, diagnostic_tpp(dir)) &&
            push!(points, dir)
    end
    sort!(points; by=dir -> something(diagnostic_tpp(dir), Inf))
    return points
end

"Torus counting of `(1,m)`-admissible configurations for a particle cut."
function torus_admissible_count(n_orbitals::Int, n_particles::Int; m::Int=3)
    n_particles == 0 && return 1
    n_orbitals >= m * n_particles || return 0
    return (n_orbitals * binomial(n_orbitals - (m - 1) * n_particles - 1,
        n_particles - 1)) ÷ n_particles
end

function plot_diagnostic_results(; phases=[:AHC, :FCI, :CDW], samples=STUDY_GEOMETRIES,
    max_flow_curves::Int=20)
    max_flow_curves >= 1 || error("max_flow_curves must be positive; got $max_flow_curves.")
    outputs = String[]
    for phase_symbol in phases, sample in samples
        phase = String(phase_symbol)
        datadirs = existing_diagnostic_points(phase, sample)
        isempty(datadirs) && continue
        data_root = joinpath(RESULT_ROOT, "diagnostics", phase, geometry_tag(sample))
        figure_root = joinpath(FIGURE_ROOT, "diagnostics", phase, geometry_tag(sample))
        for datadir in datadirs
        outdir = datadir == data_root ? figure_root : joinpath(figure_root, basename(datadir))
        mkpath(outdir)
        tpp = diagnostic_tpp(datadir)
        title_prefix = isnothing(tpp) ? "$phase — $(geometry_tag(sample))" :
                       "$phase — $(geometry_tag(sample)), t′′ = $(format_tpp(tpp))"
        state_specs = manifold_state_specs(datadir)
        flow_isolation_min = nothing

        ground_allowed_path = joinpath(datadir, "structure_ground_allowed.csv")
        ground_dense_path = joinpath(datadir, "structure_ground_dense.csv")
        manifold_allowed_path = joinpath(datadir, "structure_manifold_allowed.csv")
        manifold_dense_path = joinpath(datadir, "structure_manifold_dense.csv")
        if all(isfile, (ground_allowed_path, ground_dense_path,
                       manifold_allowed_path, manifold_dense_path))
            ground_allowed = read_simple_csv(ground_allowed_path)
            manifold_allowed = read_simple_csv(manifold_allowed_path)
            qx = [csv_float(row.qx) for row in ground_allowed]
            qy = [csv_float(row.qy) for row in ground_allowed]
            ground_sq = [csv_float(row.S_q) for row in ground_allowed]
            manifold_sq = [csv_float(row.S_q) for row in manifold_allowed]
            ground_kx, ground_ky, ground_dense = dense_sf_arrays(ground_dense_path)
            manifold_kx, manifold_ky, manifold_dense = dense_sf_arrays(manifold_dense_path)
            lo = min(minimum(ground_sq), minimum(manifold_sq),
                minimum(ground_dense), minimum(manifold_dense))
            hi = max(maximum(ground_sq), maximum(manifold_sq),
                maximum(ground_dense), maximum(manifold_dense))

            fig_sf = Figure(size=(1120, 900))
            color_source = nothing
            for (row_index, row_title, allowed_values, dense_x, dense_y, dense_values) in [
                (1, "absolute ground state", ground_sq,
                    ground_kx, ground_ky, ground_dense),
                (2, "selected-manifold projector", manifold_sq,
                    manifold_kx, manifold_ky, manifold_dense),
            ]
                ax_allowed = Axis(fig_sf[row_index, 1]; xlabel="qₓ", ylabel="qᵧ",
                    title="$row_title — finite-torus momenta", aspect=DataAspect())
                scatter!(ax_allowed, qx, qy; color=allowed_values,
                    colorrange=(lo, hi), colormap=:viridis, markersize=30)
                draw_square_bz!(ax_allowed)
                ax_dense = Axis(fig_sf[row_index, 2]; xlabel="qₓ", ylabel="qᵧ",
                    title="$row_title — dense q grid", aspect=DataAspect())
                color_source = heatmap!(ax_dense, dense_x, dense_y, dense_values;
                    colorrange=(lo, hi), colormap=:viridis)
                draw_square_bz!(ax_dense)
            end
            Colorbar(fig_sf[1:2, 3], color_source; label="connected S(q)")
            Label(fig_sf[0, 1:2], "$title_prefix structure-factor comparison";
                fontsize=18)
            sf_path = joinpath(outdir, "structure_ground_vs_manifold.svg")
            save(sf_path, fig_sf)
            push!(outputs, sf_path)
        end

        flow_path = joinpath(datadir, "spectrum_flow.csv")
        if isfile(flow_path)
            rows = [row for row in read_simple_csv(flow_path)
                    if -1e-12 <= csv_float(row.flux_over_2pi) <= 1 + 1e-12]
            groups = Dict{Tuple{Int,Int,Int},Vector{Any}}()
            for row in rows
                key = (csv_int(row.k1), csv_int(row.k2), csv_int(row.level))
                push!(get!(groups, key, Any[]), row)
            end
            curves = collect(groups)
            sort!(curves; by=curve -> begin
                key, values = curve
                reference = values[argmin([abs(csv_float(row.flux_over_2pi)) for row in values])]
                (csv_float(reference.energy_minus_flux_ground), key...)
            end)
            resize!(curves, min(max_flow_curves, length(curves)))

            focused_states = unique((state.sector..., state.level) for state in state_specs)
            if isempty(focused_states)
                focused_states = unique(curve[1]
                    for curve in curves[1:min(3, length(curves))])
            end
            length(focused_states) <= length(FOCUS_STATE_COLORS) || error(
                "Focused manifold has $(length(focused_states)) states, but only " *
                "$(length(FOCUS_STATE_COLORS)) highlight colors are configured.")
            focused_colors = Dict(state => FOCUS_STATE_COLORS[index]
                for (index, state) in enumerate(focused_states))
            is_focused_curve(key) = haskey(focused_colors, key)

            # Draw background states first so all focused-manifold curves remain visible.
            ordered_curves = vcat(
                [curve for curve in curves if !is_focused_curve(curve[1])],
                [curve for curve in curves if is_focused_curve(curve[1])],
            )
            focused_energies = Float64[]
            all_flow_energies = Float64[]
            focused_by_flux = Dict{Float64,Vector{Float64}}()
            other_by_flux = Dict{Float64,Vector{Float64}}()
            for (key, values) in curves
                energy_by_flux = is_focused_curve(key) ? focused_by_flux : other_by_flux
                for row in values
                    flux = csv_float(row.flux_over_2pi)
                    energy = csv_float(row.energy_minus_flux_ground)
                    push!(all_flow_energies, energy)
                    push!(get!(energy_by_flux, flux, Float64[]), energy)
                    is_focused_curve(key) && push!(focused_energies, energy)
                end
            end

            fig = Figure(size=(980, 540))
            ax = Axis(fig[1, 1]; xlabel="inserted flux / 2π", ylabel="E - E₀(θ)",
                title="$title_prefix spectrum flow ($(length(curves)) lowest zero-flux states)")
            for curve in ordered_curves
                (k1, k2, level), values = curve
                sort!(values; by=row -> csv_float(row.flux_over_2pi))
                focused = is_focused_curve((k1, k2, level))
                line_color = focused ? :black : OTHER_FLOW_COLOR
                point_color = focused ? focused_colors[(k1, k2, level)] : OTHER_FLOW_COLOR
                alpha = focused ? 1.0 : 0.82
                lines!(ax, [csv_float(row.flux_over_2pi) for row in values],
                    [csv_float(row.energy_minus_flux_ground) for row in values];
                    color=line_color, linewidth=focused ? 1.6 : 1.0, alpha=alpha)
                scatter!(ax, [csv_float(row.flux_over_2pi) for row in values],
                    [csv_float(row.energy_minus_flux_ground) for row in values];
                    color=point_color, markersize=focused ? 14 : 10, alpha=alpha,
                    strokewidth=0)
            end
            xlims!(ax, 0.0, 1.0)
            if !isempty(focused_energies)
                focused_min, focused_max = extrema(focused_energies)
                focused_span = focused_max - focused_min
                direct_gaps = Float64[]
                for (flux, focus_values) in focused_by_flux
                    haskey(other_by_flux, flux) || continue
                    push!(direct_gaps,
                        minimum(other_by_flux[flux]) - maximum(focus_values))
                end
                global_flow_gap = isempty(direct_gaps) ? nothing : minimum(direct_gaps)
                fully_separated = !isnothing(global_flow_gap) && global_flow_gap > 1e-10
                plot_range = fully_separated ?
                    1.5 * (focused_span + global_flow_gap) : 2.0 * focused_span

                # A perfectly flat focused branch has zero span.  Retain a small,
                # finite window so Makie receives valid limits in this singular case.
                if plot_range <= 1e-10
                    energy_scale = isempty(all_flow_energies) ? 0.0 : maximum(all_flow_energies)
                    plot_range = max(0.1 * energy_scale, 0.05)
                end
                ylims!(ax, -0.05 * plot_range, plot_range)
            end

            legend_elements = [MarkerElement(marker=:circle, color=focused_colors[state],
                strokewidth=0, markersize=10) for state in focused_states]
            legend_labels = ["focused sector ($(state[1]), $(state[2])), level $(state[3])"
                             for state in focused_states]
            push!(legend_elements, MarkerElement(marker=:circle, color=OTHER_FLOW_COLOR,
                strokewidth=0, markersize=9))
            push!(legend_labels, "other levels")
            Legend(fig[1, 2], legend_elements, legend_labels;
                framevisible=true, labelsize=11, patchsize=(20, 12))
            colgap!(fig.layout, 10)
            path = joinpath(outdir, "spectrum_flow.svg")
            save(path, fig)
            push!(outputs, path)

            # A pump is only an invariant of a manifold that remains isolated
            # from its complement throughout the twist path.  Make that
            # precondition explicit using the instantaneous global ranks.
            manifold_size = isempty(state_specs) ? phase_spec(phase).manifold_size :
                            length(state_specs)
            gap_points = NamedTuple[]
            for flux_rows in Base.values(group_rows(rows, (:flux_over_2pi,)))
                energies = sort([csv_float(row.energy) for row in flux_rows])
                length(energies) >= manifold_size + 1 || continue
                push!(gap_points, (
                    flux=csv_float(flux_rows[1].flux_over_2pi),
                    neutral_from_ground=energies[manifold_size + 1] - energies[1],
                    rank_isolation=energies[manifold_size + 1] - energies[manifold_size],
                    manifold_width=energies[manifold_size] - energies[1],
                ))
            end
            sort!(gap_points; by=point -> point.flux)
            if !isempty(gap_points)
                flow_isolation_min = minimum(point.rank_isolation for point in gap_points)
                fig_gap = Figure(size=(760, 520))
                ax_gap = Axis(fig_gap[1, 1]; xlabel="inserted flux / 2π",
                    ylabel="energy difference",
                    title="$title_prefix manifold isolation along flux")
                for (field, label, color, style) in [
                    (:neutral_from_ground, "E$(manifold_size + 1)−E₁",
                        :royalblue3, :solid),
                    (:rank_isolation,
                        "E$(manifold_size + 1)−E$(manifold_size) (pump isolation)",
                        :darkorange2, :solid),
                    (:manifold_width, "E$(manifold_size)−E₁", :slategray3, :dash),
                ]
                    y = [getproperty(point, field) for point in gap_points]
                    lines!(ax_gap, [point.flux for point in gap_points], y;
                        color=color, linestyle=style, linewidth=2, label=label)
                    scatter!(ax_gap, [point.flux for point in gap_points], y;
                        color=color, markersize=5)
                end
                hlines!(ax_gap, [0.0]; color=:black, linestyle=:dot, linewidth=1.0)
                axislegend(ax_gap; position=:rt)
                gap_path = joinpath(outdir, "manifold_gap_flow.svg")
                save(gap_path, fig_gap)
                push!(outputs, gap_path)
            end
        end

        pump_path = joinpath(datadir, "charge_pump.csv")
        if isfile(pump_path) && charge_pump_point_is_active(phase_symbol, tpp)
            rows = read_simple_csv(pump_path)
            groups = group_rows(rows, (:branch,))
            fig = Figure(size=(760, 520))
            isolation_note = if flow_isolation_min === nothing
                ""
            elseif flow_isolation_min <= 1e-8
                "\nwarning: assumed manifold touches outside states along flow"
            else
                "\nmin flow isolation = " * @sprintf("%.3g", flow_isolation_min)
            end
            ax = Axis(fig[1, 1]; xlabel="inserted flux / 2π", ylabel="pumped charge ΔQ",
                title="$title_prefix charge pump$isolation_note")
            for (key, values) in sort(collect(groups); by=x -> csv_int(x[1][1]))
                branch = csv_int(key[1])
                sort!(values; by=row -> csv_float(row.flux_over_2pi))
                x = [csv_float(row.flux_over_2pi) for row in values]
                y = [csv_float(row.pumped_charge) for row in values]
                label = if branch <= length(state_specs)
                    state = state_specs[branch]
                    "branch $branch — sector ($(state.sector[1]), $(state.sector[2])), level $(state.level)"
                else
                    "branch $branch"
                end
                lines!(ax, x, y; color=Makie.Cycled(branch), linewidth=2, label=label)
                scatter!(ax, x, y; color=Makie.Cycled(branch), markersize=6)
            end
            axislegend(ax; position=:lt)
            path = joinpath(outdir, "charge_pump.svg")
            save(path, fig)
            push!(outputs, path)
        end

        pes_path = joinpath(datadir, "particle_entanglement_spectrum.csv")
        if phase_symbol == :FCI && isfile(pes_path)
            rows = read_simple_csv(pes_path)
            momenta = sort(unique((csv_int(row.k1), csv_int(row.k2)) for row in rows))
            index = Dict(k => i for (i, k) in enumerate(momenta))
            x = [index[(csv_int(row.k1), csv_int(row.k2))] for row in rows]
            y = [csv_float(row.entanglement_energy) for row in rows]
            # Keep momentum labels legible on geometries with many sectors.
            fig = Figure(size=(760, 520))
            ax = Axis(fig[1, 1]; xlabel="subsystem momentum sector", ylabel="ξ = -log(λ)",
                title="$title_prefix momentum-resolved PES",
                xticks=(1:length(momenta), ["($(k[1]),$(k[2]))" for k in momenta]),
                xticklabelrotation=π / 3)
            scatter!(ax, x, y; color=y, colormap=:viridis, markersize=7)
            sorted_y = sort(y)
            if length(sorted_y) >= 2
                gaps = diff(sorted_y)
                gap_index = argmax(gaps)
                gap_cut = (sorted_y[gap_index] + sorted_y[gap_index+1]) / 2
                gap_label = "largest gap: $gap_index levels below"
                summary_path = joinpath(datadir, "summary.csv")
                if phase_symbol == :FCI && isfile(summary_path)
                    summary = read_simple_csv(summary_path)[1]
                    n_a = csv_int(summary.PES_NA)
                    expected = torus_admissible_count(prod(sample), n_a; m=3)
                    gap_label *= "; (1,3) expected $expected"
                end
                hlines!(ax, [gap_cut]; color=:firebrick3, linestyle=:dash,
                    linewidth=1.8, label=gap_label)
                axislegend(ax; position=:rt)
            end
            path = joinpath(outdir, "particle_entanglement_spectrum.svg")
            save(path, fig)
            push!(outputs, path)
        end
        end
    end
    @info "Wrote diagnostic figures" count=length(outputs)
    return outputs
end

function charge_gap_rows(phase)
    root = joinpath(RESULT_ROOT, "charge_gap", String(phase))
    isdir(root) || return NamedTuple[]
    rows = NamedTuple[]
    for (dir, _, files) in walkdir(root)
        "charge_gap.csv" in files || continue
        for row in read_simple_csv(joinpath(dir, "charge_gap.csv"))
            candidate_point_is_active(phase, csv_float(row.tpp_actual)) && push!(rows, row)
        end
    end
    sort!(rows; by=row -> (csv_float(row.tpp_actual), csv_int(row.n_sites)))
    return rows
end

function plot_charge_gap_results(; phases=[:AHC, :FCI, :CDW])
    outputs = String[]
    outdir = joinpath(FIGURE_ROOT, "charge_gap")
    mkpath(outdir)
    phase_data = NamedTuple[]
    colors = [:royalblue3, :darkorange2, :seagreen4, :mediumpurple3,
        :firebrick3, :goldenrod2, :deeppink3, :turquoise3,
        :slateblue3, :olivedrab3, :orchid3, :cornflowerblue]
    color_index = 0

    for phase in phases
        rows = charge_gap_rows(phase)
        isempty(rows) && continue
        by_tpp = group_rows(rows, (:tpp_actual,))
        for point_rows in values(by_tpp)
            length(point_rows) >= 2 || continue
            sort!(point_rows; by=row -> csv_int(row.n_sites))
            x = [1 / csv_int(row.n_sites) for row in point_rows]
            y = [csv_float(row.charge_gap) for row in point_rows]
            design = hcat(ones(length(x)), x)
            fit = design \ y
            intercept, slope = fit
            residual = sqrt(mean((design * fit .- y) .^ 2))
            tpp_actual = csv_float(point_rows[1].tpp_actual)
            color_index += 1
            color = colors[mod1(color_index, length(colors))]

            push!(phase_data, (phase=String(phase), rows=point_rows, x=x, y=y,
                tpp_actual=tpp_actual, intercept=intercept,
                slope=slope, color=color))

            fit_path = joinpath(RESULT_ROOT, "charge_gap", String(phase),
                "finite_size_fit_tpp_$(tpp_tag(tpp_actual)).csv")
            ensure_parent(fit_path)
            open(fit_path, "w") do io
                println(io, "phase,tpp_numerator,tpp_actual,n_sizes,intercept_at_inverse_size_zero,slope,rms_residual")
                @printf(io, "%s,%.16g,%.16g,%d,%.16g,%.16g,%.16g\n", String(phase),
                    csv_float(point_rows[1].tpp_numerator), tpp_actual,
                    length(point_rows), intercept, slope, residual)
            end
        end
    end

    isempty(phase_data) && return outputs
    fig = Figure(size=(760, 520))
    ax = Axis(fig[1, 1]; xlabel="1 / N_sites", ylabel="charge gap Δc",
        title="Charge-gap finite-size scaling with 1/N_sites → 0 extrapolation")
    hlines!(ax, [0.0]; color=:black, linestyle=:dot, linewidth=1.2)
    vlines!(ax, [0.0]; color=:black, linestyle=:dot, linewidth=1.2)

    ymaximum = 0.0
    for data in phase_data
        line_x = collect(range(0.0, maximum(data.x); length=200))
        fit_y = data.intercept .+ data.slope .* line_x
        lines!(ax, line_x, fit_y; color=data.color, linewidth=2, linestyle=:dash)
        order = sortperm(data.x)
        label = "$(data.phase) — t′′ = $(format_tpp(data.tpp_actual))"
        lines!(ax, data.x[order], data.y[order]; color=data.color, linewidth=1.8,
            label=label)
        scatter!(ax, data.x, data.y; color=data.color, markersize=10)
        scatter!(ax, [0.0], [data.intercept]; color=data.color, marker=:diamond,
            markersize=13, strokecolor=:black, strokewidth=0.5)
        ymaximum = max(ymaximum, maximum(data.y), maximum(fit_y), data.intercept)
    end
    ylims!(ax, -0.1, 1.1 * ymaximum)
    axislegend(ax; position=:lt)

    path = joinpath(outdir, "tpp_charge_gap_finite_size_scaling.svg")
    save(path, fig)
    push!(outputs, path)
    @info "Wrote combined charge-gap scaling figure" phases=length(phase_data) path
    return outputs
end

function plot_all_results()
    return vcat(
        plot_sweep_results(),
        plot_ed_spectrum_results(),
        plot_structure_factor_results(),
        plot_diagnostic_results(),
        plot_charge_gap_results(),
    )
end
