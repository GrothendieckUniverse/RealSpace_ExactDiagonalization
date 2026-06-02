# ═══════════════════════════════════════════════════════════════════════════
# flux_spectrum_flow — scan E(θ)
# ═══════════════════════════════════════════════════════════════════════════

"""
Core Method to Compute the Spectrum Flow Under the Twisted Boundary Condition
---
by scan the twisted phases `θ_α` for a given `flux_direction-α` and for all sectors in `sector_labels`. The flux-induced twisted phases directly affect the model's `bilinear_terms` via Peierls substitution. 
    
In mode `:sector`, the symmetry group and orbit catalog stabiliser phases are then updated in place via `update_orbit_stabilizer_phases!`, avoiding re-enumeration at each flux value. In mode `:identity`, the group is the trivial identity group (no symmetry resolution).

- Args:
    - `model::Real_Space_Second_Quantized_Model`: the second quantized model to be used
    - `sector_labels`: either `:identity` (full Hilbert space) or a `Vector` of sector labels, e.g. `[(0,0), (1,0)]`
- Named Args:
    - `filling_fraction::Rational{Int}`: the number of particles per flatband
    - `flux_direction::Int=1`: the direction along which the flux is applied
    - `twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`: the list of twisted phased to scan
    - `nev::Int=3`: the number of low-lying level for each symmetrized sector
    - `fig_path::Union{Nothing,String}=nothing`: path to save a SVG/PNG figure of the spectrum flow
    - `checkpoint_path::Union{Nothing,String}=nothing`: path to JLD2 checkpoint file for recovery
- Reture:
    - `res::NamedTuple`: a named tuple containing `twisted_phases_list`, `energies`, `sector_labels`, `flux_direction`, `nev`, `is_identity`, and `fig_path`
"""
function flux_spectrum_flow(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    flux_direction::Int=1,
    twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    nev::Int=3,
    fig_path::Union{Nothing,String}=nothing,
    checkpoint_path::Union{Nothing,String}=nothing,
)::NamedTuple
    # ── Resolve sector labels ──
    is_identity = sector_labels == :identity
    labels = is_identity ? [:identity] :
             sector_labels isa Tuple ? [Tuple(Int.(sector_labels))] :
             [Tuple(Int.(l)) for l in sector_labels]

    dim = model.lattice.dim
    energies = fill(NaN, length(twisted_phases_list), length(labels), nev)

    # ── Initialise model at θ=0, build ED data once ──
    flux0 = zeros(Float64, dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_group = is_identity ? build_identity_group(model.lattice.n_site) :
                 build_translation_group(model.lattice, flux0)
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=init_group)

    for (iθ, θ_val) in enumerate(twisted_phases_list)
        flux = zeros(Float64, dim)
        flux[flux_direction] = θ_val

        # ── In-place update: bilinear terms + symmetry group + catalog ──
        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux)

        if is_identity
            ed_data.symmetry_group = build_identity_group(model.lattice.n_site)
        else
            active_group = build_translation_group(model.lattice, flux)
            ed_data.symmetry_group = active_group
            update_orbit_stabilizer_phases!(ed_data.orbit_catalog, active_group, model.statistics)
        end

        for (isector, label) in enumerate(labels)
            vals, _ = ed_scan_at_irrep_matrix!(label, ed_data; nev=nev)
            energies[iθ, isector, 1:length(vals)] .= vals
        end
    end

    # ── Plot ──
    if !isnothing(fig_path)
        mkpath(dirname(fig_path))
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1];
            xlabel="Inserted Flux [2π] along Direction-$(flux_direction)",
            ylabel="E",
            title="$(model.lattice.sample_size)-sample Spectrum Flow ($(is_identity ? "Full" : "Sector-Resolved"))"
        )
        for isector in eachindex(labels), level in 1:nev
            lbl = level == 1 ? "sector $(repr(labels[isector]))" : nothing
            lines!(ax, twisted_phases_list, energies[:, isector, level]; color=Makie.Cycled(isector), alpha=(1 - (level - 1) / nev), linewidth=2, label=lbl)
            scatter!(ax, twisted_phases_list, energies[:, isector, level]; color=Makie.Cycled(isector), alpha=(1 - (level - 1) / nev), markersize=6)
        end
        is_identity || axislegend(ax; position=:rt)
        save(fig_path, fig)
        @info "spectrum flow plot saved to `$fig_path`"
    end

    # ── Checkpoint ──
    if !isnothing(checkpoint_path)
        mkpath(dirname(checkpoint_path))
        result = (; twisted_phases_list, energies, sector_labels=labels,
            flux_direction, nev, is_identity, fig_path)
        @save checkpoint_path result
    end

    res = (; twisted_phases_list, energies, sector_labels=labels, flux_direction, nev, is_identity, fig_path) # named tuple

    return res
end