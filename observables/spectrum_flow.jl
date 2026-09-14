# ═══════════════════════════════════════════════════════════════════════════
# flux_spectrum_flow — scan E(θ)  (thin wrapper over shared ED checkpoints)
# ═══════════════════════════════════════════════════════════════════════════

"""
Core Method to Compute the Spectrum Flow Under the Twisted Boundary Condition
---
Scans `twisted_phases_over_2π_list` along `flux_direction`, reading the full
symmetry-resolved ED data from canonical per-θ checkpoints produced by
[`ed_scan!`](@ref) in flux-scan mode.  If a checkpoint is missing it is
computed on the fly.

- Args:
    - `model::Real_Space_Second_Quantized_Model`: the second quantized model
    - `sector_labels`: `:all` (default), `:identity` or a `Vector` of sector tuples, e.g. `[(0,0), (1,0)]`
- Named Args:
    - `filling_fraction::Rational{Int}`: particles per flattened vertex
    - `flux_direction::Int=1`
    - `twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`
    - `checkpoint_dir::String="checkpoints"`: directory for per-θ checkpoints
    - `nev::Int=3`: eigenvalues per sector
    - `fig_path::Union{Nothing,String}=nothing`
    - `overwrite::Bool=false`: recompute even if checkpoint exists
- Returns:
    - `res::NamedTuple`: `twisted_phases_over_2π_list`, `energies`, `sector_labels`, `fig_path`, ...
"""
function flux_spectrum_flow(
    model::Real_Space_Second_Quantized_Model,
    sector_labels=:all;
    filling_fraction::Rational{Int},
    flux_direction::Int=1,
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    checkpoint_dir::String="checkpoints",
    nev::Int=3,
    fig_path::Union{Nothing,String}=nothing,
    overwrite::Bool=false,
)::NamedTuple
    # ── Resolve sector labels ──
    is_identity = sector_labels == :identity
    labels = sector_labels == :all ? [(i,j) for i in 0:model.lattice.sample_size[1]-1 for j in 0:model.lattice.sample_size[2]-1] :
             is_identity ? [:identity] :
             sector_labels isa Tuple ? [Tuple(Int.(sector_labels))] :
             [Tuple(Int.(l)) for l in sector_labels]

    # Every momentum sector is scanned; labels only restrict displayed branches.
    # ── Ensure all per-θ ED checkpoints exist ──
    flux0 = zeros(Float64, model.lattice.dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_ed_data = build_ed_data(model; filling_fraction=filling_fraction,
        symmetry_group=(is_identity ? build_identity_group(model.lattice.n_site) : build_translation_group(model.lattice, flux0)))

    checkpoint_paths = ed_scan!(init_ed_data;
        nev=nev, mode=:matrix,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
        flux_direction=flux_direction,
        checkpoint_dir=checkpoint_dir,
        overwrite=overwrite,
        scanned_sectors=nothing,
    )

    # ── Extract eigenvalues for requested sectors from each checkpoint ──
    energies = fill(NaN, length(twisted_phases_over_2π_list), length(labels), nev)
    for (iθ, ckpt_path) in enumerate(checkpoint_paths)
        ed_data = load_checkpoint(ckpt_path)
        for (isector, label) in enumerate(labels)
            irrep_idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
            if irrep_idx !== nothing && haskey(ed_data.ed_scan_res, irrep_idx)
                vals = ed_data.ed_scan_res[irrep_idx][1]
                nv = min(length(vals), nev)
                energies[iθ, isector, 1:nv] .= vals[1:nv]
            end
        end
    end

    all_energies = [sort(vcat([r[1] for r in values(load_checkpoint(p).ed_scan_res)]...)) for p in checkpoint_paths]
    global_energies = reduce(hcat, [v[1:min(nev,length(v))] for v in all_energies])'

    # ── Plot ──
    if !isnothing(fig_path)
        mkpath(dirname(fig_path))
        fig = Figure(size=(800, 600))
        sectors_str = is_identity ? "Full Hilbert Space" : join(repr.(labels), ", ")
        ax = Axis(fig[1, 1];
            xlabel="Inserted Flux [2π] along Direction-$(flux_direction)",
            ylabel="E",
            title="$(model.lattice.sample_size)-sample Spectrum Flow — sectors: $sectors_str"
        )
        for isector in eachindex(labels), level in 1:nev
            lbl = level == 1 ? "sector $(repr(labels[isector]))" : nothing
            lines!(ax, twisted_phases_over_2π_list, energies[:, isector, level];
                color=Makie.Cycled(isector), alpha=(1 - (level - 1) / nev),
                linewidth=2, label=lbl)
            scatter!(ax, twisted_phases_over_2π_list, energies[:, isector, level];
                color=Makie.Cycled(isector), alpha=(1 - (level - 1) / nev),
                markersize=6)
        end
        is_identity || axislegend(ax; position=:lt)
        save(fig_path, fig)
        @info "spectrum flow plot saved to `$fig_path`"
    end

    res = (; twisted_phases_over_2π_list, energies, sector_labels=labels,
        flux_direction, nev, is_identity, fig_path, checkpoint_paths, global_energies)
    return res
end