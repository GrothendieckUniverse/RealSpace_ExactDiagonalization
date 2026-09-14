#!/usr/bin/env julia

using RealSpace_ExactDiagonalization
using CairoMakie

"Small bosonic FCI example: all-sector flow, transverse pump, and manifold Chern number."
function run_fci_topology(;
    output_dir=joinpath(@__DIR__, "..", "figures", "fci_topology"),
    checkpoint_dir=joinpath(@__DIR__, "..", "checkpoints", "fci_topology"),
    flux_points::Int=17, chern_grid=(5, 5))
    flux_points >= 2 || throw(ArgumentError("Need at least two flux points."))
    model = build_zero_flux_bosonic_fci_second_quantized_model(; sample_size=[2, 3])
    filling = 1 // 4
    fluxes = collect(range(0.0, 1.0; length=flux_points))
    mkpath(output_dir)
    flow = flux_spectrum_flow(model, :all; filling_fraction=filling, nev=4,
        twisted_phases_over_2π_list=fluxes, checkpoint_dir=joinpath(checkpoint_dir, "flow"))
    reference = minimum(flow.energies[1, :, :])
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1]; xlabel="inserted x flux / 2π", ylabel="E − E₀(0)",
        title="Bosonic FCI 2×3 — all momentum sectors")
    for (index, label) in enumerate(flow.sector_labels), level in 1:flow.nev
        lines!(ax, fluxes, flow.energies[:, index, level] .- reference;
            color=Cycled(index), alpha=level == 1 ? 1.0 : 0.3,
            label=level == 1 ? string(label) : nothing)
    end
    axislegend(ax; nbanks=2, labelsize=10)
    save(joinpath(output_dir, "spectrum_flow.svg"), fig)
    pump = flux_charge_pump(model, :all; filling_fraction=filling, manifold_size=2,
        twisted_phases_over_2π_list=fluxes, checkpoint_dir=joinpath(checkpoint_dir, "flow"),
        fig_path=joinpath(output_dir, "charge_pump.svg"))
    chern = many_body_chern_number(model, :all; filling_fraction=filling, manifold_size=2,
        flux_grid_size=chern_grid, checkpoint_dir=joinpath(checkpoint_dir, "chern"),
        fig_path=joinpath(output_dir, "chern_number.svg"))
    println("Manifold Chern number: ", chern.chern_number)
    println("Pump branches: ", pump.pumped_charges, "; total: ", sum(pump.pumped_charges))
    println("Minimum pump gap: ", pump.min_gap,
        "; minimum projected-position singular value: ", minimum(pump.min_position_singular_values))
    return (; flow, pump, chern)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_fci_topology()
end
