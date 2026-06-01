# ═══════════════════════════════════════════════════════════════════════════
# flux_spectrum_flow — scan E(θ)
# ═══════════════════════════════════════════════════════════════════════════

"""
    flux_spectrum_flow(model, sector_labels; kwargs...)

Scan one twisted boundary angle and plot the low-energy many-body spectral flow (E vs. θ). This follows the ExactDiagonalization.jl FCI showcase.

The model is updated **in-place** at each θ via [`update_second_quantized_model_with_twisted_phases!`](@ref), which delegates to `TightBinding.generate_bilinear_terms`.  No model copies are created inside the scan loop.

# Arguments
- `model::Real_Space_Second_Quantized_Model`: the zero-flux model (mutated in-place)
- `sector_labels`: `:identity` for the full Hilbert space, or a tuple `(k₁,k₂)` or vector of tuples `[(k₁,k₂), ...]` for sector-resolved diagonalisation.

# Keyword arguments
- `filling_fraction::Rational{Int}`: particle filling
- `flux_direction::Int=1`: which periodic direction to thread the flux through
- `nθ::Int=8`: number of θ grid points
- `θ_max_per_2π::Real=1.0`: maximum flux in units of 2π (e.g. `2.0` = 4π)
- `nev::Int=3`: number of eigenvalues per sector
- `fig_path::Union{Nothing,String}=nothing`: where to save the plot (SVG/PNG)
- `checkpoint_path::Union{Nothing,String}=nothing`: where to save raw data (.jld2)

# Flux-aware symmetry
When `sector_labels` is a momentum tuple (not `:identity`), the gauge-covariant translation group is rebuilt at each θ via [`build_translation_group`](@ref) with the current flux, and the orbit catalog's stabiliser phases are updated in-place via [`update_orbit_stabilizer_phases!`](@ref) — this avoids the expensive O(C(N,νN)) Gosper iteration at every flux point.  The momentum labels `[k₁,k₂]` remain the standard integer crystal momenta.

With `:identity`, the full Hilbert space is diagonalised at each θ.
"""
function flux_spectrum_flow(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    flux_direction::Int=1,
    nθ::Int=8,
    θ_max_per_2π::Real=1.0,
    nev::Int=3,
    fig_path::Union{Nothing,String}=nothing,
    checkpoint_path::Union{Nothing,String}=nothing,
)
    # ── Resolve sector labels ──
    is_identity = sector_labels == :identity
    labels = is_identity ? [:identity] :
             sector_labels isa Tuple ? [Tuple(Int.(sector_labels))] :
             [Tuple(Int.(l)) for l in sector_labels]

    dim = model.lattice.dim
    θs = collect(range(0.0, Float64(θ_max_per_2π); length=nθ))
    energies = fill(NaN, length(θs), length(labels), nev)

    # ── Initialise model at θ=0, build ED data once ──
    flux0 = zeros(Float64, dim)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
    init_group = is_identity ? build_identity_group(model.lattice.n_site) :
                 build_translation_group(model.lattice, flux0)
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=init_group)

    for (iθ, θ_val) in enumerate(θs)
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
        fig = Figure(size=(720, 420))
        ax = Axis(fig[1, 1];
            xlabel="inserted flux θ$(flux_direction) / 2π",
            ylabel="E",
            title="Many-body spectrum flow ($(is_identity ? "full" : "sector-resolved"))")
        colors = [:royalblue3, :firebrick3, :seagreen4, :darkorange3, :purple4]
        for isector in eachindex(labels), level in 1:nev
            color = colors[mod1(isector, length(colors))]
            lbl = level == 1 ? "sector $(repr(labels[isector]))" : nothing
            lines!(ax, θs, energies[:, isector, level]; color, linewidth=2, label=lbl)
            scatter!(ax, θs, energies[:, isector, level]; color, markersize=6)
        end
        is_identity || axislegend(ax; position=:rt)
        save(fig_path, fig)
        @info "Flux flow plot → $fig_path"
    end

    # ── Checkpoint ──
    if !isnothing(checkpoint_path)
        mkpath(dirname(checkpoint_path))
        result = (; θs, energies, sector_labels=labels,
            flux_direction, nev, is_identity, fig_path, θ_max_per_2π)
        @save checkpoint_path result
    end

    return (; θs, energies, sector_labels=labels, flux_direction, nev, is_identity,
        fig_path, θ_max_per_2π)
end

# ═══════════════════════════════════════════════════════════════════════════
# Bosonic Haldane FCI model builder
# ═══════════════════════════════════════════════════════════════════════════

"""
    build_haldane_fci_model(; sample_size=[3,4], kwargs...)

Construct the hard-core boson Haldane honeycomb model that hosts a bosonic fractional Chern insulator (FCI) phase.  This is the model studied in [D.N. Sheng et al., Phys. Rev. Lett. **107**, 146803 (2011)].

The model includes:
- nearest-neighbour hopping `t` (real, inter-sublattice)
- next-nearest-neighbour hopping `t′ e^{±iϕ}` (complex, intra-sublattice,
  time-reversal breaking Haldane flux)
- third-nearest-neighbour hopping `t′′` (real, inter-sublattice)
- optional density-density interactions `V₁`, `V₂`

Default parameters are chosen so the lower band is nearly flat with Chern number ±1.  At ν = 1/2 filling of this band the ground state is an FCI with two nearly-degenerate states on the torus.

# Keyword arguments
- `sample_size::Vector{Int}=[3,4]`: unit cells [Lx, Ly]; total sites = 2·Lx·Ly
- `t::Real=1.0`, `t′::Real=0.60`, `t′′::Real=-0.58`: hopping amplitudes
- `ϕ_over_2π::Real=0.2`: Haldane flux per 2π
- `V1::Real=0.0`, `V2::Real=0.0`: interaction strengths
"""
function build_haldane_fci_model(;
    sample_size::Vector{Int}=[3, 4],
    t::Real=1.0,
    t′::Real=0.60,
    t′′::Real=-0.58,
    ϕ_over_2π::Real=0.2,
    V1::Real=0.0,
    V2::Real=0.0,
)
    params = Dict{String,Any}(
        "t" => t, "t′" => t′, "t′′" => t′′,
        "ϕ_over_2π" => ϕ_over_2π, "V1" => V1, "V2" => V2,
    )

    lattice = TightBinding.initialize_real_space_lattice(;
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [1 / 2, sqrt(3) / 2]],
        sub_crys_list=[[0.0, 0.0], [1 / 3, 1 / 3]],
        lattice_name="Haldane_Honeycomb",
        pbc_indicator=[true, true],
    )
    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="haldane_boson_FCI")

    ϕ = 2π * ϕ_over_2π

    # Nearest-neighbour (inter-sublattice, real)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, -1], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 0], 2)) => -t; is_hermitian=true)

    # NNN (intra-sublattice, complex — Haldane flux)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′ * exp(im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′ * exp(-im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 1], 1)) => -t′ * exp(im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′ * exp(-im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′ * exp(im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([-1, 1], 2)) => -t′ * exp(-im * ϕ); is_hermitian=true)

    # Third-nearest-neighbour (inter-sublattice, real)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 1), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t′′; is_hermitian=true)

    n_filled = prod(sample_size) ÷ 2
    params["n_filled"] = n_filled
    params["filling_fraction"] = n_filled // lattice.n_site
    @info "Haldane FCI: $(sample_size[1])×$(sample_size[2]) honeycomb, " *
          "$(lattice.n_site) sites, $n_filled bosons, params=$params"

    # Build the model with θ=0 hopping terms
    model = Real_Space_Second_Quantized_Model(
        params, lattice, tb_model, Bosonic(),
        TightBinding.generate_bilinear_terms(tb_model; twisted_phases_over_2π=zeros(Float64, lattice.dim)),
        Tuple{Int,Int,ComplexF64}[],
    )
    return model
end

"""
    default_fci_sectors(sample_size)

Return the crystal-momentum sector labels `[k₁, k₂]` that host the two nearly-degenerate FCI ground states (the Kramers pair on the torus).
"""
function default_fci_sectors(sample_size::Vector{Int})
    if sample_size == [2, 3]
        return [(0, 0), (1, 0)]
    elseif sample_size == [3, 2]
        return [(0, 0), (0, 1)]
    elseif sample_size == [3, 4]
        return [(0, 0), (0, 2)]
    elseif sample_size == [4, 3]
        return [(0, 0), (2, 0)]
    end
    return [(0, 0)]
end

# ═══════════════════════════════════════════════════════════════════════════
# Self-contained test: Haldane FCI spectral flow
# ═══════════════════════════════════════════════════════════════════════════

"""
    test_haldane_fci_flux_flow(; sample_size, nθ, mode, θ_max_per_2π)

Run the boundary-twisted spectrum-flow test on the bosonic Haldane FCI model and verify the expected physics: the two nearly-degenerate FCI ground states intertwine under flux insertion (Laughlin charge pump, ΔQ≈1/2 each).

Call with `mode=:sectors` (default) for momentum-sector-resolved diagonalisation, or `mode=:identity` for full Hilbert-space diagonalisation.
"""
function test_haldane_fci_flux_flow(;
    sample_size::Vector{Int}=[2, 3],
    nθ::Int=9,
    mode::Symbol=:sectors,
    θ_max_per_2π::Float64=1.0,
    flux_direction::Int=findfirst(d -> mod(sample_size[d], 2) == 0, 1:length(sample_size)),
)
    model = build_haldane_fci_model(; sample_size=sample_size)

    @assert mode in [:identity, :sectors]
    is_identity = mode == :identity
    labels = is_identity ? :identity : default_fci_sectors(sample_size)
    tag = is_identity ? "identity" : "sectors"
    dims_str = "$(sample_size[1])x$(sample_size[2])"

    PROJECT_ROOT = dirname(@__DIR__)

    @info "The chosen twisted phase direction: $flux_direction"
    result = flux_spectrum_flow(
        model,
        labels;
        filling_fraction=model.params["filling_fraction"],
        flux_direction=flux_direction,
        nθ=nθ,
        θ_max_per_2π=θ_max_per_2π,
        nev=3,
        fig_path=joinpath(PROJECT_ROOT, "figures",
            "haldane_fci_flux_flow_$(tag)_$(dims_str).svg"),
        checkpoint_path=joinpath(PROJECT_ROOT, "checkpoints",
            "haldane_fci_flux_flow_$(tag)_$(dims_str).jld2"),
    )

    # ── Physics checks ──
    @testset "Haldane FCI flux flow ($dims_str, $mode)" begin
        @test size(result.energies) == (nθ, length(result.sector_labels), 3)
        @test all(isfinite, result.energies[:, :, 1])

        # [2,3] sector-resolved: two GS intertwine over one flux quantum
        if sample_size == [2, 3] && mode == :sectors && isapprox(θ_max_per_2π, 2.0; atol=1e-12)
            mid = (nθ + 1) ÷ 2  # θ = 1.0 (one flux quantum)

            @test isapprox(result.energies[1, 1, 1], -7.16380536342344; atol=1e-8)
            @test isapprox(result.energies[1, 2, 1], -7.16337536156982; atol=1e-8)

            # Charge pump ΔQ=1/2: after one flux quantum, sectors swap
            @test isapprox(result.energies[1, 1, 1], result.energies[mid, 2, 1]; atol=1e-8)
            @test isapprox(result.energies[1, 2, 1], result.energies[mid, 1, 1]; atol=1e-8)

            # Total ΔQ=1: after two flux quanta, each GS returns to itself
            @test isapprox(result.energies[1, :, 1], result.energies[end, :, 1]; atol=1e-8)

            if nθ == 9
                expected = [
                    -7.16380536 -7.16337536
                    -7.18648366 -7.13511926
                    -7.19891147 -7.19891147
                    -7.18648366 -7.13511926
                    -7.16380536 -7.16337536
                    -7.18648366 -7.13511926
                    -7.19891147 -7.19891147
                    -7.18648366 -7.13511926
                    -7.16380536 -7.16337536
                ]
                for iθ in 1:nθ
                    @test isapprox(sort(result.energies[iθ, :, 1]), expected[iθ, :]; atol=1e-5)
                end
            end
        end
    end

    return result
end
