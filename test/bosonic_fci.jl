# ═══════════════════════════════════════════════════════════════════════════
# Bosonic Haldane FCI model builder
# ═══════════════════════════════════════════════════════════════════════════
"parameters from Wang, Gu, Gong, and Sheng's bosonic FCI [PhysRevLett.107.146803]"
const params_DNSheng = Dict(
    "t" => 1.0, # nearest-neighbor hopping
    "t′" => 0.60, # next-nearest-neighbor hopping
    "t′′" => -0.58, # next-next-nearest-neighbor hopping
    "ϕ_over_2π" => 0.2, # flux 0.4π
    "V1" => 0.0, # nearest neighbor density-density interaction strength
    "V2" => 0.0, # next-nearest neighbor density-density interaction strength
)


"""
Construct the Second Quantized Model for the Bosonic FCI Phase on Haldane Honeycomb Lattice
---
- Named Args:
    - `sample_size::Vector{Int}=[3, 4]`: the sample size along each primitive lattice vector
    - `params::Dict{String,<:Number}=params_DNSheng`: model parameters using D. N. Sheng's parameters
"""
function build_zero_flux_bosonic_fci_second_quantized_model(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_DNSheng,
)::Real_Space_Second_Quantized_Model
    t = params["t"]
    t′ = params["t′"]
    t′′ = params["t′′"]
    ϕ_over_2π = params["ϕ_over_2π"]
    V1 = params["V1"]
    V2 = params["V2"]

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

    # Build the model with initialized ZERO twisted phases
    bilinear_terms = TightBinding.generate_bilinear_terms(tb_model; twisted_phases_over_2π=zeros(Float64, lattice.dim))
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()

    model = Real_Space_Second_Quantized_Model(
        params,
        lattice,
        tb_model,
        Bosonic(),
        bilinear_terms,
        density_terms,
    )
    return model
end

"Return the crystal-momentum sector labels `[k₁, k₂]` that host the two nearly-degenerate bosonic semion FCI ground states (verified via full scan of ED)"
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

    @warn "Unknown `sample_size=$sample_size` — defaulting set sector label to [(0,0)]. The user should manually specify the FCI ground-state sectors."
    return [(0, 0)]
end

# ═══════════════════════════════════════════════════════════════════════════
# Self-contained test: Haldane FCI spectral flow
# ═══════════════════════════════════════════════════════════════════════════
"""
Test Spectrum Flow for the Bosonic FCI on Haldane Honeycomb lattice
---
- Named Args:
    - `sample_size::Vector{Int}=[2, 3]`
    - `params::Dict{String,<:Number}=params_DNSheng`
    - `filling_fraction::Rational{Int}=1 // 2`: the filling **per flatband**. Note: the `filling_fraction` input in the ED constructor is the fractional per flattened vertices!
    - `mode::Symbol=:sectors`: Here model `:sectors` for symmetry-resolved sectors, `:identity` for full Hilbert space (any sector will be the same)
    - `flux_direction::Int=findfirst(d -> mod(sample_size[d], 2) == 0, 1:length(sample_size))`: Always set to along the direction that can be divided by the GSD of the topological ordered states (2 for semion TO here)
    - `twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`: the list of twisted phased used to scan the twisted spectrum flow.
"""
function test_bosonic_fci_spectrum_flow(;
    sample_size::Vector{Int}=[2, 3],
    params::Dict{String,<:Number}=params_DNSheng,
    filling_fraction::Rational{Int}=1 // 2, # filling per flatband
    mode::Symbol=:sectors,
    flux_direction::Int=findfirst(d -> mod(sample_size[d], 2) == 0, 1:length(sample_size)), # always along the direction that can be divided by the GSD (2 for semion TO here)
    twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
)
    @info "Twisted phase list to be computed: $twisted_phases_list"
    model = build_zero_flux_bosonic_fci_second_quantized_model(; sample_size=sample_size, params=params)

    @info "The chosen twisted phase direction: $flux_direction"

    lattice = model.lattice
    filling_fraction = filling_fraction / 2 # filling fraction per flattened vertices
    n_filled = Int(lattice.n_site * filling_fraction)
    @info "Bosonic FCI over $(sample_size[1])×$(sample_size[2]) Haldane honeycomb lattice with $(lattice.n_site) sites, $n_filled bosons (vertices filling $filling_fraction)"

    @assert mode in [:identity, :sectors]
    is_identity = mode == :identity
    labels = is_identity ? :identity : default_fci_sectors(sample_size)
    tag = is_identity ? "identity" : "sectors"
    PROJECT_ROOT = dirname(@__DIR__)

    result = flux_spectrum_flow(
        model,
        labels;
        filling_fraction=filling_fraction,
        flux_direction=flux_direction,
        twisted_phases_list=twisted_phases_list,
        nev=3,
        fig_path=joinpath(PROJECT_ROOT, "figures",
            "bosonic_FCI_spectrum_flow_$(tag)_$(sample_size).svg"),
        checkpoint_path=joinpath(PROJECT_ROOT, "checkpoints",
            "bosonic_FCI_spectrum_flow_$(tag)_$(sample_size).jld2"),
    )

    return result
end

# ═══════════════════════════════════════════════════════════════════════════
# Self-contained test: Haldane FCI charge pump
# ═══════════════════════════════════════════════════════════════════════════
"""
Test the Fractional Charge Pump for the Bosonic FCI on Haldane Honeycomb Lattice
---
For the half-filled Chern band, one flux quantum should advance each polarization
branch by approximately one half charge — the finite-size version of the bosonic
ν = 1/2 Laughlin pump.

- Named Args:
    - `sample_size::Vector{Int}=[2, 3]`
    - `params::Dict{String,<:Number}=params_DNSheng`
    - `filling_fraction::Rational{Int}=1 // 2`: filling **per flatband**
    - `mode::Symbol=:sectors`: `:sectors` for symmetry-resolved, `:identity` for full Hilbert space
    - `flux_direction::Int=findfirst(d -> mod(sample_size[d], 2) == 0, 1:length(sample_size))`: flux direction (along even-size direction for semion GSD=2)
    - `twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`: flux values to scan
    - `atol::Float64=0.08`: tolerance for charge quantisation check
"""
function test_bosonic_fci_charge_pump(;
    sample_size::Vector{Int}=[2, 3],
    params::Dict{String,<:Number}=params_DNSheng,
    filling_fraction::Rational{Int}=1 // 2, # filling per flatband
    mode::Symbol=:sectors,
    flux_direction::Int=findfirst(d -> mod(sample_size[d], 2) == 0, 1:length(sample_size)),
    polarization_direction::Int=_default_polarization_direction(length(sample_size), flux_direction),
    twisted_phases_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    atol::Float64=0.08,
)
    @info "Twisted phase list to be computed: $twisted_phases_list"
    model = build_zero_flux_bosonic_fci_second_quantized_model(; sample_size=sample_size, params=params)

    lattice = model.lattice
    filling_fraction_vertex = filling_fraction // 2 # filling fraction per flattened vertices
    n_filled = Int(lattice.n_site * filling_fraction_vertex)
    @info "Bosonic FCI over $(sample_size[1])×$(sample_size[2]) Haldane honeycomb lattice with $(lattice.n_site) sites, $n_filled bosons (vertices filling $filling_fraction_vertex)"

    @info "Charge pump flux direction: $flux_direction; polarization direction: $polarization_direction"

    @assert mode in [:identity, :sectors]
    is_identity = mode == :identity
    labels = is_identity ? :identity : default_fci_sectors(sample_size)
    tag = is_identity ? "identity" : "sectors"
    PROJECT_ROOT = dirname(@__DIR__)

    result = flux_charge_pump(
        model,
        labels;
        filling_fraction=filling_fraction_vertex,
        flux_direction=flux_direction,
        polarization_direction=polarization_direction,
        twisted_phases_list=twisted_phases_list,
        nev_per_sector=1,
        fig_path=joinpath(PROJECT_ROOT, "figures",
            "bosonic_FCI_charge_pump_$(tag)_$(sample_size).svg"),
        checkpoint_path=joinpath(PROJECT_ROOT, "checkpoints",
            "bosonic_FCI_charge_pump_$(tag)_$(sample_size).jld2"),
    )

    @testset "Bosonic FCI charge pump ($(sample_size), $mode)" begin
        @test size(result.energies) == (length(twisted_phases_list), length(result.sector_labels), 1)
        @test all(isfinite, result.energies[:, :, 1])
        @test size(result.polarizations) == (length(twisted_phases_list), length(result.sector_labels))

        if mode == :sectors && isapprox(twisted_phases_list[end], 1.0; atol=1e-12)
            @test all(isapprox.(abs.(result.pumped_charges), 0.5; atol=atol))
        end
    end

    return result
end
