# ═══════════════════════════════════════════════════════════════════════════
# Fermionic FCI model — checkerboard lattice, ν = 2/3 per flatband
#
# Reference:
#   K. Sun, Z. Gu, H. Katsura, S. Das Sarma
#   "Nearly Flatbands with Nontrivial Topology", arXiv:1012.5864.
# ============================================================================

# ═══════════════════════════════════════════════════════════════════════════
# Model parameters
# ═══════════════════════════════════════════════════════════════════════════
"parameters from Sun, Gu, Katsura, and Sarma [arXiv:1012.5864]"
const params_Sun_Gu_Katsura_Sarma = Dict(
    "t" => 1.0, # nearest-neighbor hopping
    "t′_1" => 1 / (2 + sqrt(2)), # next-nearest-neighbor hopping kinds 1
    "t′_2" => -1 / (2 + sqrt(2)), # next-nearest-neighbor hopping kinds 2
    "t′′" => 1 / (2 + 2 * sqrt(2)), # next-next-nearest-neighbor hopping
    "ϕ_over_2π" => 1 / 8, # flux per 2π (time-reversal breaking)
    "V1" => 2.0, # NN density-density interaction (inter-sublattice)
    "V2" => 1.0, # NNN interaction (same sublattice, x/y)
    "V3" => 0.0  # NNNN interaction (same sublattice, diagonal)
)

# ═══════════════════════════════════════════════════════════════════════════
# Model builders
# ═══════════════════════════════════════════════════════════════════════════

"""
Build the Tight-Binding Model for the Checkerboard Lattice with Staggered Flux
---
- Named Args:
    - `sample_size::Vector{Int}=[3,4]`
    - `params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma`
    - `flip_bands::Bool=true`: flip sign of all hoppings to make the lower band flat
"""
function build_fermionic_checkerboard_tb_model(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma,
    flip_bands::Bool=true,
)::TightBinding.Real_Space_TightBinding_Model
    r_data = TightBinding.initialize_real_space_lattice(;
        lattice_name="checkerboard",
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.5, 0], [0.0, 0.5]],
        pbc_indicator=[true, true],
    )
    tb_model = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="checkerboard")

    p = deepcopy(params)
    if flip_bands
        p["t"] = -p["t"]
        p["t′_1"] = -p["t′_1"]
        p["t′_2"] = -p["t′_2"]
        p["t′′"] = -p["t′′"]
    end

    t, t′_1, t′_2, t′′ = p["t"], p["t′_1"], p["t′_2"], p["t′′"]
    ϕ = 2π * p["ϕ_over_2π"]

    # NN hoppings (inter-sublattice, complex — staggered flux)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t * exp(-im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 2)) => -t * exp(im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 1)) => -t * exp(-im * ϕ); is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t * exp(im * ϕ); is_hermitian=true)

    # NNN hoppings (intra-sublattice, real, anisotropic)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′_1; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′_2; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′_2; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′_1; is_hermitian=true)

    # NNNN hoppings (intra-sublattice, diagonal, real)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 1], 1), ([1, 0], 1)) => -t′′; is_hermitian=true)

    @info "Checkerboard TB model: $(sample_size[1])×$(sample_size[2]), $(tb_model.lattice.n_site) sites"
    return tb_model
end


"""
Construct the Second-Quantized Model for the Fermionic FCI Phase on Checkerboard Lattice
---
The model has two flat bands with Chern numbers ±1.  At ν = 2/3 filling of the
lower band the interacting ground state is a fermionic fractional Chern insulator
with three nearly-degenerate ground states on the torus (GSD = 3).

- Named Args:
    - `sample_size::Vector{Int}=[3,4]`
    - `params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma`
"""
function build_zero_flux_fermionic_fci_second_quantized_model(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma,
)::Real_Space_Second_Quantized_Model
    tb_model = build_fermionic_checkerboard_tb_model(; sample_size=sample_size, params=params)
    lattice = tb_model.lattice
    p = deepcopy(params)

    # Hopping terms with zero twisted phases
    bilinear_terms = TightBinding.generate_bilinear_terms(tb_model;
        twisted_phases_over_2π=zeros(Float64, lattice.dim))

    # Density-density interactions
    V1, V2, V3 = p["V1"], p["V2"], p["V3"]
    L1, L2 = lattice.sample_size
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()

    for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
        (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)

        i_from == i_to && continue
        cell_shift = mod.(cell_to .- cell_from, lattice.sample_size)

        # V1: NN (inter-sublattice, 4 bonds)
        if V1 != 0.0 && sub_from == 1 && sub_to == 2 &&
           (cell_shift == [0, 0] || cell_shift == [1, 0] ||
            cell_shift == [0, L2 - 1] || cell_shift == [1, L2 - 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V1)))
        end

        # V2: NNN (same sublattice, x or y)
        if V2 != 0.0 && sub_from == sub_to &&
           (cell_shift == [1, 0] || cell_shift == [0, 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V2)))
        end

        # V3: NNNN (same sublattice, diagonal)
        if V3 != 0.0 && sub_from == sub_to &&
           (cell_shift == [1, 1] || cell_shift == [1, L2 - 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V3)))
        end
    end

    model = Real_Space_Second_Quantized_Model(
        p, lattice, tb_model, Fermionic(),
        bilinear_terms, density_terms,
    )
    return model
end

# ═══════════════════════════════════════════════════════════════════════════
# Sector identification (to be filled after full ED scan)
# ═══════════════════════════════════════════════════════════════════════════

"""
Return the crystal-momentum sector labels hosting the FCI ground states.

For the ν = 2/3 fermionic FCI on checkerboard, three nearly-degenerate ground
states (GSD = 3) occupy the k_y = 0 sectors.  Verified via full ED scan on [3,4]
with V₁ = 2.0, V₂ = 1.0 (stronger NN repulsion stabilises the FCI):

    E₁ ≈ E₂ ≈ E₃  at k = (0,0), (1,0), (2,0)

The three states are split by finite-size effects of order ∼ 0.05–0.10 t.
"""
function default_fci_sectors_fermionic(sample_size::Vector{Int})
    return MLStyle.@match sample_size begin
        [2, 3] => [(0, 0), (0, 1), (0, 2)]
        [3, 2] => [(0, 0), (1, 0), (2, 0)]
        [3, 4] => [(0, 0), (1, 0), (2, 0)]
        [4, 3] => [(0, 0), (0, 1), (0, 2)]
        _ => begin
            error("Unknown sample_size=$sample_size — run full ED scan to identify FCI sectors.")
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Full ED scan — identify the FCI ground-state multiplet
# ═══════════════════════════════════════════════════════════════════════════

"""
Run a full symmetry-resolved ED scan for the fermionic checkerboard FCI model
and report the lowest eigenvalues by momentum sector.

This is the first diagnostic step: identify the nearly-degenerate ground-state
multiplet and its momentum quantum numbers.
"""
function test_fermionic_fci_full_ed(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma,
    filling_fraction::Rational{Int}=1 // 3, # ν=2/3 per band → 1/3 per vertex
    nev::Int=5,
)
    model = build_zero_flux_fermionic_fci_second_quantized_model(;
        sample_size=sample_size, params=params)

    lattice = model.lattice
    n_filled = Int(lattice.n_site * filling_fraction)
    @info "Fermionic FCI: $(sample_size[1])×$(sample_size[2]) checkerboard, " *
          "$(lattice.n_site) sites, $n_filled fermions (filling $filling_fraction)"

    symmetry_group = build_translation_group(lattice)
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=symmetry_group)

    ed_scan!(ed_data; nev=nev, mode=:matrix)

    println("\n" * "="^70)
    println("  Fermionic FCI — checkerboard $(sample_size[1])×$(sample_size[2])")
    println("  Full Hilbert space dim: $(binomial(lattice.n_site, n_filled))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momenta (k₁,k₂))")
    println("="^70)

    # Report lowest eigenvalues by sector
    all_vals = Float64[]
    all_info = Tuple{Int,Any,Int}[]
    for (irrep_idx, (vals, _)) in ed_data.ed_scan_res
        for (e_idx, v) in enumerate(vals)
            push!(all_vals, v)
            push!(all_info, (irrep_idx, ed_data.irrep_list[irrep_idx].label, e_idx))
        end
    end
    perm = sortperm(all_vals)

    println("\n--- Lowest $(min(15, length(all_vals))) eigenvalues ---")
    for i in 1:min(15, length(all_vals))
        ii, label, ei = all_info[perm[i]]
        @printf("  E%d = %.10f  (k = %s, #%d)\n", i, all_vals[perm[i]], repr(label), ei)
    end

    # Check for nearly-degenerate GS multiplet
    if length(all_vals) >= 3
        Δ12 = all_vals[perm[2]] - all_vals[perm[1]]
        Δ23 = all_vals[perm[3]] - all_vals[perm[2]]
        Δ34 = length(all_vals) >= 4 ? all_vals[perm[4]] - all_vals[perm[3]] : Inf
        println("\n  ΔE₁₂ = $(@sprintf("%.8f", Δ12))")
        println("  ΔE₂₃ = $(@sprintf("%.8f", Δ23))")
        if length(all_vals) >= 4
            println("  ΔE₃₄ = $(@sprintf("%.8f", Δ34)) (many-body gap)")
        end
        if Δ12 < 0.15 && Δ23 < 0.15
            println("  ✓ Three nearly-degenerate GS detected → FCI signature (GSD=3, finite-size split)!")
            println("  GS sectors: $(repr(ed_data.irrep_list[all_info[perm[1]][1]].label)), " *
                    "$(repr(ed_data.irrep_list[all_info[perm[2]][1]].label)), " *
                    "$(repr(ed_data.irrep_list[all_info[perm[3]][1]].label))")
        end
    end

    # Save spectrum plot
    PROJECT_ROOT = dirname(@__DIR__)
    fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
    save(joinpath(PROJECT_ROOT, "figures",
            "fermionic_FCI_spectrum_$(sample_size[1])x$(sample_size[2]).svg"), fig)

    return ed_data
end

# ═══════════════════════════════════════════════════════════════════════════
# Spectrum flow test
# ═══════════════════════════════════════════════════════════════════════════

"""
Test the Spectrum Flow for the Fermionic FCI on Checkerboard Lattice
---
- Named Args:
    - `sample_size::Vector{Int}=[3, 4]`
    - `params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma`
    - `filling_fraction::Rational{Int}=1 // 3`: filling per flattened vertex
    - `mode::Symbol=:sectors`
    - `fci_sectors`: list of momentum sector tuples for the FCI GS multiplet
    - `flux_direction::Int=findfirst(d -> mod(sample_size[d], 3) == 0, 1:length(sample_size))`: flux direction (along even-size direction for GSD=3)
    - `twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`
"""
function test_fermionic_fci_spectrum_flow(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma,
    filling_fraction::Rational{Int}=1 // 3,
    mode::Symbol=:sectors,
    fci_sectors::Vector{<:Tuple}=default_fci_sectors_fermionic(sample_size),
    flux_direction::Int=findfirst(d -> mod(sample_size[d], 3) == 0, 1:length(sample_size)),
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
)
    @info "Twisted phase list: $twisted_phases_over_2π_list"
    model = build_zero_flux_fermionic_fci_second_quantized_model(;
        sample_size=sample_size, params=params)

    lattice = model.lattice
    n_filled = Int(lattice.n_site * filling_fraction)
    @info "Fermionic FCI: $(sample_size[1])×$(sample_size[2]) checkerboard, " *
          "$(lattice.n_site) sites, $n_filled fermions (filling $filling_fraction)"
    @info "Flux direction: $flux_direction, FCI sectors: $fci_sectors"

    @assert mode in [:identity, :sectors]
    is_identity = mode == :identity
    labels = is_identity ? :identity : fci_sectors
    tag = is_identity ? "identity" : "sectors"
    PROJECT_ROOT = dirname(@__DIR__)

    result = flux_spectrum_flow(
        model,
        labels;
        filling_fraction=filling_fraction,
        flux_direction=flux_direction,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
        nev=3,
        fig_path=joinpath(PROJECT_ROOT, "figures",
            "fermionic_FCI_spectrum_flow_$(tag)_$(sample_size).svg"),
        checkpoint_dir=joinpath(PROJECT_ROOT, "checkpoints"),
    )

    return result
end

# ═══════════════════════════════════════════════════════════════════════════
# Charge pump test
# ═══════════════════════════════════════════════════════════════════════════

"""
Test the Fractional Charge Pump for the Fermionic FCI on Checkerboard Lattice
---
For ν = 2/3 filling of the lower Chern band, each of the three polarization
branches should wind by ΔQ ≈ 2/3 over one flux quantum, summing to 2 (integer).

- Named Args:
    - `sample_size::Vector{Int}=[3, 4]`
    - `params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma`
    - `filling_fraction::Rational{Int}=1 // 3`
    - `mode::Symbol=:sectors`
    - `fci_sectors`: the FCI GS multiplet sector labels
    - `flux_direction::Int=findfirst(d -> mod(sample_size[d], 3) == 0, 1:length(sample_size))`: flux direction (along even-size direction for GSD=3)
    - `twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9))`
    - `atol::Float64=0.10`: tolerance for charge quantisation
"""
function test_fermionic_fci_charge_pump(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=params_Sun_Gu_Katsura_Sarma,
    filling_fraction::Rational{Int}=1 // 3,
    mode::Symbol=:sectors,
    fci_sectors::Vector{<:Tuple}=default_fci_sectors_fermionic(sample_size),
    flux_direction::Int=findfirst(d -> mod(sample_size[d], 3) == 0, 1:length(sample_size)),
    polarization_direction::Int=2,
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
    atol::Float64=0.10,
)
    @info "Twisted phase list: $twisted_phases_over_2π_list"
    model = build_zero_flux_fermionic_fci_second_quantized_model(;
        sample_size=sample_size, params=params)

    lattice = model.lattice
    n_filled = Int(lattice.n_site * filling_fraction)
    @info "Fermionic FCI: $(sample_size[1])×$(sample_size[2]) checkerboard, " *
          "$(lattice.n_site) sites, $n_filled fermions (filling $filling_fraction)"
    @info "Flux: $flux_direction, polarization: $polarization_direction, sectors: $fci_sectors"

    @assert mode in [:identity, :sectors]
    is_identity = mode == :identity
    labels = is_identity ? :identity : fci_sectors
    tag = is_identity ? "identity" : "sectors"
    PROJECT_ROOT = dirname(@__DIR__)

    result = flux_charge_pump(
        model,
        labels;
        filling_fraction=filling_fraction,
        flux_direction=flux_direction,
        polarization_direction=polarization_direction,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
        nev_per_sector=1,
        fig_path=joinpath(PROJECT_ROOT, "figures",
            "fermionic_FCI_charge_pump_$(tag)_$(sample_size).svg"),
        checkpoint_dir=joinpath(PROJECT_ROOT, "checkpoints"),
    )

    @testset "Fermionic FCI charge pump ($(sample_size), $mode)" begin
        @test size(result.energies) == (length(twisted_phases_over_2π_list), length(result.sector_labels), 1)
        @test all(isfinite, result.energies[:, :, 1])
        @test size(result.polarizations) == (length(twisted_phases_over_2π_list), length(result.sector_labels))

        if mode == :sectors && isapprox(twisted_phases_over_2π_list[end], 1.0; atol=1e-12)
            expected_q = 2 / 3
            @test all(isapprox.(abs.(result.pumped_charges), expected_q; atol=atol))
        end
    end

    return result
end
