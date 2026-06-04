#!/usr/bin/env julia
# ============================================================================
# Example 2: Hard-core Boson FCI on Haldane Honeycomb Lattice
#
# Model from D.N. Sheng et al., Phys. Rev. Lett. 107, 146803 (2011).
#
# H = Σ t_ij b†_i b_j  +  Σ V_ij n_i n_j
#
# Lattice: 2×3 unit cells, 2 sublattices → 12 sites
# Filling: ν = 1/2 per band → 3 hard-core bosons
#
# Usage:
#   julia --project=. examples/boson_fci_haldane.jl
# ============================================================================

using RealSpace_ExactDiagonalization
using CairoMakie
using TightBinding
using Printf

# ═══════════════════════════════════════════════════════════════════════════
# Model parameters
# ═══════════════════════════════════════════════════════════════════════════
"parameters from Sun, Gu, Katsura, and Sarma [arXiv:1012.5864]"
const params_Sun_Gu_Katsura_Sarma = Dict(
    "t" => -1.0, # nearest-neighbor hopping
    "t′_1" => -1 / (2 + sqrt(2)), # next-nearest-neighbor hopping kinds 1
    "t′_2" => 1 / (2 + sqrt(2)), # next-nearest-neighbor hopping kinds 2
    "t′′" => -1 / (2 + 2 * sqrt(2)), # next-next-nearest-neighbor hopping
    "ϕ_over_2π" => 1 / 8, # flux per 2π (time-reversal breaking)
    "V1" => 2.0, # NN density-density interaction (inter-sublattice)
    "V2" => 1.0, # NNN interaction (same sublattice, x/y)
    "V3" => 0.0,  # NNNN interaction (same sublattice, diagonal)
    "λ" => 1.0, # ratio of `H_int / H_K`
)

"我自己的 optimial parameter, Quit large gap at 1/3!"
const my_optimal_param = Dict(
    "t" => -1,
    "t′_1" => -1.4 / (2 + sqrt(2)),
    "t′_2" => 1.4 / (2 + sqrt(2)),
    "t′′" => -1.0 / (2 + 2 * sqrt(2)),
    "ϕ_over_2π" => 1 / 8,
    "V1" => 2.0,
    "V2" => 0.15,
    "V3" => 0.2,
    "λ" => 3,
)



function test_fermion_fci_checkboard(;
    sample_size::Vector{Int}=[3, 4],
    params::Dict{String,<:Number}=Dict(
        "t" => -1.0,
        "t′_1" => -1 / (2 + sqrt(2)),
        "t′_2" => 1 / (2 + sqrt(2)),
        "t′′" => -1 / (2 + 2 * sqrt(2)),
        "ϕ_over_2π" => 1 / 8,
        "V1" => 2.0,
        "V2" => 1.0,
        "V3" => 0.0,
        "λ" => 1.0,
    ),
    filling_fraction_per_band::Rational{Int}=1 // 3
)
    params = merge(Dict{String,Float64}(), params_Sun_Gu_Katsura_Sarma, params)

    r_data = TightBinding.initialize_real_space_lattice(;
        lattice_name="checkerboard",
        sample_size=sample_size,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.5, 0], [0.0, 0.5]],
        pbc_indicator=[true, true],
    )
    tb_model = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="checkerboard")


    t, t′_1, t′_2, t′′ = params["t"], params["t′_1"], params["t′_2"], params["t′′"]
    ϕ = 2π * params["ϕ_over_2π"]

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

    lattice = tb_model.lattice
    n_site = lattice.n_site
    n_sub = lattice.n_sub
    @info "Checkerboard TB model: $(sample_size) unit cells × $n_sub sublattices = $n_site sites"

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. Build second-quantised hard-core boson model
    # ═══════════════════════════════════════════════════════════════════════════

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((site_from, site_to), t) in tb_model.full_hopping_map
        i_from = lattice.site_to_index_map[site_from]
        i_to = lattice.site_to_index_map[site_to]
        push!(bilinear_terms, (i_from, i_to, ComplexF64(t)))
    end

    # Hopping terms with zero twisted phases
    bilinear_terms = TightBinding.generate_bilinear_terms(tb_model;
        twisted_phases_over_2π=zeros(Float64, lattice.dim))

    # Density-density interactions
    V1, V2, V3 = params["V1"], params["V2"], params["V3"]

    # scale the interaction strength with `λ`
    V1 *= params["λ"]
    V2 *= params["λ"]
    V2 *= params["λ"]

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

    second_quantized_model = Real_Space_Second_Quantized_Model(
        params, lattice, tb_model, Fermionic(),
        bilinear_terms, density_terms,
    )

    lattice = second_quantized_model.lattice
    filling_fraction = filling_fraction_per_band // lattice.n_sub
    @info "Flatband Filling: $filling_fraction_per_band, Flatted Graph Filling: $filling_fraction"
    n_filled = Int(lattice.n_site * filling_fraction)
    @info "Fermionic FCI: $(sample_size[1])×$(sample_size[2]) checkerboard, " *
          "$(lattice.n_site) sites, $n_filled fermions (filling $filling_fraction)"

    # ═══════════════════════════════════════════════════════════════════════════
    # 4. Translation symmetry and ED
    # ═══════════════════════════════════════════════════════════════════════════

    symmetry_group = build_translation_group(lattice)
    ed_data = build_ed_data(second_quantized_model; filling_fraction=filling_fraction, symmetry_group=symmetry_group)

    println("\n" * "="^70)
    println("  Full Hilbert space dim: $(binomial(n_site, n_filled))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momenta (k₁,k₂))")
    println("="^70)

    ed_scan!(ed_data; nev=5, mode=:matrix)


    @info "Scaled interaction strengths: V1=$V1, V2=$V2, V3=$V3"


    # ═══════════════════════════════════════════════════════════════════════════
    # 5. Results — check FCI signature: two nearly-degenerate ground states
    # ═══════════════════════════════════════════════════════════════════════════

    # Optional: plot spectrum
    print_spectrum(ed_data; shift_to_zero=true)
    fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
    save(joinpath(@__DIR__, "..", "figures", "fermion_fci_checkerboard_ED_spec_$(sample_size)_ν=$(numerator(filling_fraction_per_band))_$(denominator(filling_fraction_per_band)).svg"), fig)

    println("\nDone.")
end