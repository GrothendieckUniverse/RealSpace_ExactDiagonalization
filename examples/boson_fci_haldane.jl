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
# 1. Model parameters (D.N. Sheng, PRL 107, 146803)
# ═══════════════════════════════════════════════════════════════════════════

function test_bose_hubbard()
    SAMPLE_SIZE = [2, 3]

    PARAMS = Dict(
        "t" => 1.0,        # nearest-neighbour hopping
        "t′" => 0.60,       # next-nearest-neighbour hopping
        "t′′" => -0.58,       # next-next-nearest-neighbour hopping
        "ϕ_over_2π" => 0.2,  # flux per 2π (time-reversal breaking)
        "V1" => 0.0,         # NN density interaction
        "V2" => 0.0,         # NNN density interaction
    )

    # ═══════════════════════════════════════════════════════════════════════════
    # 2. Build TightBinding model on Haldane honeycomb lattice
    # ═══════════════════════════════════════════════════════════════════════════

    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=SAMPLE_SIZE,
        brav_vec_list=[[1.0, 0.0], [1 / 2, sqrt(3) / 2]],
        sub_crys_list=[[0.0, 0.0], [1 / 3, 1 / 3]],
        lattice_name="Haldane_Honeycomb",
        pbc_indicator=[true, true],
    )

    tb_model = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="haldane_boson_FCI")

    t = PARAMS["t"]
    t′ = PARAMS["t′"]
    t′′ = PARAMS["t′′"]
    ϕ_over_2π = PARAMS["ϕ_over_2π"]

    # Nearest-neighbour hoppings (between sublattices) — real, amplitude -t
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, -1], 2)) => -t; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 0], 2)) => -t; is_hermitian=true)

    # Next-nearest-neighbour hoppings (within each sublattice) — complex, amplitude -t' e^{±iϕ}
    # Sublattice 1: clockwise = +ϕ
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([-1, 1], 1)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    # Sublattice 2: clockwise = -ϕ
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -t′ * exp(im * 2π * ϕ_over_2π); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([-1, 1], 2)) => -t′ * exp(-im * 2π * ϕ_over_2π); is_hermitian=true)

    # Next-next-nearest-neighbour hoppings (between sublattices) — real, amplitude -t′′
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 1)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 1), ([0, 1], 2)) => -t′′; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t′′; is_hermitian=true)

    lattice = tb_model.lattice
    n_site = lattice.n_site
    n_sub = lattice.n_sub
    @info "Haldane honeycomb: $(SAMPLE_SIZE) unit cells × $n_sub sublattices = $n_site sites"

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. Build second-quantised hard-core boson model
    # ═══════════════════════════════════════════════════════════════════════════

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((site_from, site_to), t) in tb_model.full_hopping_map
        i_from = lattice.site_to_index_map[site_from]
        i_to = lattice.site_to_index_map[site_to]
        push!(bilinear_terms, (i_from, i_to, ComplexF64(t)))
    end

    # Density-density interactions
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    V1 = PARAMS["V1"]
    V2 = PARAMS["V2"]
    for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
        (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)

        cell_shift = cell_to .- cell_from
        # V1: NN (inter-sublattice, same cell or shifted by one unit)
        if V1 != 0.0 &&
           ((mod.(cell_shift, lattice.sample_size) == [0, 0] && sub_from == 1 && sub_to == 2) ||
            (mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 2 && sub_to == 1) ||
            (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 2 && sub_to == 1))
            push!(density_terms, (i_from, i_to, ComplexF64(V1)))
        end
        # V2: NNN (same sublattice, shifted by one unit)
        if V2 != 0.0 &&
           ((mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 1 && sub_to == 1) ||
            (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 1 && sub_to == 1) ||
            (mod.(cell_shift, lattice.sample_size) == [1, 0] && sub_from == 2 && sub_to == 2) ||
            (mod.(cell_shift, lattice.sample_size) == [0, 1] && sub_from == 2 && sub_to == 2))
            push!(density_terms, (i_from, i_to, ComplexF64(V2)))
        end
    end

    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        PARAMS, lattice, tb_model, Bosonic(), bilinear_terms, density_terms,
    )

    n_filled = prod(SAMPLE_SIZE) ÷ 2   # 3 particles
    filling_frac = n_filled // n_site
    @info "FCI bosons: n_site=$n_site, n_filled=$n_filled, filling=$filling_frac"

    # ═══════════════════════════════════════════════════════════════════════════
    # 4. Translation symmetry and ED
    # ═══════════════════════════════════════════════════════════════════════════

    symmetry_group = build_translation_group(lattice)
    ed_data = build_ed_data(model; filling_fraction=filling_frac, symmetry_group=symmetry_group)

    println("\n" * "="^70)
    println("  Bosonic FCI — Haldane honeycomb $(SAMPLE_SIZE[1])×$(SAMPLE_SIZE[2])")
    println("  Full Hilbert space dim: $(binomial(n_site, n_filled))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momenta (k₁,k₂))")
    println("="^70)

    ed_scan!(ed_data; nev=5, mode=:matrix)

    # ═══════════════════════════════════════════════════════════════════════════
    # 5. Results — check FCI signature: two nearly-degenerate ground states
    # ═══════════════════════════════════════════════════════════════════════════

    all_vals = Float64[]
    all_info = Tuple{Int,Any,Int}[]
    for (irrep_idx, (vals, _)) in ed_data.ed_scan_res
        for (e_idx, v) in enumerate(vals)
            push!(all_vals, v)
            push!(all_info, (irrep_idx, ed_data.irrep_list[irrep_idx].label, e_idx))
        end
    end
    perm = sortperm(all_vals)

    println("\n--- Lowest 6 eigenvalues ---")
    for i in 1:min(6, length(all_vals))
        ii, label, ei = all_info[perm[i]]
        println("  E$i = $(@sprintf("%.10f", all_vals[perm[i]]))  (k = $(repr(label)), #$ei)")
    end

    if length(all_vals) >= 2
        Δ12 = all_vals[perm[2]] - all_vals[perm[1]]
        println("\n  ΔE₁₂ = $(@sprintf("%.8f", Δ12))")
        if Δ12 < 0.01
            println("  ✓ Nearly-degenerate ground-state pair detected  →  FCI signature!")
        end
    end
    if length(all_vals) >= 3
        Δ23 = all_vals[perm[3]] - all_vals[perm[2]]
        println("  ΔE₂₃ (many-body gap) = $(@sprintf("%.8f", Δ23))")
    end

    println("\nReference values (D.N. Sheng, PRL 107, 146803):")
    println("  E₀ ≈ -7.1638,  E₁ ≈ -7.1634")

    # Optional: plot spectrum
    print_spectrum(ed_data; shift_to_zero=true)
    fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
    save(joinpath(@__DIR__, "..", "figures", "haldane_fci_2x3_spectrum.svg"), fig)

    println("\nDone.")

end