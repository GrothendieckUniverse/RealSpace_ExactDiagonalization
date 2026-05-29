#!/usr/bin/env julia
# ============================================================================
# Example 1: Spin-1/2 Heisenberg Chain
#
#   H = J Σ_{⟨i,j⟩} S_i · S_j   (PBC)
#
# via hard-core boson (Matsubara-Matsuda) mapping:
#   S^z_i = n_i − ½,   S^+_i = b†_i,   S^-_i = b_i
#
#   S_i·S_j = n_i n_j + ½(b†_i b_j + h.c.) − ½ n_i − ½ n_j + ¼
#
# Summing over bonds:  Σ(−½ n_i−½ n_j) = −Σ_i n_i = −N_e,  Σ ¼ = N/4.
# At half-filling N_e = N/2, the net constant is −J·N/4.
# This is absorbed via on-site density terms (i,i,−J/2).
# →  H_code ≡ H_spin  (no external constant-shift bookkeeping needed).
#
# Half-filling sector: N_e = N_site / 2  →  total S^z = 0
# Translation symmetry: Z_{N_site} (1D momentum sectors)
#
# Usage:
#   julia --project=. examples/spin_heisenberg_chain.jl
# ============================================================================

using RealSpace_ExactDiagonalization
using TightBinding
using Printf
using CairoMakie

# ═══════════════════════════════════════════════════════════════════════════
# Model parameters
# ═══════════════════════════════════════════════════════════════════════════
function test_Heisenberg_1D()
    N_SITE = 16          # chain length
    J = 1.0         # exchange coupling (J > 0 = antiferromagnetic)
    N_E = N_SITE ÷ 2  # half-filling  →  S^z_total = 0

    # ═══════════════════════════════════════════════════════════════════════════
    # 1. Build the 1D chain lattice via TightBinding (2D with Ly=1, PBC in x)
    # ═══════════════════════════════════════════════════════════════════════════

    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=[N_SITE, 1],
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0]],
        lattice_name="1D_Chain",
        pbc_indicator=[true, false],      # PBC only along x
    )

    tb_model = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="Heisenberg_Chain")

    # Add nearest-neighbour hopping (b†_i b_j + h.c.) with amplitude J/2
    # along x-direction only.  NOTE: add_hopping_term! automatically expands
    # using translation symmetry, so we only add ONE representative bond.
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => J / 2; is_hermitian=true)

    # ═══════════════════════════════════════════════════════════════════════════
    # 2. Build the second-quantised hard-core boson model
    # ═══════════════════════════════════════════════════════════════════════════

    lattice = tb_model.lattice
    n_site = lattice.n_site

    # Collect bilinear (hopping) terms from the TightBinding model
    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for ((site_from, site_to), t) in tb_model.full_hopping_map
        i_from = lattice.site_to_index_map[site_from]
        i_to = lattice.site_to_index_map[site_to]
        push!(bilinear_terms, (i_from, i_to, ComplexF64(t)))
    end

    # Density-density terms from the spin → hard-core boson mapping:
    #
    #   S_i·S_j = n_i n_j + ½(b†_i b_j + h.c.) − ½ n_i − ½ n_j + ¼
    #
    # Summing over bonds:  Σ_{⟨i,j⟩}(−½ n_i − ½ n_j) = −Σ_i n_i = −N_e
    #                       Σ_{⟨i,j⟩} ¼ = N/4
    # →  H_spin = J·[bond terms] + J/2·[hopping] − J·N_e + J·N/4
    #            = J·[bond terms] + J/2·[hopping] − J·N/4   (at half-filling)
    #
    # We absorb the −J·N/4 by adding on-site density terms (i,i,−J/2):
    #   Σ_i (−J/2)·n_i = −(J/2)·N_e = −J·N/4  ✓
    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    # Bond terms:  J * n_i * n_j
    for x in 0:(N_SITE-1)
        i = lattice.site_to_index_map[([x, 0], 1)]
        j = lattice.site_to_index_map[([mod(x + 1, N_SITE), 0], 1)]
        push!(density_terms, (i, j, ComplexF64(J)))
    end
    # On-site terms: −J/2 * n_i  (absorbs both −½ n_i − ½ n_j from expansion AND +¼ constant)
    for x in 0:(N_SITE-1)
        i = lattice.site_to_index_map[([x, 0], 1)]
        push!(density_terms, (i, i, ComplexF64(-J / 2)))
    end

    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("J" => J, "N_site" => N_SITE),
        lattice,
        tb_model,
        Bosonic(),
        bilinear_terms,
        density_terms,
    )

    @info "Heisenberg chain: N_site=$n_site, N_e=$N_E, filling=$(N_E//n_site)"

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. Translation symmetry: Z_{N_site}
    # ═══════════════════════════════════════════════════════════════════════════

    symmetry = build_translation_group(lattice)
    @info "Translation group: |G| = $(length(symmetry.operations))"

    # ═══════════════════════════════════════════════════════════════════════════
    # 4. Build ED data and run symmetry-resolved scan
    # ═══════════════════════════════════════════════════════════════════════════

    ed_data = build_ed_data(model; filling_fraction=N_E // n_site, symmetry=symmetry)

    println("\n" * "="^70)
    println("  Heisenberg Chain ED — N=$N_SITE, S^z=0 sector")
    println("  Full Hilbert space dim: $(binomial(n_site, N_E))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momentum sectors k = 0,…,N-1)")
    println("="^70)

    # Choose matrix-free mode (better for large systems)
    ed_scan!(ed_data; nev=3, mode=:matrixfree)

    # ═══════════════════════════════════════════════════════════════════════════
    # 5. Print results
    # ═══════════════════════════════════════════════════════════════════════════

    spec = print_spectrum(ed_data; shift_to_zero=true)

    println("\n--- Ground state ---")
    # Find the global minimum across all sectors
    all_vals = Float64[]
    all_info = Tuple{Int,Any,Int}[]
    for (irrep_idx, (vals, _)) in ed_data.ed_scan_res
        for (e_idx, v) in enumerate(vals)
            push!(all_vals, v)
            push!(all_info, (irrep_idx, ed_data.irrep_list[irrep_idx].label, e_idx))
        end
    end
    perm = sortperm(all_vals)
    for i in 1:min(5, length(all_vals))
        ii, label, ei = all_info[perm[i]]
        println("  E$i = $(@sprintf("%.10f", all_vals[perm[i]]))  (k = $(repr(label)), #$ei)")
    end

    # Known result for Heisenberg chain: E₀/N ≈ -ln(2) + 1/4 ≈ -0.443147... (Bethe ansatz)
    # The on-site −J·n_i terms are already included in the density terms,
    # so H_code ≡ H_spin (up to an irrelevant uniform constant).
    E0 = minimum(all_vals)
    e0_per_site = E0 / N_SITE
    println("\n  E₀ / N = $(@sprintf("%.8f", e0_per_site))")
    println("  Bethe-ansatz thermodynamic limit: -ln(2) + 1/4 ≈ -0.443147")

    # Optional: plot spectrum
    fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
    save(joinpath(@__DIR__, "..", "figures", "heisenberg_N$(N_SITE)_spectrum.svg"), fig)

    println("\nDone.")


end