#!/usr/bin/env julia
# ============================================================================
# Example 3: Spin-½ Fermi-Hubbard Model on Square Lattice
#
#   H = -t Σ_{⟨i,j⟩,σ} (c†_{iσ} c_{jσ} + h.c.)  +  U Σ_i n_{i↑} n_{i↓}
#
# Spinful fermions handled by flattening spin dof → interleaved graph:
#   site i↑ → vertex 2i-1,   site i↓ → vertex 2i
#
# Lattice: 2×3 (6 spatial sites → 12 graph vertices)
# Filling: half-filling per spin → N_up = N_down = 3, total N_e = 6
#
# Translation symmetry: Z₂ × Z₃ (6 momentum sectors)
#
# Usage:
#   julia --project=. examples/fermion_hubbard_square.jl
# ============================================================================

using RealSpace_ExactDiagonalization
using TightBinding
using Printf

# ═══════════════════════════════════════════════════════════════════════════
# 1. Model parameters
# ═══════════════════════════════════════════════════════════════════════════


function test_fermi_hubbard()
    SAMPLE_SIZE = [2, 3]   # Lx × Ly
    t = 1.0      # hopping
    U = 8.0      # on-site Hubbard repulsion

    # ═══════════════════════════════════════════════════════════════════════════
    # 2. Build spatial lattice (spinless, one sublattice per unit cell)
    # ═══════════════════════════════════════════════════════════════════════════

    r_data = TightBinding.initialize_real_space_lattice(;
        sample_size=SAMPLE_SIZE,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0]],
        lattice_name="Square_Hubbard",
        pbc_indicator=[true, true],
    )

    # Build a spinless TB model first, then duplicate for spin
    tb_spinless = TightBinding.initialize_real_space_tightbinding_model(r_data; model_name="hubbard_spinless")

    Lx, Ly = SAMPLE_SIZE
    # NOTE: add_hopping_term! already expands via translation symmetry,
    # so we only add ONE representative bond per direction.
    # x-neighbour (PBC): cell (0,0) → cell (1,0)
    add_hopping_term!(tb_spinless, (([0, 0], 1), ([1, 0], 1)) => ComplexF64(-t); is_hermitian=true)
    # y-neighbour (PBC): cell (0,0) → cell (0,1)
    add_hopping_term!(tb_spinless, (([0, 0], 1), ([0, 1], 1)) => ComplexF64(-t); is_hermitian=true)

    n_spatial = r_data.n_site   # Lx * Ly

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. Build "spinful" lattice: 2 sublattices (↑, ↓) per unit cell
    #    Vertex ordering: spatial i  →  vertex_up = 2i-1,  vertex_dn = 2i
    # ═══════════════════════════════════════════════════════════════════════════

    r_data_spinful = TightBinding.initialize_real_space_lattice(;
        sample_size=SAMPLE_SIZE,
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.0, 0.0], [0.0, 0.0]],   # two "sublattices" = ↑, ↓
        lattice_name="Square_Hubbard_Spinful",
        pbc_indicator=[true, true],
    )
    lattice = r_data_spinful
    n_site = lattice.n_site   # 2 * Lx * Ly

    # Map: spatial site s → up vertex = index of (cell_s, 1), down vertex = index of (cell_s, 2)
    function spatial_to_vertex(lattice, x::Int, y::Int)
        up = lattice.site_to_index_map[([x, y], 1)]
        dn = lattice.site_to_index_map[([x, y], 2)]
        return up, dn
    end

    # ═══════════════════════════════════════════════════════════════════════════
    # 4. Build hopping terms: duplicate spinless hoppings for each spin species
    # ═══════════════════════════════════════════════════════════════════════════

    bilinear_terms = Vector{Tuple{Int,Int,ComplexF64}}()

    for ((site_from, site_to), tamp) in tb_spinless.full_hopping_map
        x1, y1 = site_from[1]
        x2, y2 = site_to[1]
        up1, dn1 = spatial_to_vertex(lattice, x1, y1)
        up2, dn2 = spatial_to_vertex(lattice, x2, y2)
        # Hopping for ↑ spins: vertex up1 ↔ up2
        push!(bilinear_terms, (up1, up2, ComplexF64(tamp)))
        # Hopping for ↓ spins: vertex dn1 ↔ dn2
        push!(bilinear_terms, (dn1, dn2, ComplexF64(tamp)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    # 5. Density-density term: U Σ_i n_{i↑} n_{i↓}
    # ═══════════════════════════════════════════════════════════════════════════

    density_terms = Vector{Tuple{Int,Int,ComplexF64}}()
    for y in 0:(Ly-1), x in 0:(Lx-1)
        up, dn = spatial_to_vertex(lattice, x, y)
        push!(density_terms, (up, dn, ComplexF64(U)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    # 6. Assemble the second-quantised model (Fermionic!)
    # ═══════════════════════════════════════════════════════════════════════════

    tb_model = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="Hubbard")

    model = RealSpace_ExactDiagonalization.ShortRange_Real_Space_Second_Quantized_Model(
        Dict("t" => t, "U" => U),
        lattice,
        tb_model,
        Fermionic(),                # ← fermionic statistics!
        bilinear_terms,
        density_terms,
    )

    n_filled = n_spatial          # half-filling: one particle per spatial site (= N_up + N_down = n_spatial)
    filling_frac = n_filled // n_site  # = 1/2

    @info "Fermi-Hubbard: $(SAMPLE_SIZE) → $n_spatial spatial sites, $n_site graph vertices"
    @info "  N_e=$n_filled (N↑=N↓=$(n_spatial÷2)), filling=$filling_frac, t=$t, U=$U"

    # ═══════════════════════════════════════════════════════════════════════════
    # 7. Translation symmetry and ED
    # ═══════════════════════════════════════════════════════════════════════════

    # Translation group acts on spatial coords; for spinful lattice we must
    # build a custom permutation that moves ↑ and ↓ together.
    L1, L2 = SAMPLE_SIZE
    ops = Vector{Symmetry_Operation{Tuple{Int,Int}}}()
    for (δx, δy) in Iterators.product(0:(L1-1), 0:(L2-1))
        perm = Vector{Int}(undef, n_site)
        for (cell_int, isub) in lattice.site_list
            i_spin = isub       # 1 = ↑, 2 = ↓
            shifted = (mod.(cell_int .+ [δx, δy], SAMPLE_SIZE), i_spin)
            i_old = lattice.site_to_index_map[(cell_int, isub)]
            i_new = lattice.site_to_index_map[shifted]
            perm[i_old] = i_new
        end
        push!(ops, Symmetry_Operation((δx, δy), perm))
    end
    symmetry = Finite_Symmetry_Group("translations", ops; identity_idx=1)
    @info "Translation group: |G| = $(length(symmetry.operations))"

    ed_data = build_ed_data(model; filling_fraction=filling_frac, symmetry=symmetry)

    println("\n" * "="^70)
    println("  Fermi-Hubbard — square lattice $(SAMPLE_SIZE[1])×$(SAMPLE_SIZE[2])")
    println("  Full Hilbert space dim: $(binomial(n_site, n_filled))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momenta (k₁,k₂))")
    println("="^70)

    ed_scan!(ed_data; nev=5, mode=:matrix)

    # ═══════════════════════════════════════════════════════════════════════════
    # 8. Results
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

    println("\n--- Lowest 8 eigenvalues ---")
    for i in 1:min(8, length(all_vals))
        ii, label, ei = all_info[perm[i]]
        println("  E$i = $(@sprintf("%.10f", all_vals[perm[i]]))  (k = $(repr(label)), #$ei)")
    end

    # Optional: plot spectrum
    print_spectrum(ed_data; shift_to_zero=true)
    fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
    save(joinpath(@__DIR__, "..", "figures", "fermion_hubbard_square_2x3_spectrum.svg"), fig)

    println("\nDone.")


end