# ═══════════════════════════════════════════════════════════════════════════
# off_diagonal_long_range_order — one-body density matrix in k-space
#
# ρ_{ij} = ⟨a_i^† a_j⟩                               (one-body density matrix)
# ρ(k)   = (1/N) Σ_{i,j} e^{ik·(r_i-r_j)} ⟨a_i^† a_j⟩  (momentum distribution)
#
# This is the primary diagnostic for superfluidity / Bose condensation in
# translation-invariant PBC systems.  A superfluid shows macroscopic
# occupation (ρ(k) ∼ O(N)) at the condensing momentum k*; an FCI or CDW
# insulator shows ρ(k) ∼ O(1) with no single dominant eigenvalue.
#
# Reuses `_bootstrap_ed_sector`, `_precompute_phases` from
# `observables/static_structure_factor.jl` (loaded earlier).
#
# See `doc/observables.ipynb` for the full physical discussion.
#
# ═══════════════════════════════════════════════════════════════════════════
# PERFORMANCE NOTE
#
# ρ[i,j] = ⟨a_i^† a_j⟩ is computed by an exact orbit-basis contraction:
# for each projected-basis representative we loop over its symmetry orbit,
# apply the local one-body operator, and project the scattered mask back to
# the same sector.  This avoids materializing the full Fock-space wavefunction
# as a huge Dict.  For a [4,4] Haldane system:
#
#   Orbit-basis O(N_orbits × |G| × n_filled × N)
#       ≈ 657k × 16 × 8 × 32 ≈ 2.7 B local hops, but no full-amplitude Dict
#
# vs the Dict-expansion approach:
#   Expand:       10.5 M Dict insertions  (hundreds of MB)
#   Loop:         10.5 M × 8 × 32 ≈ 2.7 B iterations + hash lookups
# ═══════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════
# Core: one-body density matrix ⟨a_i^† a_j⟩
# ═══════════════════════════════════════════════════════════════════════════

"""
    off_diagonal_long_range_order(
        model, sector_label;
        target_eigval_idx = 1,
        filling_fraction,
        ed_mode = :matrix,
        ed_data = nothing,
    ) -> Matrix{ComplexF64}

Compute the one-body density matrix `ρ[i,j] = ⟨ψ| a_i^† a_j |ψ⟩` in the
chosen symmetry sector.

`sector_label` can be `:identity` (full real-space ED) or a momentum tuple
`(k₁, k₂)` for translation-resolved ED.

Returns an `N × N` complex matrix where `ρ[i,j]` is the expectation value.
"""
function off_diagonal_long_range_order(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::Matrix{ComplexF64}
    lattice = model.lattice
    n_site = lattice.n_site
    stats = model.particle_statistics

    ed_data, basis, c = _bootstrap_ed_sector(
        model, sector_label, filling_fraction, target_eigval_idx, ed_mode, ed_data)

    odlro = zeros(ComplexF64, n_site, n_site)

    G = basis.symmetry_group
    nG = group_order(G)
    inv_nG = 1.0 / nG
    inv_sqrt_nG = 1.0 / sqrt(nG)
    cmap = CanonicalMap(ed_data.symmetry_group, stats, Dict{Mask,Tuple{Mask,Int,ComplexF64}}())

    @inbounds for (col, repr_mask) in enumerate(basis.representative_mask_list)
        c_col = c[col]
        c_col == 0 && continue

        # In the normalized projected basis
        # |[s];χ⟩ = 1/sqrt(|G||Stab(s)|) Σ_g χ(g)^* U_g |s⟩.
        #
        # For diagonal n_j, stabilizer-related copies of the same raw Fock
        # state add coherently, giving the orbit-average weight |c_col|²/|G|
        # for every group image of the representative.
        diag_weight = abs2(c_col) * inv_nG
        ket_norm = inv_sqrt_nG / sqrt(basis.stabilizer_order_list[col])

        for (gidx, op) in enumerate(G.operations)
            ket_mask, ket_group_phase = apply_operation_to_mask(repr_mask, op, stats)
            ket_basis_amp = ket_norm * conj(basis.irrep.values[gidx]) * ket_group_phase
            ket_amp = c_col * ket_basis_amp

            tmp_j = ket_mask
            while tmp_j != 0
                lsb_j = tmp_j & -tmp_j
                j = trailing_zeros(lsb_j) + 1
                tmp_j ⊻= lsb_j

                odlro[j, j] += diag_weight

                for i in 1:n_site
                    i == j && continue
                    is_site_occupied(ket_mask, i) && continue

                    new_mask = empty_site_for_mask(ket_mask, j)
                    new_mask = occupy_site_for_mask(new_mask, i)

                    proj = project_to_sector(new_mask, basis, cmap)
                    proj === nothing && continue
                    row, proj_coeff = proj

                    # If project_to_sector(new_mask) = (row, p), then the
                    # coefficient of |new_mask⟩ in the normalized row basis is
                    # sqrt(|Stab(row)|/|G|) * conj(p).
                    bra_basis_amp = sqrt(basis.stabilizer_order_list[row]) *
                                    inv_sqrt_nG * conj(proj_coeff)
                    bra_amp = c[row] * bra_basis_amp

                    phase = hopping_phase_for_stats(stats, ket_mask, j, i)
                    odlro[i, j] += conj(bra_amp) * ket_amp * phase
                end
            end
        end
    end

    return odlro
end

# ═══════════════════════════════════════════════════════════════════════════
# BZ heatmap: ρ(k) on a dense k-grid
# ═══════════════════════════════════════════════════════════════════════════

"""
    compute_odlro_map(
        model, sector_label;
        target_eigval_idx = 1,
        filling_fraction,
        k_resolution = 61,
        ed_mode = :matrix,
        ed_data = nothing,
    ) -> (kx::Vector{Float64}, ky::Vector{Float64}, odlro_map::Matrix{Float64})

Compute the momentum-space ODLRO

```math
\\rho(\\bm k) = \\frac{1}{N} \\sum_{i,j}
    e^{i\\bm k\\cdot(\\bm r_i-\\bm r_j)} \\langle a_i^\\dagger a_j\\rangle
```

on a dense regular grid spanning `kx ∈ [-1.5π, 1.5π]`, `ky ∈ [-1.5π, 1.5π]`.

This is the primary visualisation tool for diagnosing superfluidity:
a sharp macroscopic peak at `k*` signals Bose condensation at that momentum;
a broad, low-amplitude distribution signals an insulating phase (FCI or CDW).
"""
function compute_odlro_map(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    k_resolution::Int=61,
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::Tuple{Vector{Float64},Vector{Float64},Matrix{Float64}}
    lattice = model.lattice
    n_site = lattice.n_site

    odlro = off_diagonal_long_range_order(
        model, sector_label;
        target_eigval_idx=target_eigval_idx,
        filling_fraction=filling_fraction,
        ed_mode=ed_mode,
        ed_data=ed_data,
    )

    # ── Dense Cartesian k-grid over [-1.5π, 1.5π]² ──
    span = 3π
    kx = collect(range(-span / 2, span / 2; length=k_resolution))
    ky = collect(range(-span / 2, span / 2; length=k_resolution))
    q_points = [[x, y] for y in ky for x in kx]

    phases = _precompute_phases(q_points, lattice.site_cart_list)

    # Fourier transform: ρ(k) = (1/N) Σ_{i,j} e^{ik·(r_i-r_j)} ρ_{ij}
    n_q = length(q_points)
    odlro_flat = zeros(Float64, n_q)
    @inbounds for q_idx in 1:n_q
        acc = zero(ComplexF64)
        for j in 1:n_site, i in 1:n_site
            acc += conj(phases[q_idx, i]) * odlro[i, j] * phases[q_idx, j]
        end
        odlro_flat[q_idx] = real(acc) / n_site
    end

    odlro_map = reshape(odlro_flat, k_resolution, k_resolution)
    return kx, ky, odlro_map
end

"""
    plot_odlro_map(model, sector_label; kwargs...) -> Figure

Compute and plot the dense momentum-space ODLRO map `ρ(k)` with the first
Brillouin-zone boundary overlaid as solid black lines.  Keyword arguments are
forwarded to [`compute_odlro_map`](@ref), with two additional optional plotting
keywords:

- `fig_path::Union{Nothing,String}=nothing`
- `title::String="Momentum distribution ρ(k)"`
"""
function plot_odlro_map(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    k_resolution::Int=61,
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
    fig_path::Union{Nothing,String}=nothing,
    title::String="Momentum distribution ρ(k)",
)
    kx, ky, odlro_map = compute_odlro_map(
        model, sector_label;
        target_eigval_idx=target_eigval_idx,
        filling_fraction=filling_fraction,
        k_resolution=k_resolution,
        ed_mode=ed_mode,
        ed_data=ed_data,
    )

    fig = Figure(size=(650, 560))
    ax = Axis(fig[1, 1]; xlabel="k_x", ylabel="k_y", title=title, aspect=DataAspect())
    hm = heatmap!(ax, kx, ky, odlro_map; colormap=:viridis)
    _draw_bz_boundary!(ax, model.lattice)
    Colorbar(fig[1, 2], hm; label="ρ(k)")

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end

"""
    plot_odlro_map_panels(maps; fig_path=nothing, title="Momentum distribution ρ(k)")

Plot several precomputed ODLRO maps in one row with a unified color scale and
first-BZ boundary overlays.

Each entry of `maps` is a named tuple with fields:
`kx`, `ky`, `values`, `lattice`, and `title`.
"""
function plot_odlro_map_panels(
    maps::AbstractVector;
    fig_path::Union{Nothing,String}=nothing,
    title::String="Momentum distribution ρ(k)",
)
    isempty(maps) && error("maps must not be empty.")

    vmin = minimum(minimum(m.values) for m in maps)
    vmax = maximum(maximum(m.values) for m in maps)

    fig = Figure(size=(420 * length(maps) + 110, 420))
    hm_ref = nothing
    for (idx, m) in enumerate(maps)
        ax = Axis(fig[1, idx]; xlabel="k_x", ylabel="k_y", title=m.title, aspect=DataAspect())
        hm = heatmap!(ax, m.kx, m.ky, m.values;
            colormap=:viridis, colorrange=(vmin, vmax))
        _draw_bz_boundary!(ax, m.lattice)
        hm_ref === nothing && (hm_ref = hm)
    end
    Label(fig[0, 1:length(maps)], title; fontsize=18)
    Colorbar(fig[1, length(maps)+1], hm_ref; label="ρ(k)")

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end
