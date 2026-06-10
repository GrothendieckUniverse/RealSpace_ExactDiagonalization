# ═══════════════════════════════════════════════════════════════════════════
# static_structure_factor — connected density-density correlation in k-space
#
# S^{αβ}(q) = (1/N) Σ_{i,j} e^{iq·(r_i-r_j)} (⟨n_i^α n_j^β⟩ − ⟨n_i^α⟩⟨n_j^β⟩)
#
# This is the primary diagnostic for charge order in translation-invariant
# PBC systems.  A CDW/Wigner crystal shows sharp Bragg peaks at the ordering
# wavevector Q; an FCI or superfluid shows a smooth, featureless S(q).
#
# See `doc/observables.ipynb` for the full physical discussion.
#
# ═══════════════════════════════════════════════════════════════════════════
# PERFORMANCE NOTE
#
# Both n̂_i (diagonal) and a_i^† a_j (off-diagonal) are computed DIRECTLY in
# the symmetry-sector (orbit) basis to avoid expanding the eigenvector into
# the full Fock space.  For a [4,4] Haldane system (32 sites, 8 bosons):
#
#   Full Fock space      : C(32,8) ≈ 10.5 M entries  (⇐ EXPENSIVE Dict!)
#   Orbit basis (|G|=16) :  ≈ 657 k representatives   (⇐ tractable!)
#
# The expansion step _expand_sector_state_to_fock_amplitudes would create
# a Dict with ~10.5 M entries (~hundreds of MB) and then loop over it with
# O(2.7 B) operations.  The orbit-basis approach avoids this entirely by
# iterating over the group orbit of each representative on the fly.
# ═══════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════
# Internal: build the connected density-density correlation C[i,j]
#
# Computes ⟨n_i n_j⟩ from the symmetry-sector eigenvector without expanding
# to the full Fock space.  For each representative mask r with amplitude
# |c_col|², we sum over all |G| orbit members π_g(r):
#
#   ⟨n_i n_j⟩ = Σ_col |c_col|² · (1/|G|) Σ_g n_i(π_g(r_col)) n_j(π_g(r_col))
#
# Each distinct orbit member receives total weight |Stab(r)|·|c_col|²/|G|,
# which is the correct projected-basis probability.
# ═══════════════════════════════════════════════════════════════════════════

function _build_density_correlation_matrix!(
    C::Matrix{Float64},
    density_a::Vector{Float64},
    density_b::Vector{Float64},
    c::Vector{ComplexF64},
    basis::Symmetry_Sector_Basis,
    n_site::Int,
    flavor_a::Function,
    flavor_b::Function,
    particle_statistics::Particle_Statistics,
)
    fill!(C, 0.0)
    fill!(density_a, 0.0)
    fill!(density_b, 0.0)

    G = basis.symmetry_group
    nG = group_order(G)
    nG_inv = 1.0 / nG

    @inbounds for (col, repr) in enumerate(basis.representative_mask_list)
        w = abs2(c[col])
        w == 0.0 && continue
        w_eff = w * nG_inv

        for (gidx, op) in enumerate(G.operations)
            shifted, _ = apply_operation_to_mask(repr, op, particle_statistics)

            # Collect occupied vertices of this orbit member
            tmp = shifted
            while tmp != 0
                lsb = tmp & -tmp
                i = trailing_zeros(lsb) + 1
                if flavor_a(i)
                    density_a[i] += w_eff
                end
                if flavor_b(i)
                    density_b[i] += w_eff
                end
                tmp ⊻= lsb
            end

            # Two-point contribution: n_i n_j  for i,j both occupied
            tmp2 = shifted
            while tmp2 != 0
                lsb2 = tmp2 & -tmp2
                i = trailing_zeros(lsb2) + 1
                tmp2 ⊻= lsb2
                if !flavor_a(i)
                    continue
                end
                tmp3 = shifted
                while tmp3 != 0
                    lsb3 = tmp3 & -tmp3
                    j = trailing_zeros(lsb3) + 1
                    tmp3 ⊻= lsb3
                    if flavor_b(j)
                        C[i, j] += w_eff
                    end
                end
            end
        end
    end

    # ── Subtract disconnected part: ⟨n_i n_j⟩_c = ⟨n_i n_j⟩ − ⟨n_i⟩⟨n_j⟩ ──
    @inbounds for j in 1:n_site, i in 1:n_site
        C[i, j] -= density_a[i] * density_b[j]
    end
    return C
end

# ═══════════════════════════════════════════════════════════════════════════
# Fourier transform helpers (reuse TightBinding.Uniform_Grids)
# ═══════════════════════════════════════════════════════════════════════════

"Precompute phases e^{-i q·r_i} for all (q-point, site) pairs."
function _precompute_phases(kgrid::TightBinding.Uniform_Grids, positions::Vector{<:Vector{Float64}})
    n_q = kgrid.nsite
    n_site = length(positions)
    phases = Matrix{ComplexF64}(undef, n_q, n_site)
    @inbounds for (q_idx, q) in enumerate(kgrid.site_cart_list)
        for i in 1:n_site
            phases[q_idx, i] = cis(-dot(q, positions[i]))
        end
    end
    return phases
end

function _precompute_phases(q_points::Vector{Vector{Float64}}, positions::Vector{<:Vector{Float64}})
    n_q = length(q_points)
    n_site = length(positions)
    phases = Matrix{ComplexF64}(undef, n_q, n_site)
    @inbounds for (q_idx, q) in enumerate(q_points)
        for i in 1:n_site
            phases[q_idx, i] = cis(-dot(q, positions[i]))
        end
    end
    return phases
end

"Fourier-transform C[i,j] → S(q) at all q in kgrid using precomputed phases."
function _fourier_transform!(S_q::AbstractVector{Float64}, phases::Matrix{ComplexF64},
    C::Matrix{Float64}, n_site::Int)
    n_q = size(phases, 1)
    @inbounds for q_idx in 1:n_q
        acc = 0.0
        for j in 1:n_site, i in 1:n_site
            acc += phases[q_idx, i] * C[i, j] * conj(phases[q_idx, j])
        end
        S_q[q_idx] = real(acc) / n_site
    end
    return S_q
end

# ═══════════════════════════════════════════════════════════════════════════
# Brillouin-zone boundary helpers for dense k-space visualisations
# ═══════════════════════════════════════════════════════════════════════════

function _reciprocal_basis_vectors(lattice::TightBinding.Real_Space_Lattice)
    lattice.dim == 2 || error("BZ boundary plotting is implemented only for 2D lattices.")
    a1 = Float64.(lattice.brav_vec_list[1])
    a2 = Float64.(lattice.brav_vec_list[2])
    A = [a1[1] a2[1]; a1[2] a2[2]]
    B = 2π * inv(A)'
    return [B[:, 1], B[:, 2]]
end

function _clip_polygon_by_halfplane(poly::Vector{Vector{Float64}},
    normal::Vector{Float64}, offset::Float64; atol::Float64=1e-12)
    isempty(poly) && return poly
    clipped = Vector{Float64}[]
    n = length(poly)
    @inbounds for idx in 1:n
        p = poly[idx]
        q = poly[mod1(idx + 1, n)]
        fp = dot(normal, p) - offset
        fq = dot(normal, q) - offset
        p_in = fp <= atol
        q_in = fq <= atol

        if p_in && q_in
            push!(clipped, q)
        elseif p_in && !q_in
            t = fp / (fp - fq)
            push!(clipped, p .+ t .* (q .- p))
        elseif !p_in && q_in
            t = fp / (fp - fq)
            push!(clipped, p .+ t .* (q .- p))
            push!(clipped, q)
        end
    end
    return clipped
end

function _first_bz_polygon(lattice::TightBinding.Real_Space_Lattice)
    b1, b2 = _reciprocal_basis_vectors(lattice)
    extent = 2.5 * max(norm(b1), norm(b2))
    poly = [[-extent, -extent], [extent, -extent], [extent, extent], [-extent, extent]]

    Gs = Vector{Float64}[]
    for n1 in -2:2, n2 in -2:2
        n1 == 0 && n2 == 0 && continue
        push!(Gs, n1 .* b1 .+ n2 .* b2)
    end
    sort!(Gs; by=norm)

    for G in Gs
        poly = _clip_polygon_by_halfplane(poly, G, dot(G, G) / 2)
        isempty(poly) && break
    end
    return poly
end

function _draw_bz_boundary!(ax, lattice::TightBinding.Real_Space_Lattice)
    poly = _first_bz_polygon(lattice)
    isempty(poly) && return ax
    xs = [p[1] for p in poly]
    ys = [p[2] for p in poly]
    push!(xs, poly[1][1])
    push!(ys, poly[1][2])
    lines!(ax, xs, ys; color=:black, linewidth=2.5)
    return ax
end

# ═══════════════════════════════════════════════════════════════════════════
# Shared ED bootstrap (build / resolve sector / return eigenvector & basis)
# ═══════════════════════════════════════════════════════════════════════════

function _bootstrap_ed_sector(
    model::Real_Space_Second_Quantized_Model,
    sector_label,
    filling_fraction::Rational{Int},
    target_eigval_idx::Int,
    ed_mode::Symbol,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing},
)::Tuple{Symmetry_Resolved_ED_Data,Symmetry_Sector_Basis,Vector{ComplexF64}}
    lattice = model.lattice
    n_site = lattice.n_site

    if ed_data === nothing
        flux0 = zeros(Float64, lattice.dim)
        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)

        G = (sector_label == :identity) ?
            build_identity_group(n_site) :
            build_translation_group(lattice, flux0)

        ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=G)

        scanned = (sector_label == :identity) ? nothing :
                  [Tuple{Int,Int}(sector_label)]

        ed_scan!(ed_data; nev=max(target_eigval_idx, 2), mode=ed_mode,
            scanned_sectors=scanned)
    end

    # ── Resolve sector ──
    if sector_label == :identity
        irrep_idx = findfirst(irrep -> irrep.label == :identity, ed_data.irrep_list)
    else
        irrep_idx = findfirst(irrep -> irrep.label == sector_label, ed_data.irrep_list)
    end
    irrep_idx === nothing &&
        error("Sector $(repr(sector_label)) not found in ed_data.irrep_list.")
    haskey(ed_data.ed_scan_res, irrep_idx) ||
        error("Sector $(repr(sector_label)) was not scanned.")

    _, vecs = ed_data.ed_scan_res[irrep_idx]
    target_eigval_idx <= size(vecs, 2) ||
        error("target_eigval_idx=$target_eigval_idx exceeds $(size(vecs,2)) eigenvectors.")

    basis = build_symmetry_sector_basis(ed_data.orbit_catalog, ed_data.irrep_list[irrep_idx])
    c = vecs[:, target_eigval_idx]

    return ed_data, basis, c
end

# ═══════════════════════════════════════════════════════════════════════════
# Public API — static structure factor
# ═══════════════════════════════════════════════════════════════════════════

"""
    static_structure_factor(
        model, sector_label;
        target_eigval_idx = 1,
        filling_fraction,
        flavor_a = i -> true,
        flavor_b = i -> true,
        ed_mode = :matrix,
        ed_data = nothing,
    ) -> (q_points::Vector{Vector{Float64}}, S_q::Vector{Float64})

Compute the connected static structure factor

```math
S^{\\alpha\\beta}(\\bm q) = \\frac{1}{N} \\sum_{i,j}
    e^{i\\bm q\\cdot(\\bm r_i-\\bm r_j)}
    \\big(\\langle n_i^\\alpha n_j^\\beta\\rangle
         - \\langle n_i^\\alpha\\rangle \\langle n_j^\\beta\\rangle\\big)
```

at every allowed crystal-momentum vector `q` in the first Brillouin zone.
The q-point mesh is built from `TightBinding.initialize_uniform_grids(lattice)`.

`sector_label` can be `:identity` (full real-space ED) or a momentum tuple
`(k₁, k₂)`.  The connected correlator is computed from the ground-state
eigenvector in the chosen symmetry sector.

`flavor_a` / `flavor_b` are boolean site filters selecting the flavour
components (e.g., sublattice A vs B).  The default `i -> true` selects all
sites, giving the total density structure factor.
"""
function static_structure_factor(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    flavor_a::Function=i -> true,
    flavor_b::Function=i -> true,
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::Tuple{Vector{Vector{Float64}},Vector{Float64}}
    lattice = model.lattice
    n_site = lattice.n_site
    stats = model.particle_statistics

    ed_data, basis, c = _bootstrap_ed_sector(
        model, sector_label, filling_fraction, target_eigval_idx, ed_mode, ed_data)

    # Build connected density-correlation matrix (orbit-basis, no full-Fock expansion)
    C = zeros(Float64, n_site, n_site)
    density_a = zeros(Float64, n_site)
    density_b = zeros(Float64, n_site)
    _build_density_correlation_matrix!(
        C, density_a, density_b, c, basis, n_site, flavor_a, flavor_b, stats)

    # Fourier transform on the BZ mesh from TightBinding
    kgrid = TightBinding.initialize_uniform_grids(lattice)
    phases = _precompute_phases(kgrid, lattice.site_cart_list)
    S_q = zeros(Float64, kgrid.nsite)
    _fourier_transform!(S_q, phases, C, n_site)

    return kgrid.site_cart_list, S_q
end

# ═══════════════════════════════════════════════════════════════════════════
# BZ heatmap: S(q) on a dense k-grid
# ═══════════════════════════════════════════════════════════════════════════

"""
    compute_structure_factor_map(
        model, sector_label;
        target_eigval_idx = 1,
        filling_fraction,
        flavor_a = i -> true,
        flavor_b = i -> true,
        k_resolution = 61,
        ed_data = nothing,
    ) -> (kx::Vector{Float64}, ky::Vector{Float64}, S_map::Matrix{Float64})

Compute S(q) on a dense regular grid spanning
`kx ∈ [-1.5π, 1.5π]`, `ky ∈ [-1.5π, 1.5π]`, returning the three arrays
suitable for `CairoMakie.heatmap(kx, ky, S_map)`.

This is the primary visualisation tool for diagnosing charge order:
a sharp Bragg peak at `Q ≠ 0` signals a CDW/Wigner crystal, while a smooth
featureless map signals a fluid phase (FCI or superfluid).
"""
function compute_structure_factor_map(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    flavor_a::Function=i -> true,
    flavor_b::Function=i -> true,
    k_resolution::Int=61,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::Tuple{Vector{Float64},Vector{Float64},Matrix{Float64}}
    lattice = model.lattice
    n_site = lattice.n_site
    stats = model.particle_statistics

    ed_data, basis, c = _bootstrap_ed_sector(
        model, sector_label, filling_fraction, target_eigval_idx, :matrix, ed_data)

    # Build connected density-correlation matrix (orbit-basis, no full-Fock expansion)
    C = zeros(Float64, n_site, n_site)
    density_a = zeros(Float64, n_site)
    density_b = zeros(Float64, n_site)
    _build_density_correlation_matrix!(
        C, density_a, density_b, c, basis, n_site, flavor_a, flavor_b, stats)

    # ── Dense Cartesian k-grid over [-1.5π, 1.5π]² ──
    span = 3π
    kx = collect(range(-span / 2, span / 2; length=k_resolution))
    ky = collect(range(-span / 2, span / 2; length=k_resolution))
    q_points = [[x, y] for y in ky for x in kx]

    phases = _precompute_phases(q_points, lattice.site_cart_list)
    S_flat = zeros(Float64, length(q_points))
    _fourier_transform!(S_flat, phases, C, n_site)

    S_map = reshape(S_flat, k_resolution, k_resolution)

    return kx, ky, S_map
end

"""
    plot_structure_factor_map(model, sector_label; kwargs...) -> Figure

Compute and plot the dense connected static structure-factor map with the first
Brillouin-zone boundary overlaid as solid black lines.  Keyword arguments are
forwarded to [`compute_structure_factor_map`](@ref), with two additional
optional plotting keywords:

- `fig_path::Union{Nothing,String}=nothing`
- `title::String="Connected static structure factor S(q)"`
"""
function plot_structure_factor_map(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    flavor_a::Function=i -> true,
    flavor_b::Function=i -> true,
    k_resolution::Int=61,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
    fig_path::Union{Nothing,String}=nothing,
    title::String="Connected static structure factor S(q)",
)
    kx, ky, S_map = compute_structure_factor_map(
        model, sector_label;
        target_eigval_idx=target_eigval_idx,
        filling_fraction=filling_fraction,
        flavor_a=flavor_a,
        flavor_b=flavor_b,
        k_resolution=k_resolution,
        ed_data=ed_data,
    )

    fig = Figure(size=(650, 560))
    ax = Axis(fig[1, 1]; xlabel="k_x", ylabel="k_y", title=title, aspect=DataAspect())
    hm = heatmap!(ax, kx, ky, S_map; colormap=:viridis)
    _draw_bz_boundary!(ax, model.lattice)
    Colorbar(fig[1, 2], hm; label="S(q)")

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end

"""
    plot_structure_factor_map_panels(maps; fig_path=nothing, title="Connected static structure factor S(q)")

Plot several precomputed static structure-factor maps in one row with a unified
color scale and first-BZ boundary overlays.

Each entry of `maps` is a named tuple with fields:
`kx`, `ky`, `values`, `lattice`, and `title`.
"""
function plot_structure_factor_map_panels(
    maps::AbstractVector;
    fig_path::Union{Nothing,String}=nothing,
    title::String="Connected static structure factor S(q)",
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
    Colorbar(fig[1, length(maps)+1], hm_ref; label="S(q)")

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end
