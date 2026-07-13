# ═══════════════════════════════════════════════════════════════════════════
# Fermionic Checkerboard FCI — Phase Exploration
#
# Diagnose phases by scanning the dimensionless numerator x in
# t′′ = x / (2 + 2sqrt(2)) for the ν=1/3 fermionic checkerboard FCI model.
#
# All figures saved under `figures/phase_explore/`.
# ═══════════════════════════════════════════════════════════════════════════

const PROJECT_ROOT = dirname(@__DIR__)
const FIGURE_DIR = joinpath(PROJECT_ROOT, "figures", "phase_explore")
const CHECKERBOARD_TPP_DENOMINATOR = 2 + 2 * sqrt(2)

# ═══════════════════════════════════════════════════════════════════════════
# 1. Base parameters — fine-tuned for ν=1/3 FCI near t′′ numerator = -1
# ═══════════════════════════════════════════════════════════════════════════

"Fallback default parameters from Sun, Gu, Katsura, and Sarma [arXiv:1012.5864]"
const _default_params_fermionic_checkerboard = Dict{String,Float64}(
    "t" => -1.0,
    "t′_1" => -1 / (2 + sqrt(2)),
    "t′_2" => 1 / (2 + sqrt(2)),
    "t′′" => -1 / (2 + 2 * sqrt(2)),
    "ϕ_over_2π" => 1 / 8,
    "V1" => 2.0,
    "V2" => 1.0,
    "V3" => 0.0,
    "λ" => 1.0,
)

"Fine-tuned parameters for fermionic ν=1/3 FCI on checkerboard (λ=3 scales all V→λV)"
const my_optimal_param = Dict(
    "t" => -1.0,
    "t′_1" => -1.4 / (2 + sqrt(2)),
    "t′_2" => 1.4 / (2 + sqrt(2)),
    "t′′" => -1.0 / (2 + 2 * sqrt(2)),
    "ϕ_over_2π" => 1 / 8,
    "V1" => 2,
    "V2" => 0.45,
    "V3" => 0.2,
    "λ" => 2.8,
)


function _phase_explore_params_at_tpp(base_params::Dict{String,<:Number}, tpp_val::Float64)
    params = deepcopy(base_params)
    params["t′′"] = tpp_val / CHECKERBOARD_TPP_DENOMINATOR
    return params
end

# ═══════════════════════════════════════════════════════════════════════════
# 2. Reusable model builder — IDENTICAL to `test_fermion_fci_checkboard`
# ═══════════════════════════════════════════════════════════════════════════

"""
    _build_fermionic_checkerboard_model(params, sample_size, filling_fraction_per_band)

Build the fermionic checkerboard model exactly as `test_fermion_fci_checkboard`
and prepare symmetry-resolved ED data with the translation group.

Returns `(model, ed_data, n_filled, filling_fraction)`.
"""
function _build_fermionic_checkerboard_model(
    params::Dict{String,<:Number},
    sample_size::Vector{Int},
    filling_fraction_per_band::Rational{Int},
)::Tuple{Real_Space_Second_Quantized_Model,Symmetry_Resolved_ED_Data,Int,Rational{Int}}
    # Merge with defaults so any missing key gets the standard value
    params = merge(Dict{String,Float64}(), _default_params_fermionic_checkerboard, params)

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

    # ── Bilinear terms with zero twisted phases ──
    bilinear_terms = TightBinding.generate_bilinear_terms(tb_model;
        twisted_phases_over_2π=zeros(Float64, lattice.dim))

    # ── Density-density interactions (λ-scaled) ──
    V1, V2, V3 = params["V1"], params["V2"], params["V3"]

    # scale the interaction strength with `λ`
    V1 *= params["λ"]
    V2 *= params["λ"]
    V3 *= params["λ"]

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

    # ── Symmetry-resolved ED ──
    symmetry_group = build_translation_group(lattice)
    ed_data = build_ed_data(second_quantized_model; filling_fraction=filling_fraction, symmetry_group=symmetry_group)

    println("\n" * "="^70)
    println("  Full Hilbert space dim: $(binomial(n_site, n_filled))")
    println("  Orbits: $(length(ed_data.orbit_catalog.representative_mask_list))")
    println("  Irreps: $(length(ed_data.irrep_list))  (momenta (k₁,k₂))")
    println("="^70)

    @info "Scaled interaction strengths: V1=$V1, V2=$V2, V3=$V3"

    return second_quantized_model, ed_data, n_filled, filling_fraction
end

# ═══════════════════════════════════════════════════════════════════════════
# 3. Full ED spectrum scan — identify GS sectors at each t′′
# ═══════════════════════════════════════════════════════════════════════════

"""
    fermionic_phase_explore_full_ed(; t_values, sample_size, filling_fraction_per_band, nev)

Run full symmetry-resolved ED at each dimensionless t′′ numerator in `t_values`,
where the actual hopping is `t′′ = t_value / (2 + 2sqrt(2))`. Report lowest
eigenvalues with sector labels (both k-tuples and 0-based irrep indices),
identify ground-state sectors, and save spectrum plots under
`figures/phase_explore/`.

Returns `Dict{Float64, Vector{Tuple{Int,Int}}}` mapping t′′ numerator → GS
sector k-labels.
"""
function fermionic_phase_explore_full_ed(;
    t_values::Vector{Float64}=[-3.0, -1.0],
    base_params::Dict{String,<:Number}=my_optimal_param,
    sample_size::Vector{Int}=[3, 4],
    filling_fraction_per_band::Rational{Int}=1 // 3,
    nev::Int=5,
    gs_window::Float64=0.075,
)::Dict{Float64,Vector{Tuple{Int,Int}}}
    mkpath(FIGURE_DIR)
    gs_sectors_map = Dict{Float64,Vector{Tuple{Int,Int}}}()

    for t_val in t_values
        params = _phase_explore_params_at_tpp(base_params, t_val)

        model, ed_data, n_filled, filling_fraction = _build_fermionic_checkerboard_model(
            params, sample_size, filling_fraction_per_band)

        ed_scan!(ed_data; nev=nev, mode=:matrix)

        # ── Collect all eigenvalues ──
        all_vals = Float64[]
        all_info = Tuple{Int,Tuple{Int,Int},Int}[]  # (irrep_idx, k_label, eig_idx)
        for (irrep_idx, (vals, _)) in ed_data.ed_scan_res
            for (e_idx, v) in enumerate(vals)
                push!(all_vals, v)
                push!(all_info, (irrep_idx, ed_data.irrep_list[irrep_idx].label, e_idx))
            end
        end
        perm = sortperm(all_vals)

        # ── Print spectrum ──
        println("\n" * "="^70)
        println("  Phase Explore @ t′′ numerator = $t_val")
        println("  Actual t′′ = $(params["t′′"])")
        println("  $(sample_size[1])×$(sample_size[2]) checkerboard, $n_filled fermions, ν=$filling_fraction")
        println("  Hilbert dim = $(binomial(model.lattice.n_site, n_filled))")
        println("="^70)

        L1, L2 = sample_size
        println("  k-point indexing (k₁-major):")
        for idx in 0:min(11, length(ed_data.irrep_list) - 1)
            k = ed_data.irrep_list[idx+1].label
            print("    [$idx]→($(k[1]),$(k[2]))")
            mod(idx + 1, 4) == 0 && println()
        end
        mod(length(ed_data.irrep_list), 4) != 0 && println()
        println()

        n_show = min(12, length(all_vals))
        println("  Lowest $n_show eigenvalues:")
        for i in 1:n_show
            ii, k, ei = all_info[perm[i]]
            @printf("    E%-2d = %12.8f   k=(%d,%d)   irrep_idx=%d  (#%d)\n",
                i, all_vals[perm[i]], k[1], k[2], ii, ei)
        end

        # ── Identify GS sectors inside the requested low-energy window ──
        e0 = all_vals[perm[1]]
        gs_sectors = Tuple{Int,Int}[]
        gs_gap = Inf
        for i in 1:length(all_vals)
            if all_vals[perm[i]] - e0 < gs_window
                push!(gs_sectors, all_info[perm[i]][2])
            else
                gs_gap = all_vals[perm[i]] - e0
                break
            end
        end
        unique!(gs_sectors)
        sort!(gs_sectors)
        gs_sectors_map[t_val] = gs_sectors

        gs_irrep_indices = Int[]
        for k in gs_sectors
            idx = findfirst(irrep -> irrep.label == k, ed_data.irrep_list)
            idx !== nothing && push!(gs_irrep_indices, idx - 1)  # 0-based
        end

        println("\n  ── GS sector identification ──")
        println("  GSD = $(length(gs_sectors))")
        println("  GS k-labels   : $gs_sectors")
        println("  GS irrep idx  : $gs_irrep_indices  (0-based)")
        println("  GS window     : $gs_window")
        @printf("  GS gap        : %.6f\n", gs_gap)

        # ── Save spectrum plot ──
        fig_path = joinpath(FIGURE_DIR,
            "fermionic_phase_explore_spectrum_tpp=$(t_val)_$(sample_size).svg")
        fig, ax = plot_spectrum(ed_data; shift_to_zero=true)
        save(fig_path, fig)
        @info "  Saved spectrum: $fig_path"
    end

    return gs_sectors_map
end

# ═══════════════════════════════════════════════════════════════════════════
# 4. Spectrum flow — gap protection under twisted boundary conditions
# ═══════════════════════════════════════════════════════════════════════════

"""
    fermionic_phase_explore_spectrum_flow(; t_values, gs_sectors_map, sample_size, ...)

For each t′′ numerator, run `flux_spectrum_flow` using the identified GS
sectors. Saves figures under `figures/phase_explore/`.
"""
function fermionic_phase_explore_spectrum_flow(;
    t_values::Vector{Float64}=[-3.0, -1.0],
    gs_sectors_map::Dict{Float64,Vector{Tuple{Int,Int}}},
    base_params::Dict{String,<:Number}=my_optimal_param,
    sample_size::Vector{Int}=[3, 4],
    filling_fraction_per_band::Rational{Int}=1 // 3,
    flux_direction::Int=1,
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
)::Dict{Float64,Any}
    mkpath(FIGURE_DIR)
    results = Dict{Float64,Any}()

    for t_val in t_values
        haskey(gs_sectors_map, t_val) ||
            error("GS sectors not found for t′′ numerator=$t_val. Run full ED first.")
        sectors = gs_sectors_map[t_val]

        params = _phase_explore_params_at_tpp(base_params, t_val)

        model, _, n_filled, filling_fraction = _build_fermionic_checkerboard_model(
            params, sample_size, filling_fraction_per_band)

        @info "Spectrum flow @ t′′ numerator=$t_val: $(length(sectors)) sectors $sectors, flux dir=$flux_direction"

        fig_path = joinpath(FIGURE_DIR,
            "fermionic_phase_explore_spectrum_flow_tpp=$(t_val)_$(sample_size).svg")

        result = flux_spectrum_flow(
            model, sectors;
            filling_fraction=filling_fraction,
            flux_direction=flux_direction,
            twisted_phases_over_2π_list=twisted_phases_over_2π_list,
            nev=5,
            fig_path=fig_path,
            checkpoint_dir=joinpath(PROJECT_ROOT, "checkpoints"),
        )
        results[t_val] = result
        @info "  Saved spectrum flow: $fig_path"
    end

    return results
end

# ═══════════════════════════════════════════════════════════════════════════
# 5. Charge pump — many-body Chern number
# ═══════════════════════════════════════════════════════════════════════════

"""
    fermionic_phase_explore_charge_pump(; t_values, gs_sectors_map, sample_size, ...)

For each t′′ numerator, run `flux_charge_pump` using the identified GS sectors.
Saves figures under `figures/phase_explore/`.
"""
function fermionic_phase_explore_charge_pump(;
    t_values::Vector{Float64}=[-3.0, -1.0],
    gs_sectors_map::Dict{Float64,Vector{Tuple{Int,Int}}},
    base_params::Dict{String,<:Number}=my_optimal_param,
    sample_size::Vector{Int}=[3, 4],
    filling_fraction_per_band::Rational{Int}=1 // 3,
    flux_direction::Int=1,
    polarization_direction::Int=2,
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
)::Dict{Float64,Any}
    mkpath(FIGURE_DIR)
    results = Dict{Float64,Any}()

    for t_val in t_values
        haskey(gs_sectors_map, t_val) ||
            error("GS sectors not found for t′′ numerator=$t_val. Run full ED first.")
        sectors = gs_sectors_map[t_val]

        params = _phase_explore_params_at_tpp(base_params, t_val)

        model, _, n_filled, filling_fraction = _build_fermionic_checkerboard_model(
            params, sample_size, filling_fraction_per_band)

        @info "Charge pump @ t′′ numerator=$t_val: $(length(sectors)) sectors $sectors"
        @info "  Flux dir=$flux_direction, polarization dir=$polarization_direction"

        fig_path = joinpath(FIGURE_DIR,
            "fermionic_phase_explore_charge_pump_tpp=$(t_val)_$(sample_size).svg")

        result = flux_charge_pump(
            model, sectors;
            filling_fraction=filling_fraction,
            flux_direction=flux_direction,
            polarization_direction=polarization_direction,
            twisted_phases_over_2π_list=twisted_phases_over_2π_list,
            nev_per_sector=1,
            fig_path=fig_path,
            checkpoint_dir=joinpath(PROJECT_ROOT, "checkpoints"),
        )
        results[t_val] = result
        @info "  Saved charge pump: $fig_path"
    end

    return results
end

# ═══════════════════════════════════════════════════════════════════════════
# 6. Static structure factor — 3-panel unified colour scale
# ═══════════════════════════════════════════════════════════════════════════

"""
    _sf_peak_report(S_map, kx, ky; threshold=0.1)

Report the dominant peak(s) in a static structure factor map.
Returns `(peak_kx, peak_ky, peak_value, n_peaks_above_threshold)`.
"""
function _sf_peak_report(
    S_map::Matrix{Float64},
    kx::Vector{Float64},
    ky::Vector{Float64};
    threshold_frac::Float64=0.5,
)
    vmin, vmax = minimum(S_map), maximum(S_map)
    span = vmax - vmin
    span == 0.0 && return (0.0, 0.0, vmax, 0)

    threshold = vmin + threshold_frac * span
    n_peaks = count(x -> x >= threshold, S_map)

    idx = argmax(S_map)
    ix, iy = idx[1], idx[2]
    return (kx[ix], ky[iy], S_map[ix, iy], n_peaks)
end

function _sf_allowed_peak_report(
    qx::Vector{Float64},
    qy::Vector{Float64},
    S_q::Vector{Float64};
    threshold_frac::Float64=0.5,
)
    vmin, vmax = minimum(S_q), maximum(S_q)
    span = vmax - vmin
    span == 0.0 && return (0.0, 0.0, vmax, 0)

    threshold = vmin + threshold_frac * span
    n_peaks = count(x -> x >= threshold, S_q)

    idx = argmax(S_q)
    return (qx[idx], qy[idx], S_q[idx], n_peaks)
end

"""
    fermionic_phase_explore_structure_factor(; t_values, sample_size, k_resolution)

Compute the connected static structure factor S(q) for each t′′ numerator and
plot them in a single multi-panel figure with unified colour scale and first-BZ
boundary overlays. Saved under `figures/phase_explore/`.

- Sharp Bragg peaks at Q ≠ 0 → charge-density wave
- Smooth, featureless map → fluid phase (FCI or Fermi liquid)
"""
function fermionic_phase_explore_structure_factor(;
    t_values::Vector{Float64}=[-3.0, -1.0],
    gs_sectors_map::Union{Dict{Float64,Vector{Tuple{Int,Int}}},Nothing}=nothing,
    base_params::Dict{String,<:Number}=my_optimal_param,
    sample_size::Vector{Int}=[3, 4],
    filling_fraction_per_band::Rational{Int}=1 // 3,
    k_resolution::Int=61,
)
    mkpath(FIGURE_DIR)

    maps = NamedTuple[]
    peaks_info = NamedTuple[]

    for t_val in t_values
        params = _phase_explore_params_at_tpp(base_params, t_val)

        model, ed_data, n_filled, filling_fraction = _build_fermionic_checkerboard_model(
            params, sample_size, filling_fraction_per_band)

        sector = if gs_sectors_map !== nothing && haskey(gs_sectors_map, t_val) &&
                    !isempty(gs_sectors_map[t_val])
            gs_sectors_map[t_val][1]
        else
            (0, 0)
        end

        kx, ky, S_map = compute_structure_factor_map(
            model, sector;
            target_eigval_idx=1,
            filling_fraction=filling_fraction,
            k_resolution=k_resolution,
            ed_data=nothing,
        )

        pkx, pky, pval, n_peaks = _sf_peak_report(S_map, kx, ky)

        push!(maps, (
            kx=kx,
            ky=ky,
            values=S_map,
            lattice=model.lattice,
            title="t'' = $t_val, k = $sector",
        ))
        push!(peaks_info, (
            t=t_val,
            peak_kx=round(pkx, digits=4),
            peak_ky=round(pky, digits=4),
            peak_value=round(pval, digits=6),
            n_peaks_above_half=n_peaks,
        ))

        @info "  S(q) @ t′′ numerator=$t_val, sector=$sector: max=$(round(pval,digits=6)) at k=($(round(pkx,digits=3)),$(round(pky,digits=3)))"
    end

    # ── 3-panel unified figure ──
    fig_path = joinpath(FIGURE_DIR,
        "fermionic_phase_explore_structure_factor_$(sample_size).svg")
    plot_structure_factor_map_panels(
        maps;
        fig_path=fig_path,
        title="Fermionic Checkerboard S(q): t'' = $t_values  (ν=1/3 per band)",
    )
    @info "  Saved 3-panel S(q): $fig_path"

    return maps, peaks_info
end

"""
    fermionic_phase_explore_allowed_structure_factor(; t_values, sample_size, ...)

Compute the connected static structure factor only at the finite-size allowed
momenta for each t′′ numerator and plot the results in a unified multi-panel
scatter figure.  This is the cleanest plot for reading the ordering wavevector
on a finite torus.
"""
function fermionic_phase_explore_allowed_structure_factor(;
    t_values::Vector{Float64}=[-3.0, -1.0],
    gs_sectors_map::Union{Dict{Float64,Vector{Tuple{Int,Int}}},Nothing}=nothing,
    base_params::Dict{String,<:Number}=my_optimal_param,
    sample_size::Vector{Int}=[3, 4],
    filling_fraction_per_band::Rational{Int}=1 // 3,
    markersize::Real=28,
)
    mkpath(FIGURE_DIR)

    maps = NamedTuple[]
    peaks_info = NamedTuple[]

    for t_val in t_values
        params = _phase_explore_params_at_tpp(base_params, t_val)

        model, ed_data, n_filled, filling_fraction = _build_fermionic_checkerboard_model(
            params, sample_size, filling_fraction_per_band)

        sector = if gs_sectors_map !== nothing && haskey(gs_sectors_map, t_val) &&
                    !isempty(gs_sectors_map[t_val])
            gs_sectors_map[t_val][1]
        else
            (0, 0)
        end

        qx, qy, S_q = structure_factor_allowed_momenta(
            model, sector;
            target_eigval_idx=1,
            filling_fraction=filling_fraction,
            ed_data=nothing,
            fold_to_first_bz=true,
        )

        pkx, pky, pval, n_peaks = _sf_allowed_peak_report(qx, qy, S_q)

        push!(maps, (
            qx=qx,
            qy=qy,
            values=S_q,
            lattice=model.lattice,
            title="t'' = $t_val, k = $sector",
        ))
        push!(peaks_info, (
            t=t_val,
            peak_qx=round(pkx, digits=4),
            peak_qy=round(pky, digits=4),
            peak_value=round(pval, digits=6),
            n_peaks_above_half=n_peaks,
        ))

        @info "  Allowed S(q) @ t′′ numerator=$t_val, sector=$sector: max=$(round(pval,digits=6)) at q=($(round(pkx,digits=3)),$(round(pky,digits=3)))"
    end

    fig_path = joinpath(FIGURE_DIR,
        "fermionic_phase_explore_structure_factor_allowed_$(sample_size).svg")
    plot_structure_factor_allowed_momenta_panels(
        maps;
        fig_path=fig_path,
        title="Fermionic Checkerboard allowed S(q): t'' = $t_values  (ν=1/3 per band)",
        markersize=markersize,
    )
    @info "  Saved allowed-momentum S(q): $fig_path"

    return maps, peaks_info
end

# ═══════════════════════════════════════════════════════════════════════════
# 7. Master demo — run everything
# ═══════════════════════════════════════════════════════════════════════════

"""
    fermionic_phase_explore_demo(; sample_size, t_values, k_resolution, ...)

Master entry point for the fermionic checkerboard FCI phase exploration.
Here `t_values` are dimensionless t′′ numerators, i.e.
`params["t′′"] = t_value / (2 + 2sqrt(2))`.

Runs four diagnostics:
1. **Full ED spectrum** — identify ground-state sectors at each t′′
2. **Spectrum flow** — gap protection under twisted boundary conditions
3. **Charge pump** — many-body Chern number (fractional → FCI, trivial → CDW/FL)
4. **Static structure factor S(q)** — 3-panel with unified colour scale
   (Bragg peaks → CDW; smooth → fluid)

All figures are saved under `figures/phase_explore/`.
"""
function fermionic_phase_explore_demo(;
    sample_size::Vector{Int}=[3, 4],
    t_values::Vector{Float64}=[-3.0, -1.0],
    base_params::Dict{String,<:Number}=my_optimal_param,
    filling_fraction_per_band::Rational{Int}=1 // 3,
    k_resolution::Int=61,
    nev::Int=8,
    gs_window::Float64=0.075,
    flux_direction::Int=1,
    polarization_direction::Int=2,
    twisted_phases_over_2π_list::Vector{Float64}=collect(range(0.0, 1.0; length=9)),
)
    println("\n" * "="^72)
    println("  Fermionic Checkerboard FCI — Phase Exploration Demo")
    println("  Lattice: $(sample_size[1])×$(sample_size[2]) checkerboard")
    println("  Filling: $filling_fraction_per_band per band  (ν=1/3 per band)")
    println("  t′′ numerator values: $t_values")
    println("  Actual t′′ = value / $(CHECKERBOARD_TPP_DENOMINATOR)")
    println("  GS window: $gs_window")
    println("  Base params: λ=$(base_params["λ"]), V1=$(base_params["V1"]), V2=$(base_params["V2"]), V3=$(base_params["V3"])")
    println("  Figures → $FIGURE_DIR")
    println("="^72)

    # ── Step 1: Full ED scan ──
    println("\n── Step 1/4: Full ED spectrum scan ──")
    gs_sectors_map = fermionic_phase_explore_full_ed(;
        t_values=t_values,
        base_params=base_params,
        sample_size=sample_size,
        filling_fraction_per_band=filling_fraction_per_band,
        nev=nev,
        gs_window=gs_window,
    )

    # ── Step 2: Spectrum flow ──
    println("\n── Step 2/4: Spectrum flow ──")
    flow_results = fermionic_phase_explore_spectrum_flow(;
        t_values=t_values,
        gs_sectors_map=gs_sectors_map,
        base_params=base_params,
        sample_size=sample_size,
        filling_fraction_per_band=filling_fraction_per_band,
        flux_direction=flux_direction,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
    )

    # ── Step 3: Charge pump ──
    println("\n── Step 3/4: Charge pump ──")
    pump_results = fermionic_phase_explore_charge_pump(;
        t_values=t_values,
        gs_sectors_map=gs_sectors_map,
        base_params=base_params,
        sample_size=sample_size,
        filling_fraction_per_band=filling_fraction_per_band,
        flux_direction=flux_direction,
        polarization_direction=polarization_direction,
        twisted_phases_over_2π_list=twisted_phases_over_2π_list,
    )

    # ── Step 4: Static structure factor ──
    println("\n── Step 4/4: Static structure factor S(q) ──")
    sf_maps, sf_peaks = fermionic_phase_explore_structure_factor(;
        t_values=t_values,
        gs_sectors_map=gs_sectors_map,
        base_params=base_params,
        sample_size=sample_size,
        filling_fraction_per_band=filling_fraction_per_band,
        k_resolution=k_resolution,
    )

    # ── Summary diagnosis table ──
    println("\n" * "="^72)
    println("  PHASE DIAGNOSIS SUMMARY")
    println("="^72)

    L1, L2 = sample_size
    println("\n  k-point indexing for [$L1,$L2] (k₁-major, k₂ inner):")
    for k1 in 0:(L1-1), k2 in 0:(L2-1)
        idx = k1 * L2 + k2
        print("    [$idx]→($k1,$k2)")
        k2 == L2 - 1 && println()
    end
    println()

    println("  ┌──────────┬─────────────────────┬─────┬──────────────────┬──────────────────────┐")
    println("  │   t''    │  GS k-labels        │ GSD │  Charge pump ΣΔQ │  S(q) diagnosis      │")
    println("  ├──────────┼─────────────────────┼─────┼──────────────────┼──────────────────────┤")

    for t_val in t_values
        sectors = gs_sectors_map[t_val]
        gsd = length(sectors)
        k_str = join(["($(k[1]),$(k[2]))" for k in sectors], ", ")

        # Charge pump sum
        if haskey(pump_results, t_val)
            charges = pump_results[t_val].pumped_charges
            sum_q = round(sum(abs.(charges)), digits=4)
            q_str = join([@sprintf("%.3f", abs(c)) for c in charges], " + ")
            q_str *= " = $sum_q"
        else
            q_str = "N/A"
        end

        # S(q) diagnosis
        pk = sf_peaks[findfirst(p -> p.t == t_val, sf_peaks)]
        if pk.n_peaks_above_half <= 5 && abs(pk.peak_kx) < 0.3 && abs(pk.peak_ky) < 0.3
            sf_diag = "Γ-peak only → fluid"
        elseif pk.n_peaks_above_half > 5
            sf_diag = "multiple peaks → check"
        else
            sf_diag = "peak @ ($(pk.peak_kx),$(pk.peak_ky)) → CDW?"
        end

        @printf("  │ %8.2f │ %-19s │ %3d │ %-16s │ %-20s │\n",
            t_val, k_str, gsd, q_str, sf_diag)
    end
    println("  └──────────┴─────────────────────┴─────┴──────────────────┴──────────────────────┘")

    # ── FL vs Insulator caveat ──
    println("""

  ═══════════════════════════════════════════════════════════════════════
  DIAGNOSTIC LIMITATIONS
  ═══════════════════════════════════════════════════════════════════════

  What we CAN determine:
    • Topological order  → charge pump (fractional ΔQ → FCI)
    • Charge order       → S(q) (Bragg peaks → CDW/Wigner crystal)
    • Gap protection     → spectrum flow (protected under flux → topological)

  What we CANNOT determine with current tools:
    • Gaplessness        → need finite-size scaling of the many-body gap
    • Quasiparticle residue Z → need n(k) discontinuity
    • Drude weight / charge stiffness → not yet implemented

  If a phase is non-topological (charge pump = trivial), non-CDW (S(q) smooth),
  and the spectrum flow shows no protection, it is CONSISTENT with a Fermi
  liquid but CANNOT be rigorously distinguished from a trivial insulator
  without finite-size scaling or entanglement measures.

  ═══════════════════════════════════════════════════════════════════════
  All figures saved under: $FIGURE_DIR
  ═══════════════════════════════════════════════════════════════════════
  """)

    return gs_sectors_map, flow_results, pump_results, sf_maps, sf_peaks
end
