# ═══════════════════════════════════════════════════════════════════════════
# many_body_chern_number — flux-torus Chern number diagnostics
#
# Implements the non-Abelian Fukui-Hatsugai-Suzuki lattice Berry curvature for
# a selected low-energy manifold over the two twisted-boundary fluxes.
# ═══════════════════════════════════════════════════════════════════════════

function _normalize_flux_grid_size(flux_grid_size)
    if flux_grid_size isa Integer
        nx = Int(flux_grid_size)
        ny = Int(flux_grid_size)
    else
        nx = Int(flux_grid_size[1])
        ny = Int(flux_grid_size[2])
    end
    nx >= 2 && ny >= 2 || error("flux_grid_size must be at least 2 in both directions.")
    return nx, ny
end

function _chern_sector_labels(sector_labels)
    sector_labels == :identity && return [:identity]
    sector_labels isa Tuple && return [Tuple(Int.(sector_labels))]
    return [Tuple(Int.(label)) for label in sector_labels]
end

function _chern_params_digest(params)
    pairs = sort(collect(params); by=kv -> string(kv[1]))
    h = UInt64(0xcbf29ce484222325)
    for b in codeunits(repr(pairs))
        h = (h ⊻ UInt64(b)) * UInt64(0x00000100000001b3)
    end
    return string(h; base=16)
end

function _chern_flux_tag(flux)
    return join([replace(@sprintf("%.6f", θ), "." => "p", "-" => "m") for θ in flux], "_")
end

function _chern_checkpoint_path(model, filling_fraction, flux, checkpoint_dir::String)
    sample = join(model.lattice.sample_size, "x")
    digest = _chern_params_digest(model.params)
    name = "many_body_chern_$(model.tb_model.model_name)_$(sample)_nu=$(numerator(filling_fraction))_$(denominator(filling_fraction))_theta=$(_chern_flux_tag(flux))_p=$(digest).jld2"
    return joinpath(checkpoint_dir, name)
end

function _fock_amplitude_overlap(
    a::Dict{Mask,ComplexF64},
    b::Dict{Mask,ComplexF64},
)::ComplexF64
    acc = zero(ComplexF64)
    if length(a) <= length(b)
        for (mask, amp) in a
            acc += conj(amp) * get(b, mask, zero(ComplexF64))
        end
    else
        for (mask, amp) in b
            acc += conj(get(a, mask, zero(ComplexF64))) * amp
        end
    end
    return acc
end

function _subspace_overlap_matrix(states_a, states_b)
    n = length(states_a)
    length(states_b) == n || error("Chern link compares subspaces of different dimensions.")
    M = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n, j in 1:n
        M[i, j] = _fock_amplitude_overlap(states_a[i], states_b[j])
    end
    return M
end

function _unit_det_link(states_a, states_b; det_atol::Float64=1e-10)
    M = _subspace_overlap_matrix(states_a, states_b)
    z = det(M)
    if abs(z) <= det_atol
        error("Singular Berry link (|det|=$(abs(z))); refine the grid or select an isolated manifold.")
    end
    return z / abs(z), abs(z)
end

function _flux_point_states_for_chern(
    model::Real_Space_Second_Quantized_Model,
    labels;
    filling_fraction::Rational{Int},
    flux::Vector{Float64},
    nev_per_sector::Int,
    manifold_size::Int,
    mode::Symbol,
    checkpoint_dir::String,
    overwrite::Bool,
)
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux)
    active_group = labels == [:identity] ?
        build_identity_group(model.lattice.n_site) :
        build_translation_group(model.lattice, flux)
    ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=active_group)

    # Global low-energy selection must consider every momentum sector.
    scanned = nothing
    checkpoint_path = _chern_checkpoint_path(model, filling_fraction, flux,
        joinpath(checkpoint_dir,labels == [:identity] ? "solver_v2_identity" : "solver_v2_sectors"))
    ed_scan!(ed_data;
        nev=manifold_size + 1,
        mode=mode,
        scanned_sectors=scanned,
        checkpoint_path=checkpoint_path,
        overwrite=overwrite,
    )

    states = Dict{Mask,ComplexF64}[]
    energies = Float64[]
    sector_level_labels = Tuple[]

    candidates = sort([(energy=Float64(vals[level]), idx=idx, level=level)
        for (idx, (vals, _)) in ed_data.ed_scan_res for level in eachindex(vals)]; by=x->x.energy)
    length(candidates) > manifold_size || error("Need a state above the selected manifold to check its gap.")
    gap = candidates[manifold_size+1].energy - candidates[manifold_size].energy
    gap > 1e-9 || error("Low-energy manifold is not isolated at flux $flux (gap=$gap).")
    for state in candidates[1:manifold_size]
        irrep = ed_data.irrep_list[state.idx]
        basis = build_symmetry_sector_basis(ed_data.orbit_catalog, irrep)
        vecs = ed_data.ed_scan_res[state.idx][2]
        push!(states, _expand_sector_state_to_fock_amplitudes(
            vecs[:, state.level], basis, model.particle_statistics))
        push!(energies, state.energy)
        push!(sector_level_labels, (irrep.label, state.level))
    end

    return (; states, energies, sector_level_labels, checkpoint_path, gap)
end

function _plot_many_body_chern_curvature(result; fig_path::Union{Nothing,String}=nothing)
    fig = Figure(size=(650, 540))
    ax = Axis(fig[1, 1];
        xlabel="θ_x / 2π",
        ylabel="θ_y / 2π",
        title="Flux-torus Berry curvature, C=$(@sprintf("%.6f", result.chern_number))",
    )
    hm = heatmap!(ax, result.flux_x, result.flux_y, result.berry_curvature; colormap=:balance)
    Colorbar(fig[1, 2], hm; label="plaquette phase")

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end

"""
    many_body_chern_number(model, sector_labels; filling_fraction, flux_grid_size=(5, 5), ...)

Compute the Chern number of a selected many-body manifold on the two-dimensional
twisted-boundary flux torus.  `sector_labels` is the set of momentum sectors
used to infer the manifold dimension (`length(sector_labels)*nev_per_sector`).
At every flux all sectors are solved and the globally lowest `manifold_size`
states are selected. Use `:all` with explicit `manifold_size` for clarity.
Singular links and a closed isolation gap raise errors. Boundary-gauge Fock
states are periodic, so no extra sewing transformation is applied at the seam.

The algorithm uses the gauge-invariant non-Abelian lattice Berry curvature:
the U(1) link variable is the phase of the determinant of the overlap matrix
between neighboring multiplet subspaces.  The returned `chern_number` is
`sum(plaquette phases) / 2π`.
"""
function many_body_chern_number(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    flux_grid_size=(5, 5),
    nev_per_sector::Int=1,
    manifold_size::Union{Nothing,Int}=nothing,
    mode::Symbol=:matrix,
    checkpoint_dir::String="checkpoints",
    overwrite::Bool=false,
    fig_path::Union{Nothing,String}=nothing,
    det_atol::Float64=1e-10,
)
    model.lattice.dim == 2 || error("many_body_chern_number currently expects a 2D lattice.")
    nx, ny = _normalize_flux_grid_size(flux_grid_size)
    labels = sector_labels == :all ? [irrep.label for irrep in build_translation_irrep_list(build_translation_group(model.lattice),model.lattice)] : _chern_sector_labels(sector_labels)
    sector_labels == :all && manifold_size === nothing && error("Specify manifold_size with :all.")
    nstates = manifold_size === nothing ? length(labels)*nev_per_sector : manifold_size
    nstates > 0 || error("manifold_size must be positive.")

    flux_x = collect(0:(nx - 1)) ./ nx
    flux_y = collect(0:(ny - 1)) ./ ny
    states = Array{Vector{Dict{Mask,ComplexF64}},2}(undef, nx, ny)
    energies = Array{Vector{Float64},2}(undef, nx, ny)
    checkpoint_paths = Array{String,2}(undef, nx, ny)

    gaps = zeros(nx,ny)
    selected_states = Array{Vector{Tuple},2}(undef,nx,ny)
    mkpath(checkpoint_dir)
    for ix in 1:nx, iy in 1:ny
        flux = [Float64(flux_x[ix]), Float64(flux_y[iy])]
        @info "many-body Chern flux point" ix nx iy ny flux labels
        point = _flux_point_states_for_chern(
            model, labels;
            filling_fraction=filling_fraction,
            flux=flux,
            nev_per_sector=nev_per_sector,
            manifold_size=nstates,
            mode=mode,
            checkpoint_dir=checkpoint_dir,
            overwrite=overwrite,
        )
        gaps[ix, iy] = point.gap
        selected_states[ix, iy] = point.sector_level_labels
        states[ix, iy] = point.states
        energies[ix, iy] = point.energies
        checkpoint_paths[ix, iy] = point.checkpoint_path
    end

    Ux = Matrix{ComplexF64}(undef, nx, ny)
    Uy = Matrix{ComplexF64}(undef, nx, ny)
    det_abs_x = Matrix{Float64}(undef, nx, ny)
    det_abs_y = Matrix{Float64}(undef, nx, ny)
    for ix in 1:nx, iy in 1:ny
        ixp = mod1(ix + 1, nx)
        iyp = mod1(iy + 1, ny)
        Ux[ix, iy], det_abs_x[ix, iy] =
            _unit_det_link(states[ix, iy], states[ixp, iy]; det_atol=det_atol)
        Uy[ix, iy], det_abs_y[ix, iy] =
            _unit_det_link(states[ix, iy], states[ix, iyp]; det_atol=det_atol)
    end

    berry_curvature = Matrix{Float64}(undef, nx, ny)
    for ix in 1:nx, iy in 1:ny
        ixp = mod1(ix + 1, nx)
        iyp = mod1(iy + 1, ny)
        berry_curvature[ix, iy] = angle(
            Ux[ix, iy] * Uy[ixp, iy] * conj(Ux[ix, iyp]) * conj(Uy[ix, iy]))
    end

    chern_number = sum(berry_curvature) / (2π)
    res = (;
        chern_number,
        manifold_size=nstates,
        gaps, min_gap=minimum(gaps), selected_states,
        rounded_chern=round(Int, chern_number),
        berry_curvature,
        flux_x,
        flux_y,
        sector_labels=labels,
        nev_per_sector,
        flux_grid_size=(nx, ny),
        energies,
        min_link_det=min(minimum(det_abs_x), minimum(det_abs_y)),
        link_det_abs_x=det_abs_x,
        link_det_abs_y=det_abs_y,
        checkpoint_paths,
        fig_path,
    )

    fig_path !== nothing && _plot_many_body_chern_curvature(res; fig_path=fig_path)
    return res
end
