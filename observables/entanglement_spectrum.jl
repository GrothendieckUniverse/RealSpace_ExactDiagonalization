# ═══════════════════════════════════════════════════════════════════════════
# entanglement_spectrum — spatial orbital bipartition diagnostics
#
# The implementation below expands a symmetry-sector eigenvector into Fock
# amplitudes and then forms Schmidt blocks at fixed particle number in region A.
# This is intentionally general and transparent; for the [3,5] checkerboard
# ν=1/3 study the full Fock support is C(30,5), which is still manageable.
# ═══════════════════════════════════════════════════════════════════════════

function _resolve_sector_eigenvector(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    target_eigval_idx::Int,
    filling_fraction::Rational{Int},
    ed_mode::Symbol,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing},
)
    n_site = model.lattice.n_site

    if ed_data === nothing
        flux0 = zeros(Float64, model.lattice.dim)
        update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux0)
        G = sector_label == :identity ?
            build_identity_group(n_site) :
            build_translation_group(model.lattice, flux0)
        ed_data = build_ed_data(model; filling_fraction=filling_fraction, symmetry_group=G)
        scanned = sector_label == :identity ? nothing : [Tuple(Int.(sector_label))]
        ed_scan!(ed_data; nev=max(target_eigval_idx, 2), mode=ed_mode,
            scanned_sectors=scanned)
    end

    label = sector_label == :identity ? :identity : Tuple(Int.(sector_label))
    irrep_idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
    irrep_idx === nothing && error("Sector $(repr(label)) not found in ed_data.")
    haskey(ed_data.ed_scan_res, irrep_idx) ||
        error("Sector $(repr(label)) was not scanned in ed_data.")

    _, vecs = ed_data.ed_scan_res[irrep_idx]
    target_eigval_idx <= size(vecs, 2) ||
        error("target_eigval_idx=$target_eigval_idx exceeds $(size(vecs, 2)) eigenvectors.")

    basis = build_symmetry_sector_basis(ed_data.orbit_catalog, ed_data.irrep_list[irrep_idx])
    return ed_data, basis, vecs[:, target_eigval_idx]
end

function _expand_sector_state_to_fock_amplitudes(
    c::Vector{ComplexF64},
    basis::Symmetry_Sector_Basis,
    particle_statistics::Particle_Statistics;
    atol::Float64=1e-12,
)::Dict{Mask,ComplexF64}
    G = basis.symmetry_group
    nG = group_order(G)
    amplitudes = Dict{Mask,ComplexF64}()

    @inbounds for (col, repr) in enumerate(basis.representative_mask_list)
        coeff = c[col]
        abs(coeff) <= atol && continue
        prefactor = coeff / sqrt(nG * basis.stabilizer_order_list[col])

        for (gidx, op) in enumerate(G.operations)
            shifted, phase = apply_operation_to_mask(repr, op, particle_statistics)
            amp = prefactor * conj(basis.irrep.values[gidx]) * phase
            amplitudes[shifted] = get(amplitudes, shifted, zero(ComplexF64)) + amp
        end
    end

    for (mask, amp) in collect(amplitudes)
        abs(amp) <= atol && delete!(amplitudes, mask)
    end

    norm2 = sum(abs2, values(amplitudes))
    norm2 > atol || error("Expanded state has near-zero norm.")
    invnorm = inv(sqrt(norm2))
    for mask in keys(amplitudes)
        amplitudes[mask] *= invnorm
    end
    return amplitudes
end

function _subsystem_maps(partition_a::AbstractVector{<:Integer}, n_site::Int)
    seen = falses(n_site)
    for site in partition_a
        1 <= site <= n_site || error("partition_a contains invalid site index $site.")
        seen[site] && error("partition_a contains duplicate site index $site.")
        seen[site] = true
    end

    partition_b = [site for site in 1:n_site if !seen[site]]
    a_pos = zeros(Int, n_site)
    b_pos = zeros(Int, n_site)
    for (idx, site) in enumerate(partition_a)
        a_pos[Int(site)] = idx
    end
    for (idx, site) in enumerate(partition_b)
        b_pos[site] = idx
    end
    return Int.(partition_a), partition_b, a_pos, b_pos
end

function _split_mask_for_partition(mask::Mask, a_pos::Vector{Int}, b_pos::Vector{Int})
    mask_a = zero(Mask)
    mask_b = zero(Mask)
    n_a = 0
    n_site = length(a_pos)

    tmp = mask
    @inbounds while tmp != 0
        lsb = tmp & -tmp
        site = trailing_zeros(lsb) + 1
        if a_pos[site] != 0
            mask_a |= Mask(1) << (a_pos[site] - 1)
            n_a += 1
        else
            mask_b |= Mask(1) << (b_pos[site] - 1)
        end
        tmp ⊻= lsb
    end

    return mask_a, mask_b, n_a
end

@inline _partition_sign(::Bosonic, ::Mask, ::AbstractVector{Bool}) = COMPLEX_ONE

function _partition_sign(::Fermionic, mask::Mask, is_a_site::AbstractVector{Bool})::ComplexF64
    # Tensor convention: all A creation operators precede all B creation
    # operators, with each subsystem internally ordered by the original site
    # order. The sign is the parity needed to move occupied A sites past
    # earlier occupied B sites.
    b_before = 0
    swaps = 0
    tmp = mask
    @inbounds while tmp != 0
        lsb = tmp & -tmp
        site = trailing_zeros(lsb) + 1
        if is_a_site[site]
            swaps += b_before
        else
            b_before += 1
        end
        tmp ⊻= lsb
    end
    return isodd(swaps) ? -COMPLEX_ONE : COMPLEX_ONE
end

function _schmidt_blocks_from_amplitudes(
    amplitudes::Dict{Mask,ComplexF64},
    partition_a::AbstractVector{<:Integer},
    n_site::Int,
    particle_statistics::Particle_Statistics,
)
    part_a, _, a_pos, b_pos = _subsystem_maps(partition_a, n_site)
    is_a_site = falses(n_site)
    is_a_site[part_a] .= true

    block_entries = Dict{Int,Vector{Tuple{Mask,Mask,ComplexF64}}}()
    for (mask, amp) in amplitudes
        mask_a, mask_b, n_a = _split_mask_for_partition(mask, a_pos, b_pos)
        signed_amp = _partition_sign(particle_statistics, mask, is_a_site) * amp
        push!(get!(block_entries, n_a, Tuple{Mask,Mask,ComplexF64}[]),
            (mask_a, mask_b, signed_amp))
    end
    return block_entries
end

function _entanglement_probabilities_by_block(block_entries)
    rows = NamedTuple[]
    probabilities = Float64[]

    for n_a in sort(collect(keys(block_entries)))
        entries = block_entries[n_a]
        a_masks = sort!(unique(first.(entries)))
        b_masks = sort!(unique(getindex.(entries, 2)))
        a_index = Dict(mask => idx for (idx, mask) in enumerate(a_masks))
        b_index = Dict(mask => idx for (idx, mask) in enumerate(b_masks))
        M = zeros(ComplexF64, length(a_masks), length(b_masks))

        for (mask_a, mask_b, amp) in entries
            M[a_index[mask_a], b_index[mask_b]] += amp
        end

        λ = real.(svdvals(M)).^2
        λ = sort(λ[λ .> 1e-14]; rev=true)
        append!(probabilities, λ)
        for (level, p) in enumerate(λ)
            push!(rows, (
                n_a=n_a,
                level=level,
                probability=p,
                entanglement_energy=-log(p),
                dim_a=length(a_masks),
                dim_b=length(b_masks),
            ))
        end
    end

    sort!(rows; by=row -> row.entanglement_energy)
    return rows, sort(probabilities; rev=true)
end

"""
    entanglement_spectrum(model, sector_label; partition_a, filling_fraction, ...)

Compute the spatial orbital entanglement spectrum for one many-body eigenstate.

`partition_a` is a vector of 1-based flattened site indices assigned to
subsystem A; subsystem B is its complement. The result is particle-number
resolved in A and returns both Schmidt probabilities and entanglement energies
`ξ = -log(λ)`.
"""
function entanglement_spectrum(
    model::Real_Space_Second_Quantized_Model,
    sector_label;
    partition_a::AbstractVector{<:Integer},
    target_eigval_idx::Int=1,
    filling_fraction::Rational{Int},
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
)::NamedTuple
    ed_data, basis, c = _resolve_sector_eigenvector(
        model, sector_label;
        target_eigval_idx=target_eigval_idx,
        filling_fraction=filling_fraction,
        ed_mode=ed_mode,
        ed_data=ed_data,
    )

    amplitudes = _expand_sector_state_to_fock_amplitudes(
        c, basis, model.particle_statistics)
    blocks = _schmidt_blocks_from_amplitudes(
        amplitudes, partition_a, model.lattice.n_site, model.particle_statistics)
    levels, probabilities = _entanglement_probabilities_by_block(blocks)

    return (;
        levels,
        probabilities,
        entanglement_energies=[row.entanglement_energy for row in levels],
        partition_a=Int.(partition_a),
        partition_b=[site for site in 1:model.lattice.n_site if !(site in partition_a)],
        sector_label=(sector_label == :identity ? :identity : Tuple(Int.(sector_label))),
        target_eigval_idx,
        filling_fraction,
        norm_probability=sum(probabilities),
        ed_data,
    )
end

"""
    plot_entanglement_spectrum(result; fig_path=nothing, title="Entanglement spectrum")

Scatter plot of `ξ=-log(λ)` grouped by the particle number in subsystem A.
"""
function plot_entanglement_spectrum(
    result;
    fig_path::Union{Nothing,String}=nothing,
    title::String="Entanglement spectrum",
)
    fig = Figure(size=(720, 520))
    ax = Axis(fig[1, 1]; xlabel="N_A", ylabel="ξ = -log(λ)", title=title)
    xs = [row.n_a for row in result.levels]
    ys = [row.entanglement_energy for row in result.levels]
    scatter!(ax, xs, ys; markersize=8, color=ys, colormap=:viridis)

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end

# ═══════════════════════════════════════════════════════════════════════════
# Momentum-resolved particle entanglement spectrum
# ═══════════════════════════════════════════════════════════════════════════

function _fixed_particle_masks(n_site::Int, n_filled::Int)
    masks = Mask[]
    if n_filled == 0
        push!(masks, zero(Mask))
        return masks
    end
    x = _first_combination_mask(n_filled)
    upper = one(Mask) << n_site
    while x < upper
        push!(masks, x)
        x = _gosper_next(x)
    end
    return masks
end

function _occupied_sites(mask::Mask)
    sites = Int[]
    tmp = mask
    while tmp != 0
        lsb = tmp & -tmp
        push!(sites, trailing_zeros(lsb) + 1)
        tmp ⊻= lsb
    end
    return sites
end

function _particle_submasks(mask::Mask, n_take::Int)
    sites = _occupied_sites(mask)
    out = Mask[]
    function rec(start::Int, left::Int, acc::Mask)
        if left == 0
            push!(out, acc)
            return
        end
        max_start = length(sites) - left + 1
        for idx in start:max_start
            rec(idx + 1, left - 1, acc | (one(Mask) << (sites[idx] - 1)))
        end
    end
    rec(1, n_take, zero(Mask))
    return out
end

@inline _particle_partition_sign(::Bosonic, ::Mask, ::Mask) = COMPLEX_ONE

function _particle_partition_sign(::Fermionic, full_mask::Mask, mask_a::Mask)::ComplexF64
    b_before = 0
    swaps = 0
    tmp = full_mask
    while tmp != 0
        lsb = tmp & -tmp
        if (mask_a & lsb) != 0
            swaps += b_before
        else
            b_before += 1
        end
        tmp ⊻= lsb
    end
    return isodd(swaps) ? -COMPLEX_ONE : COMPLEX_ONE
end

function _particle_partition_matrix(
    amplitudes::Dict{Mask,ComplexF64},
    n_site::Int,
    n_particles_a::Int,
    n_particles_b::Int,
    particle_statistics::Particle_Statistics,
    basis_a_index::Dict{Mask,Int},
    basis_b_index::Dict{Mask,Int},
)
    M = zeros(ComplexF64, length(basis_a_index), length(basis_b_index))
    for (full_mask, amp) in amplitudes
        count_ones(full_mask) == n_particles_a + n_particles_b || continue
        for mask_a in _particle_submasks(full_mask, n_particles_a)
            mask_b = full_mask ⊻ mask_a
            ia = basis_a_index[mask_a]
            ib = basis_b_index[mask_b]
            M[ia, ib] += _particle_partition_sign(
                particle_statistics, full_mask, mask_a) * amp
        end
    end
    return M
end

function _resolve_many_body_manifold_states(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    filling_fraction::Rational{Int},
    target_eigval_idx::Int,
    manifold_states,
    ed_mode::Symbol,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing},
)
    explicit_states = if manifold_states === nothing
        [(sector=Tuple(Int.(label)), level=target_eigval_idx) for label in sector_labels]
    else
        manifold_states
    end
    labels, state_specs, _, _ =
        _normalize_manifold_state_specs(sector_labels, 1, explicit_states)
    states = Dict{Mask,ComplexF64}[]
    energies = Float64[]
    for state in state_specs
        label, level = state.sector, state.level
        ed_data, basis, c = _resolve_sector_eigenvector(
            model, label;
            target_eigval_idx=level,
            filling_fraction=filling_fraction,
            ed_mode=ed_mode,
            ed_data=ed_data,
        )
        push!(states, _expand_sector_state_to_fock_amplitudes(
            c, basis, model.particle_statistics))
        irrep_idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
        vals, _ = ed_data.ed_scan_res[irrep_idx]
        push!(energies, Float64(vals[level]))
    end
    return labels, state_specs, states, energies, ed_data
end

function _sector_embedding_matrix(
    full_basis_index::Dict{Mask,Int},
    sector_basis::Symmetry_Sector_Basis,
    particle_statistics::Particle_Statistics,
)
    dim_full = length(full_basis_index)
    dim_sector = length(sector_basis.representative_mask_list)
    U = zeros(ComplexF64, dim_full, dim_sector)
    for col in 1:dim_sector
        c = zeros(ComplexF64, dim_sector)
        c[col] = 1
        amps = _expand_sector_state_to_fock_amplitudes(c, sector_basis, particle_statistics)
        for (mask, amp) in amps
            U[full_basis_index[mask], col] = amp
        end
    end
    return U
end

"""
    particle_entanglement_spectrum(model, sector_labels; n_particles_a, filling_fraction, ...)

Compute the translation-symmetry-resolved particle entanglement spectrum (PES)
for a selected low-energy manifold.  The reduced density matrix is built by
tracing out `N_B=N-N_A` particles and is then block-diagonalized in the
many-body momentum sectors of subsystem `A`.
"""
function particle_entanglement_spectrum(
    model::Real_Space_Second_Quantized_Model,
    sector_labels;
    n_particles_a::Int,
    filling_fraction::Rational{Int},
    target_eigval_idx::Int=1,
    manifold_states=nothing,
    ed_mode::Symbol=:matrix,
    ed_data::Union{Symmetry_Resolved_ED_Data,Nothing}=nothing,
    probability_cutoff::Float64=1e-14,
)::NamedTuple
    n_site = model.lattice.n_site
    n_particles_total = Int(n_site * filling_fraction)
    0 <= n_particles_a <= n_particles_total ||
        error("n_particles_a must be between 0 and total particle number $n_particles_total.")
    n_particles_b = n_particles_total - n_particles_a

    labels, state_specs, states, energies, ed_data = _resolve_many_body_manifold_states(
        model, sector_labels;
        filling_fraction=filling_fraction,
        target_eigval_idx=target_eigval_idx,
        manifold_states=manifold_states,
        ed_mode=ed_mode,
        ed_data=ed_data,
    )

    basis_a = _fixed_particle_masks(n_site, n_particles_a)
    basis_b = _fixed_particle_masks(n_site, n_particles_b)
    index_a = Dict(mask => idx for (idx, mask) in enumerate(basis_a))
    index_b = Dict(mask => idx for (idx, mask) in enumerate(basis_b))

    ρ = zeros(ComplexF64, length(basis_a), length(basis_a))
    for amps in states
        M = _particle_partition_matrix(
            amps, n_site, n_particles_a, n_particles_b,
            model.particle_statistics, index_a, index_b)
        mul!(ρ, M, M', 1 / length(states), 1)
    end
    trρ = real(tr(ρ))
    trρ > probability_cutoff || error("Particle reduced density matrix has near-zero trace.")
    ρ ./= trρ

    G = build_translation_group(model.lattice)
    catalog_a = build_symmetry_orbit_catalog(;
        second_quantized_model=model,
        n_filled=n_particles_a,
        symmetry_group=G,
        particle_statistics=model.particle_statistics,
    )
    irrep_list = build_translation_irrep_list(G, model.lattice)

    levels = NamedTuple[]
    sector_summaries = NamedTuple[]
    for irrep in irrep_list
        sector_basis = build_symmetry_sector_basis(catalog_a, irrep)
        dim_sector = length(sector_basis.representative_mask_list)
        dim_sector == 0 && continue
        U = _sector_embedding_matrix(index_a, sector_basis, model.particle_statistics)
        ρk = Hermitian(U' * ρ * U)
        vals = sort(real.(eigvals(ρk)); rev=true)
        vals = vals[vals .> probability_cutoff]
        for (level, p) in enumerate(vals)
            push!(levels, (
                momentum=Tuple(Int.(irrep.label)),
                level=level,
                probability=p,
                entanglement_energy=-log(p),
                sector_dim=dim_sector,
            ))
        end
        push!(sector_summaries, (
            momentum=Tuple(Int.(irrep.label)),
            sector_dim=dim_sector,
            kept_levels=length(vals),
            weight=sum(vals),
        ))
    end
    sort!(levels; by=row -> row.entanglement_energy)

    return (;
        levels,
        sector_summaries,
        n_particles_a,
        n_particles_b,
        n_particles_total,
        sector_labels=labels,
        manifold_states=state_specs,
        manifold_energies=energies,
        norm_probability=sum(row.probability for row in levels),
        filling_fraction,
        ed_data,
    )
end

function plot_particle_entanglement_spectrum(
    result;
    fig_path::Union{Nothing,String}=nothing,
    title::String="Momentum-resolved particle entanglement spectrum",
)
    momenta = sort(unique([row.momentum for row in result.levels]))
    momentum_index = Dict(k => idx for (idx, k) in enumerate(momenta))
    xs = [momentum_index[row.momentum] for row in result.levels]
    ys = [row.entanglement_energy for row in result.levels]

    fig = Figure(size=(900, 520))
    ax = Axis(fig[1, 1];
        xlabel="subsystem momentum sector",
        ylabel="ξ = -log(λ)",
        title=title,
        xticks=(1:length(momenta), [repr(k) for k in momenta]),
        xticklabelrotation=π / 3,
    )
    scatter!(ax, xs, ys; markersize=7, color=ys, colormap=:viridis)

    if fig_path !== nothing
        mkpath(dirname(fig_path))
        save(fig_path, fig)
    end
    return fig
end
