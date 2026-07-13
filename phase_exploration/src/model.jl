function build_checkerboard_problem(sample::Tuple{Int,Int}, x::Real;
    n_particles::Union{Nothing,Int}=nothing,
    flux::Vector{Float64}=zeros(2),
)
    params = params_at_numerator(x)
    lattice_data = TightBinding.initialize_real_space_lattice(;
        lattice_name="checkerboard",
        sample_size=collect(sample),
        brav_vec_list=[[1.0, 0.0], [0.0, 1.0]],
        sub_crys_list=[[0.5, 0.0], [0.0, 0.5]],
        pbc_indicator=[true, true],
    )
    tb_model = TightBinding.initialize_real_space_tightbinding_model(
        lattice_data; model_name="checkerboard")

    t = params["t"]
    tp1 = params["t′_1"]
    tp2 = params["t′_2"]
    tpp = params["t′′"]
    phi = 2π * params["ϕ_over_2π"]

    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 0], 2)) => -t * exp(-im * phi); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 2)) => -t * exp(im * phi); is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 1)) => -t * exp(-im * phi); is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 1)) => -t * exp(im * phi); is_hermitian=true)

    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 0], 1)) => -tp1; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 1), ([0, 1], 1)) => -tp2; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 0], 2)) => -tp2; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([0, 1], 2)) => -tp1; is_hermitian=true)

    add_hopping_term!(tb_model, (([0, 0], 1), ([1, 1], 1)) => -tpp; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 0], 2), ([1, 1], 2)) => -tpp; is_hermitian=true)
    add_hopping_term!(tb_model, (([1, 0], 2), ([0, 1], 2)) => -tpp; is_hermitian=true)
    add_hopping_term!(tb_model, (([0, 1], 1), ([1, 0], 1)) => -tpp; is_hermitian=true)

    lattice = tb_model.lattice
    bilinear_terms = TightBinding.generate_bilinear_terms(
        tb_model; twisted_phases_over_2π=flux)

    V1 = params["λ"] * params["V1"]
    V2 = params["λ"] * params["V2"]
    V3 = params["λ"] * params["V3"]
    _, L2 = lattice.sample_size
    density_terms = Tuple{Int,Int,ComplexF64}[]
    for (i_from, (cell_from, sub_from)) in enumerate(lattice.site_list),
        (i_to, (cell_to, sub_to)) in enumerate(lattice.site_list)
        i_from == i_to && continue
        shift = mod.(cell_to .- cell_from, lattice.sample_size)
        if V1 != 0 && sub_from == 1 && sub_to == 2 &&
           (shift == [0, 0] || shift == [1, 0] || shift == [0, L2 - 1] || shift == [1, L2 - 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V1)))
        end
        if V2 != 0 && sub_from == sub_to && (shift == [1, 0] || shift == [0, 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V2)))
        end
        if V3 != 0 && sub_from == sub_to && (shift == [1, 1] || shift == [1, L2 - 1])
            push!(density_terms, (i_from, i_to, ComplexF64(V3)))
        end
    end

    model = Real_Space_Second_Quantized_Model(
        params, lattice, tb_model, Fermionic(), bilinear_terms, density_terms)
    lattice.twisted_phases_over_2π = copy(flux)
    np = isnothing(n_particles) ? default_particle_number(sample) : n_particles
    0 <= np <= lattice.n_site || error("Invalid particle number $np for $(lattice.n_site) sites.")
    filling = np // lattice.n_site
    group = build_translation_group(lattice, flux)
    ed_data = build_ed_data(model; filling_fraction=filling, symmetry_group=group)
    return model, ed_data, np, filling
end

function scan_with_resume!(ed_data;
    nev::Int,
    mode::Symbol,
    sectors=nothing,
    checkpoint_path::Union{Nothing,String}=nothing,
    overwrite::Bool=false,
)
    mode in (:matrix, :matrixfree) || error("Unknown ED mode $mode.")
    if checkpoint_path !== nothing && isfile(checkpoint_path) && !overwrite
        loaded = load_checkpoint(checkpoint_path)
        checkpoint_problem_matches(loaded, ed_data.second_quantized_model,
            ed_data.filling_fraction) || error(
            "Checkpoint `$checkpoint_path` belongs to a different model/filling. " *
            "Use --overwrite true or remove that checkpoint after changing config.jl.")
        ed_data = loaded
    elseif overwrite
        empty!(ed_data.ed_scan_res)
    end

    labels = sectors === nothing ? all_sector_labels(ed_data) : [Tuple(Int.(k)) for k in sectors]
    for label in labels
        idx = findfirst(irrep -> irrep.label == label, ed_data.irrep_list)
        idx === nothing && error("Sector $label is absent.")
        if haskey(ed_data.ed_scan_res, idx)
            values, vectors = ed_data.ed_scan_res[idx]
            length(values) >= nev && size(vectors, 2) >= nev && continue
            delete!(ed_data.ed_scan_res, idx)
        end
        ed_scan!(ed_data; nev=nev, mode=mode,
            use_distributed=(nprocs() > 1), scanned_sectors=[label])
        checkpoint_path !== nothing && save_checkpoint(ed_data, ensure_parent(checkpoint_path))
    end
    checkpoint_path !== nothing && save_checkpoint(ed_data, ensure_parent(checkpoint_path))
    return ed_data
end

function checkpoint_problem_matches(ed_data, model, filling; flux=nothing)
    loaded_model = ed_data.second_quantized_model
    Tuple(Int.(loaded_model.lattice.sample_size)) == Tuple(Int.(model.lattice.sample_size)) || return false
    ed_data.filling_fraction == filling || return false
    keys(loaded_model.params) == keys(model.params) || return false
    for key in keys(model.params)
        isapprox(Float64(loaded_model.params[key]), Float64(model.params[key]); atol=1e-13, rtol=1e-13) ||
            return false
    end
    if flux !== nothing
        loaded_flux = loaded_model.lattice.twisted_phases_over_2π
        length(loaded_flux) == length(flux) || return false
        all(isapprox.(loaded_flux, flux; atol=1e-13, rtol=1e-13)) || return false
    end
    return true
end

function ensure_flux_checkpoints(model, filling, flux_values;
    flux_direction::Int,
    nev::Int,
    mode::Symbol,
    checkpoint_dir::AbstractString,
    sectors=nothing,
    overwrite::Bool=false,
)
    mkpath(checkpoint_dir)
    paths = String[]
    for theta in flux_values
        flux = zeros(Float64, model.lattice.dim)
        flux[flux_direction] = Float64(theta)
        path = joinpath(checkpoint_dir, ed_scan_checkpoint_filename(model, flux, filling))
        if isfile(path) && !overwrite
            ed_theta = load_checkpoint(path)
            checkpoint_problem_matches(ed_theta, model, filling; flux=flux) || error(
                "Flux checkpoint `$path` does not match the requested model/flux. " *
                "Use --overwrite true after changing config.jl.")
        else
            update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π=flux)
            group = build_translation_group(model.lattice, flux)
            ed_theta = build_ed_data(model; filling_fraction=filling, symmetry_group=group)
        end
        ed_theta = scan_with_resume!(ed_theta; nev=nev, mode=mode,
            sectors=sectors, checkpoint_path=path, overwrite=overwrite)
        push!(paths, path)
    end
    update_second_quantized_model_with_twisted_phases!(
        model; twisted_phases_over_2π=zeros(Float64, model.lattice.dim))
    return paths
end
