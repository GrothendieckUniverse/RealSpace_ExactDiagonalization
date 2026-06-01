MLStyle.@data Statistics begin
    Bosonic()
    Fermionic()
end

abstract type Second_Quantized_Model end

"""
Struct `Real_Space_Second_Quantized_Model <: Second_Quantized_Model` for Short-Range Real-Space Second-Quantized Models
---
- Fields:
    - `params::Dict`: model parameters, which can include hopping amplitudes, interaction strengths, flux values, etc. depending on the specific model we are considering. For example, for the Haldane honeycomb model, we can have `params = Dict("t" => 1.0, "t′" => 0.6, "t′′" => -0.58, "ϕ" => 0.4π, "V1" => 0.2, "V2" => 0.1)`
    - `lattice::TightBinding.Real_Space_Lattice`
    - `tb_model::TightBinding.Real_Space_TightBinding_Model`
    - `statistics::Statistics`: the statistics of the particle (or the model), i.e., `Bosonic()` or `Fermionic`
    - `bilinear_terms::Vector{Tuple{Int,Int,T}}`: the bilinear terms in the Hamiltonian, each of which is a tuple `(i, j, t)`, representing `t * a†_i a_j` where `a` denotes canonical creation/annihilation operators for either bosons or fermions.
    - `density_density_terms::Vector{Tuple{Int,Int,T}}`: the density-density interaction terms, each of which is a tuple `(i, j, v)`, representing `v * n_i * n_j`.
"""
mutable struct Real_Space_Second_Quantized_Model{T} <: Second_Quantized_Model
    params::Dict
    lattice::TightBinding.Real_Space_Lattice
    tb_model::TightBinding.Real_Space_TightBinding_Model
    statistics::Statistics
    bilinear_terms::Vector{Tuple{Int,Int,T}} # `b†_i b_j` or `c†_i c_j` depending on statistics
    density_density_terms::Vector{Tuple{Int,Int,T}} # n_i n_j
end


"""
    update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π)

In-place update of the bilinear hopping terms with twisted boundary phases.
Delegates to `TightBinding.generate_bilinear_terms` which applies the Peierls
phase `exp(i·2π·θ_d·w_d)` to every hopping crossing a periodic boundary.

The model's `bilinear_terms` field is replaced; `params[\"twisted_phase_over_2π\"]`
is updated.  Density-density terms and all other fields are unchanged.
"""
function update_second_quantized_model_with_twisted_phases!(
    second_quantized_model::Real_Space_Second_Quantized_Model;
    twisted_phases_over_2π::Vector{Float64},
)::Real_Space_Second_Quantized_Model
    lattice = second_quantized_model.lattice
    length(twisted_phases_over_2π) == lattice.dim ||
        error("twisted_phases_over_2π must have length $(lattice.dim).")

    second_quantized_model.bilinear_terms = TightBinding.generate_bilinear_terms(
        second_quantized_model.tb_model;
        twisted_phases_over_2π=twisted_phases_over_2π,
    )
    second_quantized_model.params["twisted_phase_over_2π"] = twisted_phases_over_2π
    return second_quantized_model
end
