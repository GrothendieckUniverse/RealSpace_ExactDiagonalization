MLStyle.@data Statistics begin
    Bosonic()
    Fermionic()
end

abstract type Second_Quantized_Model end
abstract type Real_Space_Second_Quantized_Model <: Second_Quantized_Model end


"""
Struct `ShortRange_Real_Space_Second_Quantized_Model <: Real_Space_Second_Quantized_Model <: Second_Quantized_Model` for Short-Range Real-Space Second-Quantized Models
---
- Fields:
    - `params::Dict`: model parameters, which can include hopping amplitudes, interaction strengths, flux values, etc. depending on the specific model we are considering. For example, for the Haldane honeycomb model, we can have `params = Dict("t" => 1.0, "t′" => 0.6, "t′′" => -0.58, "ϕ" => 0.4π, "V1" => 0.2, "V2" => 0.1)`
    - `lattice::TightBinding.Real_Space_Lattice`
    - `tb_model::TightBinding.Real_Space_TightBinding_Model`
    - `statistics::Statistics`: the statistics of the particle (or the model), i.e., `Bosonic()` or `Fermionic`
    - `bilinear_terms::Vector{Tuple{Int,Int,T}}`: the bilinear terms in the Hamiltonian, each of which is a tuple `(i, j, t)`, representing `t * a†_i a_j` where `a` denotes canonical creation/annihilation operators for either bosons or fermions.
    - `density_density_terms::Vector{Tuple{Int,Int,T}}`: the density-density interaction terms, each of which is a tuple `(i, j, v)`, representing `v * n_i * n_j`.
"""
struct ShortRange_Real_Space_Second_Quantized_Model{T} <: Real_Space_Second_Quantized_Model
    params::Dict
    lattice::TightBinding.Real_Space_Lattice
    tb_model::TightBinding.Real_Space_TightBinding_Model
    statistics::Statistics
    bilinear_terms::Vector{Tuple{Int,Int,T}} # `b†_i b_j` or `c†_i c_j` depending on statistics
    density_density_terms::Vector{Tuple{Int,Int,T}} # n_i n_j
end



