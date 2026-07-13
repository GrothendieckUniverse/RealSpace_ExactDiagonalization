using Distributed

# Matrix-free Hamiltonian application uses Distributed.pmap. Workers launched
# with `julia -p N` must load the toolbox before receiving those closures.
@everywhere using RealSpace_ExactDiagonalization

include(joinpath(@__DIR__, "..", "src", "CheckerboardPhaseStudy.jl"))
using .CheckerboardPhaseStudy
