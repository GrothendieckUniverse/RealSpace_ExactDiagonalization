# Distributed workers, when present, are initialized by the Slurm entry point.
# Keeping that cluster-specific setup out of this bootstrap also lets every CLI
# script run normally in a single local Julia process.
using RealSpace_ExactDiagonalization

include(joinpath(@__DIR__, "..", "src", "CheckerboardPhaseStudy.jl"))
using .CheckerboardPhaseStudy
