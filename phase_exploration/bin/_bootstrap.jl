# CLI scripts run in one Julia process, locally or in the threaded Slurm allocation.
using RealSpace_ExactDiagonalization

include(joinpath(@__DIR__, "..", "src", "CheckerboardPhaseStudy.jl"))
using .CheckerboardPhaseStudy
