#!/usr/bin/env julia

using Pkg

isempty(ARGS) && error("Usage: run_slurm_job.jl TARGET_SCRIPT [TARGET_ARGS...]")
target_script = abspath(popfirst!(ARGS))
isfile(target_script) || error("Target Julia script does not exist: $target_script")
mode_index = findfirst(==("--mode"), ARGS)
mode_index === nothing && error("Target arguments must include --mode matrix or matrixfree.")
mode_index == length(ARGS) && error("Missing value after --mode.")
solver_mode = Symbol(lowercase(ARGS[mode_index + 1]))
solver_mode in (:matrix, :matrixfree) || error("Unknown solver mode: $solver_mode")

repo_dir = get(ENV, "PHASE_EXPLORATION_REPO", "")
isempty(repo_dir) && error("PHASE_EXPLORATION_REPO is not set.")
repo_dir = realpath(repo_dir)
ntasks = parse(Int, get(ENV, "SLURM_NTASKS", "1"))
cpus = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", string(Threads.nthreads())))
ntasks == 1 || error("The shared-memory ED solver needs one Slurm task; regenerate the job scripts.")
Threads.nthreads() == cpus || error(
    "JULIA_NUM_THREADS=$(Threads.nthreads()) must match SLURM_CPUS_PER_TASK=$cpus.")

# Both sparse matrix construction and the matrix-free operator use threads.
# Additional processes would duplicate the catalog and remain unused by ED.
Pkg.activate(repo_dir)
ENV["JULIA_PROJECT"] = dirname(Base.active_project())
using RealSpace_ExactDiagonalization
@info "Shared-memory Slurm Julia launcher" job_id=get(ENV, "SLURM_JOB_ID", missing) ntasks cpus julia_threads=Threads.nthreads() solver_mode repo_dir target_script

# ARGS now contains only arguments intended for the target CLI script.
include(target_script)
