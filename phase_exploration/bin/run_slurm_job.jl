#!/usr/bin/env julia

using Distributed
using Pkg
using SlurmClusterManager

isempty(ARGS) && error("Usage: run_slurm_job.jl TARGET_SCRIPT [TARGET_ARGS...]")
target_script = abspath(popfirst!(ARGS))
isfile(target_script) || error("Target Julia script does not exist: $target_script")

mode_index = findfirst(==("--mode"), ARGS)
mode_index === nothing && error("Target arguments must include `--mode matrix` or `--mode matrixfree`.")
mode_index == length(ARGS) && error("Missing value after `--mode`.")
solver_mode = Symbol(lowercase(ARGS[mode_index + 1]))
solver_mode in (:matrix, :matrixfree) || error("Unknown solver mode: $(ARGS[mode_index + 1])")

launcher_project = get(ENV, "JULIA_PROJECT", "")
depot_dir = get(ENV, "JULIA_DEPOT_PATH", "")
repo_dir = get(ENV, "PHASE_EXPLORATION_REPO", "")
isempty(launcher_project) && error("JULIA_PROJECT is not set.")
isempty(depot_dir) && error("JULIA_DEPOT_PATH is not set.")
isempty(repo_dir) && error("PHASE_EXPLORATION_REPO is not set.")
repo_dir = realpath(repo_dir)

# SlurmClusterManager is loaded from the shared launcher environment above.
# Switch the master to the repository project before loading any ED or phase
# code, and launch every worker with this same repository project.
Pkg.activate(repo_dir)
ed_project = dirname(Base.active_project())
ENV["JULIA_PROJECT"] = ed_project

ENV["SLURM_EXPORT_ENV"] = "ALL"
ENV["SRUN_EXPORT_ENV"] = "ALL"

@info "Slurm Julia launcher" job_id=get(ENV, "SLURM_JOB_ID", missing) ntasks=get(ENV, "SLURM_NTASKS", missing) cpus_per_task=get(ENV, "SLURM_CPUS_PER_TASK", missing) julia_threads=Threads.nthreads() solver_mode launcher_project ed_project depot_dir target_script

if solver_mode == :matrix
    Threads.nthreads() == 1 || error("Explicit distributed matrix mode requires JULIA_NUM_THREADS=1; got $(Threads.nthreads()).")
    addprocs(
        SlurmClusterManager.SlurmManager(;
            launch_timeout=900.0,
            srun_post_exit_sleep=2.0,
        );
        exeflags=["--project=$(ed_project)", "--startup-file=no"],
        env=[
            "JULIA_PROJECT" => ed_project,
            "JULIA_DEPOT_PATH" => depot_dir,
            "JULIA_NUM_THREADS" => "1",
            "SLURM_EXPORT_ENV" => "ALL",
            "SRUN_EXPORT_ENV" => "ALL",
        ],
    )

    @info "Launched workers for distributed matrix construction" nprocs=Distributed.nprocs() nworkers=Distributed.nworkers() workers=Distributed.workers()
    Distributed.nprocs() > 1 || error(
        "SlurmClusterManager launched no additional Julia processes for explicit matrix mode.")
    @everywhere using RealSpace_ExactDiagonalization, LinearAlgebra, Statistics
else
    # The matrix-free operator uses Threads.@threads and thread-local buffers.
    # Keep it in one process so all allocated CPUs contribute to H|psi>.
    # With only pid 1, Distributed reports `nworkers() == 1` and
    # `workers() == [1]`; `nprocs() == 1` is the correct single-process test.
    Distributed.nprocs() == 1 || error(
        "Matrix-free mode requires exactly one Julia process; got " *
        "nprocs=$(Distributed.nprocs()), workers=$(Distributed.workers()).")
    Threads.nthreads() > 1 || error("Matrix-free mode requires JULIA_NUM_THREADS > 1.")
    using RealSpace_ExactDiagonalization, LinearAlgebra, Statistics
    @info "Using threaded matrix-free solver" nprocs=Distributed.nprocs() julia_threads=Threads.nthreads()
end

# ARGS now contains only arguments intended for the target CLI script.
include(target_script)
