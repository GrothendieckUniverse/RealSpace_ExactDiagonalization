#!/usr/bin/env julia

using Distributed
using Pkg
using SlurmClusterManager

isempty(ARGS) && error("Usage: run_slurm_job.jl TARGET_SCRIPT [TARGET_ARGS...]")
target_script = abspath(popfirst!(ARGS))
isfile(target_script) || error("Target Julia script does not exist: $target_script")

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

@info "Slurm Julia launcher" job_id=get(ENV, "SLURM_JOB_ID", missing) ntasks=get(ENV, "SLURM_NTASKS", missing) launcher_project ed_project depot_dir target_script

addprocs(
    SlurmClusterManager.SlurmManager(;
        launch_timeout=900.0,
        srun_post_exit_sleep=2.0,
    );
    exeflags=["--project=$(ed_project)", "--startup-file=no"],
    env=[
        "JULIA_PROJECT" => ed_project,
        "JULIA_DEPOT_PATH" => depot_dir,
        "SLURM_EXPORT_ENV" => "ALL",
        "SRUN_EXPORT_ENV" => "ALL",
    ],
)

@info "Launched Julia workers" nworkers=Distributed.nworkers() workers=Distributed.workers()
Distributed.nworkers() > 0 || error("SlurmClusterManager launched no Julia workers.")

@everywhere using RealSpace_ExactDiagonalization, LinearAlgebra, Statistics

# ARGS now contains only arguments intended for the target CLI script.
include(target_script)
