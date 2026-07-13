#!/usr/bin/env julia

using Pkg

repo_dir = get(ENV, "PHASE_EXPLORATION_REPO", "")
isempty(repo_dir) && error("PHASE_EXPLORATION_REPO is not set.")
repo_dir = realpath(repo_dir)
isfile(joinpath(repo_dir, "Project.toml")) ||
    error("No Project.toml found under PHASE_EXPLORATION_REPO=$repo_dir")

launcher_project = Base.active_project()
@info "Preparing shared HPC launcher environment" launcher_project repo_dir DEPOT_PATH

# The shared environment owns SlurmClusterManager only. ED and phase-study
# dependencies must resolve from the repository project, where JLD2 and the
# remaining direct dependencies are declared.
Pkg.instantiate()
Pkg.precompile()
using SlurmClusterManager

Pkg.activate(repo_dir)
Pkg.instantiate()
Pkg.precompile()
using RealSpace_ExactDiagonalization

# Load the complete phase-study module now so a missing direct dependency fails
# in this single setup job rather than in every data job.
include(joinpath(repo_dir, "phase_exploration", "src", "CheckerboardPhaseStudy.jl"))
using .CheckerboardPhaseStudy

package_source = pathof(RealSpace_ExactDiagonalization)
package_root = realpath(joinpath(dirname(package_source), ".."))
package_root == repo_dir ||
    error("Shared environment resolved RealSpace_ExactDiagonalization at $package_root, expected $repo_dir")

@info "HPC Julia environments ready" launcher_project ed_project=Base.active_project() package_source slurm_manager_source=pathof(SlurmClusterManager)
