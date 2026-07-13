#!/usr/bin/env julia

using Pkg

repo_dir = get(ENV, "PHASE_EXPLORATION_REPO", "")
isempty(repo_dir) && error("PHASE_EXPLORATION_REPO is not set.")
repo_dir = realpath(repo_dir)
isfile(joinpath(repo_dir, "Project.toml")) ||
    error("No Project.toml found under PHASE_EXPLORATION_REPO=$repo_dir")

@info "Preparing shared HPC Julia environment" active_project=Base.active_project() repo_dir DEPOT_PATH

# Ensure the shared environment loads this exact checkout rather than a stale
# registered copy or an older development path. This operation is idempotent.
Pkg.develop(Pkg.PackageSpec(path=repo_dir))
Pkg.instantiate()
Pkg.precompile()

using Distributed
using RealSpace_ExactDiagonalization
using SlurmClusterManager

package_source = pathof(RealSpace_ExactDiagonalization)
package_root = realpath(joinpath(dirname(package_source), ".."))
package_root == repo_dir ||
    error("Shared environment resolved RealSpace_ExactDiagonalization at $package_root, expected $repo_dir")

@info "HPC Julia environment ready" active_project=Base.active_project() package_source slurm_manager_source=pathof(SlurmClusterManager)
