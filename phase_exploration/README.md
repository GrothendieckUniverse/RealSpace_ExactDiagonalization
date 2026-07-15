# Checkerboard `t''` phase exploration

This directory is a modular, point-parallel study of the spinless fermionic
checkerboard Hubbard model at band filling `nu=1/3` (flattened graph filling
`1/6`). It covers the AHC, FCI, and CDW working regions under

```text
t'' = x / (2 + 2 sqrt(2)).
```

Every generated CSV stores both `tpp_numerator=x` and `tpp_actual=t''`.
Crucially, **every sweep plot uses the physical `t''` as its horizontal
coordinate**; plot annotations show only this physical value, rounded to two
decimal places. The numerator is used only to define the grid and make
readable filenames.

## Central configuration

Edit [`src/config.jl`](src/config.jl) before a production campaign. It contains
the model parameters, geometries, 0.1-spaced numerator sweep, solver policy,
and the three common deep-phase diagnostic/scaling points:

| phase | numerator `x` | physical `t''` | default manifold |
|:--|--:|--:|--:|
| AHC | `-3.0` | `-0.6213203436` | three lowest eigenstates |
| FCI | `-1.0` | `-0.2071067812` | three lowest eigenstates |
| CDW | `1.5` | `0.3106601718` | three lowest eigenstates |

These defaults encode the archived small-cluster evidence; they are not hard
claims about the larger sizes. If the new sweep moves a large-gap plateau,
change the three values once in `PHASE_SPECS` and regenerate diagnostics/gaps
with `--overwrite true` (or remove the corresponding old result/checkpoint
directories). Checkpoint loading validates the model, filling, and flux and
will refuse incompatible cached states.

The default solver is explicit sparse-matrix ED for active geometries from 3x3
through 3x6. On Hyak, Hamiltonian matrix construction is distributed across
one-thread Julia workers. Sweeps and diagnostics currently stop at 3x6; 3x7 is
retained only for charge-gap scaling and uses one multithreaded matrix-free
process to avoid constructing the much larger `N0+1` sparse blocks. CLI
`--mode matrix` or `--mode matrixfree` always overrides this policy.

## Local jobs

One sweep point produces the ranked all-sector spectrum, allowed-momentum
`S(q)`, a dense 101x101 `S(q)` grid, and the two requested metrics:

```bash
julia --project=. phase_exploration/bin/run_sweep_point.jl \
  --geometry 3x5 --x -1.0 --task all
```

The normalized metric is evaluated literally as
`max(abs(S(q))) / abs(mean(S(q)))`. A second, numerically robust
`max(abs(S(q))) / mean(abs(S(q)))` column is also saved, but is not substituted
for the requested plot.

Run a common deep-phase diagnostic point with:

```bash
julia --project=. phase_exploration/bin/run_diagnostic_point.jl \
  --phase FCI --geometry 3x5
```

This performs a zero-flux all-sector scan, chooses the configured number of
globally lowest eigenstates (retaining both momentum sector and in-sector
level), and generates:

- all-momentum-sector spectrum flow over three flux quanta (49 points);
- the manifold charge pump over one flux quantum;
- a spatial/orbital cut ES for the absolute ground state;
- the Li-Haldane/Regnault-Bernevig momentum-resolved particle ES of the
  selected low-energy manifold (`N_A=2` by default).

Individual diagnostics can be split across jobs, for example
`--observables flow,pump` or `--observables spatial_es,pes`. Flux and zero-flux
checkpoints are shared and resumed sector by sector. Use `--refresh true` to
rebuild derived CSVs while retaining compatible checkpoints; `--overwrite
true` also recomputes the checkpoints.

See [`entanglement_counting_notes.md`](entanglement_counting_notes.md) for the
`(1,3)` PES derivation, geometry-by-geometry counting, and an explanation of
what the current spatial ES can and cannot establish.

For the corrected `3x6` FCI run, the three manifold states are levels 1--3
of the same momentum sector `(0,3)`.  The momentum-resolved particle PES uses
the equal-weight projector over all three states and has exactly 117 levels
below its largest entanglement gap, with the expected 6/7 even/odd-$K_2$
sector counting.  The spatial-orbital spectrum instead uses the absolute
ground state `(0,3,1)` and is resolved only by subsystem particle number.
The construction, numerical checks, and interpretation are recorded in
[Section 5 of the counting notes](entanglement_counting_notes.md#5-corrected-3x6-fci-construction-and-numerical-audit).

One finite-size charge-gap datum is:

```bash
julia --project=. phase_exploration/bin/run_charge_gap_point.jl \
  --phase FCI --geometry 3x5
```

It scans every momentum sector at `N-1`, `N`, and `N+1`, then stores
`Delta_c = E0(N+1) + E0(N-1) - 2 E0(N)`. Run it for 3x3, 3x4, 3x5, 3x6, and
3x7 in each phase. The plotter combines all phases in one panel, connects each
phase's raw sizes, draws its independent linear extrapolation as a matching
dashed line, and records the `1/N_sites -> 0` intercept and RMS residual. The
axis includes the gapless `Delta_c=0` reference explicitly.

After any subset of data exists, render all available figures without further
diagonalization:

```bash
julia --project=. phase_exploration/bin/plot_results.jl --kind all
```

Other plot kinds are `sweep`, `ed-spectra`, `structure`, `diagnostics`, and
`charge-gap`. The `ed-spectra` renderer loops over every available sweep point
and recreates its zero-twist, symmetry-resolved ED spectrum in the same style
as the package's `plot_spectrum`, grouped by geometry. The sweep renderer
creates separate spectrum-rank, `max|S|`, and normalized `max|S/mean(S)|`
figures for every geometry plus multi-geometry panels. Each sweep point gets a
two-panel finite-grid/dense-grid 2D structure-factor map. Diagnostic charge-pump
and spectrum-flow legends identify their momentum sectors. Spectrum-flow plots
follow the ten states that are lowest at zero flux (configurable with
`--max-flow-curves`, up to ten), using ten distinct colors and an external
legend so neither colors nor labels obscure the levels. The complete all-sector
flow remains available in the CSV.

## Output layout

```text
phase_exploration/
  results/
    sweep/<L1xL2>/x_<numerator>/
    diagnostics/<AHC|FCI|CDW>/<L1xL2>/
    charge_gap/<AHC|FCI|CDW>/<L1xL2>/
  figures/
    sweep/
    ed_spectra/<L1xL2>/
    structure_factor/<L1xL2>/
    diagnostics/<phase>/<L1xL2>/
    charge_gap/
  checkpoints/                 # ignored; large and resumable
```

CSV data and SVG figures are intended for long-term retention. JLD2 files are
restart checkpoints and can be removed after a completed campaign.

## Hyak / Klone

Edit the Julia path, depot, shared project, allocation, partition, and resource
table near the top of [`hpc/hyak_slurm_gen.sh`](hpc/hyak_slurm_gen.sh), then
generate (but do not yet submit) the independent Slurm scripts. The shared
project defaults to `${JULIA_DEPOT}/environments/v1.12`, matching the standard
Klone environment. The repository path is inferred from the generator's
location, so generated jobs follow the actual checkout under `/gscratch` or
`/mmfs1/gscratch`. `REPO_DIR` and `JULIA_PROJECT_DIR` can both be overridden
through environment variables.

```bash
bash phase_exploration/hpc/hyak_slurm_gen.sh
```

This writes one `.sbatch` file per `(geometry, t'')` sweep point, one per
`(geometry, phase)` diagnostic, and one per `(geometry, phase)` charge-gap
point under `hpc/generated/`. The accompanying `data_jobs.txt` is the exact
submission manifest. Review the collection, then queue everything with:

```bash
phase_exploration/hpc/generated/submit_all.sh
```

The helper first submits one revision-specific environment job. The shared
`v1.12` environment is used only to load `SlurmClusterManager`; the setup job
then activates the repository project, instantiates and precompiles its exact
manifest, and loads the complete phase-study module as a preflight. This split
is important: packages such as `JLD2` are direct dependencies of the repository
project but need not be directly loadable from the shared launcher environment.
Every newly queued data job receives an `afterok` dependency on setup. A broken
or incomplete environment therefore stops at one setup job instead of
producing the same package-load error in every data job.

`TightBinding` is a local path dependency and must normally be checked out as a
sibling of this repository:

```text
repo/
  RealSpace_ExactDiagonalization/
  TightBinding/
```

The tracked manifest uses the portable relative path `../TightBinding`. For a
different checkout layout, export `TIGHTBINDING_DIR=/absolute/path/TightBinding`
before running the generator. Setup validates this path before doing any long
precompile and re-develops it relative to the repository project.

The helper then calls `sbatch` for every manifest entry in a loop. Since
`sbatch` returns immediately, all data jobs are submitted asynchronously and
can run independently after setup. Before each submission it skips an exact
job-name match already shown by `squeue`, and skips a completed point when all
required result files are present. Successful generated jobs write persistent
markers under `hpc/completed/`. The timestamped submission CSV records setup,
submitted, and skipped jobs. The plot job depends on setup, newly submitted
jobs, and matching jobs that were already active.

The batch shell invokes Julia directly from the shared launcher environment.
Inside Julia, the master activates the repository project and
`SlurmClusterManager.SlurmManager` launches every allocated task with that same
repository project and depot explicitly propagated. The ED toolbox is then
loaded on every worker. `_bootstrap.jl` contains no cluster-launch logic and
remains safe for ordinary local CLI runs. Generated jobs print their resolved
paths, resources, Julia version, and the failing shell line/command to the
Slurm logs.

After an HPC-side failure, inspect the job state and corresponding logs with:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem
ls phase_exploration/hpc/logs/
```

The generated diagnostic jobs use a versioned protocol marker. This forces the
corrected `(sector, level)` manifold selection and denser flow output to run
once even when older CSVs already exist; subsequent submissions skip the
versioned completed result normally.
