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
and the characteristic diagnostic/scaling points (all entries are physical
`t''` values):

| working phase label | physical `t''` values | selected manifold |
|:--|:--|:--|
| AHC | `-0.50, -0.45` | three lowest eigenstates |
| FCI | `-0.30, -0.15, 0.00, 0.05, 0.10` | three lowest eigenstates |
| candidate CDW | `0.20, 0.30` | three lowest eigenstates |

The `CDW` key and output directory are retained for compatibility, but the
positive-side phase is not assumed to be a CDW. The structure-factor puzzles,
the rank-one versus FCI-projector comparison, and the conditions for a valid
charge pump are developed in
[`competing_phase_diagnostics.md`](competing_phase_diagnostics.md).
Diagnostic and charge-gap checkpoints live in parameter-specific
`x_<numerator>` directories.  Loading still validates the model, filling, and
flux; `--refresh true` rebuilds stale derived CSVs while resuming compatible
checkpoints, whereas `--overwrite true` also recomputes those checkpoints.

The default solver is explicit sparse-matrix ED for active geometries from 3x3
through 3x7. On Hyak, Hamiltonian matrix construction is distributed across
one-thread Julia workers. Sweeps and diagnostics currently stop at 3x6; 3x7 is
retained only for charge-gap scaling and requests 84 workers and 240 GiB per
job. This avoids the matrix-free implementation's sector-sized buffer per
thread, which exhausts memory on this geometry. CLI `--mode matrix` or
`--mode matrixfree` always overrides this policy.

## Local jobs

One sweep point produces the ranked all-sector spectrum and two versions of
the connected structure factor: the absolute ground state and the normalized
projector over the three selected FCI reference states. The latter is called
the **FCI-manifold projector structure factor**. It is invariant under any
permutation or unitary mixing within the same three-state subspace, including
a rearrangement of the FCI ground-state manifold, but not under a change of
the selected subspace itself. Each version has an
allowed-momentum grid, a dense 101x101 map, and peak metrics:

```bash
julia --project=. phase_exploration/bin/run_sweep_point.jl \
  --geometry 3x5 --x -1.0 --task all
```

The normalized metric is evaluated literally as
`max(abs(S(q))) / abs(mean(S(q)))`. A second, numerically robust
`max(abs(S(q))) / mean(abs(S(q)))` column is also saved, but is not substituted
for the requested plot.

Run one characteristic diagnostic point with the physical hopping explicitly:

```bash
julia --project=. phase_exploration/bin/run_diagnostic_point.jl \
  --phase FCI --geometry 3x5 --tpp -0.15
```

This performs a zero-flux all-sector scan, chooses the configured number of
globally lowest eigenstates (retaining both momentum sector and in-sector
level), and generates:

- absolute-ground-state and selected-manifold projector `S(q)` maps and
  metrics;
- all-momentum-sector spectrum flow over one flux quantum (21 points,
  i.e. 20 intervals, on the same grid as the charge pump);
- the manifold charge pump over one flux quantum;
- a spatial/orbital cut ES for the absolute ground state;
- the Li-Haldane/Regnault-Bernevig momentum-resolved particle ES of the
  selected low-energy manifold (`N_A=2` by default).

Individual diagnostics can be split across jobs, for example
`--observables structure,flow,pump` or `--observables spatial_es,pes`. Flux and zero-flux
checkpoints are shared and resumed sector by sector. Use `--refresh true` to
rebuild derived CSVs while retaining compatible checkpoints; `--overwrite
true` also recomputes the checkpoints.

The diagnostics renderer extracts `E4-E3` along the full stored flux path and
writes `manifold_gap_flow.svg`.  A pump plot is visibly marked with a warning
when the assumed three-state manifold touches outside states, because its
branch endpoints are then not a globally isolated-bundle invariant. Candidate-
CDW pump plots are restricted to the currently configured characteristic
points, so obsolete positive-side results such as `t''=0.21` are not rendered.

See [`entanglement_counting_notes.md`](entanglement_counting_notes.md) for the
`(1,3)` PES derivation, geometry-by-geometry counting, and an explanation of
what the current spatial ES can and cannot establish.

The sweep plotter also writes `zero_flux_gap_diagnostics_<geometry>.svg` as a
single-axis plot with only two global-energy differences, using
$E_0\leq E_1\leq E_2\leq E_3\leq\cdots$: $E_2-E_0$ is the three-state FCI
manifold width and $E_3-E_2$ is the roton gap.  No additional isolation or
same-sector-gap curves are mixed into this figure.

For the corrected `3x6` FCI run, the three manifold states are levels 1--3
of the same momentum sector `(0,3)`.  The momentum-resolved particle PES uses
the normalized projector over all three states and has exactly 117 levels
below its largest entanglement gap, with the expected 6/7 even/odd-$K_2$
sector counting.  The spatial-orbital spectrum instead uses the absolute
ground state `(0,3,1)` and is resolved only by subsystem particle number.
The construction, numerical checks, and interpretation are recorded in
[Section 5 of the counting notes](entanglement_counting_notes.md#5-corrected-3x6-fci-construction-and-numerical-audit).

One finite-size charge-gap datum is:

```bash
julia --project=. phase_exploration/bin/run_charge_gap_point.jl \
  --phase FCI --geometry 3x5 --tpp -0.15
```

As for diagnostics, use `--refresh true` to replace the result CSV while
retaining a compatible, parameter-specific checkpoint.  Reserve `--overwrite
true` for intentionally recomputing the ED sectors.

It scans every momentum sector at `N-1`, `N`, and `N+1`, then stores
`Delta_c = E0(N+1) + E0(N-1) - 2 E0(N)`. Run it for 3x3, 3x4, 3x5, 3x6, and
3x7 in each phase. The plotter combines all phases in one panel, connects each
phase's raw sizes, draws its independent linear extrapolation as a matching
dashed line, and records the `1/N_sites -> 0` intercept and RMS residual. The
axis includes the gapless `Delta_c=0` reference explicitly.

For the central `N0` calculation, a missing charge-gap checkpoint is seeded
automatically from a compatible sweep `zero_flux.jld2` when one exists. If both
the sweep and charge-gap checkpoints are partial, their completed momentum
sectors are merged before the scan continues. The merged cache is saved before
the next sector starts, then updated atomically after every completed sector;
an interrupted or OOM-killed job therefore resumes from its last finished
sector. The `N0-1` and `N0+1` particle-number sectors retain their own
charge-gap checkpoints and resume in the same way.

After any subset of data exists, render all available figures without further
diagonalization:

```bash
julia --project=. phase_exploration/bin/plot_results.jl --kind all
```

Other plot kinds are `sweep`, `ed-spectra`, `structure`, `diagnostics`, and
`charge-gap`. The `ed-spectra` renderer loops over every available sweep point
and recreates its zero-twist, symmetry-resolved ED spectrum in the same style
as the package's `plot_spectrum`, grouped by geometry. The sweep renderer
creates separate spectrum, `max|S|`, normalized `max|S/mean(S)|`, and peak
wavevector figures for every geometry plus multi-geometry panels. The two
metric curves compare the absolute ground state with the FCI-manifold
projector structure factor; dotted guides mark changes of the rank-one ground
state. The spectrum sweeps plot the
lowest 20 global levels by default and use one
default circular marker, with lines connecting the same momentum-sector and
in-sector-level state across adjacent sweep points. All states are dark gray
except the three fixed FCI reference states: level 1 in three sectors on 3x4
and 3x5, and levels 1--3 of `(0,3)` on 3x6. Those three states remain plotted
even when competing levels push them above the usual rank cutoff. This makes
gray roton levels that cross into the FCI manifold directly visible. Each sweep point gets a
two-panel finite-grid/dense-grid 2D structure-factor map. Diagnostic charge-pump
and spectrum-flow legends identify their momentum sectors; the accompanying
manifold-gap plot tests the pump's spectral-isolation prerequisite. Spectrum-flow
plots follow the twenty states that are lowest at zero flux (configurable with
`--max-flow-curves`). Every state in the focused manifold uses a black connecting
line and colored scatter points; every other curve is dark gray. Thus the 3x4
and 3x5 FCI manifolds highlight level 1 in three different sectors, whereas the
3x6 FCI manifold highlights levels 1--3 in its common sector. Spectrum-flow
markers are enlarged for visibility. If the focused trajectories remain directly
separated from every other plotted level throughout
the flux cycle, the vertical scale is 1.5 times their combined energy span plus
the minimum direct flow gap. If that separation closes, the scale is twice the
focused span. The lower limit is minus 5% of this scale so zero-energy branches
remain clearly visible. The complete
all-sector flow remains available in the CSV.

## Output layout

```text
phase_exploration/
  results/
    sweep/<L1xL2>/x_<numerator>/
    diagnostics/<AHC|FCI|CDW>/<L1xL2>/tpp_<physical-value>/
    charge_gap/<AHC|FCI|CDW>/<L1xL2>/tpp_<physical-value>/
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

This writes one `.sbatch` file per `(geometry, t'')` sweep point, plus diagnostic
and charge-gap jobs for all eight configured AHC, FCI, and candidate-CDW
characteristic points. The accompanying `data_jobs.txt` is the exact submission
manifest. Review the collection, then queue everything with:

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
Slurm logs. Every stdout/stderr filename includes the explicit Slurm job ID,
for example `tpp_dx7_3x5_fci_m0p15_slurm-12345678.out` and its matching
`.err` file.

After an HPC-side failure, inspect the job state and corresponding logs with:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem
ls phase_exploration/hpc/logs/
```

The generated diagnostic jobs use a versioned protocol marker and
parameter-specific result directory. This forces the new ground-versus-manifold
structure output, corrected `(sector, level)` selection, and four-level-per-
sector flow data to run once; subsequent submissions skip a fully completed
point normally.
