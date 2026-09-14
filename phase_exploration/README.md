# Checkerboard `t''` phase exploration

This directory is a modular, point-parallel study of the spinless fermionic
checkerboard Hubbard model at band filling `nu=1/3` (flattened graph filling
`1/6`). It covers the AHC, FCI, and CDW working regions under

```text
t'' = x / (2 + 2 sqrt(2)).
```

Parameter-summary CSVs store both `tpp_numerator=x` and `tpp_actual=t''`.
Crucially, **every sweep plot uses the physical `t''` as its horizontal
coordinate**; plot annotations show only this physical value, rounded to two
decimal places. The numerator defines sweep result/checkpoint directories;
rendered ED spectra use the physical `tpp_<value>.svg` tag shared by
diagnostics.

## Central configuration

Edit [`src/config.jl`](src/config.jl) before a production campaign. It contains
the model parameters, geometries, regular 0.1-spaced numerator sweep, solver
policy, and the characteristic diagnostic/scaling points (all entries are
physical `t''` values). The exact numerator corresponding to every
characteristic point is automatically merged into the ED sweep, so a
diagnostic such as `t''=-0.15` always has a matching zero-flux spectrum:

| working phase label | physical `t''` values | selected manifold |
|:--|:--|:--|
| candidate AHC | `-0.60, -0.55, -0.50, -0.45` | three lowest eigenstates |
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
through 3x7. On Hyak, Hamiltonian matrix construction uses shared-memory threads in
one Julia process. Sweeps and diagnostics currently stop at 3x6; 3x7 is
retained only for charge-gap scaling and requests 84 threads and 240 GiB per
job. This avoids the matrix-free implementation's sector-sized buffer per
thread, which exhausts memory on this geometry. CLI `--mode matrix` or
`--mode matrixfree` always overrides this policy.

## Local jobs

One sweep point produces the ranked all-sector spectrum and the total-density
connected structure factor for both the absolute ground state and the
normalized projector over the three selected FCI reference states. The latter is called
the **FCI-manifold projector structure factor**. It is invariant under any
permutation or unitary mixing within the same three-state subspace, including
a rearrangement of the FCI ground-state manifold, but not under a change of
the selected subspace itself. Each version has an allowed-momentum grid, a
dense 101x101 map, and peak metrics. The sweep also stores allowed-momentum
and peak-metric files for the absolute-ground-state sublattice components
$S^{AA}(\mathbf q)$ and $\operatorname{Re}S^{AB}(\mathbf q)$; no sublattice
decomposition is applied to the manifold projector.

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
- all-momentum-sector spectrum flow over one flux quantum (17 points,
  i.e. 16 intervals, on the same grid as the charge pump);
- the manifold charge pump over one flux quantum;
- for FCI candidates only, the Li-Haldane/Regnault-Bernevig
  momentum-resolved particle ES of the selected low-energy manifold
  (`N_A=2` by default).

Individual diagnostics can be split across jobs, for example
`--observables structure,flow,pump` or, for an FCI point,
`--observables pes`. Spatial/orbital ES is not part of this campaign. Flux and
zero-flux checkpoints are shared and resumed sector by sector. Use `--refresh true` to
rebuild derived CSVs while retaining compatible checkpoints; `--overwrite
true` also recomputes the checkpoints.

The diagnostics renderer uses a constant zero-flux energy reference
`E - E0(0)` for spectral flow. It scans the cached spectrum at every flux to
choose the union of branches that enter the lowest displayed ranks, then
plots each complete `(k1, k2, level)` branch. The upper panel shows these
branches and the lower panel resolves the low-energy exchange.

The renderer reads `spectrum_flow.csv` and `charge_pump.csv` and writes
`spectrum_flow.svg`, `manifold_gap_flow.svg`, and `charge_pump.svg`. The gap
plot uses global energy ranks at each flux. Pump figures show polarization
eigenbranches and their sum; their annotations use the gap and projected-position
singular values saved in the pump CSV.

Flux CSVs carry a schema version and direction fields. The diagnostic driver
checks the stored grid and sector/manifold coverage before reusing them.
Refreshing flux data moves the preceding CSVs into a `history/run_*` directory,
so an interrupted calculation cannot expose an earlier result as current.
Compatible ED checkpoints remain reusable. Diagnostic and charge-gap plots
include the currently configured parameter points.

See [`entanglement_counting_notes.md`](entanglement_counting_notes.md) for the
`(1,3)` PES derivation and geometry-by-geometry counting. Its discussion of
spatial cuts is retained as background explaining why that observable was
removed from the production campaign.

The sweep plotter also writes `zero_flux_gap_diagnostics_<geometry>.svg` as a
single-axis plot with only two global-energy differences, using
$E_0\leq E_1\leq E_2\leq E_3\leq\cdots$: $E_2-E_0$ is the three-state FCI
manifold width and $E_3-E_2$ is the roton gap.  No additional isolation or
same-sector-gap curves are mixed into this figure.

For the corrected `3x6` FCI run, the three manifold states are levels 1--3
of the same momentum sector `(0,3)`.  The momentum-resolved particle PES uses
the normalized projector over all three states and has exactly 117 levels
below its largest entanglement gap, with the expected 6/7 even/odd-$K_2$
sector counting.
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
3x7 in each phase. The plotter writes one parameter-resolved finite-size
scaling figure under `figures/charge_gap/<phase>/tpp_<value>/`, connects the
raw sizes, draws an independent linear extrapolation, and records the
`1/N_sites -> 0` intercept and RMS residual. Every figure includes the gapless
`Delta_c=0` reference explicitly.

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
as the package's `plot_spectrum`, grouped by geometry and named by physical
`t''` so characteristic points match diagnostic directory names. The sweep renderer
creates separate spectrum, normalized `max|S/mean(S)|`, and peak-wavevector
figures for every geometry plus multi-geometry panels. The absolute-ground-
state maximum plot separates `max|S^{AA}|` and `max|Re S^{AB}|`; dedicated
combined figures overlay all geometries for each component. In
$t''\in[-0.35,-0.20]$, parameter-matched curves compare the global-ground-state
`max|S(q)|` directly against the FCI-ground-state-manifold-projector
`max|S(q)|`, both in one figure per geometry and in a three-panel 3x4/3x5/3x6
summary. The projector is not displayed outside this window, where the fixed
FCI slots need not be the relevant low-energy manifold. Dotted guides mark
changes of the rank-one ground state.
The spectrum sweeps plot the
lowest 20 global levels by default and use one
default circular marker, with lines connecting the same momentum-sector and
in-sector-level state across adjacent sweep points. All states are dark gray
except the three fixed FCI reference states: level 1 in three sectors on 3x4
and 3x5, and levels 1--3 of `(0,3)` on 3x6. Those three states remain plotted
even when competing levels push them above the usual rank cutoff. This makes
gray roton levels that cross into the FCI manifold directly visible. Each sweep point gets a
two-panel finite-grid/dense-grid 2D structure-factor map. Diagnostic spectral
flows use a fixed zero-flux energy reference and retain the union of branches
low anywhere on the path. Pump legends identify polarization eigenbranches;
fresh pump CSVs also record the actual selected states, isolation gap, and
projected-position singular values at every flux.

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
    charge_gap/<phase>/tpp_<physical-value>/finite_size_scaling.svg
  checkpoints/                 # ignored; large and resumable
```

CSV data and SVG figures are intended for long-term retention. JLD2 files are
restart checkpoints and can be removed after a completed campaign.

## Hyak / Klone

The [HPC run guide](hpc/README.md) describes the 17-point full-cycle
flow/pump campaign. Edit the Julia path, depot, allocation, partition, and
resource table near the top of [`hpc/hyak_slurm_gen.sh`](hpc/hyak_slurm_gen.sh).
The Julia project defaults to the repository. Generate scripts in the HPC
checkout so their paths match that machine:

```bash
bash phase_exploration/hpc/hyak_slurm_gen.sh
bash phase_exploration/hpc/generated/submit_all.sh --kind diagnostics
```

This submits the 33 characteristic flow/pump points only, plus environment
setup and a dependent plot job. The generator also creates all sweep and
charge-gap jobs; `submit_all.sh --kind all` submits the complete study.
Set `DIAGNOSTIC_OBSERVABLES=all` when generating to also refresh structure
factors and FCI particle entanglement spectra. The default preserves them.

The setup job activates the repository project, instantiates and precompiles
its dependencies, and loads the phase-study module. Each new data job depends
on successful setup. The launcher uses one Julia process and the allocated
number of threads; no distributed-worker package is required.

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

The batch shell invokes one Julia process. The launcher checks that Slurm
allocated one task and that Julia's thread count equals `SLURM_CPUS_PER_TASK`.
Generated jobs print resolved paths, resources, Julia version, and the failing
shell line/command to the Slurm logs. Each log filename includes its job ID.
Versioned diagnostic completion markers distinguish this protocol for each requested grid.
A rejected pump (closed gap or singular projected position) leaves its flow
available and does not create a successful completion marker.

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
