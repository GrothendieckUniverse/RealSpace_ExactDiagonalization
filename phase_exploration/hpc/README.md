# Symmetry-Resolved Flux Calculations on Hyak

The default diagnostic jobs run spectral flow and charge pump on the same
**17 points, including both endpoints, over one full flux quantum**:

```text
flux / (2π) = 0, 1/16, 2/16, ..., 15/16, 1
```

Flux is inserted along direction 1 (x); the pump measures polarization along
direction 2 (y). Every momentum sector is solved at every flux. The pump
reselects the globally lowest three states at each point and verifies the
isolation gap and projected-position invertibility. Flow and pump share the
same solver-v2 checkpoints, so they do not duplicate the flux ED calculation.
No two-dimensional Chern-number grid is scheduled.

## Generate and submit on HPC

Copy the updated Julia solver, observables, and phase-exploration code to the
HPC checkout, not just this directory. Keep `TightBinding` beside the checkout
(or set `TIGHTBINDING_DIR`). From that checkout's root:

```bash
bash phase_exploration/hpc/hyak_slurm_gen.sh
bash phase_exploration/hpc/generated/submit_all.sh --kind diagnostics
```

Generate scripts on HPC so their absolute paths point to the HPC checkout.
The generator accepts `REPO_DIR`, `JULIA_BIN`, `JULIA_DEPOT`,
`JULIA_PROJECT_DIR`, and `TIGHTBINDING_DIR` overrides. The Julia project
defaults to the repository itself. Review the allocation/partition and paths
near the top of `hyak_slurm_gen.sh` for your account.

This submits setup, **33 diagnostic jobs** (11 hopping values × three sizes),
and a dependent plotting job. It does not submit sweep or charge-gap jobs.
The default diagnostics refresh only flow/pump; existing structure-factor and
PES results are retained. To regenerate every diagnostic as well:

```bash
DIAGNOSTIC_OBSERVABLES=all bash phase_exploration/hpc/hyak_slurm_gen.sh
bash phase_exploration/hpc/generated/submit_all.sh --kind diagnostics
```

The submission helper also accepts `--kind sweep`, `--kind charge-gap`, or
`--kind all` (the default). It skips matching active jobs and only considers
a diagnostic complete when all required files **and its protocol marker**
exist. The marker and job name distinguish the point counts, flux ranges,
directions, and requested observable set. A marker for another protocol cannot satisfy the current run. Do not submit a new protocol concurrently with an older job writing
the same result directory; let the older job finish or cancel it first.

## Resources and restart behavior

The current ED engine uses shared-memory threads. Each data job starts
**one Julia process**, with `JULIA_NUM_THREADS` equal to Slurm's
`--cpus-per-task`. The launcher and setup use the repository project; no distributed-worker package is required.

| Geometry | Julia threads | Memory | Time limit |
|---|---:|---:|---:|
| 3×4 | 12 | 12 GiB | 1 hour |
| 3×5 | 24 | 36 GiB | 2 hours |
| 3×6 | 48 | 72 GiB | 4 hours |

Generated jobs use `--refresh true` to regenerate CSVs while resuming compatible
solver-v2 checkpoints. Older solver checkpoints are not reused. A failed or
timed-out calculation retains completed sector checkpoints; resubmitting
continues from those. `--overwrite true` is reserved for deliberately discarding
compatible computed sectors.

Outputs are `spectrum_flow.csv` and `charge_pump.csv`. The pump CSV records
selected `(k1,k2,level)` states, the manifold gap, and the smallest singular
value of the projected position operator at each flux. Both CSVs include a
schema version and direction fields; reuse checks their grid and coverage.

Before a flux observable is refreshed, its preceding CSV is moved into
`history/run_*/` within that parameter's result directory. This preserves
previous data and prevents a failed run from leaving stale data at the
current filename. ED checkpoints continue to resume independently.

## Reading failed pump jobs

If the global three-state gap closes or the projected position operator is
singular, the pump stops with an explicit error and the job is not marked
complete. The flow CSV is written first and remains available for examining
the spectrum. A different physical manifold may be needed; repeatedly
resubmitting the same calculation does not resolve a closed gap.

The automatic plotting job uses `afterok`; if a pump is rejected, it waits
rather than silently presenting the campaign as successful. Review failed
jobs and their manifold choice. You can regenerate available figures manually:

```bash
julia --project=. phase_exploration/bin/plot_results.jl --kind diagnostics
```

The plotter uses the current CSV names only. Earlier runs in `history/` are
not substituted for missing outputs. Check the protocol marker and Slurm log
when assessing whether a calculation completed.

The complete-cycle pump is the intended topology diagnostic. Keep its measured
branch charges and total rather than rounding them to a Chern number by
construction. A 17-point grid should be checked against a finer grid when
branch matching, gap minima, or the pump winding are ambiguous. A direct FHS
many-body Chern number is a separate two-dimensional calculation.

## Local validation

`python3 test/hpc_workflow.py` checks generated Bash syntax, the 33-job
campaign selection, resource/thread matching, and rejection of old completion
markers using mocked Slurm commands. It never submits real jobs.

`julia --project=. --threads=2 test/phase_flux_protocol.jl` runs a small
17-point all-sector calculation in a temporary directory and verifies CSV
coverage, pump diagnostics, refresh history, reuse, and figure rendering.
