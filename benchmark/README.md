# Stored ED benchmark scans

The figures use the latest complete scan available in the Python companion's
`benchmark/` directory when imported on 2026-09-14. The source file modification
date is 2026-08-20; these are recorded timings, not a new timing campaign for
the current working tree.

- Julia: `benchmark_data/benchmark_raw_latest.csv` (36 rows), identical to
  `benchmark_raw.csv`. The companion's `make_readme_figures.py` uses that same
  Julia scan as input; no Python times have been relabelled as Julia times.
- Python: `benchmark_data/benchmark_raw_py_scan.csv` (36 rows), copied verbatim
  from `../RealSpace_ExactDiagonalization_PY/benchmark/benchmark_data/`.
  The companion's `benchmark_raw_latest_py.csv` is byte-identical to this scan.

Each scan has three models, six sizes per model, and matrix/matrix-free modes.
The drivers time one momentum sector with `nev=1`; the orbit catalog is built
before the timed section. Their default three repetitions discard the first
and average the remaining two. The CSV schema does not record machine, thread
count, or solver revision, so it cannot independently establish those settings
or a portable cross-language speedup. This preserves the available results
without inventing missing run metadata.

All 36 matched cases have equal sector dimensions. Their maximum ground-energy
difference is 2.13e-13.

| Largest stored sample | Sector dimension | Julia matrix (s) | Julia matrix-free (s) | Python matrix (s) | Python matrix-free (s) |
|---|---:|---:|---:|---:|---:|
| Heisenberg N=28 | 1,432,860 | 27.501 | 36.569 | 27.138 | 29.154 |
| Haldane_Boson 4×4 | 657,756 | 50.392 | 59.688 | 66.654 | 59.085 |
| Hubbard_Fermion 2×7 | 2,865,228 | 96.814 | 136.996 | 104.222 | 105.806 |

Regenerate the six Julia timing/scaling figures and the comparison figure:

```bash
julia --project=. benchmark/plot_benchmark.jl
```

The plots read only these in-repository CSVs. The Python companion checkout is
not needed to render them. An optional first argument selects another Julia CSV
and a second selects a Python CSV.

Run a new Julia multi-size scan when measuring the current code:

```bash
julia --project=. --threads=10 benchmark/benchmark.jl --scan
```

A default run without `--scan` measures one representative size per model and
keeps the stored full-size scan. Report thread count and hardware with any new
timing comparison.

Source SHA-256 (Python scan): `d20458cd3a66562816597e33d807095f09f2ff31f6c10234f186060e39a66977`.
