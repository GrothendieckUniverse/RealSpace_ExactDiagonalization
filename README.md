# Package `RealSpace_ExactDiagonalization.jl`
--- Symmetry-resolved exact diagonalization on _arbitrary_ real-space graphs, **supporting both spin/bosonic and fermionic systems with arbitrary spins or other internal degrees of freedom**.

A high-performance, statistics-agnostic Julia implementation similar to the design philosophy of [XDiag](https://github.com/awietek/xdiag). The package block-diagonalizes interacting quantum lattice Hamiltonians via bitmask encoding, orbit-stabilizer decomposition, and irrep-induced projection, with or without forming the full many-body matrix.

---

## Table of Contents

- [Design Innovations](#design-innovations)
- [Architecture](#architecture)
- [Theoretical Background](#theoretical-background)
- [⚠️ Understanding `filling_fraction`](#⚠️-understanding-filling_fraction)
- [Quick Start](#quick-start)
- [Examples](#examples)
- [Many-Body Topological Observables](#many-body-topological-observables)
- [Benchmarks](#benchmarks)
- [File Structure](#file-structure)
- [Dependencies](#dependencies)
- [References](#references)

---

## Design Innovations

1. **Bitmask encoding** — every Fock configuration $|n_1,\ldots,n_N\rangle$ (hard-core, $n_i\in\{0,1\}$) is a single `UInt` integer $m = \sum_i n_i 2^{i-1}$, enabling $O(1)$ bitwise operations via single-cycle CPU instructions (`popcount`, `ctz`, bitwise AND/OR/XOR).

2. **Orbit-stabilizer decomposition** — the many-body Hilbert space is partitioned into orbits under the symmetry group $G$. Each orbit is labelled by a canonical representative $|[\mathbf{s}]\rangle$ and its stabilizer subgroup data. Only representatives are stored, achieving the optimal $|G|$-fold compression.

3. **Irrep-induced projection** — 1D irreducible representations of finite abelian groups supply projectors $P_\chi$ that block-diagonalize the Hamiltonian without ever constructing the full matrix. An orbit contributes to irrep $\chi$ iff its stabilizer phases satisfy a compatibility condition.

4. **Two computational modes, one CanonicalMap** — a *matrix mode* that precomputes sparse CSC matrices for fast Arpack diagonalization (memory-intensive but fast), and a *matrix-free mode* that computes $H|\psi\rangle$ on-the-fly via multithreaded Lanczos (near-zero memory overhead, ~1.5–3× slower per sector). Both modes share the same `CanonicalMap` cache for O(1) canonical-representative lookups, providing 3–10× speedup in matrix construction over raw O(|G|) canonicalization.

5. **Unified boson/fermion treatment** — the entire pipeline is statistics-agnostic. Fermionic signs (permutation parity in symmetry actions and Jordan-Wigner strings in hopping) are injected via compile-time multiple dispatch on `Bosonic()` / `Fermionic()` singleton types, with zero runtime branching overhead.

6. **Twisted boundary conditions & spectral flow** — `update_second_quantized_model_with_twisted_phases!` applies Peierls substitution via `TightBinding.generate_bilinear_terms`. For flux scans the orbit catalog is built once and its stabiliser phases are updated in-place via `update_orbit_stabilizer_phases!` — avoiding the expensive O(C(N,νN)) Gosper iteration at every flux point. The flux-aware translation group preserves ordinary momentum labels.

---

## Architecture

```
                           ┌──────────────────────────────────────┐
                           │  Real_Space_Second_Quantized_Model   |
                           │  · lattice geometry                  │
                           │  · bilinear terms  (a†_i a_j)        │
                           │  · density terms   (n_i n_j)         │
                           │  · particle_statistics      (Bosonic/Fermionic)│
                           └──────────────┬───────────────────────┘
                                          │
              ┌───────────────────────────┼───────────────────────────┐
              │                           │                           │
    ┌─────────▼───────────┐    ┌───────────▼──────────┐    ┌──────────▼────────────┐
    │  Symmetry_Operation │    │  Gosper's Hack       │    │  Bitwise_Operations   │
    │  · perm:  π_g(i)    │    │  enumerate fixed-    │    │  · occupy/empty       │
    │  · phase: η_g(i)    │    │  weight bitmasks     │    │  · count_ones (popcnt)│
    │  · label g ∈ G      │    │  in lexicographic    │    │  · trailing_zeros     │
    └─────────┬───────────┘    │  order, O(1)/iter,   │    └───────────────────────┘
              │                │  zero allocation     │
              │                └───────────┬──────────┘
              │                           │
    ┌─────────▼──────────┐    ┌───────────▼───────────┐
    │  Finite_Symmetry_   │    │  Symmetry_Orbit_      │
    │  Group              │    │  Catalog              │
    │  · |G| ordered ops  │    │  · representative_mask │
    │  · identity index   │    │  · stabilizer orders   │
    └─────────┬──────────┘    │  · stabilizer phases    │
              │               └───────────┬───────────┘
              │                           │
              │               ┌───────────▼───────────┐
              │               │  OneDim_Irrep χ       │
              │               │  · label (e.g. k₁,k₂) │
              │               │  · χ(g₁), …, χ(g_|G|) │
              │               └───────────┬───────────┘
              │                           │
              │               ┌───────────▼───────────┐
              └───────────────┤  Symmetry_Sector_     │
                              │  Basis                │
                              │  · χ-compatible orbits │
                              │  · repr → index dict   │
                              └───────────┬───────────┘
                                          │
                         ┌────────────────┴────────────────┐
                         │                                 │
              ┌──────────▼──────────┐          ┌───────────▼──────────┐
              │  Matrix Mode         │          │  Matrix-Free Mode     │
              │                     │          │                      │
              │  CanonicalMap cache │          │  CanonicalMap cache  │
              │  build sparse CSC   │          │  populate (1×) then  │
              │    ↓                │          │  Threads.@threads H|ψ⟩│
              │  Arpack eigs        │          │  KrylovKit eigsolve  │
              │                     │          │                      │
              │  Memory: O(nnz)     │          │  Memory: O(dim)      │
              │  Speed:  1×         │          │  Speed:  ~1.5–3×     │
              │                     │          │  slower              │
              │  Parallel:          │          │  Parallel:           │
              │   pmap build        │          │   Threads.@threads   │
              │   BLAS multi-thread │          │   per-thread buffers │
              └─────────────────────┘          └──────────────────────┘
```

### Pipeline

1. **Define model** — construct the lattice (via `TightBinding`), add hopping and interaction terms with particle statistics.
2. **Build symmetry group** — generate `Symmetry_Operation`s (permutations + U(1) phases) and their 1D irreps.
3. **Enumerate configurations** — Gosper's hack iterates all bitmasks at fixed particle number in lexicographic order.
4. **Orbit-stabilizer decomposition** — partition bitmasks into $G$-orbits; record canonical representatives and stabilizer phases.
5. **Filter by irrep** — for a given character $\chi$, keep only orbits satisfying $\chi(h) = \alpha_h([\mathbf{s}])$ for all $h\in\mathrm{Stab}$.
6. **Build Hamiltonian block** — construct the sparse CSC matrix (matrix mode, accelerated by `CanonicalMap` cache) or pre-populate the `CanonicalMap` and run matrix-free Lanczos.
7. **Diagonalize** — Arpack (CSC) or KrylovKit (matrix-free) to obtain eigenvalues and eigenvectors.
8. **Post-process** — analyse spectra, compute correlators, checkpoint and resume.

### Parallelism Strategy (HPC-ready)

```
BLAS threads:  1 (default) — don't compete with Julia workers
CanonicalMap:  uniform O(1) cache across matrix / distributed / matrix-free modes
Build H:       Distributed.pmap across workers (matrix mode), each with own CanonicalMap
Diag:          BLAS.set_num_threads(nprocs()) temporarily for Arpack/KrylovKit
               … then restore BLAS = 1
Matvec H|ψ⟩:   Threads.@threads :static with per-thread accumulation buffers
               (CanonicalMap pre-populated single-threaded for thread safety)
GC:            explicit GC.gc(true) after each sector
Checkpoint:    JLD2 serialization of Symmetry_Resolved_ED_Data
```

---

## Theoretical Background

### Bitmask Encoding

With the hard-core constraint $n_i \in \{0,1\}$, each configuration is a binary string of length $N$:

$$\boxed{|\mathbf{s}\rangle \equiv |n_1,\ldots,n_N\rangle \;\longmapsto\; m = \sum_{i=1}^N n_i\,2^{\,i-1} \;\in\; \{0,1,\ldots,2^N-1\}}$$

Bit $i-1$ (0-based) corresponds to vertex $i$ (1-based). Key advantages:
- Occupancy test, creation, annihilation — all $O(1)$ bitwise operations.
- Lexicographic total order via integer comparison — natural for orbit representatives.
- `count_ones` (POPCNT), `trailing_zeros` (TZCNT), `&`, `|`, `^` — all single-cycle CPU instructions.

### Symmetry Action

For a finite group $G$ with unitary representation $U_g$ commuting with $H$:

$$U_g\,a_i^\dagger\,U_g^{-1} = \eta_g(i)\,a_{\pi_g(i)}^\dagger$$

On a bitmask $m$:

$$U_g|m\rangle = \alpha_g(m)\,|m'\rangle,\qquad m' = \sum_{i\in\mathrm{occ}(m)} 2^{\pi_g(i)-1},\qquad \alpha_g(m) = \prod_{i\in\mathrm{occ}(m)}\eta_g(i)$$

Fermionic statistics add a permutation parity factor: $\mathrm{sgn}(\pi_g|_{\mathrm{occ}(m)})$ computed inline via `count_ones(new_mask >> (π_g(i)-1))` during the bit-loop.

### Orbit-Stabilizer Theorem

$$|\mathrm{Orb}(\mathbf{s})| = \frac{|G|}{|\mathrm{Stab}(\mathbf{s})|}$$

The canonical representative is the smallest bitmask in the orbit (lexicographic gauge). Knowing only the $N_\text{orbits}$ representatives and their stabilizer data suffices to reconstruct the full Hilbert space, with compression ratio $\approx N_\text{orbits}/\binom{N}{N_e} \approx 1/|G|$.

### Irrep Projector and Compatibility

For a 1D irrep $\chi$ of a finite abelian group:

$$P_\chi = \frac{1}{|G|}\sum_{g\in G} \chi(g)^*\,U_g,\qquad P_\chi^2 = P_\chi,\qquad [P_\chi, H] = 0$$

An orbit representative $|[\mathbf{s}]\rangle$ contributes to sector $\chi$ iff:

$$\boxed{\chi(h) = \alpha_h([\mathbf{s}])\;\;\forall h\in\mathrm{Stab}(\mathbf{s})}$$

When compatible, the orthonormal projected basis state is:

$$|\widetilde{[\mathbf{s}];\chi}\rangle = \sqrt{\frac{|\mathrm{Stab}(\mathbf{s})|}{|G|}}\sum_{g\in G/\mathrm{Stab}(\mathbf{s})} \chi(g)^*\,U_g|[\mathbf{s}]\rangle$$

### Hamiltonian Matrix Elements

For a scattered configuration $|\mathbf{m}\rangle = H|[\mathbf{s}]\rangle$, projected back to the sector:

$$H_{\mathbf{s}',\mathbf{s}}^\chi = \mathsf{coeff}\cdot\sqrt{\frac{|\mathrm{Stab}(\mathbf{m})|}{|\mathrm{Stab}(\mathbf{s})|}}\;\delta_{[\mathbf{s}'],[\mathbf{m}]}$$

where $\mathsf{coeff} = \alpha_{g_{\mathbf{m}}}\,\chi(g_{\mathbf{m}})^*$ and $g_{\mathbf{m}}$ maps $|\mathbf{m}\rangle$ to its canonical representative.

### Jordan-Wigner String (Fermions)

The hopping sign for fermions depends only on occupied sites strictly between $i$ and $j$:

$$c_j^\dagger c_i\,|\mathbf{s}\rangle = (-1)^{\sum_{k=\min(i,j)+1}^{\max(i,j)-1} n_k}\;|\mathbf{s}; i \to j\rangle$$

Computed in $O(1)$ via `count_ones(m & between_mask)`.

---

## ⚠️ Important: Understanding `filling_fraction`

The keyword `filling_fraction` in `build_ed_data` means **particle number per _flattened_ graph vertex**, i.e. `n_filled / n_total_vertices`. This is NOT the same as the "filling per band" used in many condensed-matter communities. The two conventions differ:

| Community / Model | What they call "filling" | `filling_fraction` in this code |
|---|---|---|
| **Spinful Hubbard** (2×Lx×Ly vertices) | "half-filling" = 1 electron per spatial site | `filling_fraction = 1//2` (N_e = Lx·Ly out of 2·Lx·Ly) |
| **FCI / bosonic Hubbard** (2×Lx×Ly vertices) | "ν = 1/2 per band" = 1/4 of all graph sites | `filling_fraction = 1//4` (e.g. 3//12 for 2×3) |
| **Spinless fermions** (Lx×Ly vertices) | "half-filling" = N/2 particles | `filling_fraction = 1//2` |
| **Spin-½ Heisenberg chain** (N vertices) | "half-filling" = N/2 bosons | `filling_fraction = 1//2` |

> **Rule of thumb**: Always compute `filling_fraction = N_particles / N_total_graph_vertices`. Count ALL internal degrees of freedom (spin, sublattice, band, etc.) as separate graph vertices.

---

## Quick Start

```julia
using RealSpace_ExactDiagonalization
using TightBinding

# ── Build the Haldane honeycomb model (2×3 unit cells, 3 hard-core bosons) ──
# Step 1: real-space lattice
r_data = TightBinding.initialize_real_space_lattice(;
    sample_size=[2, 3],
    brav_vec_list=[[1.0, 0.0], [1/2, sqrt(3)/2]],
    sub_crys_list=[[0.0, 0.0], [1/3, 1/3]],
    lattice_name="Haldane_Honeycomb", pbc_indicator=[true, true])
lattice = r_data

# Step 2: tight-binding model with hoppings
tb = TightBinding.initialize_real_space_tightbinding_model(lattice; model_name="haldane")
t, t′, t′′, ϕ = 1.0, 0.60, -0.58, 0.2
add_hopping_term!(tb, (([0,0],1),([0,0],2)) => -t; is_hermitian=true)
add_hopping_term!(tb, (([0,0],1),([0,-1],2)) => -t; is_hermitian=true)
add_hopping_term!(tb, (([0,0],1),([-1,0],2)) => -t; is_hermitian=true)
# ... (see examples/boson_fci_haldane.jl for full set)

# Step 3: assemble second-quantized model
bilinear_terms = [(lattice.site_to_index_map[sf], lattice.site_to_index_map[st], ComplexF64(t))
                  for ((sf,st),t) in tb.full_hopping_map]
second_quantized_model = Real_Space_Second_Quantized_Model(
    Dict("t"=>t), lattice, tb, Bosonic(), bilinear_terms, Tuple{Int,Int,ComplexF64}[])

# ── Translation symmetry: |G| = 6 ──
symmetry = build_translation_group(lattice)

# ── Build ED data: 3 bosons / 12 graph vertices = filling_fraction 3//12 ──
#    (⚠️ NOT ν=1/2 per band — that would be 6 bosons!)
ed_data = build_ed_data(second_quantized_model; filling_fraction=3//12, symmetry_group=symmetry)

# ── Scan all momentum sectors (matrix-free, 8 threads) ──
ed_scan!(ed_data; nev=5, mode=:matrixfree)

# ── Inspect ──
print_spectrum(ed_data)
fig, ax = plot_spectrum(ed_data)
save("spectrum.svg", fig)
```

### Checkpoint-Resume (HPC-safe)

```julia
# Run with checkpoint — survives preemption
ed_scan!(ed_data; mode=:matrixfree, checkpoint_path="checkpoints/run1.jls")

# Resume — already-computed sectors are automatically skipped
ed_data = load_checkpoint("checkpoints/run1.jls")
ed_scan!(ed_data; mode=:matrixfree, checkpoint_path="checkpoints/run1.jls")
```

---

## Examples

Three self-contained, pedagogical example scripts in `examples/`. Each builds a lattice via `TightBinding`, constructs the second-quantized model, sets up translation symmetry, and runs symmetry-resolved ED.

### 1. Spin-½ Heisenberg Chain — `examples/spin_heisenberg_chain.jl`

$$H = J\sum_{\langle i,j\rangle} \mathbf{S}_i\cdot\mathbf{S}_j \qquad (\text{1D chain, PBC, } J>0)$$

Mapped to hard-core bosons via the **Matsubara-Matsuda** transformation:

$$S^z_i = n_i - \tfrac12,\qquad S^+_i = b^\dagger_i,\qquad S^-_i = b_i$$

$$\mathbf{S}_i\cdot\mathbf{S}_j = n_i n_j + \tfrac12(b^\dagger_i b_j + \mathrm{h.c.}) - \tfrac12 n_i - \tfrac12 n_j + \tfrac14$$

At half-filling ($N_e=N/2$), constant contributions sum to $-JN/4$ and are absorbed by on-site density terms $(i,i,-J/2)$. Translation symmetry $\mathbb Z^N$ gives $N$ momentum sectors.

```bash
julia --project=. examples/spin_heisenberg_chain.jl
```

**Benchmark:** $E_0/N = -\ln 2 + 1/4 \approx -0.443147$ (Bethe ansatz, thermodynamic limit).

### 2. Bosonic Fractional Chern Insulator — `examples/boson_fci_haldane.jl`

$$H = \sum_{\langle i,j\rangle} t_{ij}\,b^\dagger_i b_j + \sum_{\langle\langle i,j\rangle\rangle} t'_{ij}\,b^\dagger_i b_j + \sum_{\langle\langle\langle i,j\rangle\rangle\rangle} t''_{ij}\,b^\dagger_i b_j$$

Haldane honeycomb lattice, 2×3 unit cells (12 sites), 3 hard-core bosons at $\nu=1/2$ per band. Complex next-nearest-neighbour hoppings $t' = -0.60\,e^{\pm i\phi}$ ($\phi=0.4\pi$) break time-reversal symmetry. Parameters from D.N. Sheng et al., PRL **107**, 146803 (2011).

```bash
julia --project=. examples/boson_fci_haldane.jl
```

**Expected:** Two nearly degenerate ground states at $E \approx -7.1638$ (k=(0,0)) and $E \approx -7.1634$ (k=(1,0)), topological splitting $\Delta E \approx 4.3\times10^{-4}$.

### 3. Spinful Fermi-Hubbard Model — `examples/fermion_hubbard_square.jl`

$$H = -t\sum_{\langle i,j\rangle,\sigma} \big(c^\dagger_{i\sigma} c_{j\sigma} + \text{h.c.}\big) + U\sum_i n_{i\uparrow} n_{i\downarrow}$$

Spin degrees of freedom are handled by the **flattened-graph approach**: each spatial site $i$ generates two interleaved graph vertices — $i_\uparrow$ (vertex $2i-1$) and $i_\downarrow$ (vertex $2i$). This preserves correct fermionic anticommutation via the Jordan-Wigner string, with no modification to the `Symmetry_Operation` infrastructure.

2×3 spatial unit cells → 12 graph vertices, $N_\uparrow=N_\downarrow=3$ ($N_e=6$), $t=1$, $U=8$. Translation group $\mathbb Z^2\times\mathbb Z^3$.

```bash
julia --project=. examples/fermion_hubbard_square.jl
```

---

## Twisted-Boundary Observables

The package includes built-in flux-scan observables for twisted-boundary spectral flow and fractional charge pumping. This mirrors the ExactDiagonalization.jl FCI showcase: thread flux through one periodic direction and diagonalize the low-energy many-body spectrum, while keeping the calculation in flux-aware symmetry sectors.

### Flux-Aware Translation Symmetry

The central innovation is the **flux-aware translation group**.  When a flux $\theta$ is inserted, the ordinary translation $T$ does _not_ commute with $H(\theta)$.  We therefore build a gauge-covariant translation $T^\theta$ whose per-site phases satisfy $g(x)/g(Tx)$ with $g(x)=e^{i2\pi\theta\cdot x/L}$.  $T^\theta$ commutes with $H(\theta)$ while keeping the **standard momentum labels** $[k_1,k_2]$ unchanged — the entire flux physics is captured in the per-site phases of the group operations.

For a flux scan the orbit catalog is built **once** (at $\theta=0$) using [`build_ed_data`](@ref).  At each subsequent $\theta$ the catalog's stabiliser phases are updated in-place via [`update_orbit_stabilizer_phases!`](@ref), avoiding re-running the expensive Gosper enumeration at every flux point.

```julia
using RealSpace_ExactDiagonalization, TightBinding, CairoMakie

# Build the bosonic Haldane FCI model (ν=1/2 per band)
model = build_zero_flux_bosonic_fci_second_quantized_model(; sample_size=[2, 3])

# Sector-resolved flux scan: track momentum sectors (0,0) and (1,0)
result = flux_spectrum_flow(
    model,
    [(0, 0), (1, 0)];
    filling_fraction=1//4,   # 3 bosons / 12 vertices
    flux_direction=1,
    twisted_phases_over_2π_list=collect(range(0.0, 2.0; length=9)),
    nev=3,
    fig_path="figures/haldane_fci_flux_flow_sectors.svg",
    checkpoint_dir="checkpoints",
)

# Or full Hilbert-space scan (identity group)
result_full = flux_spectrum_flow(
    model,
    :identity;
    filling_fraction=1//4,
    twisted_phases_over_2π_list=collect(range(0.0, 1.0; length=9)),
    nev=3,
    fig_path="figures/haldane_fci_flux_flow_identity.svg",
)
```

### Fractional Charge Pump

`flux_charge_pump` computes a one-dimensional flux-cylinder pump, not a two-dimensional many-body Chern number.  In 2D the default convention is Laughlin's: insert flux along `flux_direction` and measure the periodic many-body polarization in the transverse direction,

$$
\hat U_\perp=\exp\!\left(\frac{2\pi i}{L_\perp}\sum_j x_{j,\perp}\hat n_j\right).
$$

At fractional filling, a single momentum-sector expectation value of $\hat U_\perp$ can vanish or miss the topological multiplet.  The implementation therefore projects $\hat U_\perp$ into the requested low-energy manifold and unwraps the phases of its eigenvalues.  The stored `polarizations` keep these raw unwrapped phases, while `pumped_charge_trajectories` subtract each branch's initial phase so plots start from zero.  For the bosonic Haldane FCI on `[2,3]`, the two polarization branches each wind by $\Delta Q = 1/2$ over one inserted flux quantum.

```julia
model = build_zero_flux_bosonic_fci_second_quantized_model(; sample_size=[2, 3])

pump = flux_charge_pump(
    model,
    default_fci_sectors([2, 3]);
    filling_fraction=1//4,
    flux_direction=1,        # θ_x flux
    polarization_direction=2, # U_y polarization; this is the 2D default
    twisted_phases_over_2π_list=collect(range(0.0, 1.0; length=9)),
    fig_path="figures/haldane_fci_charge_pump_sectors.svg",
)

pump.pumped_charges  # approximately [0.5, 0.5]
```

### Key Functions

| Function | Description |
|----------|------------|
| `build_zero_flux_bosonic_fci_second_quantized_model(; sample_size, params)` | Construct the bosonic Haldane FCI model |
| `build_zero_flux_fermionic_fci_second_quantized_model(; sample_size, params)` | Construct the fermionic checkerboard FCI model (ν=2/3) |
| `default_fci_sectors(sample_size)` | Return momentum sectors hosting the two bosonic FCI ground states |
| `default_fci_sectors_fermionic(sample_size)` | Return momentum sectors hosting the three fermionic FCI ground states |
| `ed_scan!(ed_data; kwargs...)` | Unified ED scan — conventional or flux-scan mode with checkpoint resume |
| `flux_spectrum_flow(model, labels; kwargs...)` | Scan E(θ) for given sector labels |
| `flux_charge_pump(model, labels; kwargs...)` | Compute the one-dimensional fractional charge pump |
| `many_body_position_phases(lattice, direction)` | Build the site phases for Resta's periodic position operator |
| `update_orbit_stabilizer_phases!(catalog, group, stats)` | In-place stabiliser update for flux scans |
| `build_translation_group(lattice, [θ])` | Build translation group, optionally with flux phases |
| `update_second_quantized_model_with_twisted_phases!(model; twisted_phases_over_2π)` | In-place Peierls substitution to hopping terms |
| `ed_scan_checkpoint_filename(model, θ, filling)` | Canonical self-describing checkpoint filename |

### Physics Background

On a torus with $L_x \times L_y$ unit cells and $n_{\text{filled}}$ particles:

- **Momentum shift**: In the boundary-gauged Hamiltonian the centre-of-mass momentum shifts by $2\pi\theta_x n_{\text{filled}}/L_x$.  Our gauge-covariant translation absorbs this shift into the group operations, keeping irrep labels fixed.
- **Spectral flow**: The two nearly-degenerate FCI ground states intertwine under one flux quantum ($\theta=1$), each contributing $\Delta Q \approx 1/2$ to the Laughlin charge pump.  After two flux quanta ($\theta=2$) each GS returns to itself.
- **Charge pump**: Threading flux in one torus direction pumps charge in the transverse direction.  The finite periodic diagnostic is the phase winding of projected Resta polarization eigenvalues, not the raw open-boundary centre of mass.
- **Sector identity**: For $[2,3]$ ($L_x=2, n_{\text{filled}}=3$), the GS at $[0,0]$ swaps to $[1,0]$ after $\theta=1$ because $n_{\text{filled}} \bmod L_x = 1$.  For $[3,4]$ ($L_x=3, n_{\text{filled}}=6$), the GS stays in its sector because $n_{\text{filled}} \bmod L_x = 0$.

### Test / Self-Check

```julia
# From within a Julia session:
using RealSpace_ExactDiagonalization
test_bosonic_fci_spectrum_flow(; sample_size=[2,3], twisted_phases_over_2π_list=collect(range(0.0,2.0;length=9)), mode=:sectors)
test_bosonic_fci_charge_pump(; sample_size=[2,3], twisted_phases_over_2π_list=collect(range(0.0,1.0;length=9)), mode=:sectors)
test_fermionic_fci_spectrum_flow(; sample_size=[3,4], twisted_phases_over_2π_list=collect(range(0.0,1.0;length=5)), mode=:sectors)
test_fermionic_fci_charge_pump(; sample_size=[3,4], twisted_phases_over_2π_list=collect(range(0.0,1.0;length=5)), mode=:sectors)
```

The spectrum-flow tests verify that the two FCI ground states exchange after one flux quantum and return after two.  The charge-pump test directly verifies the projected polarization winding $\Delta Q \approx 1/2$ for both branches.

### Fermionic Fractional Chern Insulator (Checkerboard Lattice)

The package includes a fermionic FCI model on the checkerboard lattice at $\nu=2/3$ filling of the lower Chern band, following Sun, Gu, Katsura, and Das Sarma (arXiv:1012.5864).  With nearest-neighbor repulsion $V_1=2.0$, $V_2=1.0$, the interacting ground state shows three nearly-degenerate states (GSD = 3) on the torus at crystal momenta $(0,0)$, $(1,0)$, $(2,0)$.

The fractional charge pump yields $\Delta Q \approx 2/3$ per branch (three branches summing to $2$).

> **⚠️ Minimum System Size:** The fermionic $\nu=2/3$ FCI topological phase is **only clearly visible for sample sizes $\geq [3,4]$**.  On smaller lattices (e.g. $[2,3]$ with only 12 sites and 4 fermions) the many-body gap is not well-formed and the pumped charge may not quantize to $2/3$.

```julia
model = build_zero_flux_fermionic_fci_second_quantized_model(; sample_size=[3,4])

# Spectrum flow — tracks the three FCI ground states
result = flux_spectrum_flow(
    model,
    default_fci_sectors_fermionic([3,4]);
    filling_fraction=1//3,   # ν=2/3 per band → 1/3 per vertex
    flux_direction=1,
    twisted_phases_over_2π_list=collect(range(0.0, 1.0; length=5)),
    nev=3,
    fig_path="figures/fermionic_FCI_spectrum_flow.svg",
    checkpoint_dir="checkpoints",
)

# Charge pump — each of the 3 branches winds by ΔQ ≈ 2/3
pump = flux_charge_pump(
    model,
    default_fci_sectors_fermionic([3,4]);
    filling_fraction=1//3,
    flux_direction=1,
    twisted_phases_over_2π_list=collect(range(0.0, 1.0; length=5)),
    fig_path="figures/fermionic_FCI_charge_pump.svg",
    checkpoint_dir="checkpoints",
)
pump.pumped_charges  # approximately [2/3, 2/3, 2/3]
```

---

## Benchmarks

### Comprehensive Multi-Model Benchmark

```bash
julia --project=. -p 4 -t 8 benchmark/benchmark.jl     # run all models (workstation)
julia --project=. benchmark/plot_benchmark.jl            # plot from CSV
```

Timing **one symmetry sector** per model and system size (after JIT warmup):

| Model | System Sizes | Symmetry Group | Max Vertices |
|-------|-------------|----------------|-------------|
| Spin-½ Heisenberg chain | $N = 18, 20, 22, 24, 26$ | $\mathbb Z^N$ | 26 |
| Bosonic Haldane FCI | $[2,3], [2,4], [2,5], [3,4]$ | $\mathbb{Z}^{L_1}\times\mathbb{Z}^{L_2}$ | 24 |
| Spinful Fermi-Hubbard | $[2,3], [2,4], [2,5], [3,4]$ | $\mathbb{Z}^{L_1}\times\mathbb{Z}^{L_2}$ | 24 |

- $N=26$ Heisenberg: full Hilbert space $\binom{26}{13} = 1.04\times10^7$, reduced to $\sim4.0\times10^5$ per sector.
- $[3,4]$ Hubbard: $\binom{24}{12} = 2.70\times10^6$, reduced to $\sim2.3\times10^5$ per sector.

The `CanonicalMap` cache provides an additional 3–10× speedup in matrix construction by caching repeated canonical-representative lookups.

Generated figures (in `benchmark/figures/`): system-size plots, log-log scaling plots, and a combined summary.

![Combined benchmark results](benchmark/figures/benchmark_combined.svg)

---

## Tests

```bash
julia --project=. test/runtests.jl
```

| Test | Target |
|------|--------|
| Haldane hard-core bosons, 2×3, 3 particles | FCI pair: $-7.16380536$, $-7.16337536$ |
| Matrix vs matrix-free, same system | Spectra agree to numerical precision |
| Spinless fermion open chain | Many-body energy = sum of filled single-particle levels |
| Spinful Hubbard, 3×4 open square, $U=8$, $N_\uparrow=N_\downarrow=6$ | $E_0 = -4.913259209075605$ |

The Hubbard reference is from [ExactDiagonalization.jl documentation](https://quantum-many-body.github.io/ExactDiagonalization.jl/dev/examples/HubbardModel/).

---

## File Structure

```
RealSpace_ExactDiagonalization/
├── examples/
│   ├── spin_heisenberg_chain.jl          Spin-½ Heisenberg chain (N=20)
│   ├── boson_fci_haldane.jl              Bosonic FCI on Haldane honeycomb
│   └── fermion_hubbard_square.jl         Spinful Fermi-Hubbard on square lattice
├── src/
│   ├── RealSpace_ExactDiagonalization.jl  Main module
│   ├── bitwise_operations.jl              Bitmask primitives (submodule)
│   ├── second_quantized_model.jl          Model data structures
│   ├── symmetry_resolved_ed.jl            Core ED engine
│   └── flux_utilities.jl                  Twisted-boundary-condition model builder
├── observables/
│   ├── spectrum_flow.jl                   Twisted-boundary spectrum-flow scan
│   └── charge_pump.jl                     Fractional charge-pump observable
├── doc/
│   ├── design.ipynb                       Complete design documentation (theory & code)
│   └── charge_pump.md                     Flux-cylinder charge-pump notes
├── benchmark/
│   ├── benchmark.jl                       Comprehensive multi-model benchmark
│   ├── plot_benchmark.jl                  Plotting from CSV results
│   ├── xdiag_compare.jl                   XDiag comparison benchmark
│   ├── benchmark_data/                    Generated CSV timing files
│   └── figures/                           Generated CairoMakie SVG plots
├── checkpoints/                           Serialized ED data (JLD2)
├── figures/                               Output plots
├── Project.toml
└── README.md
```

---

## Dependencies

- **Julia** ≥ 1.10
- [`TightBinding`](https://github.com/GrothendieckUniverse/TightBinding) — real-space lattice and tight-binding model construction
- `Arpack` — sparse eigensolver (matrix mode)
- `KrylovKit` — iterative eigensolver (matrix-free mode)
- `MLStyle` — algebraic data types for `Bosonic()` / `Fermionic()` dispatch
- `SparseArrays`, `LinearAlgebra`, `Distributed` — standard library
- `CairoMakie` — plotting
- `JLD2` — checkpoint serialization

---

## References

- A. Wietek, *XDiag* — [github.com/awietek/xdiag](https://github.com/awietek/xdiag); arXiv:2505.02901
- D.N. Sheng, Z.-C. Gu, K. Sun, L. Sheng, *Fractional Chern Insulator in a Bosonic Model with Flat Bands*, Phys. Rev. Lett. **107**, 146803 (2011)
- R. Gosper, *HAKMEM* Item 175 (Gosper's hack)
- Full design documentation: `doc/design.ipynb`
