# Diagnosing the two phases neighboring the fermionic FCI

## Scope and present conclusion

The fermionic FCI itself is taken as established and is not reconsidered here.
This notebook asks what happens on its two sides as $t''$ is varied:

- the negative-$t''$ anomalous Hall crystal (AHC), and
- the positive-$t''$ phase historically labeled `CDW` in the data layout.

The present evidence for the AHC is mutually consistent at zero flux: a
period-three translation multiplet, a common Bragg wavevector, an integer pump
for every selected branch, and a nonzero charge gap.  It is a strong working
finite-size identification.  The qualification “selected branch” matters:
the old rank-three bundle touches outside states at some twists, so the current
pump is not yet a globally isolated-bundle invariant.  The new moderate point
must retest this, in addition to scaling the Bragg weight and neutral gaps.

The positive-$t''$ phase is **not yet established as a CDW**.  The current
spectra, geometry-dependent low-energy multiplicities, and nonuniform pump
results justify treating its nature as an open question.  The name `CDW` is
retained in filenames and configuration only for backward compatibility.

The new common diagnostic points are

| working label | physical $t''$ | numerator $x=(2+2\sqrt2)t''$ |
|:--|--:|--:|
| AHC | $-0.54$ | $-2.6073506474$ |
| candidate CDW | $+0.25$ | $+1.2071067812$ |

The already-computed sweep points $x=-2.6$ and $x=1.2$, corresponding to
$t''=-0.53848$ and $0.24853$, are close enough to preview the zero-flux
physics.  The production pump, flow, entanglement, and charge-gap calculations
should use the exact values above.

## 1. What each diagnostic can and cannot establish

These observables answer different questions and must not be substituted for
one another.

| observable | question answered | what it does **not** establish by itself |
|:--|:--|:--|
| charge gap $\Delta_c$ | Is it costly to add or remove charge? | a neutral gap, CDW order, or Hall topology |
| fixed-$N$ low-energy spectrum | Is there an isolated low-energy manifold, and what are its momenta? | the order parameter or Hall response |
| connected $S(\mathbf q)$ | Is there charge-order correlation at $\mathbf Q$? | quantized Hall response |
| charge pump / twist Chern response | Does an isolated state or manifold carry Hall response? | CDW order; it is invalid without an isolated manifold |
| spectral flow | Does the selected manifold remain isolated and evolve into itself under flux? | a phase identity without the other tests |

This immediately resolves the apparent conflict between the large green charge
gap and the nearly closed low-energy spectrum.  The computed charge gap is

\[
\Delta_c(N)=E_0(N+1)+E_0(N-1)-2E_0(N),
\]

whereas every neutral gap below is a difference between states at the **same**
particle number $N$.  A system may be charge-incompressible and still have
very soft neutral collective modes.  A density-ordered insulator near a phase
transition is a standard setting in which this can happen.

## 2. Three distinct zero-flux gap questions

At every $t''$, sort all eigenvalues from all momentum sectors,

\[
E_1\le E_2\le E_3\le E_4\le\cdots .
\]

Because the confirmed FCI has three torus states, three useful instantaneous
rank quantities are

\[
\begin{aligned}
\delta_{3} &= E_3-E_1,
&&\text{width of the three globally lowest states},\\
\Delta_{\mathrm{rank}} &= E_4-E_3,
&&\text{isolation of those three states from rank 4},\\
\Delta_{\mathrm n}^{(0)} &= E_4-E_1
 =\delta_3+\Delta_{\mathrm{rank}},
&&\text{rank-4 neutral excitation measured from the ground state}.
\end{aligned}
\]

Both $\Delta_{\mathrm{rank}}$ and $\Delta_{\mathrm n}^{(0)}$ are plotted so
that the convention is explicit.  Calling only one of them “the neutral gap”
can hide a large finite-size splitting of the presumed ground-state manifold.

There is a second limitation: if the actual low-energy band has four or eight
states, cutting it after rank 3 is artificial.  In that case
$\Delta_{\mathrm{rank}}$ may vanish simply because ranks 3 and 4 belong to
the same larger multiplet.  The full rank plot and the adjacent gaps
$E_{r+1}-E_r$ must then be inspected to find the physically separated band.

### 2.1 Fixed FCI-reference tracking

Instantaneous ranks change identity at a level crossing.  To diagnose how the
FCI is lost, we also follow the **fixed states that constituted the confirmed
FCI manifold**, writing a state as $(k_1,k_2;\ell)$, where $\ell$ is its
level within that momentum sector:

| geometry | fixed FCI reference states |
|:--|:--|
| $3\times4$ | $(2,2;1),(1,2;1),(0,2;1)$ |
| $3\times5$ | $(0,0;1),(1,0;1),(2,0;1)$ |
| $3\times6$ | $(0,3;1),(0,3;2),(0,3;3)$ |

The in-sector level is essential.  On $3\times6$, “sector $(0,3)$” alone
does not specify the three FCI states.

Here “fixed FCI reference” means a fixed **symmetry slot**
$(k_1,k_2;\ell)$, not guaranteed adiabatic wavefunction identity.  At an exact
crossing between different sectors the slot identity is unambiguous.  Through
a same-sector avoided crossing, energy-ordered levels exchange character; a
wavefunction-overlap or fidelity calculation is then needed to follow the FCI
character.  The same-sector direct gap locates precisely where this caveat
matters.

Let ${\cal F}$ denote these three reference states.  Define

\[
E_{\cal F}^{\max}=\max_{a\in{\cal F}}E_a,
\qquad
E_{\cal F}^{\min}=\min_{a\in{\cal F}}E_a,
\]

and let $E_{\overline{\cal F}}^{\min}$ be the lowest energy of every state
not in ${\cal F}$.  The signed isolation gap is

\[
\boxed{\Delta_{\cal F}^{\mathrm{signed}}
=E_{\overline{\cal F}}^{\min}-E_{\cal F}^{\max}.}
\]

- $\Delta_{\cal F}^{\mathrm{signed}}>0$: the three FCI-reference states are
  the globally lowest isolated manifold.
- $\Delta_{\cal F}^{\mathrm{signed}}=0$: a competing state touches it.
- $\Delta_{\cal F}^{\mathrm{signed}}<0$: at least one competing state lies
  below the top of the FCI reference manifold.  The FCI manifold is no longer
  the finite-size ground manifold.

Its internal width is

\[
\delta_{\cal F}=E_{\cal F}^{\max}-E_{\cal F}^{\min}.
\]

Finally, for each momentum sector $s$ represented in ${\cal F}$, let
$m_s$ be the highest in-sector reference level.  The same-sector direct gap is

\[
\boxed{\Delta_{\cal F}^{\mathrm{direct}}
=\min_s\left[E(s,m_s+1)-E(s,m_s)\right].}
\]

Thus, on $3\times4$ and $3\times5$ we compare level 2 with level 1 in each
of the three FCI sectors.  On $3\times6$ we compare level 4 with level 3 in
sector $(0,3)$.

These two gaps distinguish two finite-size mechanisms:

- At a first-order crossing between different momentum sectors,
  $\Delta_{\cal F}^{\mathrm{signed}}$ can change sign while
  $\Delta_{\cal F}^{\mathrm{direct}}$ remains finite.
- If the relevant competitor has the same momentum, levels generally avoid
  crossing at finite size; the same-sector direct gap becomes the useful
  indicator of the avoided crossing and its size dependence.

No one of these curves is a thermodynamic critical gap.  Together they make
the state identities and finite-size crossing mechanism visible.

## 3. What the existing nearby spectra already say

At $x=1.2$, or $t''=0.248528$, the global-rank data are

| geometry | observed low-energy pattern | $E_3-E_1$ | $E_4-E_3$ | $E_4-E_1$ | $\Delta_{\cal F}^{\rm signed}$ |
|:--|:--|--:|--:|--:|--:|
| $3\times4$ | 8 states below $0.071$ | 0.0250 | $\simeq0$ | 0.0250 | $-0.2550$ |
| $3\times5$ | singlet, then a four-state cluster at $0.3432$ | 0.3432 | $\simeq0$ | 0.3432 | $-0.0842$ |
| $3\times6$ | several intruders; FCI level 3 is global rank 24 | 0.0768 | $\simeq0$ | 0.0768 | $-0.2683$ |

The $3\times4$ observation is therefore correct: the three-state isolation
gap has closed and a much larger low-energy band has formed.  It is reasonable
to call the minimum low neutral excitation “roton softening,” but once the
three-state manifold has lost its identity, a unique three-state “roton gap”
is no longer well defined.

The $3\times5$ singlet does not rescue a three-state interpretation.  The two
other FCI-sector states lie at $0.4274$, above a four-state cluster at
$0.3432$.  The signed FCI isolation is consequently negative.  On
$3\times6$, all three reference states are in $(0,3)$, but their global
ranks are 2, 7, and 24.  Selecting merely the three globally lowest states or
merely the lowest state of $(0,3)$ answers a different question.

This evidence establishes that the positive point is outside the finite-size
FCI manifold.  It does **not** identify the replacement phase.  Possibilities
include a different commensurate density wave, a larger low-energy manifold,
a phase-transition region broadened by finite size, or a more strongly
geometry-dependent/incommensurate order.

## 4. Why the old positive-side pump looks chaotic

A charge pump is topologically meaningful only when the state or subspace being
transported is spectrally isolated for the complete twist path.  For a
three-state subspace, one needs a smooth rank-three projector

\[
P(\theta)=\sum_{a=1}^{3}|\Psi_a(\theta)\rangle
\langle\Psi_a(\theta)|
\]

with a nonzero gap to every state outside the subspace at every inserted flux
$\theta$.  Inside a degenerate subspace, individual basis vectors may be
rotated by an arbitrary $U(3)$ transformation.  The invariant object is the
whole projector and its non-Abelian Berry response, not an arbitrary branch
label.

For the FCI, the isolated three-state bundle has total pump 1; a convenient
branch tracking can display $1/3$ per state together with permutation under
flux.  For the AHC data at the old point, every symmetry-related branch pumps
one charge, giving total pump 3.  This is consistent with a period-three
integer Hall crystal, but consistency is not the same as a certified global
many-body Chern invariant when the rank-three gap closes along the path.

At the old positive-side point $t''=0.31066$, the final branch values were

| geometry | branch pumps after one flux quantum |
|:--|:--|
| $3\times4$ | $0.3697,0.3697,0.3151$ |
| $3\times5$ | $0.3836,-0.1870,0.1165$ |
| $3\times6$ | $-1,-1,-1$ |

The sign of a pump depends on the flux/polarization orientation convention.
The important issue is not that the $3\times6$ value has the opposite sign;
it is that the selected three-state subspaces do not show a common isolated
manifold and a common response across geometries.

When the selected subspace touches other states, overlap-based branch matching
can switch which eigenvector is followed.  Fractional-looking or
integer-looking numbers can then occur without defining a phase invariant.
This is more fundamental than a rendering artifact in the plot.  The plotting
code may faithfully display a quantity whose topological precondition was not
satisfied.

For every future pump claim, check in this order:

1. identify the zero-flux manifold including sector **and** in-sector level;
2. verify a positive many-body gap to its complement for the entire twist path;
3. verify that the final projector returns to the initial projector, allowing
   permutations within the manifold;
4. inspect overlap singular values between neighboring twist points;
5. only then interpret the total pump and, secondarily, the branch-resolved
   trajectories.

The diagnostic plotter writes `manifold_gap_flow.svg` beside every pump plot.
It shows $E_4-E_3$, $E_4-E_1$, and $E_3-E_1$ at every available flux.  If the
orange $E_4-E_3$ curve touches zero, the assumed rank-three bundle is not
isolated and the branch pump must not be used as a phase invariant.  This is
still only a rank test; overlap singular values are required to certify smooth
parallel transport inside a nearly degenerate bundle.

## 5. AHC: why the evidence is coherent

The AHC is a charge-ordered integer Hall state.  At the nearby sweep point
$x=-2.6$, the three lowest states occur at

| geometry | three low-state momenta | $E_3-E_1$ | $E_4-E_3$ |
|:--|:--|--:|--:|
| $3\times4$ | $(1,2),(2,2),(0,2)$ | 0.0400 | 0.1388 |
| $3\times5$ | $(1,0),(2,0),(0,0)$ | 0.1442 | 0.0779 |
| $3\times6$ | $(0,3),(2,3),(1,3)$ | 0.0151 | 0.1112 |

In every geometry, $k_1=0,1,2$ appears at a common $k_2$.  This is exactly
the momentum pattern expected when translation along direction 1 is broken to
a period-three subgroup.  The connected structure-factor maximum is at
$(-2\pi/3,0)$, equivalent modulo a reciprocal vector to
$(4\pi/3,0)$, in all three geometries.  The old diagnostic point additionally
gave pump $1,1,1$.

These are not three unrelated observations: the momenta and Bragg wavevector
describe the same period-three order, while the pump distinguishes an
anomalous Hall crystal from a topologically trivial period-three CDW.  The
charge gap establishes charged incompressibility.

There is one important warning from the newly plotted old flux data.  The
minimum instantaneous rank gap $E_4-E_3$ is at numerical zero—between
$0$ and about $4\times10^{-15}$—on all three AHC geometries.  At, for example,
$3\times4$ and one quarter flux, a state outside the zero-flux AHC triplet
enters below two of its branches.  Translation symmetry prevents different
momentum sectors from hybridizing in the current calculation, so one may still
follow the selected branches and obtain the exact integers.  But a generic
symmetry-allowed definition of a globally gapped rank-three bundle is absent
at those twists.

There are two clean ways to strengthen the Hall diagnosis:

1. At the new $t''=-0.54$ point, verify that the rank-three isolation remains
   positive throughout the twist path and that the projector returns to
   itself.
2. For a symmetry-broken crystal, add an infinitesimal commensurate pinning
   field to select one CDW domain, then compute the pump/Chern response of that
   unique pinned ground state while checking its global gap.  The pinning field
   should be taken to zero only after the thermodynamic limit.

The remaining thermodynamic checks are:

- The code uses

  \[
  S(\mathbf q)=\frac1{N_s}\sum_{ij}e^{i\mathbf q\cdot
  (\mathbf r_i-\mathbf r_j)}
  \left(\langle n_i n_j\rangle-\langle n_i\rangle\langle n_j\rangle\right).
  \]

  Long-range charge order therefore requires $S(\mathbf Q)\propto N_s$, or
  a nonzero extrapolation of $S(\mathbf Q)/N_s$, at fixed compatible aspect
  ratio.
- The translation-multiplet width should shrink with size while the gap above
  that multiplet remains finite.
- The integer pump should survive the move from $t''=-0.6213$ to the new
  common point $t''=-0.54$ and remain well defined throughout the twist.

Thus “strong, coherent working diagnosis” is justified.  Calling the present
pump itself “robustly quantized” without its isolation qualification, or
calling the phase thermodynamically proven from three elongated clusters,
would be stronger than the data.

## 6. Positive side: what would establish a CDW?

A genuine finite-size CDW diagnosis requires one coherent package:

1. A common ordering wavevector $\mathbf Q$ across compatible geometries.
2. $S(\mathbf Q)/N_s$ extrapolating to a nonzero value.
3. A low-energy tower whose number and momenta follow the translations broken
   by $\mathbf Q$.
4. A shrinking splitting inside that tower and a finite neutral gap above it.
5. A pump computed for that actual isolated tower, not for a forced
   three-state selection inherited from the FCI.

At $x=1.2$, the structure-factor maxima are respectively
$(2\pi/3,-\pi/2)$, $(-2\pi/3,2\pi/5)$, and
$(2\pi/3,-2\pi/3)$ on $3\times4$, $3\times5$, and $3\times6$.
Together with the different low-energy multiplicities, this is not yet a
common CDW pattern.

Moving the diagnostic point from $t''=0.31066$ to $0.25$ is sensible:
finite-size phase boundaries can drift, and the old point is visibly far from
the FCI boundary on some clusters.  It is nevertheless not guaranteed that
$t''=0.25$ realizes the same phase on every size.  The new point is a fairer
comparison point whose phase must be inferred from the tests above, not from
its preset label.

There is also a useful filling constraint.  At one fermion per three unit
cells, a unique, fully gapped, translation-symmetric, short-range-entangled
thermodynamic ground state is obstructed by the usual flux-insertion/LSM
argument.  If the positive phase is gapped and is not topologically ordered,
it must break translation symmetry.  A finite-size singlet is compatible with
that because exact eigenstates are translation-symmetric cat states; the
partners must emerge with the predicted momenta as size grows.  The current
job is to find those partners and their order, not to infer “no CDW” from one
finite-size singlet.

## 7. New plots and how to read them

Running the sweep plotter now creates

```text
phase_exploration/figures/sweep/
  zero_flux_gap_diagnostics_3x4.svg
  zero_flux_gap_diagnostics_3x5.svg
  zero_flux_gap_diagnostics_3x6.svg
  zero_flux_gap_diagnostics_all_geometries.svg
```

Each figure has two panels:

- **Instantaneous global ranks:** $E_4-E_1$, $E_4-E_3$, and $E_3-E_1$.
- **Fixed FCI reference:** signed isolation, same-sector direct gap, and FCI
  reference width.

The vertical guides mark $t''=-0.54$ and $+0.25$.  In the right panel, the
zero line has a direct meaning: crossing below it says that other states have
entered below the top of the fixed FCI manifold.  A finite green same-sector
gap at the same crossing is evidence for an inter-sector, first-order-like
finite-size crossing rather than a same-sector avoided crossing.

The plots use the existing zero-flux sweep CSVs and require no new ED.  The new
HPC campaign is needed for exact-point charge gaps, flux response, and other
diagnostics.

## 8. Production workflow and decision table

The generated HPC jobs now rerun only AHC and candidate-CDW diagnostics and
charge gaps.  They use `--refresh true` to replace derived CSVs that previously
contained data at different $t''$.  ED checkpoints are stored in
parameter-specific `x_<numerator>` subdirectories, so an interrupted new job
resumes without loading the old point or discarding completed sectors.  Protocol
marker filenames prevent the submission helper from mistaking old result files
for the new campaign.  Confirmed FCI data are left untouched.

The spectral flow remains a **one-dimensional** boundary-twist path over three
flux quanta.  It now uses 145 points, or 48 intervals per flux quantum, so the
curves and flux-dependent isolation gaps are visually smooth.  No second twist
angle, two-dimensional twist torus, or many-body Chern integration is added.

Generate and submit with

```bash
bash phase_exploration/hpc/hyak_slurm_gen.sh
phase_exploration/hpc/generated/submit_all.sh
```

After pulling the CSVs, regenerate figures with

```bash
julia --project=. phase_exploration/bin/plot_results.jl --kind all
```

Interpret the positive-side result using the following logic:

| finite-size evidence after scaling | working interpretation |
|:--|:--|
| common $\mathbf Q$, predicted translation tower, $S(\mathbf Q)/N_s>0$, zero pump | conventional CDW |
| same charge order, isolated tower, integer pump per broken-symmetry branch | AHC or another Hall crystal |
| no common $\mathbf Q$, stable topological manifold and nontrivial total pump | topological phase, not a conventional CDW |
| inconsistent multiplicities/wavevectors and collapsing neutral gaps | transition region or strong finite-size frustration; do not assign a phase yet |
| larger stable tower with matching momenta | different enlarged-unit-cell order; count that tower rather than forcing three states |

The crucial methodological change is simple: **first determine the isolated
low-energy subspace from the spectrum and symmetry, then compute and interpret
its response**.  The phase name must follow that evidence.
