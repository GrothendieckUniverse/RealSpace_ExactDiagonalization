# Interpreting the structure-factor features around the FCI

## The puzzle

The zero-flux ED spectrum on $3\times5$ has a well-separated group of three
low states over roughly

\[
t''\in[-0.41,0.15].
\]

This was previously used as the finite-size FCI window: the FCI triplet is
lost when other levels descend into it, consistent with a roton gap closing.
The structure-factor sweeps contain substantially more structure.  On
$3\times5$, the absolute-ground-state curves have a kink near $-0.35$, a
sharp jump near $-0.25$, and another kink near $-0.06$.  None coincides with
the two spectral edges above.  In particular, the sharp jump coincides with a
reordering of the three FCI ground states while $E_3-E_0$ remains large.

These observations pose three distinct questions:

1. Is the jump near $-0.25$ only an artifact of choosing the instantaneous
   rank-one ground state from a nearly degenerate FCI manifold?
2. Are the two kinks phase transitions, changes of the maximizing momentum,
   avoided crossings, or other finite-size effects?
3. Is the absence of a visible $S(\mathbf q)$ singularity near $-0.41$ and
   $0.15$ compatible with a continuous transition driven by roton softening?

The new campaign is designed to answer these questions without assigning a
phase from any one curve.

## 1. Rank-one and manifold observables answer different questions

The old sweep evaluates the connected structure factor in the absolute ground
state,

\[
S_0(\mathbf q)=\frac{1}{N_s}\sum_{ij}e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}
\left[\langle 0|n_i n_j|0\rangle-
\langle 0|n_i|0\rangle\langle 0|n_j|0\rangle\right].
\]

At an exact crossing between different momentum sectors, the identity of
$|0\rangle$ changes discontinuously.  Its structure factor can therefore jump
even if the two states are members of the same finite-size topological
ground-state manifold.  Such a jump is a state-selection effect, not evidence
for a bulk phase transition.

The appropriate comparison inside an isolated three-state FCI manifold
$\mathcal F$ is the equal-weight density matrix

\[
\rho_{\mathcal F}=\frac{1}{3}\sum_{a\in\mathcal F}|\Psi_a\rangle
\langle\Psi_a|.
\]

The code now computes the connected correlator of this mixed state.  The
one-body densities are averaged before the disconnected contribution is
subtracted, so this is the structure factor of $\rho_{\mathcal F}$ rather than
merely the average of three already-connected curves.  The result is invariant
under a unitary rotation of a complete isolated triplet and cannot change just
because its members exchange energy rank.

The fixed reference slots used through the sweep are

| geometry | three FCI reference states $(k_1,k_2;\ell)$ |
|:--|:--|
| $3\times4$ | $(2,2;1),(1,2;1),(0,2;1)$ |
| $3\times5$ | $(0,0;1),(1,0;1),(2,0;1)$ |
| $3\times6$ | $(0,3;1),(0,3;2),(0,3;3)$ |

Here $\ell$ is the in-sector level.  Thus the $3\times6$ FCI triplet is not
spread over three momentum sectors; it consists of the first three states of
one sector.  Keeping the in-sector level is essential.

The comparison has a direct interpretation:

| observation | finite-size conclusion |
|:--|:--|
| rank-one curve jumps, projector average is smooth | reordering within the selected triplet caused the jump |
| both curves change at the same point while the triplet remains isolated | the three-dimensional low-energy subspace itself changes; investigate a genuine transition or avoided crossing |
| projector average becomes irregular after outside levels enter | the fixed FCI reference projector is no longer an isolated ground manifold, so it is only a tracking diagnostic |
| peak value kinks exactly when the maximizing $\mathbf q$ changes | kink of the `max` envelope, not necessarily a nonanalytic correlator |

This test should be read first on $3\times5$ and $3\times6$ over
$t''\in[-0.4,0]$.

## 2. Why `max(S)` and `max(S)/mean(S)` can have extra features

Even if every individual $S(\mathbf q;t'')$ is smooth, the displayed maximum

\[
M(t'')=\max_{\mathbf q}S(\mathbf q;t'')
\]

is an upper envelope.  When two allowed momenta exchange which one is largest,
$M$ can have a kink.  The new sweep output therefore retains the peak
wavevector $(q_x^{\rm peak},q_y^{\rm peak})$ and plots it beside the two peak
metrics.  A kink accompanied by a peak-wavevector switch has a simple
finite-grid explanation.

The normalized curve has another possible source of structure:

\[
R(t'')=\frac{\max_{\mathbf q}|S(\mathbf q;t'')|}
{|\operatorname{mean}_{\mathbf q}S(\mathbf q;t'')|}.
\]

A smooth but small or rapidly varying denominator can amplify a weak feature
of the numerator.  The CSVs retain the numerator, `mean_S`, and the alternative
normalization by $\operatorname{mean}|S|$ so that this can be checked directly.

Finally, a crossing entirely among excited states cannot by itself make an
exact ground-state expectation value non-smooth.  Higher levels matter only
indirectly when one approaches or hybridizes with the ground state, when a
same-sector avoided crossing rapidly changes the ground-state wavefunction, or
when an approximate state-tracking algorithm switches branches.  Therefore
the kinks near $-0.35$ and $-0.06$ should not be justified merely by pointing
to a crossing of higher levels.  The defensible audit is:

1. check whether the absolute ground-state sector changes;
2. check whether the maximizing momentum changes;
3. compare the rank-one and projector-averaged curves;
4. inspect the ground-to-excited and triplet-isolation gaps;
5. if all remain inconclusive, calculate neighboring-$t''$ ground-state
   fidelity or wavefunction overlaps to expose a same-sector avoided crossing.

The first four checks are supported by the new sweep files and plots.  The
fifth requires overlap tracking and should be added only if a kink survives
the projector and peak-wavevector tests.

## 3. Smooth $S(\mathbf q)$ at the putative spectral edges is sensible

The absence of a visible kink near $t''=-0.41$ or $0.15$ does not conflict
with a continuous transition caused by roton softening.

First, a finite Hamiltonian is generally analytic in $t''$.  States in the
same symmetry sector avoid crossing, and an order parameter or correlation
function can evolve smoothly through the rounded finite-size precursor of a
thermodynamic transition.  A sharp nonanalyticity is not required on one
finite torus.

Second, as a density-carrying neutral mode at $\mathbf Q$ softens, the FCI can
develop increasingly strong *incipient* charge correlations while the neutral
gap is still nonzero.  In the thermodynamic limit, condensation of that mode
can produce charge order continuously.  Thus the statement “the FCI smoothly
develops charge correlations on approach to a roton-driven transition” is
sensible.

It is stronger, and not yet established, to say that the same FCI phase has
long-range charge order throughout the transition.  Beyond the critical point
the resulting state could be an AHC, a conventional CDW, a topological
charge-ordered phase, or a finite-size crossover between them.  The distinction
requires simultaneous scaling of

\[
\Delta_{\rm roton},\qquad S(\mathbf Q)/N_s,\qquad
\text{low-energy momenta and multiplicity},\qquad
\text{and a valid Hall response}.
\]

The interval $[-0.41,0.15]$ should therefore be called the observed
finite-size three-state-isolation window, not a thermodynamic phase boundary
determined by ED alone.

## 4. Energy conventions and state identity

At each $t''$, order all zero-flux energies from all momentum sectors as

\[
E_0\le E_1\le E_2\le E_3\le\cdots.
\]

The two compact gap curves remain

\[
W_3=E_2-E_0,\qquad \Delta_{3\to4}=E_3-E_2.
\]

$W_3$ is the width of the three globally lowest states and
$\Delta_{3\to4}$ is their instantaneous isolation.  The ground-referenced
quantity is $E_3-E_0=W_3+\Delta_{3\to4}$.

These rank gaps do not guarantee that the same physical states are being
followed.  The spectrum sweep therefore keeps the fixed FCI slots colored and
all other levels gray.  Once a gray level descends below a colored reference
state, the fixed projector average remains useful for diagnosing ancestry, but
it must no longer be called the ground-state FCI structure factor without an
overlap calculation.

## 5. Charge pump: necessary isolation checks

A charge-pump trajectory is a topological diagnostic only for a state or
subspace that is spectrally isolated for the complete twist path.  For a
rank-three bundle one needs a smooth projector

\[
P(\theta)=\sum_{a=1}^{3}|\Psi_a(\theta)\rangle
\langle\Psi_a(\theta)|
\]

with a positive gap to every state outside it at every $\theta$.  Internal
degeneracy and permutation are allowed for a non-Abelian bundle; touching an
outside state is not.  A branch-resolved pump has the stronger requirement
that the individual branch be well defined, or that a stated symmetry protects
the crossing.

For every proposed FCI, AHC, or CDW Chern statement, apply this order:

1. identify the intended zero-flux manifold by sector and in-sector level;
2. inspect the complete spectrum flow and the minimum gap to its complement;
3. verify return of the final projector to the initial projector, allowing an
   internal permutation;
4. inspect neighboring-step overlap singular values if the bundle is nearly
   degenerate;
5. only then interpret the total pump; treat branch values as secondary.

The diagnostic renderer writes `manifold_gap_flow.svg`.  If its isolation
curve reaches numerical zero, the plotted pump is exploratory branch tracking,
not a certified many-body Chern number.  This qualification is especially
important for the old AHC and positive-side claims.

## 6. New characteristic points

The previous single FCI point near $t''=-0.21$ is replaced by three points on
the two sides of the suspicious $3\times5$ jump.  The AHC and positive-side
points are also replaced as requested.

| working label | physical $t''$ | numerator $x=(2+2\sqrt2)t''$ | purpose |
|:--|--:|--:|:--|
| AHC | $-0.50$ | $-2.4142135624$ | proposed AHC, farther left |
| AHC | $-0.45$ | $-2.1727922061$ | proposed AHC, near the left spectral edge |
| FCI | $-0.30$ | $-1.4485281374$ | before the $3\times5$ jump |
| FCI | $-0.15$ | $-0.7242640687$ | after the jump |
| FCI | $0.00$ | $0$ | near the second left-window kink |
| candidate CDW | $0.05$ | $0.2414213562$ | positive-side near-boundary point |
| candidate CDW | $0.10$ | $0.4828427125$ | positive-side comparison point |
| candidate CDW | $2.00$ | $9.6568542495$ | deep positive-side reference |

`AHC`, `FCI`, and `CDW` are working directory labels, not conclusions imposed
on the data.  In particular, the actual isolated low-energy multiplicity at
each point must be checked before interpreting a forced rank-three pump.

Each characteristic-point run now produces the zero-flux spectrum, absolute-
ground-state and selected-manifold structure factors, spectral flow, pump,
spatial entanglement spectrum, particle entanglement spectrum, and summary.
Charge gaps are generated independently for $3\times3$ through $3\times7$.
Results live in parameter-specific directories such as

```text
results/diagnostics/FCI/3x5/tpp_m0p1500/
results/charge_gap/AHC/3x6/tpp_m0p5000/
```

so different characteristic points cannot overwrite one another.

## 7. New sweep outputs and how to read them

Every sweep point retains the original rank-one files and adds

```text
structure_fci_gsd_allowed.csv
structure_fci_gsd_dense.csv
structure_fci_gsd_metrics.csv
structure_fci_gsd_state_metrics.csv
```

The metric sweeps overlay the absolute ground state and the equal-weight fixed
FCI reference projector.  Dotted vertical guides mark a change of the
rank-one ground-state symmetry slot.  Separate peak-wavevector figures show
whether a kink is caused by the `max` operation switching momentum.

For the sharp $3\times5$ feature near $-0.25$, the primary result to inspect is
therefore not another rank-one plot.  It is the difference between those two
curves:

- if only the rank-one curve jumps, classify it as FCI-manifold reordering;
- if the projector curve also jumps, check isolation and state overlap before
  assigning a transition;
- if the peak momentum switches, inspect the individual $S(\mathbf q)$ values
  rather than treating the envelope kink as a thermodynamic singularity.

The same audit is essential on $3\times6$, where all three FCI reference
states occupy one momentum sector and energy ordering inside that sector can
exchange wavefunction character through avoided crossings.

## 8. Production and interpretation

Generate the Slurm jobs with

```bash
bash phase_exploration/hpc/hyak_slurm_gen.sh
phase_exploration/hpc/generated/submit_all.sh
```

After pulling the CSVs, regenerate all available figures without new ED:

```bash
julia --project=. phase_exploration/bin/plot_results.jl --kind all
```

Use the following decision table rather than a preset phase name:

| combined observation | interpretation |
|:--|:--|
| isolated FCI triplet, smooth projector-averaged $S$, rank-one jump | topological-manifold reordering artifact |
| closing neutral mode at fixed $\mathbf Q$, growing $S(\mathbf Q)/N_s$, smooth finite-size curves | consistent with a continuous roton-driven ordering transition |
| discontinuous projector observables and an inter-sector ground-manifold crossing that sharpens with size | first-order transition candidate |
| common $\mathbf Q$, translation tower, nonzero extrapolated $S(\mathbf Q)/N_s$ | charge-ordered phase |
| charge order plus a well-isolated quantized Hall response | AHC or another Hall crystal |
| pump without an isolated twist-path manifold | no Chern conclusion |
| geometry-dependent wavevectors, multiplicities, and collapsing gaps | transition region or commensuration frustration; defer the phase label |

The central methodological rule is: determine the low-energy subspace and its
isolation first, average observables over that subspace when appropriate, and
only then interpret charge order or Hall response.
