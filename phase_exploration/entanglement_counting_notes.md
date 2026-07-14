# Entanglement-spectrum counting: a derivation-led guide

This note explains the two different entanglement spectra used in this study:

1. the **particle entanglement spectrum (PES)**, whose low-lying counting
   diagnoses bulk quasihole states; and
2. the **spatial-orbital entanglement spectrum**, whose properly
   momentum-resolved low-lying branch can diagnose boundary/edge states.

The two spectra start from the same Schmidt-decomposition idea, but they cut
the Hilbert space in different ways and therefore count different physical
objects.  Most importantly, neither diagnosis requires the microscopic FCI
Hamiltonian or its bulk wavefunction to be a conformal field theory (CFT).
CFT gives an elegant continuum construction of the Laughlin state and its
long-wavelength edge, but the PES counting also has a direct zero-mode and
combinatorial derivation.

## 1. Common language: what an entanglement spectrum is

For any bipartition of a pure state,

$$
|\Psi\rangle
  = \sum_\alpha s_\alpha
    |\alpha\rangle_A |\alpha\rangle_B,
\qquad
\sum_\alpha s_\alpha^2=1,
$$

the reduced density matrix and entanglement energies are

$$
\rho_A=\operatorname{Tr}_B|\Psi\rangle\langle\Psi|
      =\sum_\alpha \lambda_\alpha
       |\alpha\rangle\langle\alpha|,
\qquad
\lambda_\alpha=s_\alpha^2,
\qquad
\xi_\alpha=-\log\lambda_\alpha.
$$

Equivalently, define the entanglement Hamiltonian
$H_E=-\log\rho_A$, up to an additive constant.  Its eigenvalues are the
entanglement spectrum.  A useful topological statement is never merely
"there is a gap in $\xi$": one must also check that the number of levels and
their symmetry sectors below that gap agree with a theoretically expected
low-energy subspace.

## 2. Particle entanglement spectrum: bulk quasihole counting

### 2.1 What the particle cut does

The particle cut keeps $N_A$ particles and traces out $N_B=N-N_A$ particles;
it does **not** select a spatial region.  In first-quantized notation,

$$
\rho_A(Z_A,Z_A')
 = \binom{N}{N_A}
   \int dZ_B\,
   \Psi(Z_A,Z_B)\Psi^*(Z_A',Z_B),
$$

where $Z_A=(z_1,\ldots,z_{N_A})$.  In the code this integral is the finite
Fock-basis contraction $\rho_A=MM^\dagger$, followed by decomposition into
the subsystem translation sectors $(K_1,K_2)$.

On a torus the topological ground states form a multiplet.  The PES therefore
uses the basis-independent projector onto the full low-energy manifold,

$$
\rho_{\rm GS}=\frac1d\sum_{a=1}^{d}|\Psi_a\rangle\langle\Psi_a|,
\qquad
\rho_A=\operatorname{Tr}_B\rho_{\rm GS}.
$$

A unitary rotation among the $d$ ground states leaves this density matrix
unchanged.  For the $\nu=1/3$ FCI, $d=3$.

### 2.2 The Laughlin wavefunction and its clustering property

On the plane, the fermionic Laughlin state at $\nu=1/m$, with odd $m$, is

$$
\Psi_m(z_1,\ldots,z_N)
 = \prod_{i<j}(z_i-z_j)^m
   \exp\!\left[-\sum_i\frac{|z_i|^2}{4\ell_B^2}\right].
$$

For this study $m=3$.  When two fermions approach one another,

$$
\Psi_3\sim (z_i-z_j)^3.
$$

Thus relative-angular-momentum-one pairs are absent.  The state is the
densest zero-energy state of the $V_1$ Haldane pseudopotential.  This cubic
coincidence zero is the physical origin of the exclusion structure used
below.

### 2.3 What the “generalized Pauli principle” actually means

Write a Landau-level many-body basis as occupation strings of guiding-center
orbitals.  The densest or **root** pattern of the $1/3$ Laughlin state is

```text
100 100 100 100 ...
```

The corresponding $(1,3)$ admissibility rule is

$$
n_j+n_{j+1}+n_{j+2}\le 1
$$

for every three consecutive orbitals, with the indices understood cyclically
on a torus.  More generally, a $(k,r)$-admissible root pattern has at most
$k$ particles in any $r$ consecutive orbitals.

This is “generalized Pauli” because ordinary fermionic Pauli exclusion only
says $n_j\le1$ for one orbital; the $(1,3)$ rule adds a correlation hole
extending over three guiding centers.

There is an important subtlety: this does **not** say that every Slater
determinant with nonzero coefficient in the isotropic Laughlin wavefunction
obeys the literal spacing rule.  The full state also contains configurations
obtained by inward pair squeezing.  Rather, admissible occupation strings
label dominant/root configurations and index the independent quasihole
zero-mode space.  In the thin-cylinder or thin-torus limit the root string
becomes the dominant electrostatic pattern, making the rule especially
transparent.  Adiabatically returning to the isotropic torus changes the
coefficients but not the zero-mode counting.

The same structure can be expressed through Jack polynomials: the Laughlin
polynomial is generated from its root partition by squeezing, while its
clustering condition selects the admissible roots.  This is one route to the
generalized exclusion rule, but Jack-polynomial or CFT machinery is not
needed for the counting derivation below.

### 2.4 Why $(1,3)$ is the correct rule for this checkerboard system

The checkerboard lattice has two sites per unit cell.  For an
$L_1\times L_2$ torus,

$$
N_{\rm cell}=L_1L_2,
\qquad
N_{\rm site}=2L_1L_2,
\qquad
N=\frac{L_1L_2}{3}.
$$

Thus the physical site filling is $N/N_{\rm site}=1/6$, but the filling of
the relevant Chern band is

$$
\nu=\frac{N}{N_{\rm cell}}=\frac13.
$$

A Chern-one band contains one one-particle state per unit cell, so
$N_{\rm cell}$ plays the role of the Landau-level orbital number $N_\phi$.
On a torus the $1/3$ Laughlin relation is $N_\phi=3N$; there is no spherical
shift.  Consequently, the maximally dense root compatible with the filling
is precisely `100100...`, and its rule is $(1,3)$.

This correspondence does not require a uniform magnetic field.  Hybrid
Wannier states of a $C=1$ band provide a one-dimensional ordering analogous
to Landau-gauge guiding centers.  As transverse momentum winds once through
the Brillouin zone, the Wannier center shifts by one cell—the Chern pump.
That spectral winding is what lets the continuum root pattern and quasihole
counting survive on a lattice.  Lattice momenta are obtained by folding the
one-dimensional orbital labels back into the $L_1\times L_2$ Brillouin zone.

For a generic FCI the $(1,3)$ rule is therefore a **phase diagnostic**, not a
microscopic hard constraint.  If the interacting state is adiabatically
connected to Laughlin topological order, a low PES band retains the Laughlin
quasihole counting even though higher PES levels and detailed entanglement
energies are nonuniversal.

### 2.5 Non-CFT reason that the PES counts quasiholes

Hold the traced-particle coordinates $Z_B$ fixed and view
$\Psi(Z_A,Z_B)$ only as a function of the retained coordinates $Z_A$.  Any
two retained particles still have the same cubic coincidence zero.  It is
also an $N_A$-particle wavefunction at the **original** orbital number
$N_\phi$, which is larger than the minimum $3N_A$.  It therefore lies in the
$N_A$-particle quasihole zero-mode space of the Laughlin parent Hamiltonian.

As $Z_B$ varies, these functions span the row space of the particle
entanglement matrix.  Hence

$$
\operatorname{rank}\rho_A^{\rm PES}
 \le \dim \mathcal H_{\rm quasihole}(N_A,N_\phi).
$$

For the ideal Laughlin state the bound is generically saturated, including
momentum-sector resolution.  For a non-model Hamiltonian in the same phase,
extra generic levels appear, but the quasihole-like levels should remain
separated by an entanglement gap.  This is the direct bulk, zero-mode
understanding of PES counting; it invokes neither an edge nor a CFT.

### 2.6 Derivation of the cyclic $(1,m)$ counting

Consider $N_A$ occupied orbitals on a labeled ring of $L=N_\phi$ orbitals.
Let $g_i$ be the number of zeros following occupied particle $i$ before the
next occupied particle.  The $(1,m)$ rule says

$$
g_i\ge m-1,
\qquad
\sum_{i=1}^{N_A}g_i=L-N_A.
$$

Set $y_i=g_i-(m-1)\ge0$.  Then

$$
\sum_i y_i=L-mN_A.
$$

The number of ordered nonnegative solutions, after choosing a distinguished
occupied particle, is the stars-and-bars result

$$
\binom{L-(m-1)N_A-1}{N_A-1}.
$$

There are $L$ choices for the orbital of that distinguished particle, but
each physical occupation string was distinguished in $N_A$ possible ways.
Therefore

$$
\boxed{
\mathcal N_{(1,m)}(L,N_A)
 =\frac{L}{N_A}
  \binom{L-(m-1)N_A-1}{N_A-1}}
$$

for $L\ge mN_A$ and $N_A>0$.

Here $m=3$ and the calculation uses $N_A=2$, so the formula reduces to

$$
\boxed{\mathcal N_{(1,3)}(L,2)=\frac{L(L-5)}2.}
$$

There is also a simple pair-counting check.  For a particle at orbital $a$,
the second particle cannot be at $a$, $a\pm1$, or $a\pm2$, leaving $L-5$
choices.  Multiplying by $L$ and dividing by two for the unordered pair gives
the same result.

For the clusters relevant here,

| geometry | $L=N_{\rm cell}$ | expected PES levels for $N_A=2$ |
|:--|--:|--:|
| `3x4` | 12 | $12(12-5)/2=42$ |
| `3x5` | 15 | $15(15-5)/2=75$ |
| `3x6` | 18 | $18(18-5)/2=117$ |

### 2.7 Derivation of the momentum-resolved counting

Number the effective Chern-band orbitals by

$$
j=k_1+L_1k_2,
\quad
0\le k_1<L_1,
\quad
0\le k_2<L_2.
$$

Enumerate every admissible unordered pair $a<b$.  If
$a\leftrightarrow(k_1^a,k_2^a)$ and
$b\leftrightarrow(k_1^b,k_2^b)$, place that pair in

$$
(K_1,K_2)=
\bigl((k_1^a+k_1^b)\bmod L_1,
      (k_2^a+k_2^b)\bmod L_2\bigr).
$$

This direct enumeration yields:

| geometry | expected count in each momentum sector |
|:--|:--|
| `3x4` | 3 levels for even $K_2$, 4 for odd $K_2$, for each $K_1$ |
| `3x5` | 5 levels in every one of the 15 sectors |
| `3x6` | 6 levels for even $K_2$, 7 for odd $K_2$, for each $K_1$ |

Changing the Wannier-orbital origin can permute sector labels, but it cannot
change these multiplicities.  Geometry therefore affects how the universal
total is folded among lattice momenta.  It also explains why distinct
topological ground states can collapse into the same many-body momentum
sector on some commensurate clusters; momentum labels alone do not determine
the dimension of the topological manifold.

### 2.8 Comparison with the downloaded PES

| geometry | expected $(1,3)$ count | levels below largest PES gap | largest gap | verdict |
|:--|--:|--:|--:|:--|
| `3x4` | 42 | not yet available | not yet available | local/HPC generation is paused |
| `3x5` | 75 | 75 | 2.72160 | exact total and sector-by-sector match |
| `3x6` | 117 | 117 | 0.93947 | numerical match, but corrected-manifold rerun required |

For `3x5`, the five levels in every momentum sector exactly match the cyclic
admissible enumeration.  This is strong finite-size evidence for the desired
$1/3$-Laughlin FCI order.

The existing `3x6` CSV also has 117 levels below its largest gap, distributed
as 6/7 in even/odd $K_2$ sectors.  However, it was generated with the old
manifold-selection bug described in Section 5.  Its agreement is encouraging
but should not be treated as final evidence until the version-2 HPC rerun.

Blindly selecting the largest gap in non-FCI data gives 431 and 44 levels for
AHC (`3x5`, `3x6`) and 3, 10, and 33 levels for CDW (`3x4`, `3x5`, `3x6`).
These do not match the Laughlin quasihole total or momentum allocation.  A
large numerical gap without the correct counting is not a topological PES
diagnosis.

## 3. Spatial-orbital entanglement spectrum: boundary counting

### 3.1 Which cut this code actually makes

Expand the lattice state in site-occupation bases of two complementary
regions:

$$
|\Psi\rangle
=\sum_{N_A}\sum_{a\in\mathcal H_A(N_A)}
                 \sum_{b\in\mathcal H_B(N-N_A)}
 M^{(N_A)}_{ab}|a\rangle_A|b\rangle_B.
$$

For each $N_A$, the singular values of $M^{(N_A)}$ give the Schmidt
probabilities.  The code therefore computes

$$
M^{(N_A)}=U S V^\dagger,
\qquad
\lambda_\alpha^{(N_A)}=S_{\alpha\alpha}^2,
\qquad
\xi_\alpha^{(N_A)}=-\log\lambda_\alpha^{(N_A)}.
$$

The present checkerboard partition cuts a contiguous strip of lattice unit
cells along direction 2:

- `3x5`: region A has 12 sites and B has 18 sites;
- `3x6`: A and B each have 18 sites.

Because these are localized site orbitals, this is a lattice real-space cut,
also reasonably called a spatial-orbital cut.  It is not literally the
original Li-Haldane orbital cut through continuum Landau guiding-center
orbitals, although both create an effective boundary and should share a
universal low branch when the cut and state are suitable.

The current code uses the lowest single eigenstate, not the equal-weight
three-state density matrix.  On a torus, different linear combinations of
the ground multiplet can have different spatial entanglement.  The cleanest
edge spectrum is generally obtained from a minimally entangled state for the
chosen cut.  This is another reason to avoid over-interpreting raw ranks from
the current output.

### 3.2 First derive the nonuniversal rank ceiling

If A contains $M_A$ sites and B contains $M_B$ sites, then at fixed $N_A$

$$
\operatorname{rank}M^{(N_A)}
\le
\min\left\{
\binom{M_A}{N_A},
\binom{M_B}{N-N_A}
\right\}.
$$

This is ordinary Hilbert-space counting, not topological edge counting.
For the downloaded FCI spectra:

| geometry | generic maximal ranks by increasing $N_A$ | observed nonzero ranks |
|:--|:--|:--|
| `3x5` | `1, 12, 66, 153, 18, 1` | `1, 12, 66, 153, 18, 1` |
| `3x6` | `1, 18, 153, 816, 153, 18, 1` | `1, 18, 153, 676, 153, 18, 1` |

Thus the `3x5` ranks simply saturate the generic matrix-rank bound.  They are
not a universal Laughlin sequence.  The central `3x6` block has reduced rank,
which reflects structure or exact/numerical support of that finite-size
state, but `676` by itself has no standard topological interpretation.  The
reflection symmetry of the `3x6` sequence under $N_A\leftrightarrow6-N_A$
is expected from its equal A/B partition.

### 3.3 Where the edge-CFT counting comes from

Here a CFT means a gapless $(1+1)$-dimensional quantum field theory at a
scale-invariant fixed point.  Its states are organized into **primary**
sectors and oscillator **descendants** above each primary.  For a Laughlin
edge, the primary sector records the anyon/edge-charge sector, while the
descendants are propagating edge-density waves.  The integer-partition
counting comes from those descendants.

The continuum Laughlin polynomial has an exact CFT representation.  Take a
canonically normalized chiral boson $\phi$ with

$$
\langle\phi(z)\phi(w)\rangle=-\log(z-w)
$$

and electron vertex operator

$$
V_e(z)=e^{i\sqrt m\,\phi(z)}.
$$

After inserting a neutralizing background charge,

$$
\left\langle\prod_i V_e(z_i)\,\mathcal O_{\rm bg}\right\rangle
\propto\prod_{i<j}(z_i-z_j)^m.
$$

This is a wavefunction construction from a conformal block.  The same chiral
boson describes the long-wavelength Laughlin edge.  In $K$-matrix language,
it is conventional to rescale the field as
$\varphi=\phi/\sqrt m$.  Then $V_e=e^{im\varphi}$, $K=(m)$, and, for $m=3$,
the edge action is

$$
S_{\rm edge}
=\frac{3}{4\pi}\int dt\,dx\,
 \left(\partial_t\varphi\,\partial_x\varphi
       -v(\partial_x\varphi)^2\right).
$$

The electron is the local $e^{i3\varphi}$ excitation in this normalization,
while $e^{i\varphi}$ creates the charge-$e/3$ anyon.  The three values of
anyon charge modulo an electron give the three primary/topological sectors
of $U(1)_3$.

At fixed charge and fixed topological sector, edge excitations are bosonic
oscillator descendants

$$
\prod_{n=1}^{\infty}(a_{-n})^{r_n}|Q\rangle,
\qquad
r_n=0,1,2,\ldots,
$$

with excess edge momentum

$$
\Delta K=\sum_{n\ge1}n r_n.
$$

Counting descendants at a given $\Delta K$ is therefore exactly the integer
partition problem.  Its generating function is

$$
Z_{\rm edge}(q)
=\prod_{n=1}^{\infty}\left(1+q^n+q^{2n}+\cdots\right)
=\prod_{n=1}^{\infty}\frac1{1-q^n}
=1+q+2q^2+3q^3+5q^4+7q^5+11q^6+\cdots.
$$

Hence the familiar single-edge sequence

$$
1,1,2,3,5,7,11,\ldots
$$

is $p(\Delta K)$, the number of integer partitions of the momentum.  It is
not a sequence of total Schmidt ranks in successive $N_A$ blocks.

Why should a spatial entanglement spectrum see this?  A local cut through a
gapped state creates virtual boundary degrees of freedom, and
$H_E=-\log\rho_A$ is expected to be quasi-local near that boundary.  Its
long-wavelength low-energy sector can then flow to the same anomalous chiral
edge theory required by the bulk Hall topological order.  This is the
entanglement version of bulk-edge correspondence.  Microscopic velocities,
level spacings, and high entanglement levels are not universal; the symmetry
counting of a separated low branch is the useful part.

### 3.4 The two-boundary warning for this torus cut

A strip cut on a torus has **two** entanglement boundaries.  If two identical
edge oscillator towers were independent and graded by the sum of their
excitation numbers, their generating function would be

$$
Z_{\rm two\ edge}(q)=Z_{\rm edge}(q)^2
=1+2q+5q^2+10q^3+20q^4+\cdots.
$$

In an actual momentum-resolved spectrum, the two oppositely oriented
boundaries contribute momenta with opposite signs; charge/topological-sector
constraints and finite-width coupling further reorganize the branches.
Therefore even this convolution should not be applied blindly.  One must
resolve momentum parallel to the cut and identify which edge sectors are
being compared.

### 3.5 Non-CFT understanding of the spatial counting

The essential reasoning can be stated without constructing the bulk
wavefunction from a CFT:

1. A gapped, local two-dimensional ground state has entanglement dominated by
   degrees of freedom near the cut.
2. A $1/3$ Laughlin topological phase has a quantized Hall response, charge
   $e/3$ anyons, and a chiral boundary anomaly.  A purely one-dimensional
   local boundary cannot remove that anomaly while preserving charge and the
   bulk gap.
3. The minimal Abelian boundary theory is the $K=3$ chiral Luttinger liquid.
   Quantizing its density waves gives one harmonic oscillator at every
   positive momentum, hence integer partitions.

Calling the low-energy fixed point a $c=1$ CFT is an efficient description
of step 3, not an assumption that the lattice FCI bulk is conformally
invariant.  The bulk is gapped and is instead characterized at long distance
by $U(1)_3$ topological order.  CFT enters only at a gapless edge, or as an
exact mathematical construction of special trial wavefunctions.

There is also an important limitation: entanglement spectra are less robust
than quantized responses or anyon data.  The entanglement Hamiltonian can
change with the partition and microscopic deformation even while the bulk
phase stays fixed.  Edge-like counting below a stable entanglement gap is
strong supporting evidence, not a standalone definition of the phase.

### 3.6 What the present spatial CSV can and cannot establish

The current CSV resolves only $N_A$.  It does **not** store momentum
$K_\parallel$ along the cut, so it cannot organize levels by $\Delta K$ and
cannot test either the single-edge partition sequence or the appropriate
two-edge towers.  Its present uses are limited to:

- normalization and Schmidt-block consistency;
- particle-number distribution across the cut;
- A/B reflection checks for an equal bipartition;
- qualitative inspection for separated entanglement branches.

A genuine Li-Haldane-style counting test requires extending the
implementation to:

1. preserve translation parallel to the cut;
2. block-diagonalize $\rho_A$ by $(N_A,K_\parallel)$;
3. choose or construct minimally entangled ground states for that cut;
4. identify a stable low branch and compare its relative-momentum counting
   with the one- or two-edge $U(1)_3$ prediction appropriate to the geometry.

Until then, the momentum-resolved PES is the sharper counting diagnostic in
this repository.

## 4. Similarities and differences between the two countings

| question | particle ES | spatial-orbital ES |
|:--|:--|:--|
| What is subsystem A? | $N_A$ indistinguishable particles anywhere on the torus | all site orbitals in a spatial strip |
| Is a physical/virtual boundary created? | no | yes, two boundaries for the present torus strip |
| What do universal low levels count? | bulk quasihole zero modes at $(N_A,N_\phi)$ | edge oscillator/conformal towers |
| Natural quantum numbers | subsystem 2D momentum $(K_1,K_2)$ | $N_A$ and momentum $K_\parallel$ along the cut |
| Direct non-CFT derivation | clustering, parent-Hamiltonian zero modes, root patterns, cyclic combinatorics | locality near the cut, Hall anomaly, $K$-matrix/chiral density waves |
| CFT role | optional route to trial states and quasihole conformal blocks | natural low-energy language of the chiral edge |
| Main finite-size danger | wrong ground manifold or accidental PES gap | cut dependence, wrong ground-state combination, two-edge mixing, missing $K_\parallel$ |

The deep similarity is bulk-edge correspondence: for model FQH states, the
bulk quasihole counting seen by the PES and the edge counting seen by a
properly resolved orbital/real-space ES approach the same underlying
topological Hilbert-space structure in the thermodynamic limit.  They are not
the same finite-size table, because one probes particles throughout the bulk
and the other probes modes localized at a boundary.

## 5. Data-provenance warning for the `3x6` FCI

The old diagnostic driver selected the lowest state in each of three
**distinct** momentum sectors.  That is wrong when several states of the
topological manifold occupy the same sector.  For the `3x6` FCI the true
lowest manifold is

```text
(k1,k2,level) = (0,3,1), (0,3,2), (0,3,3),
```

with energy splittings `0`, `0.0160538`, and `0.0162115`, followed by a gap
of about `0.2275` to the fourth state.  The downloaded pump and PES instead
used `(0,3,1), (0,0,1), (1,3,1)`.

The toolbox now carries explicit `(sector, level)` state specifications.  The
old integer-winding `3x6` charge-pump branches are therefore a projection
artifact, not evidence against the FCI.  Likewise, only the corrected PES
rerun should be used as final evidence, even though the old PES happens to
match the 117-level admissible count.

## 6. Primary references

- Bernevig and Haldane, [*Fractional Quantum Hall States and Jack
  Polynomials*](https://arxiv.org/abs/0707.3637): root configurations,
  squeezing, clustering, and generalized exclusion principles.
- Li and Haldane, [*Entanglement Spectrum as a Generalization of Entanglement
  Entropy*](https://arxiv.org/abs/0805.0332): the orbital entanglement spectrum
  and edge-like low-level counting.
- Chandran, Hermanns, Regnault, and Bernevig,
  [*Bulk-Edge Correspondence in the Entanglement
  Spectra*](https://arxiv.org/abs/1102.2218): relation among particle,
  orbital, bulk-quasihole, and edge counting.
- Regnault and Bernevig, [*Fractional Chern
  Insulator*](https://arxiv.org/abs/1105.4867): the $(1,3)$ rule, spectral
  flow, and PES quasihole counting for a lattice $1/3$ FCI.
- Qi, [*Generic Wavefunction Description of Fractional Quantum Anomalous Hall
  States and Fractional Topological
  Insulators*](https://arxiv.org/abs/1105.4298): hybrid-Wannier mapping from a
  Chern band to Landau-level-like orbitals.
- Sterdyniak *et al.*, [*Real-Space Entanglement Spectrum of Quantum Hall
  States*](https://arxiv.org/abs/1111.2810): real-space cuts, quasihole
  counting, and edge-CFT bounds.
- Chandran, Khemani, and Sondhi, [*How Universal Is the Entanglement
  Spectrum?*](https://arxiv.org/abs/1311.2946): limitations on treating the
  detailed entanglement spectrum as a phase invariant.
