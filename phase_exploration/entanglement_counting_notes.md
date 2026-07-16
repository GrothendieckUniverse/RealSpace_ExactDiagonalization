# Entanglement spectra in this project

## What the two plots do—and do not—diagnose

This note is organized around the questions one should ask while looking at
the two entanglement-spectrum plots produced by this repository.

The short answer is:

1. The **particle entanglement spectrum (PES)** is a bulk diagnostic.  Its
   universal low band is compared with the number and momenta of Laughlin
   quasihole zero modes.
2. The present **site real-space entanglement spectrum (RSES)** makes a
   spatial strip cut.  A suitably momentum-resolved low branch can diagnose
   virtual edge physics.
3. The present RSES file is resolved only by the number of particles in the
   strip.  The number of dots in each such block is mostly ordinary
   Hilbert-space rank.  It is **not** the Li--Haldane chiral-edge counting.

Thus the physical intuition in “particle cut = quasiholes, spatial cut =
edges” is essentially correct, but it needs two qualifications.

- There are more than two named cuts in the literature.  An **orbital cut**
  divides Landau-level or hybrid-Wannier orbitals; a **real-space cut** divides
  physical positions or lattice sites; a **particle cut** divides particles.
  The first two are boundary cuts, but they are not identical constructions.
- A boundary cut is not an edge diagnostic merely because a Schmidt
  decomposition was performed.  One must resolve momentum parallel to the
  cut and identify the universal low-lying branch.  Raw Schmidt rank has no
  such meaning.

The most important distinction in this note is therefore

| quantity being counted | physical content | sequence in this project |
|:--|:--|:--|
| all nonzero spatial Schmidt values in an $N_A$ block | kinematic Fock-space rank, generally nonuniversal | e.g. `1,12,66,153,18,1` for `3x5` |
| low PES levels at fixed particle-cut size $N_A$ | bulk quasihole zero modes | 42, 75, 117 for the current `3x4`, `3x5`, `3x6` data at $N_A=2$ |
| low states versus excess edge momentum at fixed edge charge | chiral edge descendants | $1,1,2,3,5,7,11,\ldots$ for one $U(1)$ edge |

Only the last two rows can provide the topological counting tests discussed
below.

---

## 1. Common language

### 1.1 Schmidt values and entanglement energies

For a pure state and a tensor-product bipartition
$\mathcal H=\mathcal H_A\otimes\mathcal H_B$,

$$
|\Psi\rangle
 =\sum_\alpha s_\alpha
   |\alpha\rangle_A|\alpha\rangle_B,
\qquad
\lambda_\alpha=s_\alpha^2.
$$

The reduced density matrix is

$$
\rho_A=\operatorname{Tr}_B|\Psi\rangle\langle\Psi|,
$$

and the entanglement energies are

$$
\xi_\alpha=-\log\lambda_\alpha.
$$

The low part of the entanglement plot means the **largest** eigenvalues of
$\rho_A$.  Absolute entanglement energies, spacings, velocities, and the
high-energy part of the spectrum are generally nonuniversal.  The useful
fingerprint is a separated low band with the correct symmetry-resolved
counting.

### 1.2 The symbol $N_A$ means two different things in the two calculations

This is a frequent source of confusion.

- In a **PES**, one first chooses exactly $N_A$ retained particles.  Changing
  $N_A$ constructs a different reduced density matrix, with its own trace-one
  normalization and its own arbitrary additive offset in $\xi$.
- In a **spatial ES**, region A is fixed in real space.  Its particle number
  fluctuates, so $N_A=0,1,\ldots,N$ labels charge blocks inside one and the
  same $\rho_A$.

Consequently, a sweep of PES cut size should normally be shown as separate
panels, one panel per $N_A$.  By contrast, plotting all spatial $N_A$ blocks
on one axis is meaningful because their relative probability weights belong
to one density matrix.

---

## 2. Particle entanglement spectrum

### 2.1 Definition and why the full ground manifold is used

The particle cut retains $N_A$ indistinguishable particles anywhere in the
sample and traces out $N_B=N-N_A$ particles.  In first-quantized notation,

$$
\rho_A(Z_A,Z_A')
 =\binom{N}{N_A}
  \int dZ_B\,
  \Psi(Z_A,Z_B)\Psi^*(Z_A',Z_B).
$$

It creates no physical boundary.  Its natural quantum numbers on a
translation-invariant torus are the two subsystem momenta $(K_1,K_2)$.

For a topologically ordered state on a torus, the ground space has dimension
$d$.  The appropriate, basis-independent input is the normalized projector
onto that space,

$$
\rho_{\rm GS}
 =\frac1d\sum_{a=1}^d|\Psi_a\rangle\langle\Psi_a|,
\qquad
\rho_A^{\rm PES}=\operatorname{Tr}_{N_B}\rho_{\rm GS}.
$$

Any unitary rotation of the $d$ ground states leaves $\rho_{\rm GS}$
unchanged.  This matters particularly for `3x6`, where the three current FCI
manifold states are levels 1, 2, and 3 of the same momentum sector $(0,3)$.

### 2.2 What can and cannot imply the $(1,3)$ rule

The rule

$$
\boxed{\text{no more than one particle in any three consecutive orbitals}}
$$

is equivalently

$$
n_j+n_{j+1}+n_{j+2}\le 1,
$$

with cyclic indices on a torus.  Its densest occupation patterns are the
three translations of

```text
100 100 100 100 ...
```

There is an essential logical point:

> Chern number $C=1$ and filling $\nu=1/3$ do not, by themselves, prove the
> $(1,3)$ rule.

The same band and filling can support a charge-density wave, a compressible
state, or another correlated phase.  Even the long-distance $U(1)_3$
topological field theory fixes anyon charges, braiding, Hall response, and
torus degeneracy, but it does not by itself contain the finite set of
guiding-center orbitals whose 117 momentum-resolved quasihole states we count
in `3x6`.

The $(1,3)$ counting becomes the correct target after specifying the
**fermionic Laughlin zero-mode/clustering class**, or after hypothesizing
that the lattice state is adiabatically connected to that class.  The PES
then tests that hypothesis.  It does not assume the conclusion merely from
the filling.

The next three subsections give complementary routes to the rule.

### 2.3 Parent-Hamiltonian and clustering route: no CFT required

Project two-dimensional fermions into one Landau level and decompose their
interaction into relative-angular-momentum channels.  Fermi statistics
allows only odd pair angular momenta.  The $V_1$ Haldane pseudopotential is a
positive operator that penalizes the closest allowed pair channel,

$$
H_{V_1}=V_1\sum_{i<j}P_{ij}(L_{\rm rel}=1),
\qquad V_1>0.
$$

A zero mode must have no relative-angular-momentum-one component.  An
antisymmetric wavefunction already vanishes at least linearly when
$z_i\to z_j$; eliminating the linear channel forces the next possible odd
power,

$$
\Psi\sim(z_i-z_j)^3.
$$

The word “zero” now occurs in two different senses that must not be mixed up.

- A **coincidence zero** is the cubic vanishing of the first-quantized
  wavefunction as two particles meet.
- A **zero mode** is a many-body eigenstate with exactly zero energy under
  the positive-semidefinite $V_1$ Hamiltonian.
- A `0` in an occupation string merely says that one chosen guiding-center
  orbital is empty.  Such individual occupation zeros are not conserved
  under a change of torus aspect ratio.

At the largest possible density, one third, the cubic clustering condition
has a three-dimensional torus zero-mode space.  At lower density it has a
larger quasihole zero-mode space.

#### What the thin-torus limit actually says

Choose Landau-gauge orbitals, which are localized guiding-center strips, and
make one circumference very small.  After an appropriate rescaling, the
dominant terms of the projected interaction are diagonal guiding-center
repulsions.  The dominant root configurations of the densest zero modes then
have two empty guiding centers between occupied ones:

```text
...100100100...
```

The three translations of this string give the three root patterns at exact
$1/3$ filling.  With extra orbitals, every string satisfying

$$
n_j+n_{j+1}+n_{j+2}\le1
$$

is an allowed quasihole root.  Domain walls in these strings are the
thin-cylinder form of charge-$e/3$ quasiholes.

The phrase “the state follows the $(1,3)$ rule” must be interpreted carefully.

- **In the strict root/thin-torus description:** the occupation string itself
  obeys the rule.
- **At finite aspect ratio:** the exact zero mode is generally a superposition
  of the root and configurations obtained from it by pair squeezing.  Many
  Slater determinants in that superposition do not literally obey the
  three-orbital spacing rule.
- **What remains true:** the independent zero modes can be labeled by
  $(1,3)$-admissible roots, and their number and translation sectors are
  obtained by counting those roots.

Thus $(1,3)$ does not mean that the isotropic wavefunction is the single
product state $|100100\ldots\rangle$.  At exact filling that string is only
the dominant root of one thin-torus sector; its translations give the other
two roots.

#### Why the zero-mode count can remain constant

Let $H(s)$ continuously change the guiding-center metric or torus aspect
ratio, with $s=0$ thin and $s=1$ isotropic.  Assume throughout the path that

1. $H(s)\ge0$ remains the same $V_1$ parent-Hamiltonian class;
2. translations remain symmetries; and
3. a nonzero spectral gap separates the zero-energy space from positive
   energies.

Let $P_0(s)$ be the spectral projector onto zero energy.  Its rank is an
integer.  It cannot change continuously; a change would require a positive
energy level to reach zero or a zero mode to leave zero, which closes the
assumed zero-mode gap.  Similarly, if $P_K$ projects onto translation sector
$K$, then

$$
d_K(s)=\operatorname{Tr}[P_0(s)P_K]
$$

is an integer and cannot change without a symmetry-preserving level crossing
at zero.  The eigenvectors and their occupation coefficients can vary
strongly, but the zero-mode dimension $\sum_Kd_K$ and its momentum allocation
remain fixed.

This is the precise content of the adiabatic argument.  It protects a
subspace, not a binary string.

Strictly taking the circumference to zero can discard exponentially small
off-diagonal terms and may introduce accidental degeneracies in the truncated
electrostatic Hamiltonian.  The safe continuity statement starts at a small
but finite circumference with the full pseudopotential, or supplements the
thin-limit intuition with the exact clustering/zero-mode algebra.

The statement is also special to the parent-Hamiltonian path.  Under a
generic deformation within the same topological phase, the physical energy
ground state need not remain an exact $V_1$ zero mode and its coincidence
zeros need not remain exactly cubic.  What can survive is an adiabatically
identifiable Laughlin-like low PES band.  Exact polynomial zeros are stronger
microscopic structure than topological order alone.

There is also an important logical limitation: adiabaticity by itself does
not prove that the protected dimension equals the number of admissible
strings.  That correspondence comes from the additional zero-mode algebra:
clustering forces a triangular structure in the occupation basis, each
zero mode has an admissible dominant pattern, and the admissible patterns
generate independent zero modes.  Jack polynomials provide one explicit
version of this construction; purely second-quantized dominance-pattern
proofs provide another.

There are also fully second-quantized derivations: the pseudopotential can be
treated as a frustration-free guiding-center lattice Hamiltonian, and all
zero modes can be constructed algebraically from dominance patterns.  No
analytic trial polynomial or CFT is required for that formulation.

### 2.4 Jack-polynomial route, with the fermionic shift made explicit

The Jack-polynomial language is not a separate physical assumption.  It is a
special basis of symmetric polynomials in which three facts become
simultaneously visible:

1. one dominant occupation pattern labels the polynomial;
2. all other monomials lie below it in a precise squeezing order; and
3. at special negative Jack parameter, the polynomial has the clustering
   zeros required by an FQH parent Hamiltonian.

The definitions come first.

#### Partitions and occupation configurations

A partition of total degree $M$ into at most $N$ parts is an ordered list

$$
\lambda=(\lambda_1,\lambda_2,\ldots,\lambda_N),
\qquad
\lambda_1\ge\lambda_2\ge\cdots\ge\lambda_N\ge0,
\qquad
\sum_i\lambda_i=M.
$$

In a lowest-Landau-level polynomial, $\lambda_i$ is a one-particle orbital
angular momentum.  For bosons, repeated parts are allowed.  The equivalent
occupation numbers are

$$
n_m(\lambda)=\#\{i:\lambda_i=m\}.
$$

For example,

$$
\lambda=(4,2,0)
\quad\longleftrightarrow\quad
n_0n_1n_2n_3n_4=10101.
$$

The monomial symmetric polynomial $m_\lambda$ is obtained by symmetrizing
$z_1^{\lambda_1}\cdots z_N^{\lambda_N}$ over distinct permutations of the
exponents.  The collection of all $m_\lambda$ is an ordinary basis of
symmetric homogeneous polynomials.

For fermions, antisymmetry requires all occupied orbitals to be distinct, so
the corresponding partition is strict,

$$
\mu_1>\mu_2>\cdots>\mu_N\ge0,
$$

and labels a Slater determinant rather than a symmetric monomial.

#### What a Jack polynomial is

For every partition $\lambda$ and parameter $\alpha$, the monic Jack
polynomial $J_\lambda^{(\alpha)}$ is a symmetric homogeneous polynomial with
the triangular expansion

$$
J_\lambda^{(\alpha)}
=m_\lambda
 +\sum_{\kappa\prec\lambda}
   c_{\lambda\kappa}(\alpha)m_\kappa.
$$

It is uniquely fixed by this leading term together with being an eigenfunction
of the Jack/Calogero--Sutherland Laplace--Beltrami operator.  In one common
convention that operator is

$$
\mathcal D_\alpha
=\sum_i(z_i\partial_i)^2
 +\frac1\alpha\sum_{i<j}
   \frac{z_i+z_j}{z_i-z_j}
   (z_i\partial_i-z_j\partial_j).
$$

Different normalizations of $J$ or $\mathcal D_\alpha$ do not change the root
partition or the squeezing relations used here.

The symbol $\kappa\prec\lambda$ means that $\kappa$ is dominated by
$\lambda$:

$$
\sum_{i=1}^{r}\kappa_i
\le
\sum_{i=1}^{r}\lambda_i
\quad\text{for every }r,
\qquad
\sum_i\kappa_i=\sum_i\lambda_i.
$$

Equivalently, $\kappa$ can be reached through repeated inward pair squeezes.

#### What squeezing means

Choose two occupied orbitals $a<b$ and move them toward one another,

$$
(a,b)\longrightarrow(a+s,b-s),
\qquad
0<s\le\frac{b-a}{2}.
$$

This preserves particle number and total angular momentum.  For fermions the
new orbitals must remain distinct; for bosons they may coincide.  For example,
a bosonic pair in orbitals 0 and 4 can squeeze to 1 and 3 or to two particles
in orbital 2:

```text
10001  ->  01010  ->  00200.
```

The root configuration is the unique maximal or unsqueezed partition
$\lambda$ multiplying $m_\lambda$ with coefficient one.  It is not the whole
wavefunction.  The Jack is the root plus a definite linear combination of
its squeezed descendants.

#### What clustering means for the Laughlin Jack

For special negative rational values

$$
\alpha=-\frac{k+1}{r-1},
$$

well-defined Jacks with $(k,r)$-admissible roots have an FQH clustering
property.  The admissibility condition can be written

$$
\lambda_i-\lambda_{i+k}\ge r,
$$

equivalently: no more than $k$ particles occupy any $r$ consecutive orbitals.

For the bosonic Laughlin state, $k=1$, $r=2$, and $\alpha=-2$.  Its densest
root obeys

$$
\lambda_i-\lambda_{i+1}\ge2
$$

and is

```text
1010101...
```

The associated Jack vanishes quadratically when two coordinates coincide.
At the densest flux it is exactly

$$
J_{\lambda_{1/2}}^{(-2)}(z)
\propto\prod_{i<j}(z_i-z_j)^2.
$$

With additional orbitals, all $(1,2)$-admissible roots label the bosonic
Laughlin quasihole zero modes; their Jacks contain the squeezed monomials
needed to satisfy clustering.

#### Why the fermionic $1/3$ problem is written as a bosonic Jack times a
Vandermonde

Every antisymmetric polynomial is divisible by the Vandermonde determinant

$$
\Delta(z)=\prod_{i<j}(z_i-z_j)
=\det[z_i^{N-j}].
$$

After dividing by $\Delta$, the quotient is symmetric.  This is a unique
algebraic conversion: exchanging two coordinates changes the sign of both
the original polynomial and $\Delta$, so the quotient is invariant.  It is
not necessary to interpret $\Delta$ as a composite-fermion wavefunction.
It simply supplies the minimum zero required by Fermi statistics.

For the fermionic $1/3$ Laughlin polynomial,

$$
\Psi_{1/3}(z)
=\prod_{i<j}(z_i-z_j)^3
=\Delta(z)
 \underbrace{\prod_{i<j}(z_i-z_j)^2}_{
 J_{\lambda_{1/2}}^{(-2)}(z)}.
$$

Thus the symmetric quotient happens to be the bosonic $1/2$ Laughlin Jack.
This factorization is canonical because one Vandermonde is the universal
antisymmetry factor.  Refactoring the same polynomial cannot change its
intrinsic dominant partition or zero-mode dimension.  A purported
factorization that gives different counting has either changed the state or
misidentified the dominant term.

#### The staircase and why it does not disappear

Expanding the determinant shows that its orbital exponents are permutations
of

$$
\delta=(N-1,N-2,\ldots,1,0).
$$

This is called the staircase partition.  When a symmetric root $\lambda$ is
multiplied by $\Delta$, the leading fermionic Slater partition is

$$
\boxed{\mu=\lambda+\delta}
$$

component by component.  The staircase has not vanished: it converts the
bosonic weak inequalities into the strict inequalities required by Pauli
exclusion and adds one unit to every adjacent orbital spacing,

$$
\mu_i-\mu_{i+1}
=(\lambda_i-\lambda_{i+1})
 +(\delta_i-\delta_{i+1})
\ge2+1=3.
$$

For $N=3$, the complete root conversion is

$$
\lambda=(4,2,0)
\quad\longleftrightarrow\quad10101,
$$

$$
\delta=(2,1,0),
$$

$$
\mu=\lambda+\delta=(6,3,0)
\quad\longleftrightarrow\quad1001001.
$$

This is the fermionic $(1,3)$ root.  Sometimes fermionic partitions are
reported only after subtracting the universal staircase; in that convention
$\delta$ becomes invisible in the notation, but not in the physics or the
orbital spacing.

The conclusion is again about roots, not every Slater determinant.  The
fermionic Laughlin state and its quasihole zero modes contain squeezed
Slaters that can violate literal three-orbital spacing.  The independent
zero modes are indexed by the unsqueezed fermionic roots satisfying $(1,3)$.

Jack polynomials are most directly formulated on the plane, disk, or sphere,
where lowest-Landau-level states are ordinary polynomials.  A torus state is
built from theta functions and has an additional center-of-mass structure,
so it is not literally the same Jack polynomial.  What carries over is the
local clustering condition and its guiding-center root rule.  Periodic
boundary conditions turn the open $(1,3)$ rule into the cyclic rule counted
in Section 2.7, while magnetic translations determine the torus momentum
sectors.

#### How Jack ideas enter the PES argument

The three pieces now have distinct jobs.

- **Clustering** identifies the allowed vector space: after a particle cut,
  every conditional A-particle function still vanishes cubically when two A
  coordinates meet, so it lies in the Laughlin quasihole zero-mode space.
- **Admissible roots** count a basis of that zero-mode space.  For fermions
  they satisfy $(1,3)$.
- **Squeezing** explains why the actual conditional functions and PES
  eigenvectors are not single admissible strings.  Each is a superposition
  of root configurations and their squeezed descendants.

The Jack construction therefore does not assert “one PES eigenvalue equals
one product state.”  It supplies a dimension and symmetry-sector count for
the vector space in which the ideal PES has support.  Section 2.6 derives
that support/rank statement directly.

### 2.5 Why the lattice FCI inherits this as a diagnostic

For the checkerboard lattice,

$$
N_{\rm cell}=L_1L_2,
\qquad
N_{\rm site}=2L_1L_2.
$$

There are two microscopic sites per unit cell, but one isolated $C=1$ band
contains only one one-particle state per unit cell.  The Landau-level orbital
number is therefore represented by

$$
N_\phi\longleftrightarrow N_{\rm cell}=L_1L_2,
$$

not by the number of sites.  In this project,

$$
\frac{N}{N_{\rm site}}=\frac16,
\qquad
\nu_{\rm band}=\frac{N}{N_{\rm cell}}=\frac13.
$$

On a torus there is no spherical shift, so the commensurate $1/3$ ground
state has

$$
N_\phi=3N.
$$

Hybrid-Wannier orbitals provide a one-dimensional guiding-center ordering in
a $C=1$ band.  When transverse momentum winds around the Brillouin zone, the
Wannier center shifts by one unit cell; that Chern pump is the lattice
counterpart of Landau-gauge orbital flow.  Continuum FQH momentum sectors
then fold into the finite lattice Brillouin zone.

For an ideal band-projected parent Hamiltonian, $(1,3)$ is an exact zero-mode
rule.  For the unprojected, generic checkerboard Hamiltonian used here, it is
not an exact constraint on microscopic site occupations.  Instead, a low PES
band can remain adiabatically connected to the Laughlin quasihole space while
generic levels appear above it.  The entanglement gap separates these two
sets when the finite system is favorable.

### 2.6 Why a particle cut counts quasiholes

This connection is easiest to see in four steps: what the labels A and B
mean, why the flux is unchanged, why conditional A states are quasihole zero
modes, and why their span is the rank of $\rho_A$.

#### Step 1: A and B label particles, not spatial regions

Write

$$
Z_A=(z_1,\ldots,z_{N_A}),
\qquad
Z_B=(w_1,\ldots,w_{N_B}),
\qquad
N_A+N_B=N.
$$

The particles are indistinguishable; choosing the first $N_A$ coordinate
labels is only bookkeeping.  No particle is confined to a region A.  In the
reduced-density-matrix kernel,

$$
\rho_A(Z_A,Z_A')
=\binom{N}{N_A}\int dZ_B\,
 \Psi(Z_A,Z_B)\Psi^*(Z_A',Z_B),
$$

the $Z_B$ variables are integrated over, while $Z_A$ and $Z_A'$ remain as
the arguments of the operator.  “Retained coordinates” means precisely
these unintegrated $N_A$ particle variables.

#### Step 2: tracing particles does not change the one-particle flux

Before the trace, every coordinate belongs to the same one-particle
Landau-level Hilbert space

$$
\mathcal h_{N_\phi},
\qquad
\dim\mathcal h_{N_\phi}=N_\phi
$$

on a torus.  The full state lies in

$$
\bigwedge^N\mathcal h_{N_\phi}.
$$

Tracing out coordinate labels removes particles from the density matrix; it
does not change the external magnetic field, the boundary conditions, or the
number of one-particle orbitals available to the remaining labels.  The
particle-A Hilbert space is therefore

$$
\bigwedge^{N_A}\mathcal h_{N_\phi},
$$

not $\bigwedge^{N_A}\mathcal h_{3N_A}$.

For an original torus Laughlin ground state $N_\phi=3N$.  Relative to the
densest $N_A$-particle Laughlin state, which would need only $3N_A$ torus
orbitals, the retained particles see

$$
N_\phi-3N_A=3(N-N_A)=3N_B
$$

excess flux quanta.  An $N_A$-particle state at this excess flux belongs to a
many-quasihole sector.

#### Step 3: fixing $Z_B$ produces a quasihole zero mode of the A particles

For every fixed value $W=Z_B$, define the conditional A-particle function

$$
\Phi_W(Z_A)=\Psi(Z_A,W).
$$

This is not a normalized post-measurement state; it is simply a vector in the
A-particle Hilbert space whose coefficients depend on $W$.

The factorization is explicit for the planar Laughlin polynomial:

$$
\Psi_N(Z_A,W)
=\underbrace{\prod_{i<j\in A}(z_i-z_j)^3}_{\text{A--A clustering}}
 \underbrace{\prod_{p<q\in B}(w_p-w_q)^3}_{\text{constant for fixed }W}
 \underbrace{\prod_{i\in A,p\in B}(z_i-w_p)^3}_{\text{extra zeros seen by A}}.
$$

For fixed $W$, the second factor is an irrelevant scalar.  The first factor
guarantees a cubic zero whenever two retained particles meet.  The last
factor is symmetric in the A coordinates and inserts three zeros at every
$w_p$ for each retained particle.  Since one fundamental Laughlin quasihole
at position $\eta$ contributes $\prod_i(z_i-\eta)$, each fixed B coordinate
acts here like three coincident fundamental quasiholes.  There are therefore
$3N_B$ excess quasihole fluxes, exactly as found from flux counting above.

On a torus, ordinary differences are replaced by the appropriate theta
functions and a center-of-mass factor.  The conclusion is unchanged:
$\Phi_W$ has the original $N_\phi$ quasiperiodic boundary conditions and the
same cubic A--A coincidence zeros.  Hence

$$
\Phi_W\in\mathcal H_{\rm qh}(N_A,N_\phi),
$$

the $N_A$-particle zero-mode space of the $V_1$ parent Hamiltonian at the
original flux.

Notice what has and has not been shown.  We did not say that $\Phi_W$ is one
admissible occupation string.  It is generally a superposition of squeezed
configurations.  We showed that it belongs to the zero-mode vector space
whose independent basis states are labeled by admissible roots.

#### Step 4: the support of $\rho_A$ is the span of the conditional states

Using the conditional vectors, the reduced density matrix is a continuous
Gram operator,

$$
\rho_A
\propto\int dW\,|\Phi_W\rangle\langle\Phi_W|.
$$

Therefore

$$
\operatorname{supp}\rho_A
=\operatorname{span}\{\Phi_W:W\text{ varies}\}.
$$

The finite-basis version makes the rank statement elementary.  Expand

$$
|\Psi\rangle
=\sum_{a=1}^{D_A}\sum_{b=1}^{D_B}
 M_{ab}|a\rangle_A|b\rangle_B.
$$

For a fixed B-basis state $b$, the $b$th column

$$
|\Phi_b\rangle_A=\sum_aM_{ab}|a\rangle_A
$$

is the discrete version of a conditional A state.  Tracing B gives

$$
\rho_A=MM^\dagger.
$$

For any vector $v$,

$$
v^\dagger\rho_Av=\|M^\dagger v\|^2.
$$

Thus $\ker\rho_A=\ker M^\dagger$ and

$$
\boxed{
\operatorname{rank}\rho_A
=\operatorname{rank}M
=\dim\operatorname{span}\{|\Phi_b\rangle_A\}.}
$$

Since every column of the ideal Laughlin particle-entanglement matrix lies
in the quasihole space,

$$
\boxed{
\operatorname{rank}\rho_A^{\rm ideal}
\le
D_{\rm qh}(N_A,N_\phi)
=\dim\mathcal H_{\rm qh}(N_A,N_\phi).}
$$

Sector by sector, the same argument gives

$$
\operatorname{rank}\rho_A^{\rm ideal}(K)
\le D_{\rm qh}(K).
$$

Equality requires an additional spanning statement: as the B configuration
varies, the conditional states must span the entire allowed quasihole space.
This is not a consequence of the trace definition alone.  It is found for
the standard Laughlin PES in the usual $N_A\le N/2$ regime, barring accidental
linear dependencies, but proving finite-size saturation in complete
generality is subtler than proving the support bound.

For the equal mixture of $d$ torus ground states, concatenate their particle
entanglement matrices,

$$
\mathcal M=\frac1{\sqrt d}
 [M^{(1)}\;M^{(2)}\;\cdots\;M^{(d)}].
$$

Then

$$
\rho_A=\mathcal M\mathcal M^\dagger,
$$

so the support is the span of conditional states from the whole ground
manifold.  This explains why using the full torus projector is important.

#### Ideal model state versus the generic lattice FCI

For the exact parent-Hamiltonian Laughlin state, clustering gives an exact
support restriction.  When saturation holds, the number of nonzero PES
eigenvalues equals the quasihole dimension.

For a generic, unprojected lattice FCI, exact cubic clustering is absent in
the microscopic site basis.  Its reduced density matrix can have full rank.
The Laughlin descendant space is then identified not with the entire support
but with a low-entanglement-energy band separated from generic levels.  In
the current `3x6`, $N_A=2$ data,

$$
\operatorname{rank}\rho_A=\binom{36}{2}=630,
$$

while the quasihole prediction is 117.  The diagnostic is that exactly 117
levels, with the correct momenta, occur below the entanglement gap—not that
the full lattice density matrix has rank 117.

### 2.7 Total $(1,3)$ count on a torus

This section counts admissible **root occupation strings**.  The variables
$g_i$ and $y_i$ are not new physical conditions.  They are a convenient way
to rewrite the same spacing rule on a periodic ring.

To keep the notation light, write

$$
n=N_A,
\qquad
L=N_\phi.
$$

Thus we place $n$ particles in $L$ labeled orbitals numbered
$0,1,\ldots,L-1$, with orbital labels understood modulo $L$.

#### From $\mu_i-\mu_{i+1}\ge3$ to the cyclic torus rule

First consider an open line of orbitals.  Let the occupied positions in
increasing order be

$$
0\le x_1<x_2<\cdots<x_n\le L-1.
$$

The $(1,3)$ rule requires consecutive particles to be at least three orbital
labels apart:

$$
x_{i+1}-x_i\ge3.
$$

For example, particles at orbitals 2 and 5 are allowed,

```text
orbital:   2 3 4 5
           1 0 0 1
```

whereas particles at 2 and 4 are forbidden because only one empty orbital
lies between them.

The Jack-partition convention lists occupied powers in decreasing order,

$$
\mu=(\mu_1,\ldots,\mu_n)=(x_n,x_{n-1},\ldots,x_1).
$$

Therefore

$$
\mu_i-\mu_{i+1}\ge3
$$

is exactly the same condition as $x_{i+1}-x_i\ge3$ for neighboring particles
that do not cross the end of the orbital list.

On a torus, that is not yet enough.  Orbital $L-1$ is adjacent to orbital 0,
so the last particle must also be separated from the first across the
periodic boundary:

$$
x_1+L-x_n\ge3.
$$

In descending partition notation, this additional condition is

$$
L-\mu_1+\mu_n\ge3.
$$

Thus the complete torus rule is

$$
\boxed{
x_{i+1}-x_i\ge3\quad(i=1,\ldots,n-1),
\qquad
x_1+L-x_n\ge3.}
$$

Equivalently, every cyclic block of three orbitals contains at most one
particle.  The variables $g_i$ simply build the wraparound condition in from
the beginning.

#### Definition of $g_i$: actual empty orbitals between particles

Walk clockwise around the ring from one particle to the next.  Define

$$
g_i=\text{number of empty orbitals strictly between particle $i$ and
particle $i+1$}.
$$

For the increasing positions above,

$$
g_i=x_{i+1}-x_i-1,
\qquad i=1,\ldots,n-1,
$$

and the wraparound gap is

$$
g_n=x_1+L-x_n-1.
$$

The local picture is

```text
particle i                       particle i+1
    1   0   0   [extra zeros]         1
        <-------- g_i -------->
```

A separation of three orbital labels means two empty orbitals in between, so

$$
x_{i+1}-x_i\ge3
\quad\Longleftrightarrow\quad
g_i\ge2.
$$

This includes $g_n$, and is therefore the full cyclic $(1,3)$ rule.

The $n$ occupied orbitals and the $n$ gaps partition the entire ring.  There
are $n$ occupied orbitals and $L-n$ empty ones, hence

$$
\boxed{\sum_{i=1}^{n}g_i=L-n.}
$$

#### Definition of $y_i$: remove the two mandatory zeros

Every gap must contain at least two zeros.  Separate those mandatory zeros
from any additional zeros by writing

$$
g_i=2+y_i,
\qquad
y_i\ge0.
$$

Thus $y_i$ is simply the number of **extra** empty orbitals in gap $i$ beyond
the two required by $(1,3)$:

```text
1  0 0  [y_i additional zeros]  1.
```

Substitute $g_i=2+y_i$ into the sum of all gaps:

$$
\sum_i(2+y_i)=L-n.
$$

Since there are $n$ gaps,

$$
2n+\sum_i y_i=L-n,
$$

and therefore

$$
\boxed{\sum_{i=1}^{n}y_i=L-3n.}
$$

The quantity

$$
Q=L-3n
$$

is the number of orbitals left over after assigning each particle one
occupied orbital plus its two mandatory following zeros.  At exact Laughlin
filling $L=3n$, $Q=0$.  With quasihole flux, $L>3n$, and the $Q$ extra zeros
can be distributed among the $n$ cyclic gaps.

#### What the binomial coefficient counts

We must count the ordered nonnegative integer solutions of

$$
y_1+y_2+\cdots+y_n=Q.
$$

This is the stars-and-bars problem.  Represent the $Q$ extra zeros by $Q$
identical stars and use $n-1$ bars to divide them among the $n$ gaps.  For
example,

```text
*** | * | **
```

represents $(y_1,y_2,y_3)=(3,1,2)$.  There are $Q+n-1$ symbols in total, and
choosing the positions of the $n-1$ bars gives

$$
\#\{(y_1,\ldots,y_n)\}
=\binom{Q+n-1}{n-1}.
$$

Using $Q=L-3n$,

$$
\boxed{
\#\{\text{ordered gap patterns after marking a starting particle}\}
=\binom{L-2n-1}{n-1}.}
$$

This binomial does **not yet** count unmarked occupation strings.  It counts
the ordered distributions of the extra zeros after one particle has been
distinguished as the starting particle.

#### Why multiply by $L$ and divide by $n$

Choose the orbital $s$ of the distinguished starting particle.  There are
$L$ choices.  Given $s$ and an ordered tuple $(y_1,\ldots,y_n)$, walk around
the ring, putting two mandatory zeros plus $y_i$ extra zeros after particle
$i$.  This uniquely constructs a cyclic occupation string with one marked
particle.

Hence the number of marked admissible strings is

$$
L\binom{L-2n-1}{n-1}.
$$

Every physical string contains $n$ particles, and any one of them could have
been marked as the starting particle.  Each unmarked string has therefore
been counted exactly $n$ times.  Dividing by $n$ gives

$$
\boxed{
\mathcal N_{(1,3)}^{\rm torus}(L,n)
=\frac{L}{n}
 \binom{L-2n-1}{n-1}.}
$$

Restoring $n=N_A$,

$$
\boxed{
\mathcal N_{(1,3)}^{\rm torus}(L,N_A)
=\frac{L}{N_A}
 \binom{L-2N_A-1}{N_A-1}.}
$$

This holds for $N_A>0$ and $L\ge3N_A$.  Define the $N_A=0$ count to be one.

#### Check 1: exact filling gives three translated roots

If $L=3n$, then $Q=0$ and the only gap tuple is

$$
(y_1,\ldots,y_n)=(0,\ldots,0).
$$

Every gap contains exactly two zeros.  The formula gives

$$
\mathcal N_{(1,3)}^{\rm torus}(3n,n)
=\frac{3n}{n}\binom{n-1}{n-1}=3.
$$

These are exactly the three translations

```text
100100100...
010010010...
001001001...
```

on the ring.

#### Check 2: `3x4`, $L=12$, $N_A=2$

Here

$$
Q=L-3N_A=12-6=6,
\qquad
y_1+y_2=6.
$$

The seven ordered solutions are

$$
(0,6),(1,5),(2,4),(3,3),(4,2),(5,1),(6,0),
$$

which agrees with

$$
\binom{6+2-1}{2-1}=\binom71=7.
$$

If the marked particle is at orbital 0, these solutions place the second
particle at orbitals $3,4,5,6,7,8,9$, respectively.  Orbitals 1, 2, 10, and
11 are forbidden because they lie within two steps of orbital 0 around the
ring.

There are 12 choices for the marked particle's orbital and two possible
marked particles in every unmarked pair, so

$$
\mathcal N_{(1,3)}^{\rm torus}(12,2)
=\frac{12\times7}{2}=42.
$$

#### Check 3: `3x6`, $L=18$

For $N_A=2$,

$$
Q=18-6=12,
\qquad
y_1+y_2=12.
$$

There are $\binom{13}{1}=13$ marked gap patterns, and therefore

$$
\mathcal N_{(1,3)}^{\rm torus}(18,2)
=\frac{18}{2}\binom{13}{1}
=117.
$$

For $N_A=3$,

$$
Q=18-9=9,
\qquad
y_1+y_2+y_3=9.
$$

Stars and bars gives

$$
\binom{9+3-1}{3-1}=\binom{11}{2}=55
$$

marked gap patterns, and hence

$$
\mathcal N_{(1,3)}^{\rm torus}(18,3)
=\frac{18}{3}\times55
=330.
$$

#### What this count means for the PES

For the $V_1$ Laughlin parent problem, the admissible-root construction gives

$$
D_{\rm qh}(N_A,L)
=\dim\mathcal H_{\rm qh}(N_A,L)
=\mathcal N_{(1,3)}^{\rm torus}(L,N_A).
$$

This is what the combinatorial formula counts: the dimension of the ideal
quasihole zero-mode space.  Combining it with Section 2.6 gives the precise
chain of statements

$$
\underbrace{\mathcal N_{(1,3)}}_{\text{admissible roots}}
=\underbrace{D_{\rm qh}}_{\text{ideal zero-mode dimension}}
\ge\underbrace{\operatorname{rank}\rho_A^{\rm ideal}}_{
\text{conditional states that are actually spanned}}.
$$

When the ideal particle cut saturates the quasihole space, the inequality is
an equality.  In a generic FCI, $\mathcal N_{(1,3)}$ instead predicts the
number of levels in the low PES band; the total microscopic rank can be much
larger.

For $N_A=2$,

$$
\mathcal N_{(1,3)}^{\rm torus}(L,2)=\frac{L(L-5)}2.
$$

This also follows by fixing one particle: the other cannot occupy the same
orbital or either of the two orbitals on each side, leaving $L-5$ choices.
Divide the ordered count by two.

For comparison, on a sphere or open orbital chain there is no cyclic
wraparound.  The number of length-$L$ strings with $N_A$ ones separated by at
least two zeros is

$$
\mathcal N_{(1,3)}^{\rm open}(L,N_A)
=\binom{L-2(N_A-1)}{N_A}.
$$

On the sphere one must additionally use $L=N_\phi+1$ and the Laughlin shift
$N_\phi=3(N-1)$.  The torus formula, not the open formula, is the one relevant
to this repository.

### 2.8 From geometry and filling to a predicted count

For a torus calculation in this checkerboard model:

1. Compute $L=L_1L_2$, the number of unit cells and of states in one Chern
   band.
2. Check the band filling $N/L$.  At the commensurate fermionic Laughlin
   point, $L=3N$.
3. Choose a particle cut, normally $1\le N_A\le\lfloor N/2\rfloor$.
4. Insert $L$ and $N_A$ into the cyclic formula above.
5. Enumerate the admissible strings to resolve the total among lattice
   momenta.

For the study geometries, the clean half-cut range gives

| geometry | $L=N_{\rm cell}$ | $N$ | predicted low PES count(s) |
|:--|--:|--:|:--|
| `3x4` | 12 | 4 | $N_A=1:12$; $N_A=2:42$ |
| `3x5` | 15 | 5 | $N_A=1:15$; $N_A=2:75$ |
| `3x6` | 18 | 6 | $N_A=1:18$; $N_A=2:117$; $N_A=3:330$ |

The PES implementation works in the full two-site lattice basis.  Its raw
$N_A=2$ dimensions are therefore $\binom{24}{2}=276$,
$\binom{30}{2}=435$, and $\binom{36}{2}=630$.  The much smaller numbers 42,
75, and 117 predict the low quasihole-like band, not the entire microscopic
PES.

### 2.9 Momentum-resolved counting for the present convention

For the geometries in this project, label the effective Chern-band orbitals
by

$$
j=k_1+L_1k_2,
\qquad
0\le k_1<L_1,
\qquad
0\le k_2<L_2.
$$

For every cyclic $(1,3)$-admissible set $S$ of $N_A$ orbital labels, assign

$$
K_1=\sum_{j\in S}(j\bmod L_1)\pmod{L_1},
$$

$$
K_2=\sum_{j\in S}\left\lfloor\frac{j}{L_1}\right\rfloor
    \pmod{L_2}.
$$

Increment the count in that $(K_1,K_2)$ sector.  In pseudocode,

```text
counts[K1,K2] = 0
for each N_A-element subset S of {0,...,L-1}:
    if every cyclic interval of length 3 contains at most one member of S:
        K1 = sum(j mod L1 for j in S) mod L1
        K2 = sum(floor(j/L1) for j in S) mod L2
        counts[K1,K2] += 1
```

This enumeration reproduces the momentum folding of the current data.
Changing the hybrid-Wannier origin can permute momentum labels without
changing the multiset of sector counts.  For a different lattice convention,
Chern number, or aspect-ratio folding, the momentum map must be rederived;
the one-dimensional admissibility count alone does not determine label names.

At $N_A=2$ the expected multiplicities are

| geometry | low levels in each subsystem momentum sector |
|:--|:--|
| `3x4` | 3 for even $K_2$, 4 for odd $K_2$, for every $K_1$ |
| `3x5` | 5 in all 15 sectors |
| `3x6` | 6 for even $K_2$, 7 for odd $K_2$, for every $K_1$ |

Thus a correct `3x6`, $N_A=2$ plot has 18 momentum columns.  Below one
approximately horizontal entanglement gap, each even-$K_2$ column contains
six dots and each odd-$K_2$ column contains seven.  The dots need not be
degenerate or form a perfectly flat band.  Their count and momenta, rather
than their precise heights, are the fingerprint.

### 2.10 What an $N_A$ sweep should look like

Use separate momentum-resolved panels for each particle-cut size.  For
`3x6`, the predicted low-band sector counts are:

#### $N_A=1$

There are 18 quasihole-like states: one in every $(K_1,K_2)$ sector.
The full microscopic one-particle density matrix can have 36 levels, so one
expects a target low set of 18 and, in a favorable finite system, a separated
nonuniversal set above it.

#### $N_A=2$

The target is the observed 117-level band:

| $K_1$ | $K_2=0$ | 1 | 2 | 3 | 4 | 5 |
|--:|--:|--:|--:|--:|--:|--:|
| 0 | 6 | 7 | 6 | 7 | 6 | 7 |
| 1 | 6 | 7 | 6 | 7 | 6 | 7 |
| 2 | 6 | 7 | 6 | 7 | 6 | 7 |

#### $N_A=3$

The target grows to 330 levels:

| $K_1$ | $K_2=0$ | 1 | 2 | 3 | 4 | 5 |
|--:|--:|--:|--:|--:|--:|--:|
| 0 | 21 | 18 | 18 | 21 | 18 | 18 |
| 1 | 18 | 18 | 18 | 18 | 18 | 18 |
| 2 | 18 | 18 | 18 | 18 | 18 | 18 |

The microscopic $N_A=3$ lattice space has dimension
$\binom{36}{3}=7140$, so this is a substantially denser plot.  Finite-size
mixing is stronger and the entanglement gap may be less visually clean even
if the phase is unchanged.

Across the sweep, one should compare:

- the number of levels below a candidate common gap in each panel;
- their complete $(K_1,K_2)$ distribution;
- the stability of that separation with system size and Hamiltonian
  parameters.

One should **not** compare absolute $\xi$ values between different $N_A$
panels.  Each panel comes from a separately normalized density matrix.

It is usually best to stop at $N_A\le N/2$.  For a single pure state, the
nonzero singular values of complementary particle cuts are related.  For the
equal mixture of several torus ground states that simple one-state symmetry
is modified, and for $N_A>N/2$ the traced subsystem can impose an additional
rank bottleneck.  The formal zero-mode count may then exceed the PES rank
that can actually be exposed.  The standard half-cut range avoids these
ambiguities and supplies the strongest diagnostic.

### 2.11 What the current PES data show

The new FCI files use the following ground-manifold states:

| geometry | states $(K_1,K_2,\text{in-sector level})$ |
|:--|:--|
| `3x4` | $(2,2,1)$, $(1,2,1)$, $(0,2,1)$ |
| `3x5` | $(0,0,1)$, $(1,0,1)$, $(2,0,1)$ |
| `3x6` | $(0,3,1)$, $(0,3,2)$, $(0,3,3)$ |

At $N_A=2$:

| geometry | all retained PES levels | predicted low count | observed below matching gap | gap $\Delta_\xi$ | weight below gap |
|:--|--:|--:|--:|--:|--:|
| `3x4` | 276 | 42 | 42 | 2.370970 | 0.969507 |
| `3x5` | 435 | 75 | 75 | 2.721598 | 0.970693 |
| `3x6` | 630 | 117 | 117 | 2.301527 | 0.967218 |

Every geometry also matches the sector-by-sector multiplicities in Section
2.9.  This is much stronger than observing a large gap somewhere in a sorted
list.  An unrelated state can have a numerically large gap at the wrong level
count.  The correct procedure is to identify the theoretically predicted
low-band dimension and momenta first, and only then check whether a stable
gap isolates them.

For `3x6`, the gap is between

$$
\xi_{117}=5.3985263512,
\qquad
\xi_{118}=7.7000530006.
$$

The exact total and momentum-resolved agreement across all three geometries
is strong finite-size evidence for fermionic $1/3$ Laughlin-type FCI order.
It is supporting evidence to be combined with the ground-state manifold,
spectral flow, charge pump, and charge gap—not a logical proof from one plot.

---

## 3. Site real-space entanglement spectrum

### 3.1 Orbital ES, real-space ES, and the cut used here

The original Li--Haldane construction is an **orbital entanglement spectrum**:
one divides a set of Landau-level guiding-center orbitals.  An orbital is
spatially localized in one direction but extended in the other, so this is
not a sharp real-space cut.

A **real-space entanglement spectrum** instead assigns localized degrees of
freedom according to their physical position.  It retains locality along the
cut and can also expose edge structure.

The present code assigns checkerboard **site orbitals** to a contiguous strip
of unit cells.  It is therefore most precisely called a site real-space ES,
not an orbital ES in the original Landau-level sense.  The strip is cut along
direction 2:

| geometry | cells retained in direction 2 | sites in A | sites in B |
|:--|--:|--:|--:|
| `3x4` | 2 of 4 | 12 | 12 |
| `3x5` | 2 of 5 | 12 | 18 |
| `3x6` | 3 of 6 | 18 | 18 |

Because the strip contains all cells along direction 1, translation along
direction 1 remains a symmetry of the bipartition.  Momentum $K_1$ is
therefore available in principle as momentum parallel to the entanglement
boundaries.  The current implementation does not retain it in the output.

The current RSES also uses one pure absolute ground state,

$$
\rho_A^{\rm RSES}
=\operatorname{Tr}_B|\Psi_0\rangle\langle\Psi_0|,
$$

rather than the mixed ground-space projector used by the PES.  This is a
deliberate distinction: for a boundary spectrum, one generally wants a pure
minimally entangled state for the chosen cylinder cut.  On a finite torus,
the raw numerical energy eigenvector is not guaranteed to be that state,
especially when multiple ground states share one momentum sector.

Concretely, one searches within the ground space,

$$
|\Psi(c)\rangle=\sum_{a=1}^{d}c_a|\Psi_a\rangle,
\qquad \sum_a|c_a|^2=1,
$$

for combinations that minimize the entanglement entropy (or an equivalent
Rényi entropy) for the chosen cylinder cut.  In an Abelian topological phase,
these minimally entangled states approximately have definite anyon flux
through the cylinder.  Their spectra organize cleanly into one topological
sector; an arbitrary superposition can combine several sectors and obscure
the boundary towers.

### 3.2 What the binomial calculation actually proves

Expand the state in occupation bases on the two sides,

$$
|\Psi\rangle
=\sum_{N_A}\sum_{a\in\mathcal H_A(N_A)}
                \sum_{b\in\mathcal H_B(N-N_A)}
M^{(N_A)}_{ab}|a\rangle_A|b\rangle_B.
$$

If A contains $M_A$ fermionic sites and B contains $M_B$ sites, then

$$
\dim\mathcal H_A(N_A)=\binom{M_A}{N_A},
\qquad
\dim\mathcal H_B(N-N_A)=\binom{M_B}{N-N_A}.
$$

The Schmidt matrix $M^{(N_A)}$ is rectangular, so

$$
\boxed{
\operatorname{rank}M^{(N_A)}
\le
R_{\max}(N_A)
=\min\left\{
\binom{M_A}{N_A},
\binom{M_B}{N-N_A}
\right\}.}
$$

This is the entire content of the binomial derivation.  It uses neither a
Chern number nor fractionalization nor an edge theory.  A generic fermionic
state can saturate it.

The current ranks are:

| geometry | observed ranks from $N_A=0$ through $N$ | kinematic ceilings |
|:--|:--|:--|
| `3x4` | `1,12,66,12,1` | `1,12,66,12,1` |
| `3x5` | `1,12,66,153,18,1` | `1,12,66,153,18,1` |
| `3x6` | `1,18,153,676,153,18,1` | `1,18,153,816,153,18,1` |

The `3x4` and `3x5` blocks saturate the generic ceiling.  In `3x6`, only the
central $N_A=3$ block falls below it.  The implementation discards
probabilities at or below $10^{-14}$, and the smallest retained central value
is already about $1.04\times10^{-14}$.  The reported 676 is therefore
cutoff-sensitive.  It has no established universal Laughlin meaning.

These ranks and the reflection symmetry of the equal `3x4` and `3x6` cuts
are useful numerical checks.  They do **not** diagnose an FCI edge.

### 3.3 Where the single-edge sequence comes from

Suppose first that there is one isolated boundary of fermionic
$\nu=1/3$ topological order.  The long-wavelength bulk response is described
by the Abelian $K$ matrix

$$
K=(3),
\qquad t=(1).
$$

Gauge invariance of the Chern--Simons bulk in the presence of a boundary
requires a chiral boundary degree of freedom.  Quantizing its charge-density
waves gives one bosonic oscillator $a_{-n}$ for every positive boundary
momentum $n$.  At fixed edge charge and fixed topological sector, a state is

$$
\prod_{n\ge1}(a_{-n})^{r_n}|Q,a\rangle,
\qquad r_n=0,1,2,\ldots,
$$

with excess momentum

$$
\Delta K=\sum_{n\ge1}nr_n.
$$

One mode of momentum $n$ contributes

$$
1+q^n+q^{2n}+\cdots=\frac1{1-q^n}.
$$

Multiplying over independent modes gives

$$
Z_{\rm edge}(q)
=\prod_{n=1}^{\infty}\frac1{1-q^n}
=1+q+2q^2+3q^3+5q^4+7q^5+11q^6+\cdots.
$$

The coefficient of $q^{\Delta K}$ is the integer-partition number
$p(\Delta K)$.  For example, at $\Delta K=3$ the possibilities are

$$
3,\qquad2+1,\qquad1+1+1,
$$

corresponding to
$a_{-3}$, $a_{-2}a_{-1}$, and $(a_{-1})^3$.  Hence one isolated chiral edge
has

| $\Delta K$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|--:|--:|--:|--:|--:|--:|--:|--:|
| descendants | 1 | 1 | 2 | 3 | 5 | 7 | 11 |

Calling this a $U(1)_3$ chiral-boson CFT is useful but does not assume that
the **bulk** FCI is a CFT.  The bulk is gapped.  The chiral theory follows
from its Hall anomaly at a boundary; CFT is the low-energy language of that
gapless boundary.

There is also a direct wavefunction version.  Multiplying a Laughlin droplet
by symmetric power sums

$$
p_n=\sum_i z_i^n
$$

creates edge-density deformations while preserving bulk clustering.  Products
$\prod_np_n^{r_n}$ at fixed added angular momentum are again counted by
integer partitions.  This route reaches the same sequence without using a
CFT correlator to construct the bulk wavefunction.

### 3.4 Precisely where that sequence appears in an ES plot

The sequence $1,1,2,3,5,\ldots$ is **not** obtained by moving through the
spatial blocks $N_A=0,1,2,\ldots$.

For one isolated edge, one should:

1. fix the subsystem charge $N_A$ and the topological/primary sector;
2. resolve $\rho_A$ by momentum $K_\parallel$ parallel to the cut;
3. identify the edge-vacuum momentum $K_0$;
4. define the excess momentum $\Delta K=K_\parallel-K_0$ with the appropriate
   finite-size wrapping convention;
5. count only the low-lying branch below a stable entanglement gap.

The counting $1,1,2,3,5,\ldots$ runs **horizontally in momentum within one
fixed charge tower**, not vertically through total ranks of different charge
blocks.

At finite particle number the sequence can truncate or bend because edge
polynomials develop relations, momenta wrap around the Brillouin zone, and
the two sides of a narrow entanglement region interact.  Agreement is sought
first at small $\Delta K$ and then tested for stability with increasing
system size.

### 3.5 The present torus strip has two entanglement boundaries

A cylinder-shaped region cut from a torus has two boundary circles.  Their
orientations are opposite, so the subsystem contains a right-moving tower at
one boundary and a left-moving tower at the other.  Ignoring zero modes for a
moment, the oscillator content is schematically

$$
Z_{\rm strip}(q,\bar q)
=\prod_{n\ge1}\frac1{(1-q^n)(1-\bar q^n)}.
$$

If states were graded only by the **sum** of the two oscillator levels, the
counts would be the convolution

$$
\sum_{r=0}^{\ell}p(r)p(\ell-r)
=1,2,5,10,20,\ldots.
$$

That sequence is not automatically what a momentum-resolved strip plot will
show.  Parallel momentum is instead the difference of the two chiral
momenta, while entanglement energy depends roughly on their positive sum and
on generally unequal edge velocities.  Fixed total $N_A$ can also be
distributed between the two boundaries in different ways, producing several
primary towers, and finite strip width couples the boundaries.

Therefore one must not compare the current torus-strip spectrum blindly with
either $1,1,2,3,5,\ldots$ or $1,2,5,10,20,\ldots$.  The correct two-edge
character and charge constraints have to be matched to the actual cut.

If the goal is the cleanest possible single-chiral-edge sequence, an
infinite-cylinder half cut, or a disk/sphere construction with one boundary
and a conserved angular momentum, is conceptually cleaner than a finite
torus strip.  The torus strip is still valid, but its universal low spectrum
contains both entanglement edges.

### 3.6 Why a real-space cut can diagnose an edge

For a local gapped ground state, entanglement across a sharp spatial cut is
generated predominantly by degrees of freedom within a correlation length
of the cut.  The low entanglement Hamiltonian can therefore behave like a
quasi-local boundary Hamiltonian.  In a chiral topological phase, that
boundary must reproduce the anomaly of the bulk Hall response.  This is the
cut-and-glue or entanglement form of bulk--edge correspondence.

For ideal FQH model states, maps between particle, orbital, and real-space
entanglement matrices make this statement sharper: the symmetry-resolved
RSES counting is bounded by quasihole/edge counting and is expected to
saturate in the appropriate thermodynamic regime.

This does not make every detail of the spectrum a phase invariant.  The
entanglement Hamiltonian can change substantially under microscopic
deformations that leave the bulk phase unchanged.  Edge-like counting below
a stable gap is supporting evidence; raw level positions or ranks are not.

### 3.7 What the current spatial plot establishes

The current CSV contains

```text
N_A, level, probability, entanglement_energy, dim_A, dim_B
```

but no $K_\parallel$.  It can establish:

- normalization of $\rho_A$;
- probability weight in each spatial charge sector;
- Schmidt-rank and cutoff checks;
- $A\leftrightarrow B$ reflection for equal bipartitions;
- qualitative separation of groups of entanglement levels.

It cannot presently establish Li--Haldane edge counting, because all
$K_1$ sectors have been combined into each vertical $N_A$ column.  Counting
all dots in such a column simply reproduces or approaches the binomial rank
ceiling.

For `3x6`, the current block weights are

| $N_A$ | retained rank | kinematic ceiling | $\operatorname{Tr}\rho_A(N_A)$ |
|--:|--:|--:|--:|
| 0 | 1 | 1 | $1.02711\times10^{-7}$ |
| 1 | 18 | 18 | $9.10311\times10^{-4}$ |
| 2 | 153 | 153 | $0.2028837561$ |
| 3 | 676 | 816 | $0.5924116595$ |
| 4 | 153 | 153 | $0.2028837561$ |
| 5 | 18 | 18 | $9.10311\times10^{-4}$ |
| 6 | 1 | 1 | $1.02711\times10^{-7}$ |

The reflection and total trace are excellent implementation checks.  They
are not a topological edge-counting result.

### 3.8 What a genuine edge-diagnostic plot requires here

For the present strip, the next calculation should:

1. retain translation $K_1$ parallel to the cut when constructing the
   Schmidt blocks;
2. diagonalize blocks labeled by $(N_A,K_1)$;
3. choose or optimize a pure minimally entangled ground state for this cut;
4. display separate fixed-$N_A$ panels of $\xi$ versus $K_1$;
5. identify low towers and compare them with the two-edge $U(1)_3$ character,
   including allowed charge/topological sectors;
6. repeat for nearby cuts, parameters, and larger circumference to test that
   the low counting is stable while nonuniversal entanglement energies move.

Only after that extension should the spatial plot be advertised as a direct
edge-counting diagnostic.  Until then, the momentum-resolved PES is the
sharper entanglement-spectrum result in this repository.

---

## 4. Side-by-side reading guide

| question | particle ES | current site real-space ES |
|:--|:--|:--|
| What is subsystem A? | exactly $N_A$ particles anywhere on the torus | all microscopic site orbitals in a strip |
| Does $N_A$ define a new density matrix? | yes | no; it labels blocks of one $\rho_A$ |
| Is a virtual boundary created? | no | yes; two boundaries for the torus strip |
| Natural symmetry labels | $(K_1,K_2)$ | $(N_A,K_\parallel)$, although $K_\parallel$ is not yet output |
| Universal low states | bulk quasihole zero modes | boundary oscillator/primary towers |
| Counting target here | cyclic $(1,3)$ admissible roots | two-edge $U(1)_3$ character after momentum resolution |
| Meaning of total raw rank | nonuniversal microscopic rank | nonuniversal binomial-limited rank |
| Current verdict | exact total and momentum counting on `3x4`, `3x5`, `3x6` | valid Schmidt calculation, but not yet an edge-counting test |

The two successful universal countings are related by bulk--edge
correspondence, but they are not the same finite-size table.  The PES counts
quasiholes that can move throughout the bulk at fixed $N_A$ and $N_\phi$.
The spatial ES counts modes living near one or more boundaries and graded by
boundary momentum.

---

## 5. A practical checklist for future plots

### Particle ES

- Use the full low-energy torus manifold, including multiple in-sector levels
  when necessary.
- State $N_A$, $N_\phi=N_{\rm cell}$, and the microscopic PES dimension.
- Predict the total admissible count before examining numerical gaps.
- Predict its complete momentum-sector distribution.
- Mark a common low/high separation and verify the count below it.
- Prefer $N_A\le N/2$ and use separate panels for an $N_A$ sweep.
- Test stability with geometry and Hamiltonian parameters.

### Site real-space ES

- State the real-space region and number of boundary components.
- Use a pure minimally entangled ground state appropriate to the cut.
- Resolve charge and momentum parallel to the cut.
- Compare low towers within a fixed charge/topological sector, not total
  $N_A$-block ranks.
- Use the one-edge or two-edge character appropriate to the geometry.
- Treat binomial ranks only as implementation bounds.

---

## 6. Primary references

- Bernevig and Haldane,
  [*Fractional Quantum Hall States and Jack Polynomials*](https://arxiv.org/abs/0707.3637):
  dominant roots, squeezing, clustering, and generalized Pauli principles.
- Bergholtz and Karlhede,
  [*Quantum Hall system in Tao--Thouless limit*](https://arxiv.org/abs/0712.1927):
  the one-dimensional thin-torus limit, fractional domain walls, and
  continuity to bulk hierarchy states.
- Mazaheri, Ortiz, Nussinov, and Seidel,
  [*Zero modes, Bosonization and Topological Quantum Order: The Laughlin State
  in Second Quantization*](https://arxiv.org/abs/1409.3577): a polynomial-free,
  guiding-center algebraic construction of Laughlin-type zero modes.
- Li and Haldane,
  [*Entanglement Spectrum as a Generalization of Entanglement Entropy*](https://arxiv.org/abs/0805.0332):
  the original orbital-ES fingerprint and low-lying edge structure.
- Chandran, Hermanns, Regnault, and Bernevig,
  [*Bulk-Edge Correspondence in the Entanglement Spectra*](https://arxiv.org/abs/1102.2218):
  the relation between PES quasihole counting and orbital-edge counting.
- Majidzadeh Garjani, Estienne, and Ardonne,
  [*On the particle entanglement spectrum of the Laughlin states*](https://arxiv.org/abs/1501.04016):
  the polynomial-ideal formulation of PES rank and why finite-size rank
  saturation is subtler than the quasihole support bound.
- Sterdyniak, Chandran, Regnault, Bernevig, and Bonderson,
  [*Real-Space Entanglement Spectrum of Quantum Hall States*](https://arxiv.org/abs/1111.2810):
  sharp real-space cuts and their relation to quasihole and edge counting.
- Regnault and Bernevig,
  [*Fractional Chern Insulator*](https://arxiv.org/abs/1105.4867): the
  fermionic `1/3` FCI PES and $(1,3)$ quasihole fingerprint.
- Bernevig and Regnault,
  [*Emergent Many-Body Translational Symmetries of Abelian and Non-Abelian
  Fractionally Filled Topological Insulators*](https://arxiv.org/abs/1110.4488):
  folding continuum quasihole momenta into FCI Brillouin-zone sectors.
- Qi,
  [*Generic Wavefunction Description of Fractional Quantum Anomalous Hall
  States and Fractional Topological Insulators*](https://arxiv.org/abs/1105.4298):
  the hybrid-Wannier bridge between Chern bands and Landau-level orbitals.
- Chandran, Khemani, and Sondhi,
  [*How Universal Is the Entanglement Spectrum?*](https://arxiv.org/abs/1311.2946):
  why detailed entanglement spectra must not be treated as immutable phase
  invariants.
