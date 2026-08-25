# C185: fixed-point inversion and an infinite-dimensional Euler growth estimate

**Date:** 2026-08-23

**Status:** exact architecture map and a positive infinite-dimensional
linearized-Euler operator-norm estimate obtained by combining C159 with
Shvydkoy's abstract-level spectral inclusion; C189 holds the stronger
essential-spectral-radius citation pending a paper-body Theorem 4.1 check;
no unforced viscous scale-return profile, Navier--Stokes spectral gap,
nonlinear stage, or singularity

**Checker:**
[checks/essential_pde_growth_c185.py](../checks/essential_pde_growth_c185.py)

## 0. Result and boundary

The modular formulation of C136 has repeatedly moved its missing estimate
from BAFL to LBRG, RIGM, PPRG, FFCC, and CFFC.  This note changes the proof
architecture.  The object to construct is the exact stage-renormalization
operator

\[
        {\cal R}={\cal C}_{\rm exit}\circ S_T^{NS}             \tag{0.1}
\]

on a complete active-plus-wake state.  Here \(S_T^{NS}\) is exact unforced
Navier--Stokes evolution and \({\cal C}_{\rm exit}\) performs the specified
translation, rotation, dilation, amplitude normalization, phase fixing, and
exit chart.  The retained wake and all modulation coordinates are part of
the state; C125 and C175 show that deleting them makes the proposed map
non-Markovian.

The one component with no analogue in the three completed constructions
audited below is the **unforced viscous physical-velocity scale return
(UVSR)**:

> Construct a nonzero smooth boundaryless three-dimensional unforced
> Navier--Stokes orbit or approximate profile for which (0.1) returns to a
> smaller-copy structured state and the physical velocity normalization
> realizes C135's \(bq^{3/2}\) concentration.

The abstract fixed-point technology is not missing.  What is missing is the
leading unforced viscous object to which that technology could be applied.

This inversion also exposes a positive theorem already latent in the landed
corpus.  Let \(U\) be C152's smooth stationary \(A_2\) Euler field on
\(\mathbb T^3\), let \(G_t\) be the velocity-form linearized Euler group on
\(L^2_{\rm div}(\mathbb T^3)\), and let \(T\) be the certified zero-drift
orbit period.  Then

\[
 \boxed{
   \|G_{nT}\|_{L^2\to L^2}\ge e^{n/5}>
                      \left({6\over5}\right)^n
   \quad(n\ge1).}                                      \tag{0.2}
\]

Until the C189 source check is closed, the stronger displayed form
\(r_{\rm ess}(G_T)\ge e^{1/5}\) is citation-held and is not a downstream
premise.

Thus the program now has an explicit positive cocycle-growth estimate for
an actual infinite-dimensional PDE operator, not merely for a Galerkin
matrix or a frozen symbol.  It is an inviscid linear estimate.  It does not
provide UVSR, a viscous spectral gap, finite-frequency endpoint coherence,
or nonlinear closure.

## 1. The renormalized state and residual

Write a structured state as

\[
 X=(u_{\rm act},u_{\rm wake},\lambda,R,x_0,\vartheta,\zeta), \tag{1.1}
\]

where \(\lambda,R,x_0\) are scale, rotation, and translation data,
\(\vartheta\) fixes the temporal/phase gauge, and \(\zeta\) denotes the
finite endpoint coordinates used by the chart.  This list is schematic;
choosing a Banach space and a nonredundant modulation slice remains part of
UVSR.  Its purpose is to state the Markov requirement: every component that
can return at the next stage belongs to \(X\).

For an approximate state \(X_*\), the load-bearing object is the single
residual

\[
                 {\cal E}(X_*)={\cal R}(X_*)-X_*.        \tag{1.2}
\]

The former active leakage, wake leakage, pressure return, connector defect,
and chart error are projections of (1.2).  A fixed-point or trapping proof
may use projections, but it must establish one full-operator estimate of the
form

\[
 \|D{\cal R}(X_*)h\|_{X_s}\le \kappa\|h\|_{X_s}
 \quad(\kappa<1)                                      \tag{1.3}
\]

on the stable complement, together with a finite expanding/modulation
block and one explicit nonlinear remainder estimate.  Equations (1.2)--
(1.3) replace the practice of promoting each failed projection estimate to
a new architecture gate.  No value of \(\kappa\) is claimed here.

## 2. Map to completed constructions

### 2.1 Albritton--Bru\'e--Colombo: viscosity and full closure, but forced profile creation

Albritton, Bru\'e, and Colombo write forced Navier--Stokes in autonomous
similarity variables and make a smooth compactly supported vector field
\(\bar U\) a steady similarity profile.  Their force is defined from
\(\bar U\); this freedom is precisely what the unforced problem does not
have.

The importable chain is:

1. Theorem 3.1 and Corollary 3.2 transfer a genuine unstable eigenvalue to
   the infinite-dimensional self-similar Navier--Stokes linearization.
2. Theorem 4.1 selects the maximally unstable spectral block.
3. Lemma 4.4 supplies the full parabolic semigroup growth and smoothing
   estimate.
4. Propositions 4.5 and 4.7 close the entire nonlinear remainder by one
   Duhamel contraction on the unstable manifold.

This is the completed analogue of replacing C125's scalar projected
identity and all named complement channels by a spectral projection and one
Banach-space contraction.  It does not supply UVSR because the force, rather
than the unforced profile equation, creates the leading state.

Primary source: D. Albritton, E. Bru\'e, and M. Colombo,
[*Non-uniqueness of Leray solutions of the forced Navier--Stokes
equations*](https://annals.math.princeton.edu/2022/196-1/p03), especially
Sections 1.1, 3, and 4; arXiv:2112.03116.

### 2.2 Elgindi: unforced profile, coercivity, and compactness, but inviscid low regularity

Elgindi constructs an exact unforced self-similar Euler profile in an
axisymmetric no-swirl symmetry class.  The relevant architecture is:

1. symmetry and a small parameter isolate an explicitly solvable
   fundamental model;
2. the full profile is written as the model plus a correction;
3. Proposition 5.5 gives the explicit weighted coercivity

   \[
      (L_\Gamma f\,w,fw)_{L^2}\ge {1\over4}\|fw\|_2^2;   \tag{2.1}
   \]

4. Corollary 6.11 and Proposition 6.13 propagate coercivity through the
   angular and higher-order hierarchy;
5. Theorem 2 controls the nonlocal Biot--Savart correction uniformly;
6. Proposition 9.6 and Corollary 9.7 close the nonlinear error; and
7. Section 9.5 evolves an artificial-time problem, obtains uniform bounds,
   and passes to a stationary profile by compactness.

This is the completed analogue of using the \(A_2\) structure to define one
invariant symmetry class and one coercive full-operator norm.  It does not
supply UVSR: the equation is inviscid, the profile has \(C^{1,\alpha}\)
velocity, and its leading Biot--Savart mechanism is tied to the
axisymmetric small-\(\alpha\) regime.

Primary source: T. Elgindi,
[*Finite-Time Singularity Formation for \(C^{1,\alpha}\) Solutions to the
Incompressible Euler Equations on \(\mathbb R^3\)*](https://arxiv.org/abs/1904.04795),
especially Sections 5--7 and 9; Annals of Mathematics 194 (2021).

### 2.3 Chen--Hou: a validated trapping tube, but boundary Euler

Chen and Hou demonstrate that an exact fixed profile is not mandatory.  In
dynamic similarity variables they prove that the exact solution stays for
all rescaled time in an explicitly certified neighborhood of a numerical
approximate steady profile.  Their proof combines:

1. an outgoing-flow inequality which turns weighted transport into damping;
2. weighted \(L^\infty\), H\"older, and nonlocal estimates;
3. a splitting \(L=L_0+K\), with \(L_0\) coercive and \(K\) finite rank;
4. validated space-time responses for every direction in
   \(\operatorname{ran}K\); and
5. one nonlinear stability lemma containing the full profile residual and
   nonlinear error.

Part II bounds the approximate profile, residual, regular velocity, and
finite-rank responses with interval arithmetic.  This is the closest
completed model for

\[
 \text{coercive infinite-dimensional bulk}
 +\text{validated finite-dimensional bad block}
 +\text{one trapping inequality}.                       \tag{2.2}
\]

Their leading object is not UVSR: it is an inviscid axisymmetric Euler
profile on a cylinder with a physical boundary.  The already-landed C33
audit also shows that its physical focusing scale is diffusion-dominated
for fixed positive viscosity.

Primary sources: J. Chen and T. Hou,
[*Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler
equations with smooth data I: Analysis*](https://arxiv.org/abs/2210.07191),
especially Sections 2 and 5 and Lemma A.2, and
[*Part II: Rigorous Numerics*](https://arxiv.org/abs/2305.05660), especially
the validated profile and finite-rank response estimates.

### 2.4 Exact correspondence with the old modules

| Repository object | Completed analogue | Import |
|---|---|---|
| C136 one-cell return | similarity steady state or trapping tube | define \({\cal R}\) and its normalization first |
| C125 projected Duhamel gate | ABC Section 4 | full semigroup plus one contraction |
| C175 invariant graph | ABC unstable manifold; Chen--Hou finite bad block | one spectral/finite-rank splitting |
| C176 residence and chart | Elgindi symmetry coordinates; Chen--Hou modulation weights | build them into the state norm |
| C177--C180 reservoirs | approximate-profile ingredients | measure the full residual, not isolated edges |
| C178 prescribed preparation | whole-trajectory construction | absorb preparation into the initial state/fixed point |
| C181--C184 exclusions | no direct palette analogue | use them only to constrain UVSR profile search |

## 3. C185: the positive infinite-dimensional estimate

### Theorem 3.1 (linearized-Euler operator-norm growth)

Let

\[
 f(a,b)=\cos a+\cos b+{4\over5}\cos(a+b),\qquad
 U=N\times\nabla f-\sqrt2 fN                         \tag{3.1}
\]

be the smooth real periodic Beltrami field of C152, in the normalization of
C159.  Let \(G_t\) be the \(C_0\)-group generated on
\(L^2_{\rm div}(\mathbb T^3)\) by the velocity-form Euler equation
linearized about \(U\).  Let \(T\) be the period of C152's zero-drift orbit.
Then (0.2) holds.

### Proof

C159 supplies a periodic bicharacteristic \((x(t),k(t))\) and a periodic
oriented transverse frame.  In unit period time \(s=t/T\), its amplitude
system is

\[
                         z'=B(s)z.                       \tag{3.2}
\]

The outward certificate proves that \(B\) is Metzler and, for

\[
                         w=(1,3/20)^T,                   \tag{3.3}
\]

that

\[
                         B(s)w>{1\over5}w                \tag{3.4}
\]

componentwise for the entire orbit.  Cooperative comparison gives for the
one-period monodromy \(M\)

\[
                         Mw\ge e^{1/5}w.                 \tag{3.5}
\]

The positive cone is invariant, so iteration yields
\(M^nw\ge e^{n/5}w\).  The frame returns at each period and its two columns
are orthogonal.  Hence componentwise enlargement of the positive
coefficients implies enlargement in the physical Euclidean fiber norm.
The maximal Lyapunov exponent \(\mu_{\max}\) of the velocity-amplitude
bicharacteristic system therefore satisfies

\[
                         \mu_{\max}\ge {1\over5T}.        \tag{3.6}
\]

The abstract-level spectral inclusion in Shvydkoy puts the exponential of
the relevant bicharacteristic-amplitude spectral point into the spectrum of
the velocity-form linearized Euler group.  Consequently

\[
 r(G_{nT})\ge e^{n/5},\qquad
 \|G_{nT}\|_{L^2\to L^2}\ge r(G_{nT})\ge e^{n/5}.
                                                               \tag{3.7}
\]

C189's independent source-scope check confirms that this abstract
inclusion covers the velocity-form group and the C159
bicharacteristic-amplitude system at the level needed for (3.7).  The
further identification

\[
                         r_{\rm ess}(G_t)=e^{\mu_{\max}t}
                                                               \tag{3.7a}
\]

is the paper-body Theorem 4.1 step held by C189 for source verification;
it is not needed for the operator-norm conclusion (3.7).  Finally

\[
 e^{1/5}>1+{1\over5}+{1\over2}\left({1\over5}\right)^2
          ={61\over50}>{6\over5},                       \tag{3.8}
\]

which proves (0.2).  \(\square\)

Load-bearing source scope: the abstract-level spectral inclusion in
R. Shvydkoy,
[*The essential spectrum of advective equations*](https://arxiv.org/abs/math-ph/0412019),
as confirmed by C189.  The paper-body Theorem 4.1, Propositions 4.2--4.3,
equation (41), and inequality (81) remain the pending source check for the
essential-radius refinement (3.7a).

## 4. What C185 changes and what it does not

C185 discharges one strategic deficit: the certified C159 exponent is now
a theorem about the actual infinite-dimensional linearized Euler operator.
No separate finite-frequency packet construction is required merely to
prove operator growth; geometric optics is already encoded in the cited
spectral inclusion.  The essential-radius refinement remains source-held
as stated above.

The following remain open and are not consequences of (0.2):

1. a useful coherent \(q^2\)-to-\(q^3\) endpoint rather than an operator-
   norm witness;
2. persistence through positive viscosity on the stage time and frequency
   scales;
3. a spectral splitting or coercive complement estimate for
   \(D{\cal R}(X_*)\);
4. nonlinear trapping of the active-plus-wake state; and
5. UVSR itself.

In particular, essential growth about a steady Euler field is not a
Navier--Stokes scale-return theorem.  Its proper role is to provide a
literature-certified positive expanding direction once a candidate UVSR
profile and its viscous state space have been constructed.

## 5. Verification boundary

The checker verifies with exact rational arithmetic that C159's four
rounded cone inequalities imply (3.4), that the fixed positive vector gives
physical norm growth in an orthogonal returning frame, and that the
explicit lower bound (3.8) is stronger than \(6/5\).  It also verifies the
iteration and exponent ledger on finite exact test instances.

The checker does not reproduce C159's interval-Taylor orbit certificate or
prove Shvydkoy's pseudodifferential essential-spectrum theorem.  Those are
separately identified premises: the former is the landed C159 artifact and
the latter is the cited primary-source theorem.  It does not verify UVSR,
viscosity, nonlinear closure, or a Navier--Stokes singularity.
