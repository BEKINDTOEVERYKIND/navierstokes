# Countable-stagnation Liouville theorem for sub-parabolic Euler profiles

Date: 2026-08-01

## Claim boundary

This note records a derived Liouville theorem for the globally self-similar
incompressible Euler profile equation.  It is **not** a Navier--Stokes
singularity construction and is not claimed as a literature theorem.  The
argument has passed two independent internal audits, including the
nonnormal/Jordan and center--unstable cases, but still requires external
expert review before any novelty claim.

The source identities are from Constantin--Ignatova--Vicol (CIV),
[arXiv:2602.17570v3](https://arxiv.org/abs/2602.17570).

## Theorem

Let (0<gamma<1/2), and let
((U,P)in C^2(mathbb R^3)	imes C^1(mathbb R^3)) solve

[
 (1-gamma)U+gamma(ycdot
abla)U+(Ucdot
abla)U+
abla P=0,
 qquad 
ablacdot U=0
 	ag{1}
]

on (mathbb R^3), with

[
 U(0)=0,qquad U(y)=o(|y|)quad (|y|	oinfty).
 	ag{2}
]

Set

[
 V(y)=gamma y+U(y),qquad Omega=
abla	imes U,qquad
 N_V={y:V(y)=0}.
 	ag{3}
]

If (N_V) is countable, then

[
 Omegaequiv0,qquad Uequiv0.
 	ag{4}
]

Consequently, every nontrivial sub-parabolic profile satisfying the stronger
CIV far-field bound (3.8) must have an **uncountable** stagnation set.

## Proof

### 1. Backward precompactness

The vector field (V) has at most linear growth, so its flow
(Phi_t) is complete.  From (2), for sufficiently large (|y|),

[
 V(y)cdot yge rac{gamma}{2}|y|^2.
 	ag{5}
]

Thus every negative semiorbit enters a fixed ball and cannot leave that ball
as time decreases.  Each negative semiorbit is precompact.

### 2. Every negative orbit converges to one stagnation point

CIV's self-similar Bernoulli function satisfies

[
 Vcdot
abla H=(2gamma-1)|V|^2.
 	ag{6}
]

Since (2gamma-1<0), (H) is a strict Lyapunov function away from (N_V).
On a precompact negative orbit, the alpha-limit set is nonempty, compact,
connected, and invariant.  Equation (6) (or its integrated form plus uniform
continuity of (V(Phi_ta))) places that alpha-limit set inside (N_V).

A connected countable subset of (mathbb R^3) is a singleton.  Therefore,
for every (ainmathbb R^3), there is a (p=p(a)in N_V) such that

[
 Phi_t(a)longrightarrow pqquad(t	o-infty).
 	ag{7}
]

The space is partitioned by the backward basins (B^-(p)).

### 3. Low spectral-abscissa basins carry no vorticity

Fix (ain B^-(p)), put

[
 A_p=DV(p),qquad
 alpha_p=max{mathop{m Re}lambda:lambdainsigma(A_p)}.
 	ag{8}
]

Assume first that (alpha_p<gamma+1), and choose

[
 alpha_p<q_1<q<gamma+1.
 	ag{9}
]

The shifted self-similar Cauchy formula is

[
 Omega(Phi_{t_0}a)
 =e^{-(gamma+1)(t_0-t)}
 DPhi_{t_0-t}(Phi_ta),Omega(Phi_ta).
 	ag{10}
]

Because (DV(Phi_ta)	o A_p), the standard asymptotically autonomous
fundamental-matrix estimate gives, for (tll t_0),

[
 |DPhi_{t_0-t}(Phi_ta)|
 le C_a e^{q(t_0-t)}.
 	ag{11}
]

For completeness, Jordan blocks and nonnormality are handled by first using
(|e^{A_ps}|le M e^{q_1s}), then taking (DV(Phi_ta)-A_p) uniformly
small and applying Duhamel--Gronwall.  No eigenbasis or Euclidean logarithmic
norm is being assumed.

The negative orbit is compact, so (Omega(Phi_ta)) is bounded.  Combining
(10)--(11) and sending (t	o-infty) gives

[
 Omega(Phi_{t_0}a)=0.
 	ag{12}
]

In particular, (Omega(a)=0) on every basin for which
(alpha_p<gamma+1).

### 4. The remaining basins have measure zero

Since (
ablacdot V=3gamma),

[
 mathop{m tr}A_p=3gamma<gamma+1.
 	ag{13}
]

If (alpha_pgegamma+1), (13) forces at least one eigenvalue of (A_p)
to have strictly negative real part.  The local center--unstable manifold
(W^{cu}_{m loc}(p)), defined using all generalized eigenspaces with
nonnegative real part, therefore has dimension at most two.

The center--unstable containment theorem says that every full negative
semiorbit which remains near (p) lies in (W^{cu}_{m loc}(p)).  By
(7), every orbit in (B^-(p)) eventually has this property, hence

[
 B^-(p)subseteq
 igcup_{nge0}Phi_n(W^{cu}_{m loc}(p)).
 	ag{14}
]

The right-hand side is a countable union of at-most-two-dimensional immersed
(C^1) manifolds and has three-dimensional Lebesgue measure zero.

There are only countably many (pin N_V).  Thus the union of all basins in
the high-spectral-abscissa case is null, while (12) gives (Omega=0) on
every other basin.  Hence (Omega=0) almost everywhere, and continuity
gives (Omegaequiv0).

Finally, (
ablacdot U=0) and (
abla	imes U=0) imply
(Delta U=0).  An entire sublinear harmonic vector field is constant, and
(U(0)=0), so (Uequiv0).

## Safe normally-hyperbolic extension

The same conclusion holds if (N_V) is a finite pairwise-disjoint union of
compact connected (C^2) normally hyperbolic equilibrium submanifolds.

For a component (S) of dimension (k), normal hyperbolicity supplies
asymptotic phase and an unstable lamination.

- If (S) has a stable normal direction, its backward basin lies in
  (W^{cu}(S)), whose dimension is at most two, and is null.
- If (S) is normally repelling, its tangent eigenvalues are zero and all
  normal real parts are positive with sum (3gamma).  Its spectral
  abscissa is therefore at most (3gamma<gamma+1), and the Cauchy estimate
  kills vorticity throughout its possibly full-dimensional basin.

A finite union completes the same almost-everywhere argument.

This extension should **not** be stated for an arbitrary Whitney-stratified
or general uncountable equilibrium set.  Without normal hyperbolicity,
asymptotic phase can fail; a trajectory can drift tangentially toward a
continuum, and uncountably many lower-dimensional basins can collectively
have positive measure.

## Consequence for the proposed Euler-to-Navier--Stokes route

CIV leave the formal window (2/5legamma<1/2) only for genuinely
nonaxisymmetric profiles without their outgoing property.  The theorem above
narrows that window further: a surviving profile needs an uncountable,
dynamically degenerate stagnation set which is not a finite union of compact
normally hyperbolic components.

This also blocks a perturbative transplant of Gavrilov's compact steady
tori.  In similarity variables, (6) forbids nonstationary recurrence, while

[
 int_{partial D}Vcdot n
 =int_D
ablacdot V
 =3gamma|D|>0
 	ag{15}
]

forbids (V) from being tangent to the boundary of any bounded invariant
solid torus.

The result is a useful no-go, not a positive blow-up mechanism.  It demotes
the sub-parabolic self-similar route relative to the physical-variable
Palasek--Gavrilov three-scale rectifier.
