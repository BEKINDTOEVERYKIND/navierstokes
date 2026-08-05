# Countable-stagnation Liouville theorem for sub-parabolic Euler profiles

**Date:** 2026-08-01

**Status:** corrected source rendering and independently rechecked argument.
This is a derived theorem, not a Navier--Stokes singularity construction and
not a literature-novelty claim.

The source identities are from Constantin--Ignatova--Vicol,
[arXiv:2602.17570v3](https://arxiv.org/abs/2602.17570).

## Theorem

Let \(0<\gamma<1/2\), and let

\[
 (U,P)\in C^2(\mathbb R^3)\times C^1(\mathbb R^3)
\]

solve the globally self-similar incompressible Euler profile equation

\[
 (1-\gamma)U+\gamma(y\cdot\nabla)U
 +(U\cdot\nabla)U+\nabla P=0,
 \qquad \nabla\cdot U=0                                      \tag{1}
\]

on \(\mathbb R^3\), with

\[
 U(0)=0,\qquad U(y)=o(|y|)\quad (|y|\to\infty).              \tag{2}
\]

Set

\[
 V(y)=\gamma y+U(y),\qquad
 \Omega=\nabla\times U,\qquad
 N_V=\{y:V(y)=0\}.                                           \tag{3}
\]

If \(N_V\) is countable, then

\[
                         \Omega\equiv0,\qquad U\equiv0.       \tag{4}
\]

Consequently, every nontrivial sub-parabolic profile satisfying the stronger
CIV far-field bound (3.8) must have an uncountable stagnation set.

## Proof

### 1. Backward precompactness

The vector field \(V\) has at most linear growth, so its flow \(\Phi_t\) is
complete.  From (2), for sufficiently large \(|y|\),

\[
                         V(y)\cdot y\ge {\gamma\over2}|y|^2.  \tag{5}
\]

Thus every negative semiorbit enters a fixed ball and cannot leave that ball
as time decreases.  Every negative semiorbit is precompact.

### 2. Every negative orbit converges to one stagnation point

CIV's self-similar Bernoulli function satisfies

\[
                      V\cdot\nabla H=(2\gamma-1)|V|^2.       \tag{6}
\]

Because \(2\gamma-1<0\), \(H\) is a strict Lyapunov function away from
\(N_V\).  On a precompact negative orbit, the alpha-limit set is nonempty,
compact, connected, and invariant.  Applying the integrated form of (6) on
that limit set places it inside \(N_V\).

A connected countable subset of \(\mathbb R^3\) is a singleton.  Therefore,
for every \(a\in\mathbb R^3\), there is a \(p=p(a)\in N_V\) such that

\[
                    \Phi_t(a)\longrightarrow p
                    \qquad (t\to-\infty).                    \tag{7}
\]

The space is partitioned by the backward basins \(B^-(p)\).

### 3. Low-spectral-abscissa basins carry no vorticity

Fix \(a\in B^-(p)\), and put

\[
 A_p=DV(p),\qquad
 \alpha_p=\max\{\operatorname{Re}\lambda:
                   \lambda\in\sigma(A_p)\}.                  \tag{8}
\]

Assume first that \(\alpha_p<\gamma+1\), and choose

\[
                  \alpha_p<q_1<q<\gamma+1.                  \tag{9}
\]

The shifted self-similar Cauchy formula is

\[
 \Omega(\Phi_{t_0}a)
 =e^{-(\gamma+1)(t_0-t)}
  D\Phi_{t_0-t}(\Phi_ta)\,\Omega(\Phi_ta).                   \tag{10}
\]

Because \(DV(\Phi_ta)\to A_p\), the standard asymptotically autonomous
fundamental-matrix estimate gives, for \(t\ll t_0\),

\[
              |D\Phi_{t_0-t}(\Phi_ta)|
              \le C_a e^{q(t_0-t)}.                          \tag{11}
\]

For completeness, Jordan blocks and nonnormality are handled by first using
\(|e^{A_ps}|\le M e^{q_1s}\), then taking
\(DV(\Phi_ta)-A_p\) uniformly small and applying Duhamel--Gronwall.  No
eigenbasis or Euclidean logarithmic norm is assumed.

The negative orbit is compact, so \(\Omega(\Phi_ta)\) is bounded.  Combining
(10)--(11) and sending \(t\to-\infty\) gives

\[
                         \Omega(\Phi_{t_0}a)=0.              \tag{12}
\]

In particular, \(\Omega(a)=0\) on every basin with
\(\alpha_p<\gamma+1\).

### 4. The remaining basins have measure zero

Since

\[
                 \operatorname{tr}A_p=\nabla\cdot V=3\gamma
                 <\gamma+1,                                 \tag{13}
\]

the condition \(\alpha_p\ge\gamma+1\) forces at least one eigenvalue of
\(A_p\) to have strictly negative real part.  The local center--unstable
manifold \(W^{cu}_{\rm loc}(p)\), formed from all generalized eigenspaces
with nonnegative real part, therefore has dimension at most two.

The center--unstable containment theorem says that every full negative
semiorbit which remains near \(p\) lies in \(W^{cu}_{\rm loc}(p)\).
By (7), every orbit in \(B^-(p)\) eventually has this property, hence

\[
 B^-(p)\subseteq
 \bigcup_{n\ge0}\Phi_n\!\left(W^{cu}_{\rm loc}(p)\right).    \tag{14}
\]

The right side is a countable union of at-most-two-dimensional immersed
\(C^1\) manifolds and has three-dimensional Lebesgue measure zero.

There are only countably many \(p\in N_V\).  Thus the union of all
high-spectral-abscissa basins is null, while (12) gives \(\Omega=0\) on
every other basin.  Hence \(\Omega=0\) almost everywhere, and continuity
gives \(\Omega\equiv0\).

Finally, \(\nabla\cdot U=0\) and \(\nabla\times U=0\) imply
\(\Delta U=0\).  An entire sublinear harmonic vector field is constant, and
\(U(0)=0\), so \(U\equiv0\).

## Safe normally-hyperbolic extension

The same conclusion holds if \(N_V\) is a finite pairwise-disjoint union of
compact connected \(C^2\) normally hyperbolic equilibrium submanifolds.

For a component \(S\) of dimension \(k\), normal hyperbolicity supplies
asymptotic phase and an unstable lamination.

- If \(S\) has a stable normal direction, its backward basin lies in
  \(W^{cu}(S)\), whose dimension is at most two, and is null.
- If \(S\) is normally repelling, its tangent eigenvalues are zero and all
  normal real parts are positive with sum \(3\gamma\).  Its spectral
  abscissa is therefore at most \(3\gamma<\gamma+1\), and the Cauchy estimate
  kills vorticity throughout its possibly full-dimensional basin.

A finite union completes the same almost-everywhere argument.

This extension must not be stated for an arbitrary Whitney-stratified or
general uncountable equilibrium set.  Without normal hyperbolicity,
asymptotic phase can fail; a trajectory can drift tangentially toward a
continuum, and uncountably many lower-dimensional basins can collectively
have positive measure.

## Consequence for the proposed Euler-to-Navier--Stokes route

CIV leave the formal window \(2/5\le\gamma<1/2\) only for genuinely
nonaxisymmetric profiles without their outgoing property.  The theorem above
narrows that window further: a surviving profile needs an uncountable,
dynamically degenerate stagnation set which is not a finite union of compact
normally hyperbolic components.

It also excludes an exact bounded invariant solid torus for the similarity
field \(V\), because

\[
 \int_{\partial D}V\cdot n
 =\int_D\nabla\cdot V
 =3\gamma|D|>0.                                            \tag{15}
\]

Together with (6), this rules out the literal transplantation of a recurrent
Gavrilov torus as an invariant similarity-flow cell.  It does not rule out
all perturbative constructions: a perturbation may destroy torus invariance
or create an uncountable degenerate stagnation set.

The result is a useful no-go, not a positive blow-up mechanism.  It demotes
the sub-parabolic self-similar route relative to the physical-variable
Palasek--Gavrilov transition program.
