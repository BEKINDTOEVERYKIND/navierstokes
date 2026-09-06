# C182: conditional fixed-time no-focus under a scale-adapted entrance norm

**Date:** 2026-08-14

**Status:** exact linearized-Euler pressure reduction and conditional
fixed-time finite-\(p\) no-focus theorem under explicit unnormalized
entrance-norm hypotheses; localization/tails, nonlinear depletion,
logarithmic time, C125, RIGM, BAFL, and the one-cell stage remain open

**Checker:**
[checks/holder_terminal_no_focus_c182.py](../checks/holder_terminal_no_focus_c182.py)

## 0. Claim boundary

C180 constructs \(Q=q+O(1)\) same-shell Beltrami reservoir modes whose
one-edge Fourier row is an almost-isometric \(q\)-star. Its leading edge is
transport, and the missing step was the **finite-frequency composition
closure (FFCC)**.

Two tempting endpoint proofs are invalid here. After rescaling the parent
frequency, the spatial domain is the large torus

\[
             \mathbb T_L^3=(\mathbb R/L\mathbb Z)^3,
             \qquad L=2\pi K\ge 2\pi .                  \tag{0.1}
\]

A generic inhomogeneous Hölder/Riesz constant need not be uniform at the
lowest frequencies as \(L\to\infty\). A generic homogeneous
\(B^0_{\infty,1}\) replacement also does not close the stretching product:
high--high interactions can land in arbitrarily low shells. Neither route
is used below.

Instead use the inhomogeneous Bessel-potential norm
\(W^{1/4,16}(\mathbb T_L^3)\), with **unnormalized** Lebesgue measure. Let
\(U_q\) be a smooth divergence-free periodic base on one fixed normalized
terminal interval \(0\le t\le T_0\), with

\[
       \sup_q\int_0^{T_0}\|U_q(t)\|_{C^2}\,dt\le M<\infty .
                                                               \tag{0.2}
\]

Assume that the entrance packet obeys the three global bounds

\[
 \boxed{
 \begin{aligned}
  \|v_{q,0}\|_{L^2(\mathbb T_L^3)}&\le C_2 b,\\
  \|v_{q,0}\|_{L^\infty(\mathbb T_L^3)}&\le C_\infty bq,\\
  \|v_{q,0}\|_{W^{1/4,16}(\mathbb T_L^3)}
      &\le C_Bq^{1/4}\|v_{q,0}\|_{L^{16}(\mathbb T_L^3)}.
 \end{aligned}}                                             \tag{0.3}
\]

The last line is a scale-\(q\) spectral-tail hypothesis. Exact support in a
\(Cq\) Fourier ball implies it by Bernstein, but a localized packet is not
exactly band limited; proving the corresponding tail estimate is open.

For the complete linearized Euler solution about \(U_q\), including every
daughter, reverse edge, repeated sideband, and Leray-pressure cancellation,
(0.2)--(0.3) imply

\[
 \boxed{
   \sup_{0\le t\le T_0}\|v_q(t)\|_\infty
       \le C_{M,T_0,C_2,C_\infty,C_B}\,bq^{9/8}.}        \tag{0.4}
\]

Relative to C161's raw coherent target \(bq^{3/2}\), this is

\[
                         Cq^{-3/8}\longrightarrow0.      \tag{0.5}
\]

C180's narrowed spatial bounding box carries an honest \(J^2\)
point-normalization tax. Even relative to the reduced target
\(bq^{3/2}/J^2\), the ratio is

\[
                         CJ^2q^{-3/8}\longrightarrow0    \tag{0.6}
\]

for \(J=O(\log q)\). Thus a fixed-time linear C180 splitter satisfying
(0.3) cannot provide the missing focus.

This is a **conditional** no-focus theorem, not an unconditional FFCC
closure. C176/C180 do not construct a localized packet satisfying the
unnormalized bounds (0.3). The theorem also excludes logarithmic terminal
intervals, nonlinear/depleting leading trajectories, and C179's
\(qK\)-frequency planar reservoir. In particular, C177's same-curl
reservoir shares the pump's heat scalar and cannot simply be switched on for
only the last fixed interval. If it is preloaded throughout C159's
\(O(\log q)\) gain window, the entrance state at the last interval need
not satisfy (0.3).

## 1. The exact pressure reduction has no derivative loss

Let \(v\) solve the periodic linearized Euler equation

\[
 \partial_tv+U\cdot\nabla v+v\cdot\nabla U+\nabla p=0,
 \qquad \nabla\cdot U=\nabla\cdot v=0.                    \tag{1.1}
\]

Taking the divergence and using incompressibility gives

\[
 -\Delta p
   =2\,\partial_iU_j\,\partial_jv_i
   =2\,\partial_i\bigl((\partial_jU_i)v_j\bigr).         \tag{1.2}
\]

Indeed, the nominal second derivatives in

\[
 \partial_i(U_j\partial_jv_i)
 +\partial_i(v_j\partial_jU_i)
\]

vanish after differentiating the two divergence constraints, while the
two remaining products agree after swapping \(i,j\). Hence

\[
 \boxed{
 \nabla p
   =2\,\nabla(-\Delta)^{-1}\partial_i
       \bigl((\partial_jU_i)v_j\bigr).}                   \tag{1.3}
\]

The pressure-gradient zero mode is set to zero. The multiplier in (1.3)
is order zero. Its \(L^{16}\) and \(W^{1/4,16}\) bounds are finite-\(p\)
periodic Mikhlin estimates whose constants are uniform for \(L\ge2\pi\).

## 2. Scale-uniform finite-\(p\) propagation

For \(0<s<1\) and \(1<p<\infty\), the fractional transport commutator,
finite-\(p\) product estimate, and (1.3) give

\[
 {d\over dt}\|v(t)\|_{W^{s,p}}
 \le C_{s,p}\|U(t)\|_{C^2}\|v(t)\|_{W^{s,p}}.           \tag{2.1}
\]

The deliberately strong \(C^2\) base norm makes every product and
commutator coefficient harmless. The constants are uniform on
\(\mathbb T_L^3\): equivalently, one may use the scale-uniform periodic
singular-integral and difference-quotient proofs. Gronwall with
\(s=1/4,p=16\) yields

\[
 \|v(t)\|_{W^{1/4,16}}
 \le \exp\!\left(C\int_0^t\|U(r)\|_{C^2}\,dr\right)
       \|v(0)\|_{W^{1/4,16}}.                            \tag{2.2}
\]

Since \(1/4>3/16\), the scale-uniform Sobolev embedding gives

\[
                 \|v\|_\infty\le C\|v\|_{W^{1/4,16}}.  \tag{2.3}
\]

For the inhomogeneous norm, its \(L^{16}\) part controls the mean with no
large-domain loss:

\[
 |\langle v\rangle|
       \le L^{-3/16}\|v\|_{16}\le\|v\|_{16}
       \qquad(L\ge2\pi).                                 \tag{2.4}
\]

The same conclusion for linearized Navier--Stokes follows from the
standard scale-uniform transport-diffusion \(W^{s,p}\) estimate, uniformly
for \(\nu\ge0\). It is not inferred merely from heat dissipativity, and it
does not construct the localized viscous stage.

## 3. The exact C180 base meets the uniform hypothesis

Rescale the parent frequency \(K\) to one and factor the parent velocity
scale as in C180. The old \(A_2\) pump consists of six unit-frequency modes
with fixed coefficients, hence has \(C^2\) norm \(O(1)\).

C180's reservoir has \(Q=q+O(1)\) unit-frequency modes, each with relative
coefficient

\[
                      \rho_*\asymp{1\over q\sqrt Q}.      \tag{3.1}
\]

For every fixed derivative order, the triangle inequality gives

\[
 \|U_{\rm res}\|_{C^2}
      \le C Q\rho_*
      \le C{\sqrt Q\over q}
      =O(q^{-1/2}).                                      \tag{3.2}
\]

The unforced common-shell heat factor is at most one. Therefore the exact
C121 pump plus exact same-curl C180 reservoir satisfies (0.2) on every
fixed normalized interval. This uses the same-shell fact: the reservoir
frequency is one in parent coordinates. C179's planar palette instead has
normalized frequency \(q\), so (3.2) does not apply to PPRG.

## 4. The entrance interpolation and its unresolved localization

The unnormalized entrance bounds in (0.3) give exactly

\[
 \|v_{q,0}\|_{16}
 \le \|v_{q,0}\|_2^{1/8}\|v_{q,0}\|_\infty^{7/8}
 \le Cbq^{7/8}.                                        \tag{4.1}
\]

The scale-\(q\) tail hypothesis then gives

\[
                 \|v_{q,0}\|_{W^{1/4,16}}
                         \le Cbq^{9/8}.                 \tag{4.2}
\]

Equations (2.2)--(2.3) prove (0.4).

The three assumptions in (0.3) are not consequences of C180's arithmetic.
They describe one localized parent packet after the \(y=Kx\) rescaling.
Exact compact support and exact Fourier bandlimiting are incompatible, so
the last line requires a genuine tail theorem. Moreover, normalized-Haar
control on the original torus cannot be substituted silently: conversion
to the unnormalized norm on \(\mathbb T_L^3\) introduces powers of \(K\),
including \(K^{3/16}\) at the \(L^{16}\) level.

The former heuristic that filled the C180 bounding volume
\(J^4q^{-3}\) at point scale \(bq\) is explicitly rejected: its resulting
\(L^2\) upper scale \(bJ^2q^{-1/2}\) is much smaller than the scheduled
\(b\), so it cannot represent the intended packet.

On the schedule \(q=n^8\) and \(J=O(\log n)\),

\[
 q^{-3/8}=n^{-3},
 \qquad J^2q^{-3/8}=O((\log n)^2n^{-3}),                \tag{4.3}
\]

and both tend to zero.

## 5. What survives

C182 closes only the **conditional fixed-time linear** same-shell branch.
It proves that a complete C180 propagator satisfying the explicit entrance
norms cannot turn its bright one-edge tight frame into the required physical
focus. It does not remove the autonomous preloaded C180 route, whose full
\(O(\log q)\) history remains FFCC.

The surviving possibilities on the same one-cell geometry are:

1. prove that the physical C180 packet violates (0.3) in a controlled way
   and use that concentration or spectral tail as part of the leading stage;
2. use a nonlinear/depleting leading orbit rather than the frozen
   linearization;
3. use a terminal interval longer than \(O(1)\), paying its semigroup,
   residence, heat, and BAFL costs; or
4. close C179/PPRG, whose broad planar reservoir remains a genuinely
   time-ordered polarization/pressure problem.

For an \(O(\log q)\) interval, (2.2) may contribute a power of \(q\);
no no-focus conclusion is asserted. C161 puts its abstract splitter on one
normalized interval, but the exact C177 reservoir supplies no just-in-time
switch separating that splitter from the preceding C159 evolution. A
long-window proof must conjugate the complete sideband system into a smooth
C159 Floquet frame and control its off-diagonal pressure/polarization
connection. No one-cell stage or singularity is claimed.

## 6. Verification boundary

The dependency-free checker verifies the exact pressure-divergence
identity on nontrivial incompressible Fourier pairs, including the factor
two and order-zero pressure multiplier; the C180 mode-count/base-norm
ledger; the \(L^2\)--\(L^\infty\) interpolation exponent; the fractional
scale factor; the raw and \(J^2\)-taxed target ratios; and the
fixed-time/logarithmic-time claim boundary. It does not prove the standard
finite-\(p\) Mikhlin, transport-commutator, Sobolev, or
transport-diffusion theorems, establish any entrance hypothesis in (0.3),
construct localization/tails, or verify nonlinear depletion, C125, RIGM,
BAFL, or a one-cell stage.
