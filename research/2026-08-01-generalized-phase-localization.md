# Exact curl localization for a two-phase child writer

Date: 2026-08-01

## Claim boundary

This note proves exact algebraic identities for a compactly localizable
oscillatory block.  It does **not** construct a Navier--Stokes singularity.
Its purpose is narrower: decide whether localizing the nonlinear phase of a
high--high-to-low child writer necessarily creates an error as large as the
desired low-frequency output.

The answer is asymmetric.  Compact transverse phase collars can create an
`O(Q)` **low-frequency buffer wake**, but every high-frequency and
self-interaction remainder displayed below retains at least one inverse
carrier factor.  Compact localization in the carrier direction forces a
nonzero `O(Q^2/Lambda)` high remainder somewhere; it does not force an
`O(Q)` high remainder.

All identities in Sections 1--5 are direct calculations in the constant
amplitude core and are independently checked by
`checks/generalized_phase_localization.py`.

## 1. An exactly divergence-free two-phase block

Let `t,w` be constant vectors with

\[
 |t|=1,\qquad t\cdot w=0,
\]

and put

\[
 c=w\times t,\qquad \varepsilon=\frac1{2\Lambda}.
\]

For a real phase difference `Psi`, define

\[
 q=\nabla\Psi,\qquad q_t=t\cdot q,\qquad
 d=w\cdot q,\qquad v=q_t w-dt,
\tag{1.1}
\]

\[
 \phi_\sigma=\Lambda t\cdot x+\frac\sigma2\Psi,
 \qquad \sigma\in\{+1,-1\}.
\tag{1.2}
\]

For any smooth scalar envelope `A`, set

\[
 U_\sigma=\frac1{i\Lambda}\nabla\times
 \left(Ac\,e^{i\phi_\sigma}\right).
\tag{1.3}
\]

Then `div U_sigma=0` exactly and

\[
 U_\sigma=e^{i\phi_\sigma}
 \left[Aa_\sigma+\frac{\nabla A\times c}{i\Lambda}\right],
 \qquad
 a_\sigma=w+\sigma\varepsilon v.
\tag{1.4}
\]

On the core `A=1`,

\[
 a_\sigma=\frac{\nabla\phi_\sigma}{\Lambda}\times c,
 \qquad
 a_\sigma\cdot\nabla\phi_\sigma=0,
 \qquad
 a_\sigma\cdot q=d.
\tag{1.5}
\]

Thus phase localization is incorporated into the polarization without an
approximate Leray projection.

## 2. Exact low, cross-high, and self interactions

Write `B(u,v)=(u dot grad)v`.  On `A=1`, direct differentiation gives

\[
\begin{aligned}
 \mathcal L
 &: =B(U_+,\overline{U_-})+B(\overline{U_-},U_+)\\
 &=e^{i\Psi}\left(2id\,w-
       \frac{(v\cdot\nabla)v}{2\Lambda^2}\right),\\[1ex]
 \mathcal H
 &: =B(U_+,U_-)+B(U_-,U_+)\\
 &=e^{2i\Lambda t\cdot x}\left(
       \frac{id}{\Lambda}v-
       \frac{(v\cdot\nabla)v}{2\Lambda^2}\right),\\[1ex]
 \mathcal S_\sigma
 &: =B(U_\sigma,U_\sigma)\\
 &=e^{2i\phi_\sigma}\left(
       \frac{\sigma}{2\Lambda}(w\cdot\nabla)v+
       \frac{(v\cdot\nabla)v}{4\Lambda^2}\right).
\end{aligned}
\tag{2.1}
\]

These formulas separate the desired low beat from the two types of unwanted
high output before any asymptotic estimate is made.

## 3. Explicit pressure decomposition

For the cross-high interaction take

\[
 \pi_H=-\frac{d^2}{2\Lambda^2}e^{2i\Lambda t\cdot x}.
\]

Then

\[
 \mathcal H-\nabla\pi_H
 =e^{2i\Lambda t\cdot x}\left[
   \frac{idq_t}{\Lambda}w+
   \frac{\nabla(d^2)-(v\cdot\nabla)v}{2\Lambda^2}
 \right].
\tag{3.1}
\]

Hence localization along `t` creates the generically sharp term

\[
 \frac{idq_t}{\Lambda}w=O(Q^2/\Lambda).
\tag{3.2}
\]

Relative to an `O(Q)` low daughter, this is `O(Q/Lambda)`.  If `q_t=0`,
the entire `1/Lambda` term vanishes and the cross-high remainder is

\[
 O\!\left(\frac{Q\|\nabla^2\Psi\|}{\Lambda^2}\right).
\tag{3.3}
\]

For self-interaction define

\[
 m=D_wD_t\Psi=D_td,\qquad n=D_w^2\Psi=D_wd,
\]

and take

\[
 \pi_\sigma=\frac{\sigma i n}{4\Lambda^2}e^{2i\phi_\sigma}.
\]

Then

\[
 \mathcal S_\sigma-\nabla\pi_\sigma
 =e^{2i\phi_\sigma}\left[
   \frac{\sigma m}{2\Lambda}w+
   \frac{(v\cdot\nabla)v+nq-\sigma i\nabla n}{4\Lambda^2}
 \right].
\tag{3.4}
\]

In particular,

\[
 \|\mathbb P\mathcal S_\sigma\|
 \lesssim \frac{|D_wD_t\Psi|}{\Lambda}
 +\frac{Q\|\nabla^2\Psi\|+\|\nabla^3\Psi\|}{\Lambda^2}.
\tag{3.5}
\]

Again, `q_t=0` removes the `1/Lambda` term.

## 4. The quadratic affine-strain core

Let `r,h,t` be orthonormal, set `s=r dot x`, and choose

\[
 w=r+Hh,\qquad \Psi=as^2.
\tag{4.1}
\]

Then

\[
 q=2asr,\qquad d=2as,qquad q_t=0,qquad v=-2ast.
\tag{4.2}
\]

The low interaction becomes

\[
 \mathcal L=4ias\,w\,e^{ias^2}
 =2\nabla e^{ias^2}+4iaHs\,h\,e^{ias^2}.
\tag{4.3}
\]

Thus the exact post-pressure daughter is the affine ridge

\[
 4iaHs\,h\,e^{ias^2}.
\tag{4.4}
\]

After subtracting the pressures from Section 3, the cross-high and self
remainders are respectively

\[
 \frac{4a^2s}{\Lambda^2}r\,e^{2i\Lambda t\cdot x},
 \qquad
 \frac{a^2s}{\Lambda^2}r\,
 e^{i(2\Lambda t\cdot x+\sigma as^2)}.
\tag{4.5}
\]

The earlier one-dimensional chirp calculation is therefore recovered
exactly, including its `psi' psi''/Lambda^2` obstruction.

## 5. Exact high-darkness no-go under carrier localization

The residual (3.2) is not merely an artifact of the curl polarization.
At a point, let arbitrary principal polarizations satisfy

\[
 b_\pm\perp(\Lambda t\pm q/2),
 \qquad \alpha_\pm=b_\pm\cdot q.
\]

The principal low and cross-high vectors are

\[
 L_0=\alpha_+b_-+\alpha_-b_+,
 \qquad
 H_0=-\alpha_+b_-+\alpha_-b_+.
\tag{5.1}
\]

One always has `H_0 dot q=0`.  If `q_t` is nonzero and the high term is pure
pressure, then `H_0` must be parallel to `t`; hence `H_0=0`.  But `H_0=0`
forces `L_0=0`: if both coefficients are nonzero, the polarizations are
parallel, therefore perpendicular to both carrier covectors and to `q`, a
contradiction; the cases with a zero coefficient are immediate.  Nonzero
daughter production and exact high-pressure cancellation therefore cannot
coexist where `q_t` is nonzero.

Compact localization in the `t` direction forces such a region.  Indeed, if
`d q_t=0` everywhere, then on the open set `{d != 0}` one has `q_t=0` and

\[
 D_td=D_tD_w\Psi=D_wD_t\Psi=D_wq_t=0.
\]

Along each `t`-line, `d` is constant on every component of `{d != 0}`.
Continuity and `d=0` for large `|t dot x|` rule out every nonzero component.
Thus `d=0`, contradicting the affine core.  The `O(Q^2/Lambda)` high leakage
is qualitatively unavoidable, although this argument gives no quantitative
lower bound: a long collar can reduce its coefficient.

## 6. Transverse buffer routing

For `Psi=as^2 chi`,

\[
 d=D_w\Psi=2as\chi+as^2D_w\chi.
\tag{6.1}
\]

A generic transverse cutoff creates the unsuppressed second term.  There is,
however, an exact characteristic direction.  With

\[
 \zeta=(h-Hr)\cdot x,\qquad \chi=\chi(t\cdot x,\zeta),
\]

one has `D_w chi=0`, so `d=2as chi` exactly.  The direction of `q` still
changes in this collar and can rotate the projected low daughter by `O(Q)`.
This is a real low-frequency buffer wake, not high-frequency leakage.  It
must be routed into a harmless region or included in the next-stage
bookkeeping.

The carrier envelope can be tapered only after `q` has returned to zero.  In
that separated region the exact curl correction in (1.4) removes the naive
first-order divergence defect; the remaining cutoff interactions carry
inverse powers of `Lambda`.  A full Sobolev estimate for this separated
envelope is still required.

## 7. Consequence for the cascade program

This calculation eliminates one feared obstruction and sharpens the real
one:

1. compact phase localization does not force fatal `O(Q)` high/self error;
2. its unavoidable carrier-direction high error is perturbative when
   `Q/Lambda -> 0`;
3. transverse collars create `O(Q)` low wakes, so geometry and support
   routing remain load-bearing;
4. the explicit sources in (3.1) and (3.4) are candidates for a one-order
   WKB corrector;
5. none of this supplies the small effective parent--child coupling or the
   finite-curvature Floquet block needed for a complete transition.

## Primary-source context

- A. V. Gavrilov, *A steady Euler flow with compact support*:
  https://arxiv.org/abs/1810.08020
- A. Cheskidov, M. Dai, and S. Palasek, *Instantaneous Type I blow-up and
  non-uniqueness of smooth solutions of the Navier--Stokes equations*:
  https://arxiv.org/abs/2511.09556
- S. Palasek, *Finite-time blow-up in an elementary model of the 3D
  Navier--Stokes equations*: https://arxiv.org/abs/2605.13827

