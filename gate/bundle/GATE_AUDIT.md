# Navier–Stokes breakthrough gate audit

## Bottom line

The earlier direct cascade computations established a real but one-step
critical transfer. They did not reproduce at the next shell, so that mechanism
has been retired rather than polished further.

The highest-upside route now is an embedding of Palasek's super-exponentially
separated Obukhov blow-up into the true Navier–Stokes nonlinearity. The audit
below reduces that program to one local construction:

> Build a divergence-free, time-dependent high-frequency cell whose resonant
> quadratic stress is a prescribed lower-frequency parent mode, while its
> nonresonant self-interaction is second order in the two scale ratios.

This is not yet a proof. It is, however, a falsifiable gate with an explicit
nonempty exponent window.

## 1. Target shell dynamics

The target ODE is

\[
X_k'=-\nu N_k^2X_k+
N_{k-1}^{\alpha}X_{k-1}X_k-
N_k^{\alpha}X_{k+1}^2+f_k,
\qquad N_k=N_0^{b^k}.
\]

Palasek proves forced viscous blow-up when

\[
b>1,\qquad 2b<\beta<\alpha,\qquad \alpha>2,
\]

with terminal amplitudes \(X_k(0)\asymp N_k^{\beta-\alpha}\).
The forcing can be smooth and vanish near the singular time. The open step is
to realize the two quadratic interactions inside the genuine Euler bilinear
form while controlling every other interaction.

## 2. Why the previous branches fail

### Direct geometric-shell recurrence

The validated direct \(Q_3\) crossing is a one-shell event. Optimizing the next
shell sacrifices the first crossing, and the two-shell recurrence does not
close. That is evidence against repeating one fixed geometric packet, not
against every possible singularity mechanism.

### Compact steady point packets

For a compactly supported steady Euler field,

\[
\int u_i u_j\,dx=-\delta_{ij}\int p\,dx.
\]

Its total Reynolds stress is isotropic. More directly, for any divergence-free
parent \(U\),

\[
\langle u,(u\cdot\nabla)U\rangle
=-\langle U,(u\cdot\nabla)u\rangle=0.
\]

Exact steadiness therefore removes the desired amplifier together with the
unwanted self-error.

### Steady Mikado tubes

A steady tube can localize in two transverse directions, giving at most the
critical value \(\alpha=2\). Viscous Palasek blow-up requires \(\alpha>2\).
Cutting the tube into a finite segment creates an endpoint/self-advection
error. If the child points into a contracting direction of the parent shear,
it crosses the parent's thin core too quickly; the resulting inequalities are
incompatible.

### Exact finite Fourier networks

Finite-mode Euler solutions have been classified: the only possibilities are
stationary 2D-like or Beltrami flows. Thus an exact finite Fourier cascade cell
cannot be the missing mechanism.

### Backward self-similar axisymmetric shortcut

A recent preprint claiming a torus singularity from one lifted scalar profile
does not pass the equation-level audit. The true axisymmetric-with-swirl system
has two dynamical quantities (swirl and azimuthal vorticity) coupled through a
Biot–Savart/streamfunction relation. Replacing the second evolution equation
by an elliptic recovery from the square of the first variable is not an exact
reduction. In addition, localized backward Leray self-similar profiles in the
claimed decay classes are covered by established nonexistence results.

## 3. Lagged intermittency

The useful adjustment is to distinguish a mode's carrier frequency from its
localization frequency.

Let shell \(k\) oscillate at \(N_k\), but let its envelope vary at
\(N_{k-1}\). It then occupies volume

\[
V_k\asymp N_{k-1}^{-3}=N_k^{-3/b}.
\]

Consequently

\[
\frac{\|u_k\|_\infty}{\|u_k\|_2}
\asymp N_{k-1}^{3/2}
=N_k^{3/(2b)}.
\]

The parent strain coefficient is

\[
N_{k-1}N_{k-2}^{3/2}
=N_{k-1}^{\,1+3/(2b)}.
\]

Thus the effective Obukhov exponent is

\[
\boxed{\alpha=1+\frac{3}{2b}}.
\]

This exceeds two when \(b<3/2\), and it is compatible with the viscous trapping
condition \(2b<\alpha\) precisely for

\[
1<b<\frac{1+\sqrt{13}}4\approx1.15139.
\]

## 4. An exact rational parameter certificate

Choose

\[
b=\frac{11}{10},\qquad
\alpha=\frac{26}{11},\qquad
\beta=\frac94.
\]

Then

\[
\beta-2b=\frac1{20}>0,\qquad
\alpha-\beta=\frac5{44}>0.
\]

The uncorrected same-scale envelope leakage is faster than the intended parent
amplification by the exponent

\[
(b-1)(\beta-1)=\frac18
\]

when measured in the parent frequency.

If a nonlinear WKB/normal-form construction cancels the nonresonant leakage to
second order, the two available small ratios leave positive margins:

\[
2(b-1)-(b-1)(\beta-1)=\frac3{40}>0
\]

for the carrier ratio, and

\[
2\left(1-\frac1b\right)-(b-1)(\beta-1)
=\frac5{88}>0
\]

for the envelope-hierarchy ratio.

First-order cancellation is insufficient. Second order is the sharp concrete
gate produced by this audit.

## 5. The missing cell lemma

Write \(q=N_{k-2}\ll p=N_{k-1}\ll N=N_k\). The desired cell consists of a
parent \(U_p\) and a child \(W_N\) such that:

1. \(U_p\) and \(W_N\) are exactly divergence-free.
2. \(W_N\) has carrier \(N\), envelope bandwidth \(p\), and
   \(L^\infty/L^2\asymp p^{3/2}\).
3. The resonant part of
   \(\mathbb P\operatorname{div}(W_N\otimes W_N)\) is a signed multiple of
   \(U_p\), with coefficient \(p\,q^{3/2}\).
4. By energy duality, the corresponding linearized parent-child interaction
   amplifies the child with the opposite coefficient.
5. All nonresonant same/high-frequency output is reduced by
   \(O((p/N)^2+(q/p)^2)\), after incompressibility and nonlinear WKB
   correctors.
6. Transport of the time-dependent envelope does not destroy the stress
   alignment over the activation interval.

Items 3 and 4 are the high-high-low/low-high energy pair. Item 5 is where a
second-order Beltrami or Kelvin-wave normal form must be derived. Item 6 is the
remaining Floquet/transport issue and is why an instantaneous coefficient
alone is not a proof.

## 6. Decisive computation

The companion GPU program searches directly for a finite-resolution analogue
of Item 5. For a normalized divergence-free child \(W\), it decomposes

\[
F(W)=\mathbb P\operatorname{div}(W\otimes W)
\]

into a prescribed parent band and its orthogonal leakage, and minimizes

\[
\frac{\|F_{\rm leak}\|_2}{\|F_{\rm parent}\|_2}
\]

while keeping the parent interaction strength nontrivial.

The ratio sweep uses carrier/parent separations \(4,5,8\). The make-or-break
diagnostic is

\[
\rho^2
\frac{\|F_{\rm leak}\|_2}{\|F_{\rm parent}\|_2},
\qquad \rho=N/p.
\]

If that quantity remains bounded (or decreases) while the normalized parent
strength stays comparable, the second-order cell deserves an analytic
construction. If it grows, or the optimizer obtains a small ratio only by
collapsing the parent interaction, this branch should be stopped.

## 7. What a positive result would and would not mean

A positive sweep would not solve Navier–Stokes. It would establish numerical
plausibility for the only cell order that fits a fully explicit viscous
blow-up parameter window. The next analytic tasks would be:

1. derive the second-order divergence and nonlinear correctors;
2. prove a time-dependent transport/Floquet estimate for the child cell;
3. place the cell estimates into a perturbed version of Palasek's trapping
   region;
4. prove that the residual defines a smooth force vanishing near blow-up;
5. pass from Galerkin truncations to a classical solution up to the singular
   time.

A negative sweep would be equally useful: it would close the current
Palasek-embedding branch before more large simulations are attempted.

## References

- C. Fefferman, *Existence and smoothness of the Navier–Stokes equation*
  (official Clay problem statement):
  <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
- S. Palasek, *Finite-time blow-up in an elementary model of the 3D
  Navier–Stokes equations* (2026): <https://arxiv.org/abs/2605.13827>
- M. P. Coiculescu and S. Palasek, *Non-uniqueness of smooth solutions of the
  Navier–Stokes equations from critical data*:
  <https://arxiv.org/abs/2503.14699>
- T. Buckmaster and V. Vicol, *Nonuniqueness of weak solutions to the
  Navier–Stokes equation* (intermittent Beltrami waves):
  <https://arxiv.org/abs/1709.10033>
- A. V. Gavrilov, *A steady Euler flow with compact support*:
  <https://arxiv.org/abs/1810.08020>
- P. Constantin and V. Vicol, *Remarks on high Reynolds numbers
  hydrodynamics and the Euler equations* (compact steady-flow stress
  identity): <https://arxiv.org/abs/1903.11699>
- N. Kishimoto and T. Yoneda, *Characterization of three-dimensional Euler
  flows supported on finitely many Fourier modes*:
  <https://arxiv.org/abs/2110.08039>
- B. R. Fabijonas and D. D. Holm, *Multi-frequency Craik–Criminale solutions
  of the Navier–Stokes equations*: <https://arxiv.org/abs/nlin/0304049>
- T. Y. Hou and C. Li, *Dynamic stability of the 3D axisymmetric
  Navier–Stokes equations with swirl* (the exact coupled axisymmetric
  equations): <https://arxiv.org/abs/math/0608295>
- Q. Jiu, Y. Wang, and W. Wei, *Leray's backward self-similar solutions to
  the 3D Navier–Stokes equations in Morrey spaces*:
  <https://arxiv.org/abs/2006.15776>
