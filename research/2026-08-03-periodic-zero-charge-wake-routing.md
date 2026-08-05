# Periodic zero-charge routing and the terminal wake condition

**Date:** 2026-08-03

**Status:** exact periodic symmetric anti-divergence, Gevrey-2 multiplier
bound, and terminal heat-wake condition.  The nonlinear endpoint map remains
open.

**Scope:** the zero-charge channel isolated in C85.  This note uses the
periodic branch of the Navier--Stokes problem.  It does not construct a
transition, erase the terminal wake, or prove singularity formation.

## 1. Outcome

On the three-torus the zero charge has a simpler compatibility theory than
a compactly supported Euclidean cell.  There is no angular-momentum
cokernel and no exterior-decay condition.  A periodic vector field is the
divergence of a periodic symmetric tensor if and only if its spatial mean
vanishes.

The inverse is explicit, has order minus one, and is bounded on the same
Gevrey-2 Fourier spaces used by the charge majorant.  Consequently the
three-phase recursion now has the following exact split.

* Every nonzero material charge is inverted by its phase ellipticity, as in
  C85.
* Every zero charge which is a projected divergence, Laplacian, or time
  derivative of a fixed-mean field is routed through a global periodic
  symmetric stress without a small divisor or factorial loss.

This does **not** make the wake disappear.  The corresponding forced heat
equation has a generally nonzero terminal state.  Requiring that state to
vanish imposes one weighted time-moment condition for every nonzero spatial
Fourier mode.  The exact remaining target is therefore a nonlinear
terminal wake/carry-forward theorem, rather than a zero-charge elliptic
inverse.

## 2. Exact periodic symmetric anti-divergence

Work on
\[
                 \mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3
\]
and use
\((\operatorname{div}R)_i=\partial_jR_{ij}\).  Let \(f\) be a smooth
periodic vector field with
\[
                         \int_{\mathbb T^3}f\,dx=0.             \tag{2.1}
\]
Let \(\Delta^{-1}\) denote the zero-mean inverse of the periodic
Laplacian and define
\[
 (\mathcal R f)_{ij}
 =\partial_i\Delta^{-1}f_j+\partial_j\Delta^{-1}f_i
  -\delta_{ij}\partial_k\Delta^{-1}f_k.                       \tag{2.2}
\]
Then \(\mathcal Rf\) is symmetric and
\[
\begin{aligned}
 \partial_j(\mathcal Rf)_{ij}
 &=\partial_i\partial_j\Delta^{-1}f_j+f_i
   -\partial_i\partial_k\Delta^{-1}f_k\\
 &=f_i.                                                        \tag{2.3}
\end{aligned}
\]

Conversely, the integral of a periodic divergence is zero, so (2.1) is
necessary.  Thus zero mean is the complete compatibility condition.

For \(k\ne0\), the Fourier symbol is
\[
 \widehat{\mathcal Rf}_{ij}(k)
 =-{i\over|k|^2}
   \left(k_i\widehat f_j+k_j\widehat f_i
          -\delta_{ij}k\cdot\widehat f\right).                 \tag{2.4}
\]
It obeys
\[
 i k_j\widehat{\mathcal Rf}_{ij}=\widehat f_i,
 \qquad
 \|\widehat{\mathcal Rf}(k)\|_{\rm op}
 \le {3\over|k|}|\widehat f(k)|.                              \tag{2.5}
\]
If \(f\) is divergence free, then \(\mathcal Rf\) is trace free.  The
formula commutes with translations, orthogonal lattice symmetries, and
time differentiation, so the reflection and helicity symmetries already
used by the transition ledger can be retained.

## 3. Gevrey-2 tameness

For \(\tau>0\) and \(s\in\mathbb R\), define
\[
 \|f\|_{\mathcal G^2_{\tau,s}}^2
 =\sum_{k\in\mathbb Z^3}
   e^{2\tau|k|^{1/2}}\langle k\rangle^{2s}|\widehat f(k)|^2.   \tag{3.1}
\]
Equations (2.4)--(2.5) give
\[
 \|\mathcal Rf\|_{\mathcal G^2_{\tau,s+1}}
 \le C\|f\|_{\mathcal G^2_{\tau,s}}.                         \tag{3.2}
\]
The Leray projector is bounded on the same space.  For \(s>3/2\), this is
an algebra because
\[
 |k|^{1/2}\le|\ell|^{1/2}+|k-\ell|^{1/2}.                     \tag{3.3}
\]
Thus a bilinear zero-charge source has the same Gevrey convolution bound
as the nonzero charges.  Applying \(\mathcal R\) gains one spatial
derivative and introduces no additional factorial.

There is also an exact frequency split.  If \(P_{\ge N}\) is the spatial
Fourier projector onto \(|k|\ge N\), then
\[
 \|\mathcal RP_{\ge N}f\|_{\mathcal G^2_{\tau,s}}
 \le {C\over N}\|P_{\ge N}f\|_{\mathcal G^2_{\tau,s}}.        \tag{3.4}
\]
High spatial frequencies therefore receive the expected scale gain.  Low
spatial modes do not.  They are precisely the global periodic wake and
must be retained rather than estimated as a localized order-minus-one
error.

## 4. Why every formal zero-charge residual is compatible

Let \(v(t)\) be periodic, divergence free, and of fixed spatial mean.  Let
\(Q_0(t,x)=Q_0^T\) be the zero material charge produced by the three-phase
covariance and lower-order interactions.  The projected zero-charge
Navier--Stokes residual has the form
\[
 F_0
 =\partial_tv+\mathbb P\operatorname{div}(v\otimes v+Q_0)
  -\nu\Delta v.                                                \tag{4.1}
\]
Every term on the right has zero mean.  It is also divergence free.
Consequently
\[
 R_0=\mathcal RF_0                                             \tag{4.2}
\]
is a global periodic, symmetric, trace-free tensor satisfying
\[
                         \operatorname{div}R_0=F_0.            \tag{4.3}
\]

The trace-free property is useful: an isotropic pressure gauge can be
added independently.  If the leading covariance has a strict positive
margin, every sufficiently small higher-order \(R_0\) is absorbed into
that margin without changing its divergence or its work against an
incompressible parent.

The exact work identity remains
\[
 \int_{\mathbb T^3}R_0:\operatorname{sym}\nabla U\,dx
 =-\int_{\mathbb T^3}F_0\cdot U\,dx.                           \tag{4.4}
\]
Thus the periodic inverse does not create a spurious energy degree of
freedom.  It routes the zero charge while preserving the mandatory work
ledger.

Equations (3.2)--(4.3), together with C85, remove a formal all-order
solvability obstruction: neither the zero charge nor any nonzero charge
causes a Gevrey-2 small divisor or factorial proliferation.  This is a
statement about the Reynolds/WKB hierarchy, not yet an exact
Navier--Stokes transition.

## 5. Exact global wake and its terminal condition

The velocity version of the global routing is the forced heat wake.  With
zero incoming wake, set
\[
 z(t)=\int_{t_-}^t e^{\nu(t-s)\Delta}F_0(s)\,ds.                \tag{5.1}
\]
Then
\[
 \partial_tz-\nu\Delta z=F_0,
 \qquad \operatorname{div}z=0,
 \qquad \int z\,dx=0.                                        \tag{5.2}
\]
The heat semigroup is a contraction in (3.1), so
\[
 \|z(t)\|_{\mathcal G^2_{\tau,s}}
 \le\int_{t_-}^t\|F_0(s)\|_{\mathcal G^2_{\tau,s}}\,ds.      \tag{5.3}
\]
If \(F_0\) is flat at \(t_-\), then (5.2) recursively shows that \(z\) is
flat there as well.

At the outgoing endpoint, however,
\[
 \widehat z(k,t_+)
 =\int_{t_-}^{t_+}
   e^{-\nu|k|^2(t_+-s)}\widehat F_0(k,s)\,ds,
 \qquad k\ne0.                                                 \tag{5.4}
\]
Therefore a bare zero-wake endpoint requires the infinite family
\[
 \boxed{
 \int_{t_-}^{t_+}
   e^{-\nu|k|^2(t_+-s)}\widehat F_0(k,s)\,ds=0
 \quad\text{for every }k\ne0.}                                \tag{5.5}
\]
There is no reason for (5.5) to hold in a generic localized transition.
This recovers the earlier endpoint-jet obstruction in an exact periodic
Fourier form.

The viable endpoint state is instead
\[
                         z_+=z(t_+),                            \tag{5.6}
\]
carried into the next stage.  A construction must either propagate this
state through the complete nonlinear endpoint map or supply enough
time-dependent controls to meet (5.5).  The first option is consistent
with the mandatory work-carrying wake already identified in C36 and C61.

## 6. Remaining theorem

For the periodic route, the all-order phase problem is no longer an
elliptic existence problem.  The remaining load-bearing statement is:

> Couple the strain-amplified three-phase packet to the global zero-charge
> mild equation, carry its nonzero terminal state into the next scale, and
> prove a uniform Gevrey-2 endpoint map whose derivative contains the
> rank-five affine capture chart.

The nonlinear terms involving the wake are of turnover size, not
automatically perturbative.  The note therefore does not claim that a
standard short-time contraction is uniform across the cascade.  It only
removes the periodic compatibility, small-divisor, and factorial-growth
parts of the zero-charge obstruction and states the exact remaining
terminal condition.
