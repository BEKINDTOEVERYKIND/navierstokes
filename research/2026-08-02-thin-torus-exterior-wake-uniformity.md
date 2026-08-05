# Aspect-uniform extraction of the thin Gavrilov second-jet wake

**Date:** 2026-08-02
**Status:** self-derived; exact scaling and tensor algebra checked; not
cross-audited
**Scope:** the second Navier--Stokes jet only.  The result is uniform for
the exterior wake after it is split from the singular local thin-layer
piece.  It is not an all-order endpoint theorem.

## 1. Result

Let \(\varepsilon\) denote the geometric minor-to-major radius ratio of a
thin Gavrilov torus.  Thus the localization parameter in
\(\chi(p/\delta)\), with \(p\sim r^2/2\) near the seed circle, is
\(\delta=\varepsilon^2\).  Normalize the localized velocity by its actual
pointwise size \(v\), let \(R\) be its major radius, and let
\(\tau=T R/v\).

The thin profile has transverse derivative scale
\((\varepsilon R)^{-1}\), but its second-jet pressure tensor satisfies the
aspect-uniform estimate

\[
 \left\|\nu\left(
 U\otimes\Delta U+\Delta U\otimes U\right)\right\|_{L^1}
 \le C\nu v^2R.                                           \tag{1.1}
\]

Consequently, at distance \(d\ge C_0R\),

\[
 |\nabla^m(\tau^2u_{2}^{\rm ext})(x)|
 \le C_m\nu T^2R^3d^{-4-m},                              \tag{1.2}
\]

with constants independent of \(\varepsilon\).  Moreover the actual
exterior field admits a divergence-free extension \(Z_2\), equal to the
second velocity jet outside a fixed enlargement of the torus ball, such
that

\[
 \|\nabla^m(\tau^2Z_2)\|_{L^2}
 \le C_m\nu T^2R^{1/2-m}.                                \tag{1.3}
\]

Thus the aspect-ratio caveat in the separated second-wake estimate is
removable: no moment cancellation is needed for a uniform \(d^{-4}\)
bound or for a uniform exterior \(H^m\) bound.

There is still a real nonuniform local object.  The raw pressure-gradient
piece inside the enclosing ball has the estimate

\[
 \|\nabla^m\nabla p_1\|_2
 \le C_m\nu v^2R^{-3/2-m}\varepsilon^{-2-m},              \tag{1.4}
\]

and its leading flat-cylinder profile is nonzero.  It must be recombined
with the compact nonlinear part of the second jet before doing thin-limit
local estimates.  Formula (1.3) applies to the extracted exterior wake,
not to the raw pressure term on the core.

For a packed stage of \(K^3\) microtori of major radius
\(\delta=\ell/K\), common actual speed \(v\), and duration
\(\tau=T\ell/v\), absolute summation gives the aspect-uniform outer-stage
bound

\[
 |\nabla^m Z_{2,\rm stage}(x)|
 \le C_m\nu T^2\ell^3K^2d^{-4-m}.                        \tag{1.5}
\]

At a core in the same packed stage the nearest-neighbor lattice instead
costs

\[
 |Z_{2,\rm same}(x)|\le C\nu T^2{K^3\over\ell}.           \tag{1.6}
\]

For \(v_j=\ell_j^{-\gamma}K_j^\gamma\),
\(1<\gamma<3/2\), geometric \(\ell_j\), polynomial \(K_j\), and bounded
or polynomial \(T_j\), both ratios to \(v_j\),

\[
 \nu T_j^2\ell_j^{\gamma-1}K_j^{2-\gamma},
 \qquad
 \nu T_j^2\ell_j^{\gamma-1}K_j^{3-\gamma},               \tag{1.7}
\]

tend to zero.  The packed exterior second wake therefore remains
summable even when \(\varepsilon_j\to0\).  The same-stage \(K^3\) cost,
which is one power worse than the previously recorded inter-stage
\(K^2\) cost, has to be included in a nonlinear wake--core induction.

## 2. Uniform thin-profile scaling

Let \(W_\varepsilon^{\rm nat}\) be Gavrilov's natural localization

\[
 W_\varepsilon^{\rm nat}
 =\chi(p/\varepsilon^2)u
\]

around the unit seed circle.  The analytic expansion used in the
second-jet obstruction is

\[
 u=O(r),\qquad p={r^2\over2}+O(r^3),
\]

where \(r\) is distance to the circle.  In fixed rescaled tubular
coordinates \(y=r/\varepsilon\), repeated chain rule gives, for every
fixed \(m\),

\[
 \begin{aligned}
 \left|\operatorname{supp}W_\varepsilon^{\rm nat}\right|
     &\le C\varepsilon^2,\\
 \|\nabla^mW_\varepsilon^{\rm nat}\|_\infty
     &\le C_m\varepsilon^{1-m}.                          \tag{2.1}
 \end{aligned}
\]

Define the unit-speed profile

\[
 \widehat W_\varepsilon=\varepsilon^{-1}
 W_\varepsilon^{\rm nat}
\]

and the physical torus

\[
 U(x)=vQ\widehat W_\varepsilon
 \left(Q^T{x-c\over R}\right).
\]

Then

\[
 \|\nabla^mU\|_\infty
 \le C_mvR^{-m}\varepsilon^{-m},\qquad
 |\operatorname{supp}U|\le CR^3\varepsilon^2.            \tag{2.2}
\]

If

\[
 S=\nu\left(U\otimes\Delta U+\Delta U\otimes U\right),
\]

then (2.2) gives the more general estimate

\[
 \|\nabla^mS\|_1
 \le C_m\nu v^2R^{1-m}\varepsilon^{-m}.                  \tag{2.3}
\]

The case \(m=0\) is (1.1).  It is the cancellation of the
\(\varepsilon^{-2}\) transverse Laplacian with the
\(\varepsilon^2\) tube volume.  This is a scaling cancellation, not a
pressure-moment cancellation.

This normalization is independent of notation for amplitude.  If the
natural multiplier is \(A\), then the actual speed is \(v\asymp
A\varepsilon\); substitution into (1.1) gives
\(\|S\|_1=O(\nu A^2R\varepsilon^2)\).  If instead energy per major-radius
cell is held fixed, then \(v\) grows like \(\varepsilon^{-1}\).  The
instantaneous pressure moment grows like \(\varepsilon^{-2}\), but the
major-turnover time \(R/v\) shrinks by \(\varepsilon\), so the normalized
coefficient \(\tau^2\|S\|_1\) is still uniform.  In general the only
duration factor is

\[
 T={\tau v\over R},\qquad
 \tau^2\|S\|_1\le C\nu T^2R^3.                           \tag{2.4}
\]

An aspect loss occurs only if the transition is deliberately made
\(\varepsilon^{-1}\) major turnovers long, not from the thin geometry
itself.

## 3. Exact second jet and the local/exterior split

At a steady Euler endpoint \(u(0)=U\), the first two formal
Navier--Stokes velocity jets are

\[
 u_1=\nu\Delta U,
\]

\[
 u_2=\nu^2\Delta^2U
 -\nu\left((U\cdot\nabla)\Delta U
          +(\Delta U\cdot\nabla)U\right)
 -\nabla p_1,                                            \tag{3.1}
\]

where

\[
 -\Delta p_1=\partial_i\partial_jS_{ij}.                 \tag{3.2}
\]

Outside the support of \(U\), all local terms in (3.1) vanish and

\[
 u_2=-\nabla p_1.                                        \tag{3.3}
\]

The global \(L^2\) estimate for the raw pressure piece follows from the
order-one Fourier multiplier in
\(\nabla\partial_i\partial_j(-\Delta)^{-1}\):

\[
 \|\nabla^m\nabla p_1\|_2
 \le C_m\|\nabla^{m+1}S\|_2
 \le C_m\nu v^2R^{-3/2-m}\varepsilon^{-2-m},             \tag{3.4}
\]

which is (1.4).

The loss in (3.4) is not merely a crude support estimate.  The leading
flat-cylinder profile from Gavrilov's Taylor expansion is, in transverse
polar coordinates,

\[
 V_\perp=-r\chi(s)e_\theta,\qquad
 V_\parallel={r\chi(s)\over\sqrt2},\qquad
 s={r^2\over2}.
\]

Writing

\[
 A(s)=4\chi'(s)+2s\chi''(s),
\]

one computes exactly

\[
 (V\cdot\nabla)\Delta V+(\Delta V\cdot\nabla)V
 =-2r\chi(s)A(s)e_r.                                    \tag{3.5}
\]

This is a nonzero radial gradient for every nonzero smooth compactly
supported \(\chi\).  Indeed \(A=0\) would imply
\((s^2\chi')'=0\), which is incompatible with nonzero compact support.
The leading local pressure gradient cancels (3.5) in the projected
velocity jet, but the pressure piece by itself has the
\(\varepsilon^{-2}\) \(L^2\) scale in (3.4).  Separating that raw piece
globally is therefore the wrong thin-uniform decomposition.

## 4. Uniform far field and a solenoidal exterior extension

Let \(\Gamma(x)=(4\pi|x|)^{-1}\).  From (3.2),

\[
 p_1=\partial_i\partial_j\Gamma*S_{ij}.
\]

Because the entire torus lies in a ball of radius \(C R\), (1.1) and the
Newton-kernel derivative estimates give, for
\(d=|x-c|\ge C_0R\),

\[
 |\nabla^m(-\nabla p_1)(x)|
 \le C_m\nu v^2R\,d^{-4-m}.                             \tag{4.1}
\]

Every constant here is independent of \(\varepsilon\).

There is also an aspect-uniform global representative of precisely this
exterior tail.  Choose a radial cutoff \(\zeta\) that is zero on a ball
containing the torus and one outside a fixed larger ball.  Put

\[
 E=-\nabla p_1,\qquad f=\nabla\zeta\cdot E.
\]

The exterior field is divergence free, and its flux through every
enclosing sphere is zero because the right side of (3.2) has zero total
mass.  Hence \(f\) has mean zero on the fixed transition annulus.  A
compact-support Bogovskii solution \(b\) on that annulus satisfies

\[
 \operatorname{div}b=f
\]

with the scale-covariant estimates of the standard annular operator.
After dilation by \(R\), this is one fixed annulus; its operator constants
therefore have no \(\varepsilon\)-dependence.
Therefore

\[
 Z_2=\zeta E-b                                             \tag{4.2}
\]

is divergence free, is zero in a smaller enclosing ball, and equals the
true velocity jet \(u_2\) outside the larger ball.  Integrating (4.1) and
using the annular Bogovskii estimate yields

\[
 \|\nabla^mZ_2\|_2
 \le C_m\nu v^2R^{-3/2-m}.                               \tag{4.3}
\]

Multiplying by \(\tau^2=T^2R^2/v^2\) proves (1.3).
The complementary field \(u_2-Z_2\) is compactly supported in the
enlarged torus ball.  Its thin norms are not uniform; it belongs to the
local carrier problem rather than the separated wake problem.

On a fixed flat three-torus, write the periodic Green function in the
coordinate chart as \(\Gamma_{\mathbb T^3}=\Gamma+H\), with \(H\)
smooth.  The image contribution obeys

\[
 \|\nabla^m(\tau^2\nabla\partial_i\partial_jH*S_{ij})
   \|_\infty
 \le C_m\nu T^2R^3,                                     \tag{4.4}
\]

again uniformly in \(\varepsilon\).  It is a smooth global centre field,
not a decaying exterior tail, and must be retained in a periodic endpoint
map.

## 5. The exact leading moment and optional cancellations

The zeroth tensor moment is

\[
 M_{ij}:=\int S_{ij}\,dx
 =-2\nu{\mathsf G}(U)_{ij},\qquad
 {\mathsf G}(U)_{ij}
 =\int\partial_kU_i\,\partial_kU_j\,dx.                  \tag{5.1}
\]

For the thin normalized profile,

\[
 {\mathsf G}(U)
 =v^2R\varepsilon^{-2}
 Q{\mathsf G}(W_\varepsilon^{\rm nat})Q^T.              \tag{5.2}
\]

The calculation in the second-jet obstruction gives

\[
 {\mathsf G}(U)^{\rm TF}
 =v^2R\left[
 D_\chi\left(n\otimes n-\frac13I\right)+o(1)\right],
\qquad
 D_\chi=\pi^2\int_0^\infty\chi(s)^2\,ds>0.               \tag{5.3}
\]

Thus the aspect-uniform leading tail is generally nonzero.  Its pressure
term is

\[
 p_1(x)
 =-\frac{\nu}{2\pi d^3}
 \left(3\widehat d\otimes\widehat d-I\right):
 {\mathsf G}(U)
 +O(\nu v^2R^2d^{-4}).                                  \tag{5.4}
\]

No cancellation is required for (4.1).  To improve the velocity decay
from \(d^{-4}\) to \(d^{-5}\), the necessary zeroth-order condition is

\[
 \left(\sum_\alpha{\mathsf G}(U_\alpha)\right)^{\rm TF}=0.
\tag{5.5}
\]

For equal axisymmetric tori, three rotations with mutually orthogonal
axes satisfy (5.5) exactly for every \(\varepsilon\), not only in the
limit.  For sufficiently thin members, where the common anisotropy is
nonzero and has one sign, two positive copies cannot do so.

To remove the next \(d^{-5}\) velocity coefficient as well, the first
tensor moments of \(S\) must vanish modulo the contractions killed by
harmonicity.  A simple sufficient construction is to take, for each of
the three orthogonal orientations, a centrally inverted pair about the
cluster centre.  Central inversion leaves its zeroth tensor moment
unchanged and reverses every first moment.  The resulting six-torus
cluster has an isotropic total zeroth moment and zero full first moment.
Its exterior velocity therefore begins at order

\[
 |u_2^{\rm cluster}(x)|
 \le C\nu v^2R^3d^{-6}.                                 \tag{5.6}
\]

Making the exterior pressure identically zero would require the full
infinite family

\[
 \sum_\alpha\int S_\alpha:\nabla^2q\,dx=0
\]

for every harmonic polynomial \(q\).  Finite quadrupole cancellation
does not supply those all-order moment identities.

## 6. Packed-stage ledger

Take \(N\lesssim K^3\) microtori, each of major radius
\(\delta\asymp\ell/K\), in a macro region of diameter \(\ell\).  Use
actual speed \(v\), a fixed support/separation margin between the
enclosing microballs, and stage duration

\[
 \tau=T{\ell\over v}.
\]

For one microtorus, (4.1) after multiplication by \(\tau^2\) has moment
coefficient

\[
 \tau^2\nu v^2\delta
 =\nu T^2\ell^2\delta.                                  \tag{6.1}
\]

Absolute summation proves

\[
 |\nabla^mZ_{2,\rm stage}(x)|
 \le C_m\nu T^2\ell^3K^2d^{-4-m},                       \tag{6.2}
\]

which is (1.5).  The triangle estimate for the extracted exterior fields
is likewise

\[
 \|\tau^2Z_{2,\rm stage}\|_2
 \le C\nu T^2\ell^{1/2}K^{9/2}.                         \tag{6.3}
\]

Both estimates are independent of the thin aspect ratio.

There are two distinct interaction distances:

1. For a different stage at distance comparable to \(\ell\), (6.2)
   gives size \(C\nu T^2K^2/\ell\).
2. At a core in the same packed stage, summing the
   three-dimensional lattice kernel
   \(\sum_{q\ge1}q^2(q\delta)^{-4}\) gives
   \(C\nu T^2K^3/\ell\), proving (1.6).

The second cost is not part of the separated-annulus estimate and must
not be silently identified with it.  Nevertheless (1.7) shows that both
are perturbative for the geometric/polynomial cascade.

Optional local clustering can reduce the outer polynomial loss.  Group
the microtori into six-torus clusters of diameter \(O(\delta)\) satisfying
the zeroth- and first-moment cancellations in Section 5.  The first
surviving moment is then \(O(\nu v^2\delta^3)\).  At macro distance
\(d\asymp\ell\), summing \(O(K^3)\) clusters after the stage time gives

\[
 |Z_{2,\rm clustered}(x)|
 \le C{\nu T^2\over\ell},                               \tag{6.4}
\]

removing both powers of \(K\) from the outer wake.  This does not improve
the field seen by partners inside the same cluster, so it is an optional
far-field optimization rather than an all-order closure.

The minor-scale viscous parameter is

\[
 \Theta={\nu\tau\over(\varepsilon\ell/K)^2}
 ={\nu T K^2\over v\ell\varepsilon^2}.                  \tag{6.5}
\]

The outer and same-stage wake ratios can be written

\[
 {\nu T^2K^2\over v\ell}=T\varepsilon^2\Theta,\qquad
 {\nu T^2K^3\over v\ell}=TK\varepsilon^2\Theta.          \tag{6.6}
\]

Hence the thin limit does not enlarge the exterior second wake.  It does
make the local viscous/carrier problem harder through
\(\varepsilon^{-2}\), which is a separate obligation.

## 7. Claim boundary

The second-jet separated wake now has an aspect-uniform formulation:

* normalize by actual velocity and record the duration \(T\);
* extract only the solenoidal exterior field by (4.2);
* keep the singular local pressure/nonlinear cancellation inside the
  carrier block;
* budget \(K^2\) for inter-stage and \(K^3\) for same-stage absolute
  interactions, or use six-copy clusters to remove the outer \(K^2\);
* do not infer compact pressure from finitely many moment cancellations.

What remains open is the decisive higher-jet statement.  At later orders
the local thin derivatives feed new pressure sources, and (1.1) has not
been generalized to a uniform Gevrey majorant through
\(M_j\asymp j^2/\log j\).  The exterior split prevents the already-known
second jet from blocking the cascade, but it does not construct the
transition.

## Primary source

* A. V. Gavrilov,
  [A steady Euler flow with compact support](https://arxiv.org/abs/1810.08020),
  especially the local expansion near the seed circle and the final
  pressure localization.
