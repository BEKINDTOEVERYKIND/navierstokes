# The third thin-torus Navier--Stokes jet: the first aspect loss

**Date:** 2026-08-02
**Status:** exact recursion, moment algebra, profile integrals, and exponents
independently cross-audited; the one-viscosity tubular projection lemma is
conditional
**Scope:** one steady Gavrilov endpoint and the exterior part of its third
time jet.  No all-order or nonlinear transition theorem is claimed.

## 1. Result

Let \(U\) be a steady Euler bubble localized to a torus of major radius
\(R\), geometric minor-to-major aspect
\(\varepsilon\), and actual velocity size \(v\).  Let

\[
 u_n=\partial_t^nu(0),\qquad p_n=\partial_t^np(0)
\]

denote the formal unforced Navier--Stokes jets at \(u_0=U\), and use the
duration

\[
 \tau=T{R\over v}.
\]

The third velocity jet is

\[
u_3=\nu\Delta u_2
-\mathbb P\left(
 (U\cdot\nabla)u_2
 +2(u_1\cdot\nabla)u_1
 +(u_2\cdot\nabla)U\right),                              \tag{1.1}
\]

where

\[
 u_1=\nu\Delta U,
\]

\[
 u_2=\nu^2\Delta^2U-\nu\mathbb PF,\qquad
 F=(U\cdot\nabla)\Delta U+(\Delta U\cdot\nabla)U.        \tag{1.2}
\]

The second pressure derivative satisfies

\[
 -\Delta p_2=\partial_i\partial_j({\mathcal T}_2)_{ij},
\]

\[
 {\mathcal T}_2
 =U\otimes u_2+2u_1\otimes u_1+u_2\otimes U.            \tag{1.3}
\]

There are two distinct pressure-source blocks:

\[
 {\mathcal T}_2
 =\nu^2{\mathcal A}_2-\nu{\mathcal C}_2,                \tag{1.4}
\]

\[
 \begin{aligned}
 {\mathcal A}_2
 &=U\otimes\Delta^2U+\Delta^2U\otimes U
   +2\Delta U\otimes\Delta U,\\
 {\mathcal C}_2
 &=U\otimes\mathbb PF+\mathbb PF\otimes U.
 \end{aligned}                                         \tag{1.5}
\]

Conditional on the global tubular projected-remainder lemma (3.3), the
one-viscosity block retains an aspect-uniform \(L^1\) estimate after using
the flat-cylinder gradient cancellation:

\[
 \|{\mathcal C}_2\|_1\le C v^3.                         \tag{1.6}
\]

The two-viscosity block has a genuine loss:

\[
 \|{\mathcal A}_2\|_1
 \le C v^2R^{-1}\varepsilon^{-2}.                       \tag{1.7}
\]

More strongly, its zeroth moment is exactly

\[
 \boxed{\displaystyle
 \int_{\mathbb R^3}{\mathcal A}_2\,dx
 =4\int_{\mathbb R^3}\Delta U\otimes\Delta U\,dx.}       \tag{1.8}
\]

Thus the \(\varepsilon^{-2}\) factor in (1.7) is not an artifact of
absolute values.  For a generic thin Gavrilov cutoff, the trace-free part
of (1.8) is nonzero and forces an
\(\varepsilon^{-2}d^{-4}\) third-velocity tail.

The exterior field admits the same solenoidal cutoff construction as the
second jet, now with

\[
 \begin{aligned}
 \|\nabla^m(\tau^3Z_3)\|_2
 \le C_m\big(
 &\nu T^3R^{1/2-m}\\
 &+\nu^2T^3v^{-1}R^{-1/2-m}\varepsilon^{-2}\big).
                                                               \tag{1.9}
 \end{aligned}
\]

At distance \(d\asymp R\), relative to the main velocity \(v\), this is

\[
 T^3\mu+T^3\mu^2\varepsilon^{-2},\qquad
 \mu={\nu\over vR}.                                     \tag{1.10}
\]

If

\[
 \Theta={\nu\tau\over(\varepsilon R)^2}
 ={T\mu\over\varepsilon^2}                              \tag{1.11}
\]

is the minor-scale viscous parameter, then (1.10) becomes

\[
 T^2\varepsilon^2\Theta
 +T\varepsilon^2\Theta^2.                               \tag{1.12}
\]

Therefore the third exterior wake remains perturbative in the intended
\(\Theta\ll1\) window, but an aspect-uniform coefficient theorem with
only the macro Reynolds number is false beginning at this jet.

## 2. Exact third-jet recursion

Differentiate

\[
 \partial_tu+(u\cdot\nabla)u+\nabla p-\nu\Delta u=0,
 \qquad \operatorname{div}u=0
\]

twice at \(t=0\).  The binomial coefficient in the middle interaction
gives

\[
 \begin{aligned}
 u_3
 &+(u_2\cdot\nabla)U
 +2(u_1\cdot\nabla)u_1
 +(U\cdot\nabla)u_2\\
 &+\nabla p_2-\nu\Delta u_2=0.                           \tag{2.1}
 \end{aligned}
\]

Projection gives (1.1), while divergence gives (1.3).  Substitution of
(1.2) gives (1.4)--(1.5).

Although \(u_2\) is global, \({\mathcal T}_2\) is compactly supported:
every occurrence of \(u_2\) in (1.3) is multiplied by \(U\), and \(u_1\)
is compactly supported.  Moreover, outside the endpoint support,

\[
 u_2=-\nabla p_1,\qquad \Delta u_2=0.
\]

All local terms in (2.1) consequently vanish there, so

\[
 u_3=-\nabla p_2                                        \tag{2.2}
\]

outside the bubble.  The third exterior velocity is again exactly a
harmonic pressure wake.

## 3. Why the one-viscosity block remains uniform

The crude bound

\[
 \|\mathbb PF\|_2\le\|F\|_2
 \le Cv^2R^{-3/2}\varepsilon^{-2}
\]

would give
\(\|{\mathcal C}_2\|_1\le Cv^3\varepsilon^{-1}\).
That loses one aspect power.  It misses a cancellation special to the
thin Gavrilov profile.

In flat transverse polar coordinates, the leading rescaled field is

\[
 V_\perp=-r\chi(s)e_\theta,\qquad
 V_\parallel={r\chi(s)\over\sqrt2},\qquad
 s={r^2\over2}.
\]

Put

\[
 A(s)=4\chi'(s)+2s\chi''(s).
\]

The leading transverse vector Laplacian is

\[
 \Delta_yV_\perp=-rA(s)e_\theta.
\]

All tangential scalar profiles are radial, so their advection vanishes at
this order.  Direct polar differentiation gives

\[
 (V\cdot\nabla_y)\Delta_yV
 +(\Delta_yV\cdot\nabla_y)V
 =-2r\chi(s)A(s)e_r
 =\nabla_y\Phi_0.                                      \tag{3.1}
\]

The leading \(\varepsilon^{-3}\) vector is therefore a pure gradient.
It is killed by \(\mathbb P\).

More explicitly, the analytic tubular chart and its metric have uniform
expansions on the fixed support of \(\chi\):

\[
 \begin{aligned}
 U&=v\left(V(y)+O(\varepsilon)\right),\\
 \nabla&=R^{-1}\left(\varepsilon^{-1}\nabla_y+O(1)\right),\\
 \Delta&=R^{-2}\left(\varepsilon^{-2}\Delta_y
                  +O(\varepsilon^{-1}\nabla_y)+O(1)\right).
                                                               \tag{3.2}
 \end{aligned}
\]

The local expansion formally gives

\[
 F={v^2\over R^3}\left[
 \varepsilon^{-3}\nabla_y\Phi_0
 +\varepsilon^{-2}{\mathcal R}_\varepsilon\right],
\]

where the desired \({\mathcal R}_\varepsilon\) is uniformly bounded in
rescaled tubular variables.  To convert this local cancellation into a
theorem, one still must prove that the tubular expansion and its support
seams supply a global scalar \(\Phi_\varepsilon\), constant off the thin
shell, for which

\[
 \left\|F-\nabla\Phi_\varepsilon\right\|_2
 \le Cv^2R^{-3/2}\varepsilon^{-1}.                      \tag{3.3}
\]

The derivative counting is consistent with (3.3): every term in the
remainder has replaced at least one transverse
derivative \((\varepsilon R)^{-1}\) by a curvature or tangential
derivative \(R^{-1}\).  This does not by itself establish the global
Leray/seam estimate.  **Assuming (3.3)** and using
\(\mathbb P\nabla\Phi_\varepsilon=0\),

\[
 \|\mathbb PF\|_2
 \le Cv^2R^{-3/2}\varepsilon^{-1}.                      \tag{3.4}
\]

Together with

\[
 \|U\|_2\le CvR^{3/2}\varepsilon,
\]

Cauchy--Schwarz proves (1.6) conditionally.  Thus the volume factor would
cancel the remaining projected thin loss exactly once (3.3) is proved.

This is the only place where the leading flat-cylinder structure is
needed.  Without (3.1), the best immediate estimate loses
\(\varepsilon^{-1}\).

## 4. The two-viscosity moment cannot use that cancellation

The derivative bounds and tube volume give

\[
 \|\Delta^2U\|_\infty
 \le CvR^{-4}\varepsilon^{-4},\qquad
 |\operatorname{supp}U|\le CR^3\varepsilon^2.
\]

Every term of \({\mathcal A}_2\) therefore has the \(L^1\) scale in
(1.7).  Integration by parts gives the stronger identity

\[
 \begin{aligned}
 \int U_i\Delta^2U_j\,dx
 &=\int\Delta U_i\Delta U_j\,dx,\\
 \int\Delta^2U_iU_j\,dx
 &=\int\Delta U_i\Delta U_j\,dx.
 \end{aligned}
\]

Adding the explicit heat--heat term proves (1.8).  This is a positive
semidefinite covariance, so no cancellation between the three
two-viscosity tensors is possible.

There can still be harmonic cancellation if the covariance happens to
be isotropic.  The thin-profile limit makes this condition explicit.
With

\[
 k(r)=-rA(s),\qquad
 l(r)={1\over\sqrt2}\left({\chi(s)\over r}+rA(s)\right),
\]

the axial and either Cartesian transverse components of the limiting
covariance are

\[
 \begin{aligned}
 H_{\parallel,0}
 &=2\pi^2\int_0^\infty k(r)^2r\,dr,\\
 H_{\perp,0}
 &=\pi^2\int_0^\infty
       \left(k(r)^2+2l(r)^2\right)r\,dr.                 \tag{4.1}
 \end{aligned}
\]

Their difference is

\[
 \boxed{\displaystyle
 H_{\parallel,0}-H_{\perp,0}
 =\pi^2\left(
 4\int_0^\infty s\chi'(s)^2\,ds
 -{1\over2}\int_0^\infty{\chi(s)^2\over s}\,ds
 \right).}                                              \tag{4.2}
\]

Correspondingly,

\[
 \int\Delta U\otimes\Delta U\,dx
 ={v^2\over R\varepsilon^2}
 \left(H_0+O(\varepsilon)\right).                       \tag{4.3}
\]

The scalar in (4.2) is generically nonzero.  A simple sufficient
condition for positivity is obtained by setting
\(f(t)=\chi(e^t)\).  If the logarithmic support length is \(L\), the
Dirichlet Poincare inequality gives

\[
 H_{\parallel,0}-H_{\perp,0}
 \ge\pi^2\left({4\pi^2\over L^2}-{1\over2}\right)
 \int_{\mathbb R}f(t)^2\,dt.                            \tag{4.4}
\]

Thus every nonzero cutoff with
\(L<2\sqrt2\pi\) gives a strictly nonisotropic leading tensor.  An
exceptional cutoff satisfying equality in (4.2) removes this particular
trace-free coefficient but does not improve the \(L^1\) scale (1.7).

## 5. Exterior third-wake estimate

Let

\[
 M_2=\int{\mathcal T}_2\,dx.
\]

Equations (1.6)--(1.8) give

\[
 \|{\mathcal T}_2\|_1
 \le C\left(
 \nu v^3+\nu^2v^2R^{-1}\varepsilon^{-2}\right),          \tag{5.1}
\]

and

\[
 M_2
 =4\nu^2\int\Delta U\otimes\Delta U\,dx
 -\nu\int{\mathcal C}_2\,dx.                            \tag{5.2}
\]

For the generic cutoffs in Section 4,

\[
 M_2^{\rm TF}
 ={4\nu^2v^2\over R\varepsilon^2}H_0^{\rm TF}
 +O\left({\nu^2v^2\over R\varepsilon}+\nu v^3\right).   \tag{5.3}
\]

At \(d=|x-c|\ge C_0R\), the Newton kernel and the moment bounds imply

\[
 |\nabla^mu_3(x)|
 \le C_m\left(
 \nu v^3+\nu^2v^2R^{-1}\varepsilon^{-2}\right)
 d^{-4-m}.                                             \tag{5.4}
\]

The leading term is

\[
 p_2(x)=\partial_i\partial_j\Gamma(x-c)(M_2)_{ij}
 +O\left(R\|{\mathcal T}_2\|_1d^{-4}\right).            \tag{5.5}
\]

As for the second jet, cut off \(-\nabla p_2\) on a fixed annulus at
major-radius scale and solve the scalar divergence error with a
compact-support Bogovskii operator.  The annulus is fixed after
\(x\mapsto x/R\), so its constants are independent of
\(\varepsilon\).  The resulting divergence-free \(Z_3\) equals \(u_3\)
outside the enlarged ball and obeys

\[
 \|\nabla^mZ_3\|_2
 \le C_m\left(
 \nu v^3R^{-5/2-m}
 +\nu^2v^2R^{-7/2-m}\varepsilon^{-2}\right).             \tag{5.6}
\]

Multiplication by \(\tau^3=T^3R^3/v^3\) proves (1.9).

Three orthogonal rotations of the isolated self tensor in (1.8)
isotropize its zeroth moment because it is axisymmetric.  For a
multi-bubble endpoint, however, \(u_2\) from one bubble is already
nonzero on every other core.  These cross-wake terms also enter
\({\mathcal T}_2\).  Exact aggregate cancellation therefore requires a
fully symmetric cluster or a direct moment check; rotating only the
three self tensors is not by itself sufficient.

## 6. Packed-stage consequence

For \(K^3\) separated microtori of major radius
\(\delta=\ell/K\), actual speed \(v\), common aspect
\(\varepsilon\), and macro duration

\[
 \tau=T{\ell\over v},
\]

the self-source bounds summed in absolute value give, at macro distance
\(d\asymp\ell\),

\[
 {|\tau^3u_{3,\rm stage}^{\rm ext}|\over v}
 \le C\left[
 {\nu T^3K^3\over v\ell}
 +{\nu^2T^3K^4\over v^2\ell^2\varepsilon^2}
 \right].                                               \tag{6.1}
\]

The cross-wake contribution to \(U_m\otimes u_{2,n}\) is no worse:
at one core the three-dimensional neighbor sum is
\(O(\nu v^2\delta^{-3})\); multiplication by the core volume
\(O(\delta^3\varepsilon^2)\) and by \(v\) gives
\(O(\nu v^3\varepsilon^2)\) per core.  It is bounded by the first
self-source term in (6.1).

For

\[
 v_j=\ell_j^{-\gamma}K_j^\gamma,
\qquad 1<\gamma<3/2,
\]

the two ratios in (6.1) are

\[
 \nu T_j^3\ell_j^{\gamma-1}K_j^{3-\gamma},
\]

\[
 \nu^2T_j^3\ell_j^{2\gamma-2}
 K_j^{4-2\gamma}\varepsilon_j^{-2}.                    \tag{6.2}
\]

Both are summable for geometric \(\ell_j\) when
\(K_j\), \(T_j\), and \(\varepsilon_j^{-1}\) grow polynomially.  If
\(\varepsilon_j=\ell_j^\beta\), the two-viscosity term instead requires

\[
 \beta<\gamma-1.                                       \tag{6.3}
\]

This is the restriction from the third exterior coefficient alone.  The
minor-scale high-Reynolds condition is stronger.  From (1.11), for
polynomial \(K_j,T_j\),

\[
 \Theta_j
 =\nu T_j\ell_j^{\gamma-1-2\beta}K_j^{2-\gamma}\to0
 \quad\Longrightarrow\quad
 \boxed{\beta<{\gamma-1\over2}}.                     \tag{6.4}
\]

Thus (6.4), not (6.3), is the compatible design gate when the transition
also requires perturbative viscosity on the minor radius.

## 7. Claim boundary

The second-jet aspect-uniform exterior split does extend one step, but
with a new coefficient:

* conditionally on the global tubular projection lemma (3.3), the projected
  one-viscosity source is uniform because its flat \(\varepsilon^{-3}\)
  part is a gradient;
* the two-viscosity source has exact moment
  \(4\nu^2\int\Delta U\otimes\Delta U\) and loses
  \(\varepsilon^{-2}\);
* this loss is generally visible in the exterior harmonic quadrupole;
* the correct perturbative quantity is the minor-scale parameter
  \(\Theta\), not the macro Reynolds number alone; and
* polynomial thinning remains compatible with the current
  geometric/polynomial cascade, while geometric thinning
  \(\varepsilon_j=\ell_j^\beta\) must satisfy the stronger existing
  minor-viscosity gate \(\beta<(\gamma-1)/2\).

At the fourth and later jets, products involving the compact local
remainder and accumulated global wakes proliferate.  Nothing here proves
an order-\(M_j\) Gevrey bound or closes those cross interactions.

## Primary source

* A. V. Gavrilov,
  [A steady Euler flow with compact support](https://arxiv.org/abs/1810.08020),
  for the analytic toroidal seed and its pressure localization.
