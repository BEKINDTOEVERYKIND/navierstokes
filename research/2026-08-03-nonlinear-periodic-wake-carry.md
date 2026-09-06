# Nonlinear periodic wake carry: the Reynolds no-go and the scale-shift gate

**Date:** 2026-08-03

**Status:** exact nondimensionalization, low-mode contraction no-go, and
stage-dilation ledger; the one-stage wake propagator and its weighted-space
closure remain open.

**Scope:** the nonlinear terminal-state problem left by C86--C88.  This note
does not construct the three-phase transition, prove a wake fixed point, or
prove singularity formation.  Every new conclusion below is self-audited.

## 1. Outcome

There is no Reynolds-number contraction for the periodic zero-charge wake.
After one turnover, envelope-scale transport, stretching, wake self-advection,
and an active high--high zero-charge source are all order one.  Viscosity is
only

\[
                  \mu_j={\nu\over a_j\ell_j}=\operatorname{Re}_j^{-1},
\]

and the retained carrier hierarchy deliberately satisfies
\(\mu_j(M_jK_j)^2\to0\).  Heat therefore approaches the identity on every
retained mode.  In fact, even with no background and no forcing, the
derivative of the terminal map is \(e^{\mu S\Delta}\), whose norm on the
lowest mean-zero torus mode is \(e^{-\mu S}\to1\).  High Reynolds number
cannot yield a viscosity-uniform strict contraction.

There is, however, a different possible small parameter.  Passing a terminal
wake from scale \(\ell_j\) to \(\ell_{j+1}=\ell_j/r\) applies the exact
renormalization

\[
 (\mathcal D_jz)(y)
 =q_j^{-\gamma}z(y/r),\qquad
 q_j=r{K_{j+1}\over K_j},                                    \tag{1.1}
\]

for the polynomial-carrier amplitudes
\(a_j=\ell_j^{-\gamma}K_j^\gamma\).  On an annular velocity space with
weight \(|y|^\alpha\),

\[
 \|\mathcal D_j\|_{X_\alpha^m\to X_\alpha^m}
 =r^\alpha q_j^{-\gamma}
 \le r^{\alpha-\gamma}.                                      \tag{1.2}
\]

Thus:

* the physical critical wake weight \(\alpha=\gamma\) is asymptotically
  neutral, not contracting;
* every slightly weaker weight \(0<\alpha<\gamma\) gains the genuine scale
  contraction \(r^{-(\gamma-\alpha)}\); and
* if the **one-stage nonlinear wake endpoint map** has a uniform Lipschitz
  constant \(L\), the carried map is a contraction whenever
  \(Lr^{\alpha-\gamma}<1\).

This isolates the correct analytic target.  One must prove a wake-specific
uniform endpoint bound (or one-way spatial export).  The full carrier
semigroup estimate over the actual \(G_j\asymp j^2\) gain window is not
enough: its available upper bound is exponential in \(G_j\), which no fixed
geometric scale jump can beat.

## 2. Exact zero-charge equation

Use one parent-envelope length and velocity as units:

\[
 x=x_j+\ell_jy,\qquad
 t=t_j+{\ell_j\over a_j}\tau,\qquad
 u=a_jV.                                                       \tag{2.1}
\]

Then projected Navier--Stokes is

\[
 \partial_\tau V+\mathbb P\operatorname{div}(V\otimes V)
 =\mu_j\Delta V,
 \qquad \mu_j={\nu\over a_j\ell_j}.                          \tag{2.2}
\]

Decompose the complete phase-resolved field into

\[
                         V=B+Z+W.                              \tag{2.3}
\]

Here \(B\) is the prescribed zero-charge parent/child path, \(Z\) is the
retained zero-charge wake, and \(W\) contains all nonzero material charges.
Let \(Q\) denote the complete zero-charge part of \(W\otimes W\), including
all opposite-charge pairs and lower-order correctors.  The zero-charge
projection of (2.2) is

\[
 \partial_\tau(B+Z)
 +\mathbb P\operatorname{div}
 \big((B+Z)\otimes(B+Z)+Q\big)
 =\mu_j\Delta(B+Z).                                           \tag{2.4}
\]

Define the designed-path defect

\[
 F_j:=\mu_j\Delta B-\partial_\tau B
      -\mathbb P\operatorname{div}(B\otimes B+Q).             \tag{2.5}
\]

The exact nonlinear wake equation is therefore

\[
 \boxed{
 \partial_\tau Z-\mu_j\Delta Z
 +\mathbb P\operatorname{div}
  (B\otimes Z+Z\otimes B+Z\otimes Z)=F_j.}                   \tag{2.6}
\]

All fields in (2.6) have zero spatial mean.  The mean is conserved, so C86
supplies the periodic symmetric routing at every time without adding a
compatibility condition.

The Stokes mild form is

\[
\begin{aligned}
 Z(\tau)&=e^{\mu_j\tau\Delta}Z_-\\
 &\quad+\int_0^\tau e^{\mu_j(\tau-s)\Delta}
 \left[F_j(s)-\mathbb P\operatorname{div}
 \big(B\otimes Z+Z\otimes B+Z\otimes Z\big)(s)\right]ds.
                                                                    \tag{2.7}
\end{aligned}
\]

Equivalently, if \(\mathcal U_{B,\mu}(\tau,s)\) is the evolution family for
the two terms linear in \(Z\),

\[
 Z(\tau)=\mathcal U_{B,\mu}(\tau,0)Z_-
 +\int_0^\tau\mathcal U_{B,\mu}(\tau,s)
 \left[F_j-\mathbb P\operatorname{div}(Z\otimes Z)\right](s)ds.
                                                                    \tag{2.8}
\]

Equation (2.8), not the uncoupled heat formula in C86, is the actual
carry-forward map.

## 3. Complete one-turnover scaling

Let the dimensionless wake amplitude be \(b\), its frequency be \(N\), and
take \(B\), \(\nabla B\), and the active covariance \(Q\) to have unit size.
Before cancellations, the sizes in (2.6), integrated for one turnover
\(\Delta\tau\asymp1\), are

\[
\begin{array}{c|c}
\text{term}&\text{one-turnover size}\\ \hline
\partial_\tau Z&b\\
B\cdot\nabla Z&Nb\\
Z\cdot\nabla B&b\\
Z\cdot\nabla Z&Nb^2\\
\mu_j\Delta Z&\mu_jN^2b\\
\mathbb P\operatorname{div}Q&1\\
\text{exact-collar defect from C87}&\mu_j.
\end{array}                                                     \tag{3.1}
\]

The large \(Nb\) transport is not an amplitude growth rate: for
divergence-free \(B\), its principal part is skew in \(L^2\), and material
phases absorb it.  The energy-relevant linear coefficient is
\(\|\operatorname{sym}\nabla B\|_\infty=O(1)\).  Similarly,
\(Z\cdot\nabla Z\) cancels in the \(L^2\) energy identity, although it is
order one in higher norms and in the endpoint derivative.  Consequently an
energy estimate has the schematic form

\[
 {d\over d\tau}\|Z\|_{H^s}
 \lesssim
 (\|\nabla B\|_\infty+\|\nabla Z\|_\infty)\|Z\|_{H^s}
 +\|F_j\|_{H^s},                                               \tag{3.2}
\]

with no favourable factor of \(\mu_j\).

There are three relevant frequency regimes:

\[
\begin{array}{c|c}
N&\text{heat action per turnover}\\ \hline
1&\mu_j\\
K_j&\theta_j:=\mu_jK_j^2\\
M_jK_j&\mu_j(M_jK_j)^2.
\end{array}                                                     \tag{3.3}
\]

The current construction requires the last quantity to tend to zero.
Hence heat is asymptotically negligible not only on the low wake but on
every retained WKB harmonic.

The distinction between the two wake sources is exact and important.

* The C87 collar has \(F_j=O(\mu_j)\), so its generated wake is
  \(Z=O(\mu_j)\) for a fixed turnover.  Then the self-interaction is
  \(O(\mu_j^2)\), although its linear transport by \(B\) remains part of an
  order-one propagator.
* An active high--high transfer has \(Q=O(1)\) near capture.  Its designed
  part may be cancelled by the chosen evolution of \(B\), but every
  uncancelled component of \(F_j\) has order-one Duhamel scale.  Reynolds
  number does not make that component a perturbation.

If the growing packet has amplitude
\(A(\tau)=e^{\lambda(\tau-G)}\) on \(0\le\tau\le G\), its quadratic source
does **not** accumulate a factor \(G\):

\[
 \int_0^GA(\tau)^2d\tau
 ={1-e^{-2\lambda G}\over2\lambda}\le{1\over2\lambda}.       \tag{3.4}
\]

Thus a terminally concentrated transfer can create an order-one rather than
an order-\(G\) new wake.  This does not control an incoming wake, which is
present throughout the whole gain interval.

## 4. A decisive Reynolds-number no-go

Consider the easiest special case of (2.6): \(B=0\), \(F_j=0\), and
linearize at \(Z=0\) on \((\mathbb R/2\pi\mathbb Z)^3\).  For a fixed
turnover length \(S>0\), the derivative of the terminal map is

\[
                  D\Phi_{\mu,S}(0)=e^{\mu S\Delta}.            \tag{4.1}
\]

On mean-zero \(H^s\),

\[
             \|e^{\mu S\Delta}\|_{H^s_0\to H^s_0}
             =e^{-\mu S}.                                     \tag{4.2}
\]

It follows that

\[
             \lim_{\mu\downarrow0}
             \|D\Phi_{\mu,S}(0)\|=1.                         \tag{4.3}
\]

Therefore no argument using only \(\operatorname{Re}\to\infty\) can give a
terminal-map Lipschitz constant bounded by one fixed \(c<1\).  On the
expanding normalized torus associated with a fixed physical periodic box,
the lowest normalized frequency tends to zero and the conclusion is even
stronger.

The same obstruction holds at the active carrier: over \(S\) turnovers its
linear heat factor is \(e^{-S\theta_j}\).  For the actual gain window
\(G_j\asymp j^2\), the geometric decay of \(\mu_j\) implies

\[
 \mu_jG_j\to0,
 \qquad
 \mu_j(M_jK_j)^2G_j\to0                                      \tag{4.4}
\]

after including any fixed additional polynomial demanded by the existing
heat gate.  Viscosity still approaches the identity.

Conversely, a generic energy bound for an incoming wake over that interval
is

\[
 \|Z(G_j)\|\lesssim
 e^{C G_j}\left(\|Z_-\|+\int_0^{G_j}\|F_j(s)\|ds\right).      \tag{4.5}
\]

The strain cap in C82 replaces \(C\) by the sharp edge rate; it does not
make the exponential disappear.  Since \(G_j\asymp j^2\), no fixed factor
coming from \(\ell_{j+1}/\ell_j=1/r\) can beat (4.5).  A successful proof
therefore needs a wake-specific power-bounded propagator, spatial export
that makes the old wake invisible to the active core, or a selected
invariant/stable graph.  The unrestricted carrier semigroup estimate cannot
close the carry chain.

This is a no-go for an **unstructured high-Reynolds contraction**, not for
the cascade itself.

## 5. Exact terminal-state renormalization

Use the current polynomial-carrier schedule

\[
 \ell_j=r^{-j},\qquad K_j=j^A,\qquad
 a_j=\ell_j^{-\gamma}K_j^\gamma,
 \qquad 1<\gamma<\frac32.                                    \tag{5.1}
\]

For clarity take nested stages with the same center.  If
\(z(t_j^+,x)=a_jZ_j^+((x-x_j)/\ell_j)\), then in the next stage variables

\[
 Z_{j+1}^-(y)
 ={a_j\over a_{j+1}}Z_j^+(y/r)
 =q_j^{-\gamma}Z_j^+(y/r),
 \quad q_j=r{K_{j+1}\over K_j}.                              \tag{5.2}
\]

A bounded center displacement only composes (5.2) with a bounded
translation and does not change its exponents.

For an annular wake define the homogeneous seminorm

\[
 \|f\|_{X_\alpha^m}
 :=\max_{0\le k\le m}\sup_{y\ne0}
       |y|^{\alpha+k}|\nabla^kf(y)|.                          \tag{5.3}
\]

For components supported away from the new core, differentiation of (5.2)
gives the exact identity

\[
 \|\mathcal D_jf\|_{X_\alpha^m}
 =r^\alpha q_j^{-\gamma}\|f\|_{X_\alpha^m}.                 \tag{5.4}
\]

The derivative order cancels: \(r^{-k}\) from differentiation is paired
with \(r^k\) from the weight.  Because \(K_{j+1}\ge K_j\),

\[
 r^\alpha q_j^{-\gamma}
 =r^{\alpha-\gamma}
   \left({K_j\over K_{j+1}}\right)^\gamma
 \le r^{\alpha-\gamma}.                                     \tag{5.5}
\]

At the critical physical tail \(\alpha=\gamma\), the factor is

\[
 \left({j\over j+1}\right)^{A\gamma}\longrightarrow1,       \tag{5.6}
\]

and from stage \(i\) to stage \(j\) it telescopes only to

\[
                         \left({i\over j}\right)^{A\gamma}.   \tag{5.7}
\]

Thus the geometric part of the critical \(|y|^{-\gamma}\) wake is exactly
neutral.  This agrees with the earlier observation that the renormalized
outer wake is not in unweighted \(L^2\): for a shell expanded by \(R\),

\[
 \|R^{-\gamma}f(\cdot/R)\|_2
 =R^{3/2-\gamma}\|f\|_2,                                    \tag{5.8}
\]

which grows because \(\gamma<3/2\).

On every subcritical wake weight \(0<\alpha<\gamma\), however, (5.5)
gives a uniform scale contraction.  This is not viscous damping; it is the
amplitude increase of the child normalization outrunning the weaker annular
weight.

## 6. The exact contraction criterion for the nonlinear carry map

Let \(\Phi_j\) denote the complete endpoint map of (2.6), including the
prescribed phase-resolved source, on a ball in \(X_\alpha^m\).  Suppose --
this is the missing analytic input -- that

\[
 \|\Phi_j(z)-\Phi_j(\widetilde z)\|_{X_\alpha^m}
 \le L_j\|z-\widetilde z\|_{X_\alpha^m}.                     \tag{6.1}
\]

The next-stage carry map is \(\Psi_j=\mathcal D_j\Phi_j\).  Equations
(5.4)--(5.5) give the exact bound

\[
 \operatorname{Lip}(\Psi_j)
 \le r^{\alpha-\gamma}
    \left({K_j\over K_{j+1}}\right)^\gamma L_j.              \tag{6.2}
\]

Consequently, if for some \(0<\alpha<\gamma\)

\[
                    \boxed{
 \sup_j L_j\,r^{\alpha-\gamma}<1,}                           \tag{6.3}
\]

then terminal-wake differences contract geometrically from stage to stage.
If also \(\sup_j\|\Psi_j(0)\|<\infty\), the affine estimate

\[
 \|z_{j+1}\|\le\rho\|z_j\|+C,\qquad \rho<1,                \tag{6.4}
\]

gives a bounded carried wake and exponential forgetting of its initial
state.  The usual backward-limit argument then gives a unique bounded
pullback sequence whenever the stage maps are specified for all earlier
indices.

At \(\alpha=\gamma\), (6.3) cannot follow from scale separation: the
dilation is asymptotically an isometry.  One must retain the critical shell
coefficient as part of a genuine invariant graph.  In an explicit shell
model, the propagation is the unilateral shift

\[
 (c_0,c_1,c_2,\ldots)
 \longmapsto(h_{j+1},\mathcal D_jc_0,
                  \mathcal D_jc_1,\ldots),                   \tag{6.5}
\]

where \(h_{j+1}\) is the newly generated wake.  It is neutral in the
critical weighted \(\ell^\infty\) norm and contracting in every
\(\alpha<\gamma\) norm.

For the long spectral gain window, the presently available full-space
estimate permits \(L_j\sim e^{Cj^2}\), so (6.3) fails for every fixed
\(r\).  The new theorem target is therefore precise:

> Prove that the endpoint derivative restricted to incoming annular
> zero-charge wakes is uniformly bounded in a slightly subcritical weighted
> Gevrey-2 space, or prove a one-way core/wake estimate with the same effect.
> Then choose \(r\) so that (6.3) holds and solve the nonautonomous carried
> wake chain.  At the critical weight, construct the corresponding invariant
> graph rather than claiming contraction.

This reduction separates three facts that had been conflated: the collar
wake is Reynolds-small, the active work channel has order-one scaling before
designed cancellations, and the only scalar contraction exposed by the
present ledger is stage renormalization plus a wake-specific endpoint bound.

## 7. Claim boundary

What is exact here:

1. the dimensionless nonlinear zero-charge equation and mild forms
   (2.6)--(2.8);
2. the one-turnover scaling and the distinction between collar and active
   wakes;
3. the heat-semigroup no-go (4.1)--(4.3);
4. the terminal dilation and weighted norm factors (5.2)--(5.7); and
5. the conditional Lipschitz criterion (6.2)--(6.3).

What remains open:

1. closure of the periodic Leray/pressure tails in the proposed annular
   Gevrey space, including wrap-around on the fixed physical torus;
2. a uniform wake-restricted endpoint constant \(L_j\) through the
   \(G_j\asymp j^2\) carrier interval;
3. compatibility of that map with the rank-five affine capture chart and
   the prescribed C87--C88 terminal collar jet; and
4. the nonlinear invariant graph at the critical \(|y|^{-\gamma}\) tail.

Accordingly, this note is a rigorous reduction and a no-go for one tempting
argument, not a solution of Navier--Stokes regularity.
