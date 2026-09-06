# An edge-matched carrier in the locked-pitch Gavrilov family

**Date:** 2026-08-02

**Status:** exact straight-column BAS/profile theorem and dependency-free
rational-interval certificate; compactification argument self-derived;
Albritton--Ożański proof transfer self-audited but not independently
cross-audited

**Scope:** a pressure-modulated Gavrilov straight limit and its linearized
Euler spectrum.  This note does not bend the mode into a finite-curvature
torus, prove a frequency-uniform PDE semigroup bound, construct the nonlinear
transition, or prove a Navier--Stokes singularity.

## 1. Outcome

The base-profile mismatch in C67 has a repair inside the actual standard
Gavrilov straight family.  It is not necessary to realize the Batchelor
profile.

There is an explicit smooth hollow locked-pitch column

\[
 V(r)=r\Omega(r),\qquad W(r)={V(r)\over\sqrt2},                 \tag{1.1}
\]

with a compact annular version such that:

1. one resonant ring is the unique maximizer of the **full** straight-column
   BAS exponent, allowing every radius, helical ratio, and radial covector;
2. the same ring is a strict local fixed-helical-ratio maximum of the AO
   coefficient `b` and the unique negative minimum of the helical phase
   `Lambda`;
3. its squared edge is the exact algebraic number
   \[
   \lambda_*^2={5+\sqrt{22}\over3}\,\Omega_0^2;
                                                               \tag{1.2}
   \]
4. the uniform physical BAS propagator argument of C66 applies with no
   polynomial prefactor; and
5. the quotient `q` in the literal statement of Albritton--Ożański
   Assumption A necessarily has a pole for every nonzero compact
   locked-pitch bump, but that quotient cancels exactly from the physical
   Rayleigh coefficients.

A source-level audit of the proof of Albritton--Ożański Theorem 1.1 finds
that the gluing argument uses only the regularized Rayleigh coefficients,
the local trapping conditions, their behavior at the axis, and isolation of
the `Lambda` level.  It never subsequently uses positivity or boundedness of
the artificial quotient `q`.  This gives a direct q-free proof-transfer
lemma and, subject to independent audit, exact high-frequency ring modes for
the compact locked-pitch column.  The weaker local-Weber-quasimode route
survives even if that proof transfer is not accepted.

The principal sources are:

- D. Albritton and W. Ożański, *Linear and nonlinear instability of vortex
  columns*, arXiv:2310.20674v3:
  <https://arxiv.org/abs/2310.20674>;
- A. V. Gavrilov, *A steady Euler flow with compact support*,
  arXiv:1810.08020: <https://arxiv.org/abs/1810.08020>; and
- P. Constantin and V. La, *Remarks on a paper by Gavrilov: Grad--Shafranov
  equations, steady solutions of the three-dimensional incompressible Euler
  equations with compactly supported velocities, and applications*,
  arXiv:1903.11699: <https://arxiv.org/abs/1903.11699>.

## 2. Exact locked-pitch BAS reduction

Let

\[
 V=r\Omega,\qquad W=\sigma V,\qquad \sigma^2={1\over2},
 \qquad h={r\Omega'\over\Omega}.                    \tag{2.1}
\]

In the notation of the full BAS reduction C63,

\[
 S=r\Omega'=\Omega h,\qquad
 T=W'=\sigma\Omega(1+h).                            \tag{2.2}
\]

For the helical covector

\[
 m=-{n\over r},\qquad k=\beta n,
\]

the radial covector drift is

\[
 D=mS+kT=n(\beta W'-\Omega').                       \tag{2.3}
\]

Every nonresonant trajectory `D != 0` has zero Lyapunov exponent by the
exact symmetrization in C63.  At resonance, putting `y=beta r`,

\[
                       y={h\over\sigma(1+h)}.        \tag{2.4}
\]

For radial covector `ell` and
`K_h^2=m^2+k^2`, direct substitution into C63 gives

\[
 \lambda^2
 =\Omega^2 f(h){K_h^2\over\ell^2+K_h^2},
 \qquad
 f(h)=-2h\,{3h^2+6h+1\over3h^2+2h+1}.              \tag{2.5}
\]

The denominator in (2.5) is everywhere positive.  Hence

\[
 f(h)>0
 \quad\Longleftrightarrow\quad
 h<-1-\sqrt{2\over3}
 \quad\hbox{or}\quad
 -1+\sqrt{2\over3}<h<0.                             \tag{2.6}
\]

There is an invariant form which also covers the endpoint hidden by the
helical ratio.  For normalized horizontal covector
`(p_h,q_h)=(m,k)/K_h`, exact resonance is
`p_h S+q_h T=0`.  If `S^2+T^2>0`, substitution of the two resonant unit
directions gives

\[
 \mu_{\rm res}
 =-2\Omega S\,{S^2+T^2+2\Omega S\over S^2+T^2}.    \tag{2.7a}
\]

For the locked-pitch profile, (2.7a) reduces to \(\Omega^2f(h)\), including
the `m=0`, `h=-1` endpoint where it equals \(-2\Omega^2\).  If `S=T=0`
but `Omega` is nonzero, direct substitution gives
\(\mu=-4\Omega^2q_h^2\le0\); if `Omega=S=T=0`, it gives zero.  None of
these endpoint cases adds positive growth.

The unrestricted positive BAS edge is therefore

\[
 \lambda_{\rm BAS}^2
 =\sup_{r>0}[\Omega(r)^2f(h(r))]_+,                 \tag{2.7}
\]

and positive equality can occur only at `ell=0`.

For later use, before imposing resonance the AO coefficient is

\[
 {b_\beta\over\Omega^2}
 =-2y\,{\sigma(1+h)+y(h+2)\over1+y^2}.             \tag{2.8}
\]

At a resonant ring, transverse stationarity in `y` is equivalent to

\[
                    3h^2+10h+1=0.                  \tag{2.9}
\]

At either root of (2.9),

\[
 f(h)=-h,\qquad {f'(h)\over f(h)}=-{3h+1\over4h}.  \tag{2.10}
\]

Writing a dot for `d/d(log r)`, stationarity in the radial direction then
requires

\[
                   \dot h={8h^2\over3h+1}.          \tag{2.11}
\]

The descending-flank root is

\[
 h_*={-5-\sqrt{22}\over3}=-3.2301385866\ldots,
 \qquad
 c_*={8h_*^2\over3h_*+1}=-9.6048756028\ldots .      \tag{2.12}
\]

## 3. An explicit isolated full edge

Put `t=log r`, fix `Omega_0>0`, and define

\[
 \Omega_\infty(r)
 =-\Omega_0\exp\left(h_*t+{c_*\over2}t^2\right),
 \qquad
 V_\infty=r\Omega_\infty,
 \qquad
 W_\infty={V_\infty\over\sqrt2}.                  \tag{3.1}
\]

Then

\[
 h(t)=h_*+c_*t.                                     \tag{3.2}
\]

The negative quadratic in (3.1) dominates every linear power at both ends.
Thus the profile extends smoothly and flatly through the axis and decays
faster than every power at infinity.

Choose

\[
 \beta_*={\sqrt2h_*\over1+h_*}
          =2.0483506384\ldots .                     \tag{3.3}
\]

At `r_*=1`, equations (2.9)--(2.12) give

\[
 b_*=-h_*\Omega_0^2,
 \qquad
 \lambda_*=\sqrt{-h_*}\,\Omega_0
            =1.7972586310\ldots\Omega_0.            \tag{3.4}
\]

### 3.1 Global BAS maximum

Use `h` as the global coordinate.  Equations (3.1)--(3.2) give

\[
 \Omega_\infty^2
 =\Omega_0^2\exp\left({h^2-h_*^2\over c_*}\right),
\]

so the positive resonant growth is

\[
 \mathcal B(h)=\Omega_0^2
 \exp\left({h^2-h_*^2\over c_*}\right)f(h).         \tag{3.5}
\]

On the steep positive lobe set `x=-h`,

\[
 A=3x^2-6x+1,\qquad B=3x^2-2x+1,
 \qquad f(-x)={2xA\over B}.
\]

After putting `(log f)''` over the positive denominator
`x^2 A^2 B^2`, its numerator is `-R(x)`, where

\[
\begin{aligned}
R(x)={}&81x^8-216x^7+612x^6-1152x^5+1110x^4\\
      &\quad-504x^3+132x^2-16x+1.                  \tag{3.6}
\end{aligned}
\]

The translated polynomial has strictly positive coefficients:

\[
\begin{aligned}
 R(1+z)={}&48+224z+672z^2+1632z^3+2640z^4\\
          &+2520z^5+1368z^6+432z^7+81z^8.          \tag{3.7}
\end{aligned}
\]

Thus `(log f)''<0` throughout the steep lobe, where `x>1`.  Since
`2/c_*<0`, \(\log\mathcal B\) is strictly concave there.  Equations
(2.10)--(2.12) make its derivative vanish at `h_*`, so this is its unique
steep-lobe maximum.

On the mild positive lobe, put again `x=-h`.  Then

\[
 0<x<1-\sqrt{2\over3}<{1\over5},
 \qquad 0<f(-x)<2x<{2\over5}.                      \tag{3.8}
\]

Moreover

\[
 \exp\left({h^2-h_*^2\over c_*}\right)
 \le \exp\left({4+\sqrt{22}\over8}\right)<4.      \tag{3.9}
\]

Consequently the entire mild lobe has squared growth below
\(8\Omega_0^2/5\), whereas (3.4) is greater than \(3\Omega_0^2\).
The stable intervals have no positive exponent.  This proves that `r=1`
is the unique full BAS edge.

### 3.2 Local AO geometry and the isolated phase level

For the fixed ratio (3.3),

\[
 \Lambda=\beta_*W-\Omega.
\]

At the selected ring,

\[
 \Lambda_*={\Omega_0\over1+h_*}
            =-0.4484026266\ldots\Omega_0,           \tag{3.10}
\]

and

\[
 \Lambda''(1)
 =\Omega(1)\left(h_*-{c_*\over1+h_*}\right)
 =7.5369900354\ldots\Omega_0>0.                    \tag{3.11}
\]

To see global isolation, write

\[
 z(t)={\beta_*r\over\sqrt2}=z_*e^t,
 \qquad z_*={h_*\over1+h_*}>1.
\]

Then

\[
 {d\Lambda\over dt}=\Omega[z(1+h)-h].              \tag{3.12}
\]

Critical points solve `z=h/(1+h)`.  Since `h` decreases strictly while
`z` increases strictly, there is one critical point on `h>0`, none on
`-1<h<0`, and exactly one on `h<-1`, namely `t=0`.  The first is a positive
maximum of `Lambda`; the designed point is its unique negative minimum.
In particular no other radius has the level `Lambda_*`.

The rational-interval checker additionally gives, with `t=log r` and
`Omega_0=1`,

\[
 -193<{d^2\over dt^2}b_{\beta_*}(1)<-192,           \tag{3.13}
\]

\[
 -207<{d^2\over dt^2}\mathcal B(1)<-205,           \tag{3.14}
\]

and the transverse second derivative of (2.8) is negative.  Thus the ring
has the complete strict local trapping geometry used by the Weber
construction.

## 4. Compactification and Gavrilov realization

Let `chi_L(t)` be a smooth nonnegative annular cutoff equal to one on a
large interval containing zero, rising monotonically on its inner
transition and falling monotonically on its outer transition.  Put

\[
                    \Omega_L=\chi_L\Omega_\infty.  \tag{4.1}
\]

All logarithmic-slope formulas in this paragraph are understood on the open
set where `chi_L>0`.  On the inner transition, `h` is large and positive and
`d(log chi_L)/dt >= 0`; hence the modified logarithmic slope stays positive
and creates no BAS growth.  On the outer transition,

\[
 f(h)\lesssim1+|h|
\]

and

\[
 (\chi_L\Omega_\infty)^2
 \left|h+{d\over dt}\log\chi_L\right|
 \lesssim \Omega_\infty^2
 \left(\chi_L^2(1+|h|)+|\chi_L\dot\chi_L|\right).  \tag{4.2}
\]

Translating a fixed cutoff outward makes (4.2) uniformly arbitrarily small
because \(\Omega_\infty^2\) has Gaussian decay in `t`.  Choosing the transitions
far enough away therefore preserves the exact jet at `r=1` and the unique
full edge.

At a flat cutoff endpoint, `h` itself is undefined.  The invariant physical
formula (2.7a) extends there: `Omega_L=S_L=T_L=0`, hence the resonant growth
is zero.  Estimate (4.2) proves that the positive expression from the
interior tends to this value.

Both velocity components acquire the same cutoff, so

\[
                     \Lambda_L=\chi_L\Lambda_\infty. \tag{4.3}
\]

Multiplication by a number in `[0,1]` cannot lower a value below the unique
negative minimum `Lambda_*`; (4.3) preserves its isolation exactly.

The standard pressure-modulated Gavrilov straight limit is

\[
 V_G(s)=\kappa_0s g(s),\qquad
 W_G(s)={\kappa_0s g(s)\over\sqrt2}.                \tag{4.4}
\]

The hollow radial function in (4.1) is smooth and is zero in neighborhoods
of both `s=0` and infinity.  Choosing \(g=\Omega_L/\kappa_0\) in (4.4) therefore
realizes exactly this compact profile at leading straight-tube order.  In
Gavrilov's construction this is the common smooth pressure modulation.  The
sign of `g` is harmless because the modulated Euler pressure is obtained
from `dP_g=g(p)^2 dp`.

This closes C67's **principal-symbol base-profile compatibility** by its
second advertised option: a new spectral certificate in the realizable
Gavrilov family.  It does not yet quantify convergence of the curved torus
coefficients to (4.4), nor the residual of a bent eigenfunction.

## 5. Why literal AO Assumption A fails

Albritton--Ożański define

\[
 q_{\rm AO}=-{(rV)'\over W'}.
\]

For a locked-pitch column this becomes

\[
 q_{\rm AO}
 =-{V+rV'\over\sigma V'}
 =-{r(h+2)\over\sigma(h+1)}.                       \tag{5.1}
\]

Every nonzero smooth compact swirl profile with `V(0)=0` attains a nonzero
interior extremum of `|V|`.  At that point

\[
 V'=W'=0,\qquad (rV)'=V\ne0,                       \tag{5.2}
\]

so (5.1) has a genuine pole.  Therefore no nonzero smooth compact
locked-pitch bump satisfies the literal requirement
`0<inf q<sup q<infinity`.  For (3.1), the pole is at `h=-1`, separated from
the unstable ring `h_*<-3`.  The selected outer-flank ring also has
`q_AO<0` locally.  In addition, `q_AO=0` at `h=-2`, and it is the
indeterminate quotient `0/0` on the hollow core and exterior zero regions.

This is an obstruction to citing Assumption A verbatim, not a physical
singularity of the linearized equation.

## 6. The q-free Rayleigh equation

Let

\[
                    \Gamma=rV.
\]

On every open set where AO's quotient expressions are defined, their two
coefficients simplify algebraically to

\[
 a(r)=r{d\over dr}\left[
 {\beta r^2W'-\Gamma'\over r(1+\beta^2r^2)}
 \right],                                          \tag{6.1}
\]

\[
 b(r)=-{2\beta V(W'+\beta\Gamma')
              \over1+\beta^2r^2}.                  \tag{6.2}
\]

Indeed, \(qW'=-\Gamma'\) in `a`, while

\[
 \Phi={1\over r^3}(\Gamma^2)'={2V\Gamma'\over r^2}
\]

cancels both the numerator and denominator occurrences of `q` in `b`.
This cancellation alone does not define the equation on an open zero region
or at `q=0`; those are `0/0` in the quotient notation.

Instead, directly eliminate `u_theta,u_z` and the pressure from the Fourier
linearized Euler system AO (1.7), using the divergence equation but never
dividing by `W'` or `Gamma'`.  The result is (6.3) with (6.1)--(6.2).
Thus these formulas are the primary coefficients, not merely continuous
extensions of the quotient notation.  They are smooth when `W'=0`; for the
compact hollow profile they vanish identically near the axis and outside a
bounded annulus.  Hence `a,b=O(r^N)` at the axis for every `N`.

The scalar Rayleigh equation

\[
 {d\over dr}\left({r\over1+\beta^2r^2}
 {d\over dr}(ru_r)\right)
 -n^2\left(1+{a\over n\gamma}+{b\over\gamma^2}\right)u_r=0,
 \qquad
 \gamma=n\Lambda-\omega,                           \tag{6.3}
\]

is therefore globally meaningful with no quotient `q_AO`.  More explicitly,
put `D=1+beta^2 r^2` and `H_r=r u_r`.  The direct elimination gives

\[
 P=-{i\gamma r\over n^2D}H_r'
   +{i(\beta r^2W'-\Gamma')\over nrD}H_r.          \tag{6.3a}
\]

Conversely,

\[
 u_\theta={i\Gamma'\over r\gamma}u_r
           +{nP\over r\gamma},\qquad
 u_z={iW'\over\gamma}u_r-{\alpha P\over\gamma}.   \tag{6.3b}
\]

For `r>0`, these formulas are nonsingular when `Im omega>0`, because
`Im gamma=-Im omega`.  At the axis and infinity, the scalar equation
supplies the missing endpoint regularity.  With

\[
 \phi=\left({r^3\over D}\right)^{1/2}u_r,
\]

the hollow-axis indicial exponents are `1/2 +/- n`; the `H_0^1` branch has
`phi=O(r^(n+1/2))` and `u_r=O(r^(n-1))`.  Equations (6.3a)--(6.3b) then give
`P=O(r^n)`, `u_theta=-i u_r+O(r^(n+1))`, and `u_z=O(r^n)`, the smooth
Cartesian branch.  In the zero exterior the positive `beta^2 n^2` mass
gives exponential decay.  This completes the direct equivalence for the
hollow profile rather than relying on continuity across an undefined
quotient.

### 6.1 Proof-transfer audit of AO Theorem 1.1

The v3 source was checked through the complete linear gluing proof.  After
equations (1.9)--(1.10), Section 2 uses only

\[
 p={1+\beta^2r^2\over r^2},\qquad d,\qquad
 a,\qquad b,\qquad\Lambda,                          \tag{6.4}
\]

and their derivatives.  The only explicit later invocation of Assumption A
in the potential estimate, Lemma 2.9, is a **uniform** gap from the level
`Lambda(r_*)` outside a neighborhood of `r_*`.  The axis estimate before
Lemma 2.11 explicitly uses only `a,b=O(r^2)`.  The inner expansion uses
bounded local `C^3` jets together with `b'=Lambda'=0`, `b_*>0`, and
`Lambda''>0`.  The remaining outer estimates use bounded smooth
coefficients, the endpoint potential behavior, and decay.  In the far zone,
the estimate needed is `k=pn^2(1+O(n^-2))`; the dimensionally inconsistent
printed absolute estimate in Lemma 2.9 is not used in stronger form.

Consequently the same proof gives the following specialized q-free version,
which is sufficient here.

> **Compact-hollow q-free AO lemma (self-derived).**  Let `beta>0` and
> `V,W` be real smooth functions compactly supported in an annulus of
> `(0,infinity)`, and put `Omega=V/r`.  Define `Gamma,Lambda,a,b` by
> (6.1)--(6.2).  Suppose that at `r_*>0`
> \[
> b'(r_*)=\Lambda'(r_*)=0,\quad b_*>0,
> \quad\Lambda''(r_*)>0,\quad\Lambda_*\ne0,         \tag{6.4a}
> \]
> and `Lambda^(-1)(Lambda_*)={r_*}`.  Then the existence, asymptotics,
> uniqueness, and Gaussian localization conclusions of AO Theorem 1.1 hold
> for (6.3), without an assumption on `q_AO`; (6.3a)--(6.3b) reconstruct
> smooth unstable Euler modes.

Compact support makes `a,b` smooth, bounded, flat at the axis, and zero in
the exterior.  Since `Lambda=0` near both endpoints and `Lambda_*!=0`, the
unique-level condition implies the quantitative gap

\[
 \inf_{|r-r_*|\ge\rho}|\Lambda(r)-\Lambda_*|>0
 \quad\hbox{for every fixed }\rho>0.               \tag{6.4b}
\]

Thus every hypothesis holds for the compact profile of Section 4.  The
repository protocol still requires independent cross-audit before upgrading
the claim status, but the proof-transfer lemma itself gives exact modes with

\[
\begin{aligned}
 \omega_{m,n}={}&n\Lambda_*+i\lambda_*\\
 &+(1-i)(2m-1)n^{-1/2}
 \left({\lambda_*\Lambda''_*\over8p_*}\right)^{1/2}
 +O(n^{-1+\delta}),                                \tag{6.5}
\end{aligned}
\]

where `p_*=1+beta_*^2`.  Since Section 3 identified `lambda_*` with the
full BAS edge,

\[
 {\operatorname{Im}\omega_{m,n}\over\lambda_{\rm BAS}^{\rm full}}
 \longrightarrow1.                                 \tag{6.6}
\]

This proof-transfer statement is deliberately not labeled cross-verified.
If a hidden use of `q` is found, the local cutoff Weber packet still gives a
q-free quasimode; a complete velocity-reconstruction and residual-norm
ledger would then be the remaining local task.

## 7. Uniform BAS envelope for this profile

The physical two-dimensional BAS theorem of C66 is not specific to the
Batchelor formula.  Its proof needs:

1. bounded `Omega,S,T` and hence bounded normalized coefficients;
2. compactness of the sector where `mu` is within a fixed fraction of the
   positive full resonant edge; and
3. a uniform transverse derivative of the resonance function
   `delta=pS+qT` in that sector.

The compact hollow profile has the first two properties.  At a resonance,
the tangential derivative of `(p,q) -> delta` on the unit circle has
magnitude `sqrt(S^2+T^2)`.  If this vanished, then `S=T=0`, and the BAS
coefficients would give

\[
 A=2\Omega q,\qquad c=-2\Omega q,
 \qquad \mu=Ac=-4\Omega^2q^2\le0.                  \tag{7.1}
\]

Thus it cannot vanish in a compact sector bounded away from zero positive
growth.  The projection-to-resonance estimate in C66 follows, as do the
excess-action and degenerate-sector bounds.  Therefore

\[
 \|\Phi(t;r,\ell,m,k)\|
 \le C e^{\lambda_*|t|}                             \tag{7.2}
\]

uniformly over the full straight-column BAS phase space.  The principal
symbol again has polynomial prefactor `d=0`.

Equation (7.2) is not a frequency-uniform Fourier-block PDE semigroup
estimate.  AO explicitly allow their modal semigroup constant to depend on
`(n,alpha)` after their equation (3.32); that separate gate remains open.

## 8. Updated load-bearing chain

This profile changes the spectral chain as follows.

1. **Closed at principal-symbol level:** the standard pressure-modulated
   Gavrilov family contains a compact profile with a unique, isolated,
   full-edge-matched BAS ring.
2. **Strong self-audited advance:** the AO eigenmode proof appears to
   transfer exactly after replacing the quotient definitions by
   (6.1)--(6.2).  Independent audit is mandatory before treating (6.5) as
   established.
3. **Still open:** make the Fourier-sector semigroup bound uniform over the
   polynomial gain window, including lower-order pressure reconstruction.
4. **Still open:** bend the real mode into the finite-curvature Gavrilov
   torus and prove its normalized residual, expected at first order in the
   aspect ratio.
5. **Still open:** couple the amplified real carrier to the multi-colour
   child-writing transition and the all-order global wake/Gevrey induction.

The exact and interval calculations are reproduced by
`checks/locked_pitch_gavrilov_edge.py`.
