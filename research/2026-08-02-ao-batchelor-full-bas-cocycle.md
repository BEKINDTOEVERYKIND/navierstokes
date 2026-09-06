# Full BAS reduction for the explicit AO Batchelor profile

**Date:** 2026-08-02
**Status:** exact self-derived cocycle reduction and rationally certified
counterexample; not independently cross-audited
**Scope:** the velocity-form bicharacteristic-amplitude system of the
straight Euler vortex column.  This is not a uniform PDE semigroup estimate
and does not address bending into a vortex ring.

## 1. Result

For the exact profile in
`2026-08-02-ao-batchelor-global-bas-certificate.md`, every non-resonant BAS
trajectory has Lyapunov exponent zero.  In the single helical sector with
the selected ratio `beta`, the only positive resonant exponent is at `r_0`
and is at most `sqrt(b_0)`.

The corresponding statement for the **full** BAS cocycle, allowing the
helical ratio to vary, is false.  An explicit second resonant radius

\[
             r_1^2=\log(7/4)                         \tag{1.1}
\]

has squared exponent greater than `0.21034`, whereas the selected AO edge
is less than `0.21020`.  Thus

\[
 {\lambda_{\rm full}\over\sqrt{b_0}}>1.0003.         \tag{1.2}
\]

The gap is small but strict and theorem-grade.  Therefore global
maximization of `b(r)` at fixed `beta` does not by itself identify the full
BAS edge.

The BAS equations used below are the standard velocity-amplitude equations

\[
 \dot x=u(x),\qquad \dot\xi=-(\nabla u)^T\xi,
 \qquad
 \dot a=-\nabla u\,a
       +2{\xi\otimes\xi\over|\xi|^2}\nabla u\,a,         \tag{1.3}
\]

as stated, for example, in R. Shvydkoy, *The essential spectrum of
advective equations*, arXiv:math-ph/0412019v2, equation (2):

<https://arxiv.org/abs/math-ph/0412019>

The vortex profile and ring-mode asymptotics come from the primary source
D. Albritton and W. Ożański, *Linear and nonlinear instability of vortex
columns*, arXiv:2310.20674v3:

<https://arxiv.org/abs/2310.20674v3>

## 2. Exact cylindrical reduction

Fix a particle trajectory at radius `r>0`.  (At the axis the Batchelor
profile has only solid rotation to first order, so the BAS exponent is
zero.)  In the orthonormal frame
`(e_r,e_theta,e_z)` rotating with angular velocity `Omega=V/r`, put

\[
 S=r\Omega'=V'-\Omega,\qquad T=W'.                   \tag{2.1}
\]

The velocity gradient and frame-rotation matrices are

\[
 A=\begin{pmatrix}
 0&-\Omega&0\\ \Omega+S&0&0\\ T&0&0
 \end{pmatrix},\qquad
 R=\begin{pmatrix}
 0&-\Omega&0\\ \Omega&0&0\\0&0&0
 \end{pmatrix}.                                     \tag{2.2}
\]

Write the covector in this frame as

\[
                 \xi=(\ell,m,k).                    \tag{2.3}
\]

Transforming (1.3) into the rotating frame gives

\[
 \dot\xi=-(A-R)^T\xi,
 \qquad
 \dot a=-(A+R)a+2\xi{\xi^TAa\over|\xi|^2}.         \tag{2.4}
\]

Consequently

\[
 \dot\ell=-D,\qquad \dot m=\dot k=0,
 \qquad D:=mS+kT.                                   \tag{2.5}
\]

If `K_h^2=m^2+k^2=0`, the divergence constraint gives a zero radial
amplitude and (2.4) leaves the two transverse components constant, so its
exponent is zero.  Hence assume `K_h>0`, put `K^2=ell^2+K_h^2`, and introduce

\[
 e_s={1\over K_h}(0,m,k),\qquad
 e_h={1\over K_h}(0,k,-m),\qquad
 p=e_r-{\ell\over K_h}e_s.                          \tag{2.6}
\]

The divergence constraint `xi dot a=0` permits the unique representation

\[
                       a=xp+ye_h.                   \tag{2.7}
\]

A direct substitution into (2.4) gives the exact two-dimensional system

\[
 \begin{aligned}
 \dot x&={2\ell D\over K^2}x
        +{2\Omega kK_h\over K^2}y,\\
 \dot y&={mT-k(S+2\Omega)\over K_h}x.
 \end{aligned}                                      \tag{2.8}
\]

Define the constants

\[
 A_0=2\Omega kK_h,qquad
 C={mT-k(S+2\Omega)\over K_h},qquad
 \mu=A_0C.                                          \tag{2.9}
\]

Since `(K^2)'=-2ell D`, setting `z=K^2x` reduces (2.8) further to

\[
                 \dot z=A_0y,qquad
                 \dot y={C\over K^2}z.             \tag{2.10}
\]

This identity is the useful global symmetrization of the column BAS.

## 3. Non-resonant trajectories have zero exponent

Suppose `D!=0`.  Then `ell(t)=ell_0-Dt`.  If `C=0`, equation (2.10) makes
`y` constant and `z` at most linear.  If `A_0=0`, `z` is constant and `y`
is bounded because `1/K^2` is integrable in time.

In the remaining case, eliminating `z` yields

\[
                  (K^2\dot y)^{\displaystyle\cdot}=\mu y.    \tag{3.1}
\]

Set

\[
 \tau={\ell_0-Dt\over K_h},\qquad
 s=\operatorname{arsinh}\tau,
 \qquad \kappa={\mu\over D^2}.                     \tag{3.2}
\]

Equation (3.1) becomes

\[
             y_{ss}+\tanh(s)y_s-\kappa y=0.         \tag{3.3}
\]

The coefficient matrix of the first-order system for `(y,y_s)` is bounded
by a constant depending only on `kappa`.  Gronwall therefore gives

\[
 |y(s)|+|y_s(s)|
 \le C_\kappa e^{C_\kappa|s-s_0|}
 \le C'_{\kappa,\ell_0,K_h}(1+|t|)^{C_\kappa}.      \tag{3.4}
\]

Equation (2.10) gives the same type of polynomial bound for `z` and `x`.
The change from `(x,y)` to the physical amplitude has only polynomial
condition number because `|p|^2=1+ell^2/K_h^2`.  Applying the same estimate
to the inverse fundamental matrix (equivalently, reversing the `s` interval)
also gives a polynomial lower bound for every nonzero solution.  Thus every
nonzero non-resonant solution has Lyapunov exponent exactly zero.

This conclusion is not uniform as `D` tends to zero: `kappa=mu/D^2`
diverges, and the polynomial power and constants in (3.4) degenerate.  On
the time scale `|t| << K_h/|D|`, a near-resonant trajectory can
mimic the constant-coefficient exponential dynamics.  This is the precise
cocycle-level reason that a frequency-uniform semigroup prefactor does not
follow merely from the Lyapunov calculation.

## 4. Exact resonance

If `D=0`, then `ell` and `K` are constant and (2.8) has eigenvalues

\[
                 \lambda_\pm=\mathord\pm\sqrt{\mu/K^2}       \tag{4.1}
\]

when `mu>0`; otherwise its Lyapunov exponents are zero.

For the helical phase `exp(i(alpha z-n theta))`, set

\[
 m=-{n\over r},\qquad k=\alpha=\beta n.             \tag{4.2}
\]

Then

\[
 D=n(\beta W'-\Omega')=n\Lambda'(r).                \tag{4.3}
\]

At resonance, using `q=-(rV)'/W'`, direct algebra gives

\[
 \mu=b_\beta(r)K_h^2,
 \qquad
 \lambda_+^2=b_\beta(r){K_h^2\over\ell^2+K_h^2}.   \tag{4.4}
\]

For the selected `beta` in the companion note, `Lambda'` has the unique
positive zero `r_0`, and `b_beta(r_0)=b_0`.  Equations (3.4) and (4.4)
therefore prove the fixed-helical-sector bound

\[
            \lambda_{\rm BAS}^{(\beta)}\le\sqrt{b_0},         \tag{4.5}
\]

with equality at `r=r_0`, `ell=0`.

## 5. A strict full-cocycle counterexample

Keep the swirl parameter `Q` of the exact profile,

\[
 X=\log(9/5),\qquad
 \beta_0=(4-X)^{-1/2},\qquad
 Q={\beta_0X^2\over4/5-X}.                          \tag{5.1}
\]

For the Batchelor column, every positive radius `x=r^2` becomes resonant
for the helical ratio

\[
             \beta_{\rm res}(x)
       =Q{e^x-1-x\over x^2}.                        \tag{5.2}
\]

Choose

\[
                    Y=\log(7/4).                    \tag{5.3}
\]

Then `e^{-Y}=4/7`, and put `beta_1=beta_res(Y)`.  The squared exponent at
radial tilt `ell=0` is exactly

\[
 B_1=4\beta_1Q(1-\beta_1Q)
     {{4\over7}{3\over7}\over1+\beta_1^2Y}.         \tag{5.4}
\]

The dependency-free rational interval checker uses the positive atanh
series for both logarithms and certifies

\[
 \begin{aligned}
 0.53579&<\beta_1<0.53581,\\
 0.21034&<B_1<0.21036,\\
 b_0&<0.21020,\\
 B_1&>(1.0003)^2b_0.
 \end{aligned}                                      \tag{5.5}

Therefore the resonant orbit `(Y,beta_1,ell=0)` grows faster than the AO
ring selected in (5.1).  If a rational Fourier ratio is required, the
inequality persists by continuity and rational ratios are dense; the
resonance radius varies continuously with `beta` by AO Appendix A.1.

## 6. Consequences and repair target

1. Non-resonant trajectories are not the source of a larger asymptotic
   exponent.  They have polynomial growth only.
2. The exact profile gives a sharp bound in its invariant fixed-helical
   sector, but not for the unrestricted cocycle.
3. The failure is a transverse optimization failure: the chosen point is
   stationary in `r` at fixed `beta`, but not stationary when `beta` is
   allowed to vary along the resonance curve.
4. A nearby repair appears possible.  Along the resonance curve, simultaneous
   stationarity in `r` and `beta` adds the exact condition

   \[
        \beta Q={1-\beta^2x\over2}.                 \tag{6.1}
   \]

   Numerical scouting places the resulting Batchelor parameter near
   `x=0.59671215`, `Q=0.82785756`.  This is only a repair candidate here:
   uniqueness and global maximality along the full resonance curve have not
   been interval-certified.
5. Even after repairing the transverse maximum, the degeneration
   `mu/D^2` near resonance and the PDE resolvent constants remain separate
   uniformity problems.

Thus the full-cocycle bound requested for the exact `log(9/5)` profile is
refuted, while its fixed-sector version is proved exactly.
