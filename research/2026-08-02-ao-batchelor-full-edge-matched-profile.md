# A Batchelor AO profile matched to the full BAS edge

**Date:** 2026-08-02
**Status:** exact analytic construction with rational interval certificate;
self-derived, not independently cross-audited
**Scope:** straight Euler vortex columns and their complete velocity BAS
cocycle.  This note does not provide a frequency-uniform PDE semigroup
bound or a finite-curvature vortex-ring estimate.

## 1. Main result

There is a unique explicitly characterized Batchelor swirl parameter for
which the Albritton--Ożański ring radius is also the global maximizer over
**all** resonant radii and helical ratios.  For this profile every BAS
Lyapunov exponent is bounded by the AO leading growth rate.

Let `x_*` be the unique positive zero of

\[
 J(x):=(2x+1)e^{2x}-(7x+2)e^x+(5x+1).              \tag{1.1}
\]

Define

\[
 \begin{aligned}
 h(x)&={e^x-1-x\over x^2},\\
 g(x)&={2-e^x\over (x+1)e^x-(2x+1)},\\
 \beta_*&=\sqrt{g(x_*)},\qquad
 Q_*={\beta_*\over h(x_*)},\qquad
 r_* =\sqrt{x_*}.                                  \tag{1.2}
 \end{aligned}
\]

Then the Batchelor column

\[
 V(r)={Q_*\over r}(1-e^{-r^2}),\qquad W(r)=e^{-r^2} \tag{1.3}
\]

satisfies AO Assumption A at `(r_*,beta_*)`.  Moreover, if
`lambda_BAS` denotes any Lyapunov exponent of the full velocity-form BAS,
with arbitrary radius, covector, and helical ratio, then

\[
             \lambda_{\rm BAS}\le\sqrt{b_*},        \tag{1.4}
\]

where equality is attained at `(r_*,beta_*,ell=0)`.

The checker certifies

\[
\begin{aligned}
0.59671214&<x_*<0.59671216,\\
0.5101679&<\beta_*<0.5101682,\\
0.8278572&<Q_*<0.8278581,\\
0.2090085&<b_*<0.2090092,\\
0.4571745&<\sqrt{b_*}<0.4571753.
\end{aligned}                                       \tag{1.5}
\]

The AO source is D. Albritton and W. Ożański, *Linear and nonlinear
instability of vortex columns*, arXiv:2310.20674v3, especially Assumption A
and Appendix A.1:

<https://arxiv.org/abs/2310.20674v3>

The full BAS reduction, including the proof that every non-resonant
trajectory has exponent zero, is given in the companion note
`2026-08-02-ao-batchelor-full-bas-cocycle.md`.

## 2. Resonance curve for a fixed Batchelor profile

Write `x=r^2`.  For a column with fixed swirl `Q>0`, the helical phase
`exp(i(alpha z-n theta))`, `beta=alpha/n`, is resonant precisely when

\[
 \Lambda'(r)=0,
 \qquad
 \beta=\beta_Q(x):=Qh(x).                           \tag{2.1}
\]

The integral representation

\[
 h(x)=\int_0^1(1-s)e^{sx}\,ds                       \tag{2.2}
\]

shows that `h` is strictly increasing from `1/2` to infinity.  Thus (2.1)
parametrizes every positive-ratio resonance.  Negative ratios and ratios
at most `Q/2` have no positive-radius resonance.

At radial tilt `ell=0`, the squared resonant exponent is

\[
 R_Q(x)=4p(x)(1-p(x))
 {e^{-x}(1-e^{-x})\over1+Q^2h(x)^2x},
 \qquad p(x):=Q^2h(x).                              \tag{2.3}
\]

For nonzero radial tilt, this is multiplied by
`K_h^2/(ell^2+K_h^2)<=1`.  The function `R_Q` is positive exactly while
`p<1`; it vanishes at `x=0` and at the unique endpoint where `p=1`.

## 3. Eliminating the two stationarity conditions

For fixed `beta`, differentiating the Batchelor `b` in `x` gives a factor

\[
 H_a(x)=e^{-x}[2(1+ax)+a]-(1+ax+a),
 \qquad a=\beta^2.                                  \tag{3.1}
\]

The derivative

\[
 H_a'(x)=e^{-x}(a-2-2ax)-a<0                        \tag{3.2}
\]

for every `a>0`, so `H_a` has one zero.  Solving that zero for `a` gives

\[
                  a=g(x).                           \tag{3.3}
\]

This is AO's fixed-`beta` condition `b'(r)=0`.

At fixed `x,Q`, differentiation in `beta` shows that the sign of
`partial_beta b` is the sign of

\[
              1-2\beta Q-\beta^2x.                 \tag{3.4}
\]

Hence transverse stationarity is

\[
             \beta Q={1-\beta^2x\over2}.            \tag{3.5}
\]

Along resonance `beta=Qh`, conditions (3.3) and (3.5) become

\[
 Q^2={g(x)\over h(x)^2}
     ={1\over h(x)(2+xh(x))}.                       \tag{3.6}
\]

Equivalently,

\[
                  g(x)={h(x)\over2+xh(x)}.          \tag{3.7}
\]

Clearing the positive denominators in (3.7) yields exactly `J(x)=0` from
(1.1).

## 4. The positive root is unique

The Taylor expansion of (1.1) is

\[
 J(x)=-2x^2+
 \sum_{n=3}^{\infty}
 {2^n(n+1)-(7n+2)\over n!}x^n.                      \tag{4.1}
\]

Every coefficient in the sum is positive.  At `n=3` its numerator is `9`,
and if

\[
 a_n=2^n(n+1)-(7n+2),
\]

then

\[
             a_{n+1}-2a_n=2^{n+1}+7n-5>0.           \tag{4.2}
\]

It follows that `J(x)/x^2` is strictly increasing from `-2` to infinity on
`x>0`; hence it has exactly one zero.  Rational Taylor remainders for the
exponential give the bracket in (1.5), in particular

\[
                    0<x_*<3/5<\log2.                \tag{4.3}
\]

Here the last inequality follows, for example, from
`log(1+t)>2t/(2+t)` at `t=1`, which gives `log 2>2/3`.

## 5. Global maximality along every resonance

It remains to prove that the joint critical point is global, not merely
local.  Two monotonicities make this exact.

First, `g` is strictly decreasing on `(0,log 2)`.  One self-contained way
to see this is to use (3.1): `H_a` is strictly decreasing in `x`, and

\[
 \partial_aH_a=e^{-x}(2x+1)-(x+1)<0                \tag{5.1}
\]

for `x>0`.  The latter inequality follows because
`(x+1)e^x-(2x+1)` has positive derivative after its double zero at the
origin.  Therefore the value `a=g(x)` required to restore `H_a=0`
decreases strictly as `x` increases.  Since `h` increases, the function

\[
                      {h(x)^2\over g(x)}             \tag{5.2}
\]

is strictly increasing on `(0,log 2)`.

Second, `h(x)(2+xh(x))` is strictly increasing for all `x>0`.

Let `Q=Q_*` and follow the resonance curve `beta(x)=Q_*h(x)`.  By (3.6),
at `x=x_*` both partial derivatives `partial_x b` and `partial_beta b`
vanish.  If `0<x<x_*`, monotonicity gives

\[
 Q_*^2h(x)^2<g(x),\qquad
 Q_*^2h(x)(2+xh(x))<1.                              \tag{5.3}
\]

Equations (3.1) and (3.4) show that both partial derivatives are then
strictly positive.  Because `beta'(x)=Q_*h'(x)>0`,

\[
 {d\over dx}R_{Q_*}(x)
 =\partial_xb+\partial_\beta b\,\beta'(x)>0.         \tag{5.4}
\]

For `x_*<x<log 2`, both inequalities and both signs reverse.  For
`x>=log 2`, formula (3.1) is negative for every `a>0`, while the expression
in (3.4) stays negative by monotonicity.  Hence

\[
 {d\over dx}R_{Q_*}(x)<0                            \tag{5.5}
\]

throughout the rest of its positive interval.  Thus `x_*` is the unique
global maximum of every positive resonant exponent for the fixed profile.

The full-BAS conclusion (1.4) now follows from three facts established in
the companion cocycle note: non-resonant exponents are zero, radial tilt
only lowers (2.3), and every exact positive-ratio resonance lies on (2.1).

## 6. AO Assumption A and the spectral ratio

At the constructed point,

\[
 {\beta_*\over Q_*}=h(x_*)>{1\over2},
 \qquad
 \beta_*Q_*={1\over2+x_*h(x_*)}<{1\over2}<1.        \tag{6.1}
\]

Thus `Q_*/2<beta_*<1/Q_*`.  The Batchelor calculation in AO Appendix A.1,
or the elementary proof in the earlier profile note, gives:

* `Lambda` has a unique strict global minimum at `r_*`;
* `Lambda''(r_*)>0` and no other radius has the same `Lambda` value;
* `b_beta_*` has a unique strict global maximum at `r_*`;
* the regularity, parity, decay, and bounded-positive-`q` requirements hold.

Therefore AO Theorem 1.1 applies.  Since the full BAS edge now equals
`sqrt(b_*)`, its fixed-`m` ring eigenvalues satisfy

\[
 {\operatorname{Im}\omega_{m,n}\over\Lambda_{\rm BAS}^{\rm full}}
 =1-D_mn^{-1/2}+O(n^{-1+\delta})\longrightarrow1.   \tag{6.2}
\]

This closes the profile-selection and full-cocycle-ratio gates.  It does
not close the non-effective AO remainder, the frequency-uniform PDE
propagator prefactor, or the finite-curvature residual.
