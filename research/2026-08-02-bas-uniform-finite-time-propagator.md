# A frequency-uniform finite-time bound for the full column BAS cocycle

**Date:** 2026-08-02
**Status:** exact self-derived BAS theorem with a dependency-free arithmetic
checker; degenerate-sector proof independently cross-audited in session
**Scope:** the two-dimensional velocity-amplitude cocycle of a straight
Batchelor vortex column.  This is a principal-symbol statement, not a bound
for the Euler or Navier--Stokes PDE semigroup and not a curved-ring estimate.

## 1. Result

Let `lambda_* > 0` be the largest resonant BAS exponent of the full
Batchelor column, allowing the radius and the helical ratio to vary.  Then
the physical velocity-amplitude propagator satisfies

\[
       \|\Phi(t;r,\ell,m,k)\|_{\mathbb R^3\to\mathbb R^3}
       \le C e^{\lambda_* |t|}                       \tag{1.1}
\]

on the divergence-free amplitude plane.  The constant is independent of
`ell,m,k`, after excluding the zero covector, and is uniform in `r > 0`.
Thus the smallest nonnegative polynomial exponent in

\[
       \|\Phi(t)\|\le C(1+|t|)^d e^{\lambda_*|t|}   \tag{1.2}
\]

is

\[
                         \boxed{d=0}.               \tag{1.3}
\]

The apparent singularity `mu/D^2` in the Legendre reduction is a coordinate
singularity.  At the growing edge, its possible excess action is bounded by

\[
 {4\over3|D/K_h|\,\lambda_*^2}
 \left(\mu+{(D/K_h)^2\over4}-\lambda_*^2\right)_+^{3/2}.       \tag{1.4}
\]

For the Batchelor column, the bracket is `O(|D|/K_h)` as `D` tends to zero.
Consequently (1.4) is `O(sqrt(|D|/K_h))`, not an inverse power of `D`.

The proof applies both to the explicit `log(9/5)` profile and to the nearby
simultaneously optimized Batchelor profile constructed in
`2026-08-02-ao-batchelor-full-edge-matched-profile.md`.  It uses the actual full
resonant maximum in each case; it does not assert that the selected AO ring
of the first profile attains that full maximum.

## 2. Frequency normalization and the physical two-by-two system

Use the notation of `2026-08-02-ao-batchelor-full-bas-cocycle.md`.  For
`K_h=sqrt(m^2+k^2)>0`, put

\[
 p={m\over K_h},\qquad q={k\over K_h},\qquad
 \tau={\ell\over K_h},\qquad p^2+q^2=1,             \tag{2.1}
\]

and define

\[
 \delta=pS+qT={D\over K_h},\qquad
 A=2\Omega q,\qquad c=pT-q(S+2\Omega),\qquad
 \mu=Ac.                                            \tag{2.2}
\]

The frequency magnitude has disappeared.  Since `tau'=-delta`, the exact
system (2.8) of the companion note is

\[
 \dot x={2\tau\delta\over1+\tau^2}x
        +{A\over1+\tau^2}y,
 \qquad
 \dot y=cx.                                         \tag{2.3}
\]

The physical amplitude is `a=x p_vec+y e_h`, with
`|p_vec|^2=1+tau^2`.  Hence

\[
             u=\sqrt{1+\tau^2}\,x,
 \qquad |a|^2=u^2+y^2.                              \tag{2.4}
\]

Writing `g=1+tau^2`, (2.3) becomes

\[
 {d\over dt}\binom u y=
 \begin{pmatrix}
  \delta\tau/g&A/\sqrt g\\
  c/\sqrt g&0
 \end{pmatrix}\binom u y.                          \tag{2.5}
\]

This is the useful normalization: it is homogeneous of degree zero in the
covector and its Euclidean norm is exactly the physical amplitude norm.

If `K_h=0`, the divergence constraint and the original BAS equations give
a constant transverse amplitude, so (1.1) is immediate.

## 3. The growing sector and the exact excess-action interval

Suppose `mu>0`.  With a harmless sign in the second coordinate, conjugate
(2.5) by

\[
                  \operatorname{diag}(\sqrt{|c|},\sqrt{|A|}). \tag{3.1}
\]

The matrix becomes symmetric:

\[
 N(t)=\begin{pmatrix}a(t)&h(t)\\h(t)&0\end{pmatrix},
 \qquad
 a={\delta\tau\over g},\qquad h={\sqrt\mu\over\sqrt g}.       \tag{3.2}
\]

Reflecting `tau` if necessary, set `d=|delta|` and write `a=d tau/g`.
The upper eigenvalue is

\[
 q_+(\tau)={1\over2}\left({d\tau\over g}
       +\sqrt{{d^2\tau^2\over g^2}+{4\mu\over g}}\right).     \tag{3.3}
\]

Fix any comparison rate `lambda>0`.  Since the lower eigenvalue is
non-positive,

\[
 q_+>\lambda
 \quad\Longleftrightarrow\quad
 F(\tau):=\lambda^2(1+\tau^2)-\lambda d\tau-\mu<0.             \tag{3.4}
\]

Completing the square gives

\[
 F(\tau)=\lambda^2\left(\tau-{d\over2\lambda}\right)^2-E,
 \qquad E=\mu+{d^2\over4}-\lambda^2.                \tag{3.5}
\]

Furthermore, evaluating the characteristic polynomial at `lambda` yields

\[
 0<q_+-\lambda
 \le {-F(\tau)\over\lambda(1+\tau^2)}
 \le {-F(\tau)\over\lambda}                         \tag{3.6}
\]

on the interval where `F<0`.  Since `|dt|=|d tau|/d`, integration over the
entire real `tau` axis proves the exact majorant

\[
 \int (q_+-\lambda)_+dt
 \le {4E_+^{3/2}\over3d\lambda^2}.                  \tag{3.7}
\]

For `d=0`, the left side vanishes whenever `mu<=lambda^2`.

Equation (3.7) is sharper than Gronwall in the Legendre variable.  It is
also the complete explanation of the near-resonant limit: only the finite
parabolic interval (3.5), not an `O(1/d)` time interval, can grow faster
than the resonant comparison rate.

## 4. Why the excess is uniform for a Batchelor column

Define the full resonant edge by

\[
 \lambda_*^2=max\{\mu(r,p,q)_+:
       p^2+q^2=1,\ \delta(r,p,q)=0\}.               \tag{4.1}
\]

The maximum is positive.  For the Batchelor profile

\[
 \begin{aligned}
 \Omega&=Q{1-e^{-r^2}\over r^2},\\
 S&={2Q\over r^2}\big((r^2+1)e^{-r^2}-1\big),\\
 T&=-2re^{-r^2},\\
 S+2\Omega&=2Qe^{-r^2}.                             \tag{4.2}
 \end{aligned}
\]

Thus `A,c,delta` are bounded uniformly in `r,p,q`.  Moreover
`sup_(p,q) mu_+(r,p,q)` tends to zero both as `r` tends to zero and as
`r` tends to infinity.  The set

\[
       \mathcal K=\{(r,p,q):\mu\ge\lambda_*^2/2\}   \tag{4.3}
\]

is consequently compact and stays away from the axis and infinity.

At every fixed `r>0`, the map

\[
       (p,q)\longmapsto\delta=pS+qT                 \tag{4.4}
\]

has nonzero tangential derivative at its two zeros on the unit circle,
because `T<0`.  On the compact set (4.3), its size is uniformly bounded
away from zero.  Projecting `(p,q)` to the nearest resonant direction and
using the uniform Lipschitz bound for the quadratic function `mu=Ac`
therefore gives

\[
       \mu-\lambda_*^2\le L|\delta|
       \quad\hbox{whenever}\quad
       \mu\ge\lambda_*^2/2.                         \tag{4.5}
\]

Here the resonant value at the projected direction is at most
`lambda_*^2` by (4.1).  Points not lying in a sufficiently small resonant
tube are absorbed by increasing `L`.

Taking `lambda=lambda_*` in (3.7), (4.5) and the boundedness of `delta`
give

\[
       \sup_{r,p,q,\,\mu\ge\lambda_*^2/2}
       \int(q_+-\lambda_*)_+dt<\infty.              \tag{4.6}
\]

More sharply, the right side tends to zero like `O(sqrt(d))` as `d` tends
to zero.  Because `Ac>=lambda_*^2/2` in this region while `A,c` are bounded,
both `|A|` and `|c|` are bounded away from zero.  The symmetrizer (3.1) and
its inverse therefore have uniform condition number.  Equations (3.2) and
(4.6) prove (1.1) in the only sector that can approach the exponential
edge.

## 5. Degenerate sectors do not require a polynomial prefactor

It remains to make uniformity through `A=0`, `c=0`, and `mu<=0` explicit.
This is where a norm estimate for the symmetrized system alone is
insufficient: the diagonal conjugation degenerates.  The missing factors
are recovered by estimating its two columns separately.

### 5.1 Positive but separated from the edge

On `0<mu<lambda_*^2/2`, choose once and for all

\[
             {\lambda_*\over\sqrt2}<\lambda_0<\lambda_*.      \tag{5.1}
\]

Formula (3.7), now with `lambda=lambda_0`, is uniformly bounded on every
time subinterval.  For small `d` its bracket is negative by the strict gap
in (5.1), while for `d` bounded away from zero all coefficients lie in a
compact set.  Hence the symmetrized propagator satisfies

\[
             \|\Psi(t,s)\|\le C e^{\lambda_0|t-s|}.            \tag{5.2}
\]

Let `H_a(t,s)` be the scalar propagator for `z'=a z` and put
`h_0=sqrt(mu)`.  For the solution starting in the second coordinate, the
first coordinate is

\[
 x(t)=\int_0^t H_a(t,s){h_0\over\sqrt{g(s)}}y(s)\,ds.          \tag{5.3}
\]

For the solution starting in the first coordinate, the second coordinate
is the analogous integral of `h_0 x/sqrt(g)`.  Estimate (5.2) and the scalar
geometric estimate (5.7) below give

\[
 |\Psi_{12}(t)|+|\Psi_{21}(t)|
 \le C h_0(1+|t|)e^{\lambda_0|t|}.                  \tag{5.4}
\]

Conjugating back multiplies these entries by `sqrt(|A/c|)` and
`sqrt(|c/A|)`, respectively.  The factor `h_0` converts them exactly to
`|A|` and `|c|`, which are uniformly bounded.  Thus

\[
 \|\Phi(t)\|\le C(1+|t|)e^{\lambda_0|t|}
              \le C' e^{\lambda_*|t|}.              \tag{5.5}
\]

### 5.2 Oscillatory and triangular sectors

If `mu<0`, the same conjugation turns (2.5) into

\[
       \begin{pmatrix}a&h\\-h&0\end{pmatrix},
       \qquad h={h_0\over\sqrt g},\qquad h_0=\sqrt{-\mu}.     \tag{5.6}
\]

Its symmetric part is `diag(a,0)`.  Put
`P_[s,t]=1+d|t-s|`.  Along every monotone `tau` interval the energy identity
gives

\[
 \|\Psi(t,s)\|\le\exp\left(\int_s^t a_+d\zeta\right)
 \le P_{[s,t]}.                                      \tag{5.7}
\]

The last inequality follows from
`integral tau/(1+tau^2) d tau=(1/2)log(1+tau^2)` and the fact that
`sqrt(1+tau^2)` is one-Lipschitz.  Applying (5.3) to the second column and
integrating `y'=-h_0x/sqrt(g)` for the first column gives

\[
 |\Psi_{12}(t)|\le h_0|t|P_{[0,t]}^2,
 \qquad
 |\Psi_{21}(t)|\le h_0|t|P_{[0,t]}.                 \tag{5.8}
\]

Conjugating back again turns `h_0` into `|A|` or `|c|`.  Since `A,c,d` are
uniformly bounded, including at the axis and infinity,

\[
                    \|\Phi(t)\|\le C(1+|t|)^3.      \tag{5.9}
\]

At `mu=0`, (2.5) is triangular.  If `A=0`, then
`u=H_a u_0` and `y-y_0=c integral(u/sqrt(g))`; if `c=0`, then `y=y_0` and
`u=H_a u_0+A integral(H_a/sqrt(g))y_0`.  Estimate (5.7) again gives a
uniform polynomial.  Since `lambda_*>0`, (5.9) is bounded by a constant
times `exp(lambda_*|t|)`.

The same interval estimates apply backward in time after reversing the
`tau` path.  Combining Sections 4 and 5 proves (1.1), including uniform
passage through all coefficient degeneracies.  Since (1.1) already has no
polynomial loss, `d=0` is minimal among nonnegative exponents.

## 6. Relation to the special-function equation

When `delta!=0`, eliminating the second coordinate gives

\[
  (1+\tau^2)y_{\tau\tau}+2\tau y_\tau
       -{\mu\over\delta^2}y=0.                     \tag{6.1}
\]

This is a Legendre equation at imaginary argument.  With
`tau=sinh(s)` and `y=(cosh s)^(-1/2)v`, it becomes

\[
 v_{ss}-\left({\mu\over\delta^2}+{1\over4}
                    +{1\over4}\operatorname{sech}^2s\right)v=0. \tag{6.2}
\]

The large parameter in (6.2) explains the near-resonant exponential on a
fixed `s` interval, but it obscures the comparison with physical time.
The excess-action calculation (3.4)--(3.7) retains that comparison and
shows why no `exp(C/|delta|)` prefactor survives after the resonant rate is
subtracted.

## 7. What this does and does not close

This closes the BAS-level frequency/radius/helical-ratio prefactor problem:
the full column amplitude cocycle has a uniform bound with `d=0` at its
true full resonant edge.

It does **not** identify the essential semigroup edge of a Fourier block.
The BAS bound concerns the principal-symbol cocycle only.  Lower-order
transport, pressure reconstruction, radial boundary behavior, and the
mode-dependent constants in the Albritton--Ozanski PDE argument remain.
It also does not realize the Batchelor base profile in a compact Euler ring.
The standard Gavrilov thin-tube profile has a strict first-jet mismatch, as
recorded in `2026-08-02-batchelor-gavrilov-profile-compatibility.md`.
Base-profile realization and the subsequent curvature residual are separate
obligations in the long-gain quasimode ledger.
