# C193: the certified \(A_2\) cocycle has a contracting line and a fixed-energy polarization filter

**Date:** 2026-08-28

**Status:** exact \(A_2\) jet theorem; outward-certified stable Floquet
line on the C159 orbit; exact principal two-polarization concentration
lemma; no off-ray bundle, finite-frequency composite witness, viscosity,
nonlinear stage, UVSR, or singularity claim

**Checkers:**
[global jets](../checks/a2_global_jet_bounds_c193.py),
[stable line](../checks/c159_stable_line_c193.py), and
[filter algebra](../checks/curl_wkb_polarization_filter_c193.py)

## 0. Result and claim boundary

C191 correctly showed that a scalar multiplier cannot improve
\(\|v\|_\infty/\|v\|_2\).  The C159 monodromy is not scalar: it is
hyperbolic.  C193 certifies its contracting line and proves that the two
polarizations can, at the principal-cocycle level, exchange which spatial
profile carries a fixed amount of energy.

Let \(M\) be the C159 one-period transverse Kelvin monodromy.  In its
returning C159 coefficient frame, its expanding and contracting eigenlines
are

\[
 e_+\parallel(1,x),\qquad e_-\parallel(1,-x),\qquad
                 \boxed{\frac{13}{100}<x<\frac15}.       \tag{0.1}
\]

The C159 frame is orthogonal but not normalized.  In the corresponding
physical orthonormal frame the two slopes are \(\pm y\), where

\[
 y=x|k_0|,\qquad \frac12<y<2.                          \tag{0.2}
\]

Its eigenvalues are \(\rho,\rho^{-1}\), with

\[
                         \boxed{\rho>3000}.              \tag{0.3}
\]

The physical angle between the two lines satisfies

\[
 \sin\angle(e_+,e_-)=\frac{2y}{1+y^2}>\frac45.          \tag{0.4}
\]

Thus the two physical oblique spectral projectors have norm below \(5/4\),
and the physical eigenbasis has condition number below \(2\).

On the schedule \(q=n^8\), put

\[
 R_{\rm filt}=\left\lceil\frac38\log n\right\rceil+1
              =R_\Delta^{(192)}+1,\qquad G=\rho^{R_{\rm filt}}.\tag{0.5}
\]

Then \(G>3000n^3\).  An exact two-profile construction below has the same
\(L^2\) energy at its entrance and endpoint, never exceeds that energy at
any discrete return, and improves its concentration quotient by more than

\[
                            n^3=q^{3/8}.                 \tag{0.6}
\]

This removes C191's *finite-dimensional normalization collision*.  It
does not yet produce a finite-frequency Euler or Navier--Stokes witness.
For that, the broad profile must be put in the covariant stable bundle on
the same shrinking base/covector tube as the localized expanding profile.
C159/C192 certify only the central ray, not that bundle.  C193 therefore
does not credit (0.6) to UVSR or fire any architecture verdict.

## 1. Exact global \(A_2\) jet constants

Retain

\[
 U=N\times\nabla f-\sqrt2 fN,\qquad
 f=\cos a+\cos b+\frac45\cos(a+b).                      \tag{1.1}
\]

In the orthonormal horizontal basis

\[
 e=\frac{r_1+r_2}{\sqrt2},\qquad
 d=\frac{r_1-r_2}{\sqrt6},\qquad
 p=\frac{a+b}{2},\quad q_0=\frac{a-b}{2},               \tag{1.2}
\]

write \(P=\cos p,S=\sin p,C=\cos q_0,Q=\sin q_0\).  If
\(A=DU\), then on the horizontal input plane

\[
 A^TA=3{\cal M},\qquad {\cal M}=H^2+2gg^T,              \tag{1.3}
\]

where

\[
\begin{aligned}
{\cal M}_{11}&=\frac{139+25C^2-75P^2+240CP-160CP^3}{25},\\
{\cal M}_{12}&=\frac{8\sqrt3}{5}SQ(2P^2+1),\\
{\cal M}_{22}&=3(1+3P^2-C^2).
\end{aligned}                                           \tag{1.4}
\]

The checker proves \(12I-{\cal M}\ge0\) exactly.  The only nontrivial
sign reduces, after \(x=4/9+(5/9)t\), to a quartic whose Bernstein
coefficients are

\[
 -\frac{169675667}{6561},\quad-\frac{41568439}{972},\quad
 -\frac{1094911}{18},\quad-\frac{277743}{4},\quad-51732. \tag{1.5}
\]

All are negative.  Equality cases and the higher mode-frame sums then give
the sharp global constants

\[
 \boxed{\|DU\|_{\rm op}\le6,\qquad
        \|D^2U\|_{\rm mult}\le3\sqrt6,\qquad
        \|D^3U\|_{\rm mult}\le9.}                        \tag{1.6}
\]

On the exact zero level the first bound sharpens to

\[
                  \boxed{\sup_{f=0}\|DU\|_{\rm op}
                         =\frac{12\sqrt6}{5}.}           \tag{1.7}
\]

These are genuine physical Euclidean operator/multilinear norms, not
coordinatewise or Frobenius substitutes.  In particular
\(6<350/57\) by the exact margin \(8/57\), which will be used by the
finite-frequency bridge rather than treated as an unnamed \(O(1)\).

For completeness, the higher-jet proof is also modewise and exact.  For
each root \(q_i\) and phase \(\theta_i=q_i\cdot x\), the output brackets in
the second and third derivatives are

\[
 W_{i,2}=\sin\theta_i\,N\times q_i+\sqrt2\cos\theta_i\,N,\qquad
 W_{i,3}=\cos\theta_i\,N\times q_i-\sqrt2\sin\theta_i\,N,\tag{1.8}
\]

and both have norm \(\sqrt6\).  The weighted \(A_2\) frame has eigenvalues
\(13/5,3\), so Cauchy gives the \(3\sqrt6\) second-jet bound.  Generalized
Hölder reduces the third jet to

\[
 S(v)=\sum_i c_i|q_i\cdot v|^3\le3\sqrt{3/2}.           \tag{1.9}
\]

The checker expands both absolute-value branches, verifies their exact
polynomial factorizations, and checks equality cases.  Equality occurs for
the second jet at \(a=b=0\) with both inputs in the \(d\) direction, and for
the third jet at \((a,b)=(\pi/2,-\pi/2)\) with all inputs in that direction.

## 2. Reversibility determines the contracting line

Let \(\Phi(s,t)\) denote the C159 two-dimensional Kelvin propagator in
unit-period time and put

\[
                         A_{1/2}=\Phi(1/2,0).             \tag{2.1}
\]

The \(A_2\) root swap is an orientation-reversing physical symmetry.  The
temporal step is as follows.  If
\(\Gamma=T[3\sqrt2\beta+(2c/3)U^THU/h^2]\), then
\(\gamma'=\Gamma\).  Central inversion \(C(a,b)=(-a,-b)\) gives
\(X(s+1/2)=CX(s)\), sends \(g,U\) to \(-g,-U\), fixes \(H\), and hence
\(\Gamma\circ C=\Gamma\).  Since the certified periodic covector has
\(\gamma(1)=\gamma(0)=0\), the equal half-integrals vanish:
\(\gamma(1/2)=0\) and \(\gamma(s+1/2)=\gamma(s)\).

The reversor \(J(a,b)=(-b,-a)\) fixes the return section and gives
\(X(-s)=JX(s)\).  Physical coordinate reversal \(P\) fixes \(N\) and sends
\(g(-s)=Pg(s)\), \(H(-s)=PH(s)P\), and \(U(-s)=-PU(s)\).  Thus
\(\Gamma(X(-s))=\Gamma(X(s))\) and \(\gamma(-s)=-\gamma(s)\).
Substitution in the actual C159 frame formula makes the diagonal entries
of \(B\) odd and the off-diagonal entries even.  Together with
half-periodicity this gives the exact generator identities

\[
 B(s+1/2)=B(s),\qquad
 B(1/2-s)=-{\cal R}B(s){\cal R},\qquad
 {\cal R}=\operatorname{diag}(1,-1).                    \tag{2.2}
\]

Consequently

\[
 M=A_{1/2}^2,\qquad A_{1/2}{\cal R}A_{1/2}={\cal R}.     \tag{2.3}
\]

Strict Metzler positivity forces

\[
 A_{1/2}=\begin{pmatrix}a&b\\c&a\end{pmatrix},\qquad
 a,b,c>0,\qquad a^2-bc=1.                               \tag{2.4}
\]

No full unstable column is interval-integrated.  Instead, for the first
column write \(r=z_2/z_1\).  On each of C192's first 1024 cells, directed
scalar comparison encloses

\[
 r'=B_{21}+(B_{22}-B_{11})r-B_{12}r^2.                  \tag{2.5}
\]

The complete outward run gives

\[
 0.1396458293744483\le\frac ca\le0.1912097459836738,\qquad
                             a>29.475.                  \tag{2.6}
\]

The stable slope is

\[
 x=\sqrt{\frac cb}
   =\frac{c/a}{\sqrt{1-a^{-2}}}.                         \tag{2.7}
\]

The lower bound in (0.1) follows from \(x>c/a\).  For the upper bound the
checker uses only

\[
 \frac{(24/125)^2}{1-1/29^2}
       =\frac{20184}{546875}<\frac1{25}.                 \tag{2.8}
\]

Equations (2.4), (2.7), and C192's \(\rho>3000\) prove (0.1) and the
multiplier bound (0.3).
The stable multiplier is exactly \(\rho^{-1}<1/3000\).  Moreover, at the
return section

\[
 |k_0|^2=\frac{378}{25}+3\beta^2,\qquad
 \frac{567}{20}<|k_0|^2<\frac{12123}{400},             \tag{2.9}
\]

because the certified \(\beta\) interval lies inside
\((21/10,9/4)\).  Since \(|E_2|/|E_1|=|k_0|\), (0.1) and (2.9) give
(0.2) and (0.4).  This normalization is essential: applying the angle formula
directly to the coefficient slope \(x\) would be incorrect.

## 3. Exact fixed-energy polarization filter

Work first in the exact principal direct-sum model, and choose
\(e_+,e_-\) as physical unit eigenvectors on the two certified lines.
Let scalar profiles \(\phi,\psi\) satisfy

\[
 \|\phi\|_2=\|\psi\|_2=1.                                \tag{3.1}
\]

They and the two eigenvectors need not be orthogonal.  Put

\[
 c=\Re\langle\phi e_+,\psi e_-\rangle,\qquad |c|\le1.     \tag{3.2}
\]

For \(0\le m\le R_{\rm filt}\), define

\[
 v_m=\frac b{\sqrt D}\left(\rho^{m-R_{\rm filt}}\phi e_+
                    +\rho^{-m}\psi e_-\right),\qquad
 D=1+G^{-2}+\frac{2c}{G}.                                \tag{3.3}
\]

The cross coefficient is \(\rho^{m-R_{\rm filt}}\rho^{-m}=G^{-1}\),
independent of \(m\).  If \(x_m=\rho^{2m}\), then

\[
\begin{aligned}
 \|v_m\|_2^2/b^2
 &=D^{-1}\left(\frac{x_m}{G^2}+\frac1{x_m}+\frac{2c}{G}\right),\\
 D-\left(\frac{x_m}{G^2}+\frac1{x_m}+\frac{2c}{G}\right)
 &=\frac{(x_m-1)(G^2-x_m)}{G^2x_m}\ge0.                 \tag{3.4}
\end{aligned}
\]

Therefore

\[
 \boxed{\|v_0\|_2=\|v_{R_{\rm filt}}\|_2=b,\qquad
        \|v_m\|_2\le b\quad(0\le m\le R_{\rm filt}).}    \tag{3.5}
\]

There is no discrete-return energy overshoot.  Within-period overshoot is
not controlled by (3.5); the presently available bound is only a fixed
factor depending on one-period evolution.

For an explicit kinematic profile, on the normalized three-torus let

\[
 D_L(x)=\sum_{j=0}^{L-1}e^{2\pi ijx},\qquad
 \phi_n(x)=\frac{D_{4n^2}(x_1)D_{4n^2}(x_2)D_{4n^2}(x_3)}{(4n^2)^{3/2}},
 \qquad \psi=1.                                         \tag{3.6}
\]

Then

\[
 \|\phi_n\|_2=\|\psi\|_2=1,\qquad
 \|\phi_n\|_\infty=8n^3,\qquad\|\psi\|_\infty=1.         \tag{3.7}
\]

Since \(G>3000n^3\), the triangle and reverse-triangle inequalities give

\[
 \sqrt D\,\|v_0\|_\infty/b<1+\frac8{3000}<\frac{251}{250},\qquad
 \sqrt D\,\|v_{R_{\rm filt}}\|_\infty/b>7n^3.           \tag{3.8}
\]

The common factor \(b/\sqrt D\) cancels between the two concentration
quotients.  Hence

\[
 \boxed{\frac{{\cal C}(v_{R_{\rm filt}})}{{\cal C}(v_0)}
 >\frac{1750}{251}n^3>q^{3/8},\qquad
 {\cal C}(v)=\frac{\|v\|_\infty}{\|v\|_2}.}              \tag{3.9}
\]

Equation (3.9) is the missing normalization mechanism in the complex
principal two-fiber model: the localized expanding profile is seeded at
size \(G^{-1}\), while the broad contracting profile carries the entrance
energy; their roles reverse at the endpoint.  A real, divergence-free,
frequency-localized completion has not yet been constructed and is not
part of this lemma.

The filter clock pays one extra period.  Using C159's \(T<76/25\),

\[
 T R_{\rm filt}<\frac{57}{400}\log q+\frac{152}{25}.     \tag{3.10}
\]

The power coefficient is unchanged; any future explicit error bound must
pay the doubled additive constant.

## 4. What is proved and what is still missing

C193 proves three positive facts with explicit constants:

1. the exact global \(A_2\) flow jets through order three;
2. an outward-certified central contracting Floquet line with a controlled
   eigenbasis; and
3. an exact fixed-energy \(q^{3/8}\) concentration filter in the principal
   two-polarization model.

The third item is not yet an Euler solution.  A genuine packet uses a
three-dimensional off-ray Fourier cap.  Restricting every mode to the
exact C152 periodic plane supplies only \(O(N^2)\) modes in an \(N\)-cap
and cannot realize three-dimensional \(N^{3/2}\) concentration.  Whether
one closed phase can realize the needed action-dependent periodic-covector
family remains an additional compatibility problem; C193 does not claim
an obstruction theorem for it.

The next theorem must therefore certify, on one shrinking but genuinely
three-dimensional packet tube:

* robust expanding and covariant stable bundles with explicit projector
  derivatives;
* the pressure-resolving finite-frequency evolution and transported band;
* packet heat loss, active extraction, depletion, and wake; and
* the endpoint inequality (3.9) for the same exact solution.

Until those are proved, C193 overturns only the scalar-normalization
non-implication, not C191's verdict on the landed C185 operator norm.

## 5. Verification boundary

The three dependency-free checkers:

1. verify the global PSD/quartic/Bernstein jet proof and equality cases;
2. rerun the complete C159 path tube, C192 gain product, half-period
   Riccati comparison, and exact reversibility reduction; and
3. verify the curl/Kelvin sign algebra, rational residual coefficients,
   and fixed-energy filter identities.

They do not certify an off-ray continuity radius, stable bundle over a
localized packet, full endpoint band, viscosity, nonlinear return, UVSR,
or a Navier--Stokes singularity.
