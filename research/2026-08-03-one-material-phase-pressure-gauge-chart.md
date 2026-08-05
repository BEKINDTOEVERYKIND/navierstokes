# One material phase is enough off isotropy: a pressure-gauge boundary chart

**Date:** 2026-08-03
**Status:** exact principal algebra; dynamic realization remains open
**Scope:** pointwise covariance, affine Kelvin compatibility, and the
transition architecture. This note does not construct a Navier--Stokes
singularity or an exact localized transition.

## 1. Correction to the current strategic target

C53 proves the right statement about a covariance neighborhood of isotropy:
one incompressible phase has covariance of rank at most two and therefore
cannot cover a neighborhood of a positive definite matrix. It does **not**
prove that one phase is insufficient for the Reynolds stress modulo pressure.
The pressure quotient changes the dimension count.

Let \(\operatorname{Sym}_3\) denote the six-dimensional space of real
symmetric \(3\times3\) matrices and let

\[
 {\cal B}_2^+
 =\{Q\in\operatorname{Sym}_3:Q\geq0,\qquad
              \operatorname{rank}Q=2\}.                       \tag{1.1}
\]

This is a smooth five-dimensional manifold. A point of \({\cal B}_2^+\)
consists of a kernel direction \(n\in\mathbb S^2/\{\pm1\}\) and a positive
symmetric \(2\times2\) block on \(n^\perp\), hence \(2+3=5\) parameters.
Adding the scalar pressure gauge gives six parameters.

> **Boundary stress chart.** Let \(R_*\in\operatorname{Sym}_3\) have a
> simple largest eigenvalue. In a neighborhood of \(R_*\), every symmetric
> tensor \(R\) admits the unique decomposition
> \[
>        R=\rho(R)I-Q(R),\qquad Q(R)\in{\cal B}_2^+,          \tag{1.2}
> \]
> where
> \[
>        \rho(R)=\lambda_{\max}(R),\qquad
>        Q(R)=\lambda_{\max}(R)I-R.                           \tag{1.3}
> \]
> The maps are smooth, \(\ker Q(R)\) is the top eigendirection of \(R\),
> and \(Q(R)\) is uniformly positive on its transverse plane after the
> neighborhood is made sufficiently small.

Thus one material phase is pointwise sufficient for a full local stress chart
**modulo pressure**, provided the chart is centered off the multiple-eigenvalue
isotropic locus. The covariance lies on the smooth rank-two boundary of the
positive cone rather than in its interior.

This gives a potentially simpler alternative to the current
three-positive-phase principal target:

\[
 \boxed{\text{one material phase + a biased rank-two covariance
               + the global zero-charge wake}.}              \tag{1.4}
\]

It removes all cross-colour interactions at principal order. It does not
remove the material-integrability, localization, endpoint, or wake equations.

## 2. Spectral proof and smoothness

Write the ordered eigenvalues of \(R_*\) as

\[
 \lambda_1(R_*)>\lambda_2(R_*)\geq\lambda_3(R_*).             \tag{2.1}
\]

The simple eigenvalue \(\lambda_1\), its spectral projection \(P_1\), and a
locally oriented unit eigenvector \(n\) are smooth functions of \(R\). For
\(R\) close to \(R_*\), set \(\rho=\lambda_1(R)\). In an eigenbasis,

\[
 Q=\rho I-R
  =\operatorname{diag}
     \big(0,\lambda_1-\lambda_2,\lambda_1-\lambda_3\big).      \tag{2.2}
\]

Both nonzero eigenvalues are positive. Their lower bound is controlled by
the top spectral gap of \(R_*\). Equations (1.2)--(1.3) follow immediately.
Conversely, if \(Q\in{\cal B}_2^+\) has kernel \(n\), then \(n\) is the
simple top eigendirection of \(R=\rho I-Q\). Hence (1.2) is a genuine local
chart, not merely a dimension count.

The isotropic gauge is invisible both to the projected equation and to work:
for every trace-free strain \(S\),

\[
 {\mathbb P}\operatorname{div}(\rho I)=0,\qquad
 (\rho I):S=0.                                               \tag{2.3}
\]

Consequently \(R\) and \(-Q\) have identical projected divergence and work,
up to the convention in which the covariance is moved to the opposite side
of the Euler--Reynolds equation.

## 3. Exact tangent calculation

The spectral formula can be checked without perturbation theory. Take

\[
 n_0=e_1,\qquad Q_0=\operatorname{diag}(0,a,b),\qquad a,b>0. \tag{3.1}
\]

Use local kernel coordinates

\[
 n(x,y)=\frac{(1,x,y)}{\sqrt{1+x^2+y^2}}.                    \tag{3.2}
\]

A smooth transverse frame has first variations

\[
 v_2=e_2-xe_1+O(x^2+y^2),\qquad
 v_3=e_3-ye_1+O(x^2+y^2).                                   \tag{3.3}
\]

Parameterize the positive transverse block by its \(22,33,23\) entries.
At the base point the five covariance variations are

\[
 E_{22},\quad E_{33},\quad E_{23}+E_{32},\quad
 -a(E_{12}+E_{21}),\quad -b(E_{13}+E_{31}).                  \tag{3.4}
\]

Modulo the isotropic line these five vectors form a basis of
\(\operatorname{Sym}_3/\mathbb RI\). Adding the gauge variation \(I\) gives
a basis of all of \(\operatorname{Sym}_3\). In coordinates
\((11,22,33,12,13,23)\), the determinant has magnitude \(ab\). Therefore

\[
 (\rho,n,B)\longmapsto \rho I-Q(n,B)                         \tag{3.5}
\]

is a local diffeomorphism whenever the two transverse eigenvalues are
positive. The checker evaluates both the full determinant and the quotient
rank in exact rational arithmetic.

The two kernel rotations in (3.4) are essential. If the phase direction is
held fixed, a one-phase covariance has only its three transverse-block
coordinates and cannot provide a full trace-free chart.

Every such covariance is realized by one frozen fast phase without choosing
eigenvectors through a repeated transverse eigenvalue. Choose a smooth
orthonormal transverse frame \(E(n)\), put
\(B=E(n)^TQE(n)>0\), and let \(c_1,c_2\) be the columns of
\(E(n)B^{1/2}\). Then

\[
 Q=c_1\otimes c_1+c_2\otimes c_2,\qquad c_1,c_2\perp n.      \tag{3.6}
\]

then the mean-zero profile

\[
 w(\theta)=\sqrt2\,c_1\cos\theta+\sqrt2\,c_2\cos(2\theta)     \tag{3.7}
\]

satisfies \(n\cdot w=0\) and

\[
 \langle w\otimes w\rangle_\theta=Q.                        \tag{3.8}
\]

The two distinct harmonics also give zero averaged helicity by Fourier
orthogonality. A same-harmonic sine--cosine pair would realize the same
covariance but generally carry nonzero mean helicity.

At frozen slow coefficients (3.7) is an exact nonlinear shear:
\((w(n\cdot x)\cdot\nabla)w(n\cdot x)=0\). Slow modulation and curvature
reintroduce lower-order equations, but no second fast colour is needed for
the principal covariance. Distinct harmonics do acquire different viscous
multipliers and additional same-phase charged terms; these belong to the
one-phase hierarchy rather than a cross-colour hierarchy.

## 4. Compatibility with compression and negative stress work

The boundary chart can be centered at a covariance that simultaneously

1. uses a compressive material phase, so its wave number grows;
2. contains a growing transverse Kelvin polarization; and
3. has the sign needed to transfer energy out of the affine parent.

Use the already checked affine strain

\[
 S=\operatorname{diag}\left(-1,-\frac54,\frac94\right)       \tag{4.1}
\]

and choose

\[
 n_0=e_1,\qquad Q_0=\operatorname{diag}\left(0,1,\frac14\right).
                                                                    \tag{4.2}
\]

Then

\[
 Q_0:S=-\frac54+\frac9{16}=-\frac{11}{16}<0.                \tag{4.3}
\]

For a material phase with initial gradient \(e_1\), the affine Kelvin wave
vector is

\[
 k(t)=e^{-S^Tt}e_1=e^t e_1.                                \tag{4.4}
\]

Thus the carrier frequency grows at rate \(1\). The transverse \(e_2\)
polarization obeys the homogeneous Kelvin law \(A'=-SA\), without a pressure
correction, and grows at rate \(5/4\). The \(e_3\) component decays at rate
\(9/4\), but its positive weight \(1/4\) keeps the covariance strictly
positive on \(n_0^\perp\) and preserves the full chart.

More generally, for

\[
 S=\operatorname{diag}(-\alpha,-\beta,\alpha+\beta),\qquad
 Q_0=\operatorname{diag}(0,a,b),                            \tag{4.5}
\]

one has

\[
 Q_0:S=-a\beta+b(\alpha+\beta)<0
 \quad\Longleftrightarrow\quad
 0<b<\frac{a\beta}{\alpha+\beta}.                           \tag{4.6}
\]

This is an open parameter region with \(a,b>0\). Hence positivity,
compression, and the work sign are compatible rather than competing scalar
conditions.

The pointwise statement does not by itself give nonzero global work on a
periodic domain: a constant stress does zero total work against a periodic
velocity gradient. A localized envelope and its mandatory boundary source
remain necessary, exactly as in C38.

## 5. Why this may remove the multiphase bottleneck

Three coordinate material phases have excellent lattice properties (C85),
but different colours generically interact at size \(O(K)\). The lattice gap
prevents small divisors; it does not cancel those leading cross-colour source
terms. A one-phase bath has the exact all-chain \(K\)-cancellation already
derived in the forward multiphase parametrix: every high--high interaction
shares one fast charge and incompressibility removes the nominal factor \(K\).

The rank-two boundary chart supplies the piece that was thought to require
multiple colours. In particular:

* pressure supplies the missing scalar direction;
* the transverse positive block supplies three coordinates;
* rotations of the material kernel supply the remaining two coordinates;
* the one-phase charge algebra avoids leading cross-colour interactions;
* the single-parent/three-partner sideband map already supplies a checked
  rank-five low-strain response and can be rotated so its fast direction is
  the compressive \(e_1\) axis.

The last bullet is only a compatibility observation. It does not prove that
the covariance chart coordinates and sideband endpoint coordinates can be
controlled by one common nonlinear transition.

## 6. What remains genuinely open

The simplification is exact but it trades an interior positive chart for a
constrained boundary chart. The new load-bearing theorem is:

> **Material-aligned one-phase transition theorem (open).** Construct an
> exact localized forward Navier--Stokes transition whose **principal**
> covariance lies in \({\cal B}_2^+\) on the active interior and approaches
> zero through its closure at clean collars; whose kernel has a transported
> material integrating factor; whose endpoint map realizes the required
> rank-five collar variables through order \(M\); and whose zero-charge
> velocity wake is carried into the next affine-capture stage with a uniform
> \(C^M(M!)^2\) bound. Solve the finite-\(K\) longitudinal covariance
> correction by a uniformly controlled perturbative stress chart.

Six hazards are not resolved by the algebra above.

1. **Finite-\(K\) stress.** Exact slow divergence correction generally adds
   an \(O(K^{-1})\) longitudinal velocity and an \(O(K^{-2})\) third
   covariance eigenvalue. Thus only the principal covariance should be
   constrained to \({\cal B}_2^+\); a perturbed pressure-gauge implicit
   function theorem must absorb the full stress while the two transverse
   gaps stay uniform.
2. **Material integrability.** The top eigendirection selected pointwise by
   (1.3) need not admit an integrating factor. Locally the Frobenius
   condition is \(n\cdot\operatorname{curl}n=0\), not
   \(\operatorname{curl}n=0\). The material covector
   \(\xi=\kappa n=\nabla\Phi\) must also obey
   \(D_t\xi+(\nabla U)^T\xi=0\). On \(\mathbb T^3\), a circle-valued phase
   has fixed integral cohomology, so global kernel rotations are constrained.
   C92 subsequently realizes the two constant rotations by exact localized
   fixed-cohomology phases in a frozen core.  It does not prove the coupled
   dynamic covariance submersion.
3. **Endpoint degeneration.** As the covariance vanishes, the two
   kernel-rotation columns scale like its transverse eigenvalues and
   \(D(AA^T)|_{A=0}=0\). Vanishing physical amplitudes through order \(M\)
   forces stress flatness through order \(2M+2\). One must use the doubled
   collar splice of the endpoint-interpolation note, recur/export a
   nonvanishing carrier, or prove a dynamical interior-control theorem.
4. **Joint endpoint coupling.** The two-polarization bath creates additional
   parent--partner couplings. Kernel rotations, transverse-block controls,
   and the checked three sideband partners must yield one uniform endpoint
   inverse, not separate algebraic charts.  C95 subsequently proves that the
   pulse-onset rank survives and all new coefficients are \(O(1)\), but also
   that no finite Fourier corrector closes the charged output.  A growing or
   infinite charged hierarchy is mandatory.
5. **Global zero charge.** C86 routes the formal periodic stress, but the
   nonlinear terminal heat wake and its interaction with the next carrier
   remain.
6. **Localization.** The boundary source required by the exact work identity
   must be included in the wake; it cannot be erased by a pressure gauge.

Accordingly this note does not close the transition theorem. It does identify
a strictly simpler principal architecture that avoids the most dangerous
uncontrolled feature of the three-colour proposal.

## 7. Reproducibility

[The exact checker](../checks/one_material_phase_pressure_gauge_chart.py)
checks, using rational arithmetic:

* the six-dimensional Jacobian determinant and five-dimensional quotient
  determinant;
* the spectral decomposition \(R=\rho I-Q\) at a rational base point;
* positive transverse eigenvalue gaps;
* the work value \(-11/16\);
* affine wave-number growth and transverse polarization growth; and
* the open inequality in (4.6).
