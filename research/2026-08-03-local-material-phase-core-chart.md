# Localized one-phase material core: cohomology repair and finite-\(K\) stability

**Date:** 2026-08-03

**Status:** exact fixed-cohomology phase construction and frozen-core stress
algebra; conditional pointwise perturbative inverse.  The dynamically
coupled localized oscillatory transition remains open.

**Scope:** this note repairs two qualifications in C89 for the affine active
core. It does not solve a general spatial stress field, the endpoint splice,
or the global wake.

## 1. Outcome

The two kernel-rotation directions in the one-phase pressure-gauge chart can
be realized by exact transported phases without changing the phase's
normalized integral cohomology class on \(\mathbb T^3\).

Let \(B_0\Subset B\) be coordinate balls in the torus. There is a
two-parameter family of circle-valued phases \(\Phi_\eta\), all in the fixed
class \([d\Phi_\eta/(2\pi)]=[dx_1/(2\pi)]\), such that

\[
 d\Phi_\eta=e_1+\eta_2e_2+\eta_3e_3
 \quad\hbox{on }B_0.                                         \tag{1.1}
\]

At the initial frozen core the phase parameters supply exactly the two
kernel rotations which, together with the three transverse covariance
coordinates and one pressure gauge, make the C89 algebraic stress map a
six-dimensional local diffeomorphism.

If \(\Phi_\eta\) is transported by one common, parameter-independent coarse
flow, its covector satisfies the material equation exactly and its
cohomology class is preserved. Thus Frobenius integrability and torus
cohomology are not finite-dimensional rank obstructions for the two phase
controls in a localized frozen affine core. This statement alone does not
prove that the **dynamically evolved** velocity covariance retains the
stress submersion. The companion Kelvin-covariance note (C96) proves that
submersion for a common parameter-independent affine flow, while retaining
the conditioning and localization qualifications below.

The derivative of the chart at one point is also stable under a finite-\(K\)
stress correction. If the
complete principal-plus-corrector stress map is

\[
 {\cal F}_K(z)={\cal F}_0(z)+E_K(z)                            \tag{1.2}
\]

and \(A=D{\cal F}_0(z_0)\), then

\[
 \|A^{-1}DE_K\|\le\vartheta<1                                \tag{1.3}
\]

implies

\[
 \|(D{\cal F}_K)^{-1}\|
 \le\frac{\|A^{-1}\|}{1-\vartheta}.                           \tag{1.4}
\]

Therefore an \(O(K^{-1})\) longitudinal divergence correction does not
destroy pointwise derivative invertibility once its parameter derivative is
proved uniformly small. A right inverse onto one fixed target neighborhood
additionally needs a uniform \(C^0\) estimate, neighborhood-uniform
derivative control, and a quantitative second-derivative/Lipschitz bound.
None of those analytic WKB estimates is a result of this note.

## 2. Fixed-cohomology localized phases

Use coordinates \(x\in(\mathbb R/2\pi\mathbb Z)^3\). Choose \(B\) inside
one injective standard coordinate chart and choose a Gevrey--2 periodic
cutoff \(\chi\) with \(\operatorname{supp}\chi\Subset B\) and \(\chi=1\)
on \(B_0\). Let \(y_2,y_3\) be the corresponding single-valued lifted
coordinate functions on \(B\), and set

\[
 \psi_2=\chi y_2,\qquad \psi_3=\chi y_3,                      \tag{2.1}
\]

extended by zero outside \(B\). These are globally smooth real-valued
periodic functions. Define the circle-valued phase

\[
 \Phi_\eta(x)=x_1+\eta_2\psi_2(x)+\eta_3\psi_3(x)
 \pmod{2\pi}.                                                 \tag{2.2}
\]

Because \(d\psi_2,d\psi_3\) are exact,

\[
 [d\Phi_\eta/(2\pi)]=[dx_1/(2\pi)]
 \in H^1(\mathbb T^3;\mathbb Z)                              \tag{2.3}
\]

for every \(\eta\). On \(B_0\), (2.1) gives

\[
 d\psi_2=e_2,\qquad d\psi_3=e_3,                              \tag{2.4}
\]

and hence (1.1). If

\[
 |\eta_2|\|d\psi_2\|_\infty+
 |\eta_3|\|d\psi_3\|_\infty<1,                               \tag{2.5}
\]

then \(d\Phi_\eta\) is nonzero everywhere. Normalizing
\(n_\eta=d\Phi_\eta/|d\Phi_\eta|\), one has on \(B_0\)

\[
 \partial_{\eta_2}n_\eta|_{\eta=0}=e_2,\qquad
 \partial_{\eta_3}n_\eta|_{\eta=0}=e_3.                      \tag{2.6}
\]

Thus both tangent kernel rotations are realized inside one fixed integral
cohomology class. No arbitrary pointwise eigendirection field has been
selected.

## 3. Exact material transport

Let \(X(t,a)\) be the volume-preserving flow map of a common coarse velocity
which is independent of \(\eta\), and let \(F=D_aX\). Transport the scalar
phase:

\[
 \Phi_\eta(t,X(t,a))=\Phi_\eta(0,a).                          \tag{3.1}
\]

Then

\[
 D_t\Phi_\eta=0,\qquad
 \xi_\eta(t,X(t,a)):=\nabla_x\Phi_\eta
 =F(t,a)^{-T}\nabla_a\Phi_\eta(0,a).                          \tag{3.2}
\]

Differentiating gives the exact covector equation

\[
 D_t\xi_\eta+(\nabla U)^T\xi_\eta=0.                          \tag{3.3}
\]

The flow map is isotopic to the identity, so the integral class (2.3) is
preserved. If \(F\) and \(F^{-1}\) are bounded, nonvanishing of the initial
covector gives a uniform material phase gap. Put
\(n(t)=F^{-T}e_1/|F^{-T}e_1|\). On \(X(t,B_0)\), the exact normalized phase
tangents are

\[
 \partial_{\eta_j}n_\eta\big|_{\eta=0}
 =|F^{-T}e_1|^{-1}P_{n(t)^\perp}F^{-T}e_j,
 \qquad j=2,3.                                               \tag{3.4}
\]

They remain linearly independent because \(F^{-T}\) is invertible. If the
controlled construction instead uses an \(\eta\)-dependent coarse flow
\(X_\eta\), derivatives of \(F_\eta\) enter (3.4) and this rank conclusion
does not follow from the present calculation.

This solves phase integrability by construction for the finite-dimensional
core family. It does not show that an arbitrary spatially varying spectral
eigendirection is a material gradient.

## 4. The localized material-core stress submersion

At the initial core take

\[
 n_0=e_1,\qquad Q_0=\operatorname{diag}(0,a,b),\qquad a,b>0. \tag{4.1}
\]

Parameterize the positive covariance on \(n_\eta^\perp\) by its three
transverse block entries. The five first variations of \(Q\) are

\[
 E_{22},\quad E_{33},\quad E_{23}+E_{32},\quad
 -a(E_{12}+E_{21}),\quad -b(E_{13}+E_{31}),                  \tag{4.2}
\]

where the last two are supplied by the exact phase variations (2.6).
Adding the pressure coordinate to

\[
                       R=\rho I-Q                             \tag{4.3}
\]

negates all five covariance columns; the checker records the \(Q\)-column
sign convention, and the determinant magnitude is unchanged.  The resulting
six-by-six stress Jacobian has determinant magnitude \(ab\). Hence, for
uniform transverse gaps \(a,b\ge c_0>0\), the material-core map has a
uniform local inverse after its parameter neighborhood is fixed.

The phase-control rank in (3.4) survives at later material times. The
corresponding **stress** conclusion is not automatic: under a purely
kinematic vector pushforward an initial covariance transforms by a
congruence such as \(FQF^T\), whereas a Kelvin/WKB velocity amplitude follows
its own projected propagator. Covectors transform by \(F^{-T}\), but
covariances do not. A coupled calculation of that amplitude propagator and
its parameter derivative is still required for a control-dependent or
spatially localized coarse flow. C96 supplies the exact common-affine
calculation and finds determinant \(ab/|F^{-T}n|^8\).

This is a chart for the constant/affine stress variables in the active core,
not a right inverse for an arbitrary tensor field on the entire torus.
The cutoff region \(B\setminus B_0\) carries phase bending and can be placed
outside the principal fast-amplitude support at the level of the principal
symbol. This is not exact support separation: viscosity and pressure are
global, and their induced localization source belongs to the mandatory
boundary wake.

## 5. Stability under the longitudinal corrector

Exact slow modulation of a one-phase WKB velocity generally requires a
longitudinal amplitude of size \(O(K^{-1})\). Its covariance can add a small
third eigenvalue and perturb the five principal stress coordinates. It is
therefore unnecessarily rigid to demand that the **full** finite-\(K\)
covariance remain exactly rank two.

Let \(z\in\mathbb R^6\) denote the pressure, three transverse-block, and two
phase parameters. Write the complete stress map as (1.2). At the principal
base, \(A=D{\cal F}_0(z_0)\) is invertible by Section 4. Factor

\[
 D{\cal F}_K=A\big(I+A^{-1}DE_K\big).                         \tag{5.1}
\]

If (1.3) holds at the base point, the Neumann series gives (1.4) there. In
particular, an independently established bound

\[
 \|DE_K\|\le\frac{C}{K},\qquad
 K\ge2C\|A^{-1}\|                                             \tag{5.2}
\]

gives \(\vartheta\le1/2\) and

\[
                    \|(D{\cal F}_K)^{-1}\|
                    \le2\|A^{-1}\|.                           \tag{5.3}
\]

Thus the small third eigenvalue is not itself an algebraic obstruction to
pointwise rank. The constants used in the checker are hypotheses, not
derived corrector estimates. The real task is a uniform \(C^1\) estimate for
the exact divergence, pressure, viscous, and sideband correctors in the same
parameter norm, plus the neighborhood estimates stated after (1.4).

## 6. Remaining gates

This note closes only the fixed-cohomology/Frobenius objection for two
constant kernel-rotation controls in a localized frozen affine core. It
leaves:

1. the full spatial stress/covector equation and the control-dependent
   Kelvin covariance derivative in the transition collar;
2. a proof of the \(C^1\) corrector estimate (5.2);
3. the growing/infinite charged hierarchy forced by the C95
   two-harmonic/partner no-go;
4. the doubled-order clean endpoint or a recurrent nonzero carrier;
5. exact solenoidal localization and its work-carrying wake; and
6. the nonlinear carried-wake bound in C90.

No claim here is a Navier--Stokes singularity theorem.

## 7. Reproducibility

[The exact checker](../checks/local_material_phase_core_chart.py) verifies
the finite-dimensional cohomology/tangent ledger, the linearized covector
identity, the frozen stress Jacobian determinant, and one conditional
Neumann arithmetic example. Its inverse and corrector constants are explicit
assumptions; it does not verify any PDE corrector estimate or dynamic stress
submersion.
