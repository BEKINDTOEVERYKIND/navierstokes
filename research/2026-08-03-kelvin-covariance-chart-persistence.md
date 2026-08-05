# Kelvin transport preserves the one-phase pressure-gauge chart

**Date:** 2026-08-03

**Status:** exact finite-dimensional affine/Kelvin theorem, including a
finite-time viscous two-harmonic variant. The localized, spatially varying,
finite-\(K\) transition remains open.

**Scope:** this note resolves the dynamic covariance question isolated in
the localized material-core note when the affine coarse flow is common to
all five physical chart parameters. It does not construct a finite-energy
periodic return cell or a Navier--Stokes singularity.

## 1. Outcome

The phase covector and velocity amplitude do not obey the same transport
law. If \(F\in SL(3,\mathbb R)\) is the deformation gradient, then

\[
 k=F^{-T}k_0                                                   \tag{1.1}
\]

is the transported phase covector, while a transverse Kelvin velocity
amplitude is

\[
 {\cal L}_F(k_0)a_0
 =\frac{F(k_0\times a_0)\times k}{|k|^2},
 \qquad k_0\cdot a_0=0.                                      \tag{1.2}
\]

It is the vorticity amplitude \(k\times a\), not the velocity amplitude,
which is pushed forward by \(F\).

Despite this distinction, the desired rank does persist. Formula (1.2) is
an isomorphism

\[
 {\cal L}_F(k_0):k_0^\perp\longrightarrow k^\perp.            \tag{1.3}
\]

Consequently the induced inviscid covariance map is a diffeomorphism of the
positive rank-two boundary:

\[
 {\cal T}_F(Q)
 ={\cal L}_F(n)Q{\cal L}_F(n)^T,
 \qquad n\in\ker Q,\quad |n|=1.                              \tag{1.4}
\]

Thus

\[
 (\rho,Q)\longmapsto \rho I-{\cal T}_F(Q)                    \tag{1.5}
\]

is again a full pressure-gauge chart for every finite nonsingular affine
deformation. No small-time assumption is needed for algebraic rank,
although its conditioning can deteriorate strongly.

At

\[
 Q_0=\operatorname{diag}(0,a,b),\qquad a,b>0,\qquad n=e_1,
                                                                    \tag{1.6}
\]

put \(s=|F^{-T}n|\). In the same six controls used by the frozen chart, the
finite-time stress Jacobian has exact magnitude

\[
                         |\det J_F|=\frac{ab}{s^8}.           \tag{1.7}
\]

The checker verifies (1.7) for a non-diagonal integer shear, where
\(s^2=3\) and the determinant is \(-1/135\).

Unequal viscous damping of the two labelled harmonics also preserves rank
at every finite time. The damping makes the output depend on the chosen
harmonic factorization of the covariance, but the resulting local fibre
map is still a smooth diffeomorphism.

## 2. Exact Kelvin amplitude propagator

Let

\[
 F'=S(t)F,\qquad F(0)=I,\qquad \operatorname{tr}S=0.          \tag{2.1}
\]

The Kelvin equations along the affine trajectory are

\[
 k'=-S^Tk,\qquad
 a'=-Sa+2k\,\frac{k\cdot Sa}{|k|^2},\qquad k\cdot a=0.       \tag{2.2}
\]

Define the vorticity amplitude

\[
                         b=k\times a.                         \tag{2.3}
\]

Using (2.2) and \(\operatorname{tr}S=0\) gives

\[
                              b'=Sb.                          \tag{2.4}
\]

Therefore

\[
 k(t)=F(t)^{-T}k_0,\qquad
 b(t)=F(t)(k_0\times a_0).                                  \tag{2.5}
\]

Because \(b(t)\cdot k(t)=0\), velocity is recovered from vorticity by

\[
 a(t)=\frac{b(t)\times k(t)}{|k(t)|^2},                      \tag{2.6}
\]

which is (1.2). Conversely, given \(a\in k^\perp\), first form
\(b=k\times a\), pull it back by \(F^{-1}\), and recover

\[
 a_0=\frac{F^{-1}b\times k_0}{|k_0|^2}.                     \tag{2.7}
\]

This proves (1.3). It also proves the composition law

\[
 {\cal L}_{F_2}(F_1^{-T}k_0){\cal L}_{F_1}(k_0)
 ={\cal L}_{F_2F_1}(k_0),                                   \tag{2.8}
\]

so the inverse is the Kelvin map for \(F^{-1}\).

The transverse determinant is explicit:

\[
 \left|\det\nolimits_\perp {\cal L}_F(k_0)\right|
 =\frac{|k_0|}{|F^{-T}k_0|}.                                \tag{2.9}
\]

For unit \(k_0\), the three factors in (1.2) have two-dimensional area
scales \(1\), \(s\), and \(s^{-2}\), respectively. Their product is
\(s^{-1}\), which proves (2.9).

## 3. Covariance is transported by a boundary-cone diffeomorphism

Let \({\cal B}_2^+\) be the manifold of positive semidefinite symmetric
three-by-three matrices of rank two. For \(Q\in{\cal B}_2^+\), choose either
unit generator \(n\) of its kernel and define (1.4). Replacing \(n\) by
\(-n\) does not change \({\cal L}_F(n)\), so this definition is
unambiguous.

The image is positive, has rank two, and has kernel parallel to
\(F^{-T}n\). Equations (2.7)--(2.8) show

\[
 {\cal T}_{F^{-1}}{\cal T}_F(Q)=Q.                           \tag{3.1}
\]

Hence \({\cal T}_F\) is a smooth diffeomorphism of
\({\cal B}_2^+\). Composing it with the static spectral chart

\[
 (\rho,Q)\longmapsto \rho I-Q                               \tag{3.2}
\]

proves (1.5). In particular, all five physical covariance directions and
the scalar pressure gauge survive finite affine evolution. This is
stronger than a small-pulse Neumann argument.

The exact conditioning factor (1.7) follows from four elementary
Jacobians. The projective kernel map

\[
 n\longmapsto\frac{F^{-T}n}{|F^{-T}n|}                      \tag{3.3}
\]

has spherical Jacobian \(s^{-3}\). Congruence by the two-dimensional
Kelvin map has determinant \(s^{-3}\) on \(\operatorname{Sym}_2\). The
target transverse covariance has determinant \(ab\,s^{-2}\), which is the
Jacobian density of its two kernel-rotation columns in the ambient
symmetric tensor space. Multiplication gives

\[
          ab\,s^{-2}\,s^{-3}\,s^{-3}=ab\,s^{-8}.             \tag{3.4}
\]

Thus there is no finite-time rank loss, but a large compressed wave number
\(s\) can make the inverse expensive. A fixed compression factor per
cascade stage is compatible with a uniform chart; an unbounded
within-stage compression requires this \(s^8\) loss in the derivative
budget.

For the explicit C89 strain
\[
 S=\operatorname{diag}(-1,-5/4,9/4),\qquad n=e_1,
\]
one has \(s(t)=e^t\), so
\[
                         |\det J_{F(t)}|=ab\,e^{-8t}.        \tag{3.5}
\]
If controls are parameterized at the beginning of a gain window of length
\(G_j\), the raw inverse can therefore lose \(e^{8G_j}\).  Algebraic rank
persistence is not a tame endpoint estimate.  A viable construction must
either launch/reparameterize the controls on a bounded terminal pulse,
renormalize the control coordinates with a proved physical budget, or show
that the full dynamic endpoint map cancels this loss.

## 4. Why neither \(FQF^T\) nor \(F^{-T}QF^{-1}\) is the answer

Take

\[
 F=
 \begin{pmatrix}
 1&1&0\\
 0&1&1\\
 0&0&1
 \end{pmatrix},
 \qquad k_0=e_1.                                             \tag{4.1}
\]

Then

\[
 k=F^{-T}e_1=(1,-1,1),                                      \tag{4.2}
\]

but (1.2) gives

\[
 {\cal L}_F(e_1)e_2=\frac13(2,1,-1),\qquad
 {\cal L}_F(e_1)e_3=\frac13(-1,1,2).                        \tag{4.3}
\]

The first vector is neither \(Fe_2=(1,1,0)\) nor
\(F^{-T}e_2=(0,1,-1)\). For
\(Q_0=\operatorname{diag}(0,1,1)\), the true output is

\[
 {\cal T}_F(Q_0)
 =\frac19
 \begin{pmatrix}
 5&1&-4\\
 1&2&1\\
 -4&1&5
 \end{pmatrix},                                             \tag{4.4}
\]

whose kernel is (4.2). This supplies an exact counterexample to both naive
congruence rules while confirming the Kelvin boundary-cone theorem.

## 5. Two distinct harmonics with viscosity

Fix a local oriented transverse frame \(E(n)\), write the covariance block
as \(B>0\), and use the labelled harmonic amplitudes given by the columns
of

\[
                         E(n)B^{1/2}.                         \tag{5.1}
\]

For harmonic numbers \(m_1=1\) and \(m_2=2\), put

\[
 I(n,t)=\int_0^t|F(\tau)^{-T}n|^2\,d\tau,\qquad
 d_j=\exp\{-\nu m_j^2 I(n,t)\}.                              \tag{5.2}
\]

The viscous amplitude is \(d_j{\cal L}_F(n)c_j\). Before the common Kelvin
congruence, the transverse covariance therefore undergoes

\[
 \Psi_D(B)=B^{1/2}D^2B^{1/2},\qquad
 D=\operatorname{diag}(d_1,d_2).                             \tag{5.3}
\]

This is a diffeomorphism of the positive two-by-two cone. Indeed, if
\(X=B^{1/2}\) and \(Y=XD^2X\), then

\[
 D Y D=(DXD)^2,\qquad
 X=D^{-1}(DYD)^{1/2}D^{-1},\qquad B=X^2.                    \tag{5.4}
\]

All operations are smooth while \(d_1,d_2>0\), which holds at every finite
time. Dependence of \(D\) on \(n\) adds only an off-diagonal block to the
base--fibre derivative and cannot remove rank.

For an exact local determinant, take

\[
 B_0=\operatorname{diag}(a,b),\qquad
 \delta_j=d_j^2.
\]

The derivative of (5.3) in the \(B_{23}\) direction is multiplied by

\[
 c_D=\frac{\delta_1\sqrt a+\delta_2\sqrt b}
           {\sqrt a+\sqrt b}>0.                             \tag{5.5}
\]

The full six-dimensional stress determinant becomes

\[
 |\det J_{F,\nu}|
 =\frac{ab}{s^8}(\delta_1\delta_2)^2c_D>0.                  \tag{5.6}
\]

It reduces to (1.7) at zero viscosity. Conditioning can still degenerate
as a damping factor tends to zero, but no finite-time algebraic rank
obstruction occurs.

The unequal-damping statement is tied to the labelled factorization (5.1).
Inviscid transport is intrinsically a map of \(Q\); viscous transport of
distinct harmonics remembers which square-root column was assigned to
which harmonic.

## 6. Exact same-phase affine dynamics

The preceding calculation is not merely a formal covariance rule. If the
affine base \(U=S(t)x\) is itself an exact Euler or Navier--Stokes solution,
then a profile

\[
 w(t,x)=\sqrt2\sum_{j=1}^2
 d_j(t){\cal L}_{F(t)}(k_0)c_j
 \cos\!\big(m_j k(t)\cdot x\big)                             \tag{6.1}
\]

is an exact global affine Kelvin perturbation with the corresponding
linear pressure modes. Every amplitude is perpendicular to the common
phase covector, and hence

\[
 (w\cdot\nabla)w
 =(w\cdot k)\,\partial_\theta w=0.                           \tag{6.2}
\]

Thus the two covariance harmonics do not create a nonlinear same-phase
interaction. Viscosity acts diagonally through (5.2). The obstruction is
localization, not the affine interior evolution.

## 7. Exact claim boundary

This note closes the specific dynamic-rank qualification in the frozen
localized material core under a common parameter-independent affine flow.
It leaves five load-bearing issues.

1. A global affine field on \(\mathbb R^3\) has infinite energy and is not a
   periodic velocity. Cutting it off creates collar, pressure, and
   divergence sources.
2. If the coarse flow \(F_\eta\) depends on the phase controls, derivatives
   of \(F_\eta\) enter the chart. Small \(C^1\) feedback preserves rank by
   perturbation, but no such feedback estimate is proved here.
3. The scalar \(\rho\) is a pressure-gauge coordinate, not an independently
   generated velocity covariance. A localized transition must realize its
   pressure consistently with the full equation.
4. Slow modulation, finite-\(K\) solenoidal corrections, sidebands, and the
   global zero-charge wake are absent from the affine theorem.
5. The covariance still loses its kernel-rotation columns when it is
   forced to vanish at a clean endpoint. The doubled-order endpoint gate is
   unchanged.

Accordingly, the correct finite-dimensional conclusion is positive:
material covectors and Kelvin velocity amplitudes evolve differently, but
their coupled evolution preserves the full one-phase pressure-gauge stress
chart exactly at finite affine time. The remaining obstruction is the
localized nonlinear transition and its wake, not affine covariance rank.
