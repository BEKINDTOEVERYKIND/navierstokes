# C167: reality pairing overdetermines a static pure-normal charge star

**Date:** 2026-08-06
**Status:** exact limiting-cone full-block obstruction to the independent-edge
scalar star and scalar-return aggregate within a static autonomous pure-normal
palette on an open angular/radial sheet; time ordering, cross-daughter metrics,
and the nonlinear one-cell stage remain open
**Checker:**
[checks/reality_paired_gate_bundle_rank_no_go_c167.py](../checks/reality_paired_gate_bundle_rank_no_go_c167.py)

## 0. Claim boundary

C164 left open whether assigning a different horizontal polarization to
each pure-normal gate height could repair the C161 terminal converter.  This
note rules out the independent-edge normalized scalar star, and a particular
scalar-return full-fiber relaxation, within the static autonomous first-step
pure-normal palette.  It uses the complete two-by-two source/daughter blocks
rather than a selected compression.

There are three exact obstructions.

1. For one nonzero reality pair of heights, a non-tangential source line (or
   an edgewise common positive \(2\)-by-\(2\) source metric with separate
   daughter metrics) can be shared by the two directed returns on an open
   angular arc only when the global horizontal gate polarization is circular.
   The shared line then has slope \(2\sigma i/y\), where
   \(\sigma=\pm1\) is the helicity.
2. Two distinct positive charge magnitudes cannot share that line or metric,
   regardless of their helicities.  Hence a reality-complete multi-daughter
   bundle cannot realize the independent-edge normalized scalar star.
3. Keeping both source polarizations does not produce a scalar-return finite
   tight-frame repair.  For every nonzero reality pair, the two diagonal
   eigenvalues of the pair-summed source return have a strict, sign-definite
   gap.  Positive mixtures cannot cancel it.  Moreover, at every nondegenerate
   source angle both oscillatory eigenfrequencies vary strictly with source
   radius on the safe band \(|y|<2\), so a fixed static pulse cannot
   synchronize an open radial interval.

The coefficient normalization is retained explicitly below; all three
obstructions are homogeneous and survive the C161 \(q^{-1/2}\) scaling.
The universal tangential line is not counted as a repair because C164--C165
prove that its full charge ladder is pointwise multiplication by a phase and
cannot focus.

This is not a no-go for noncommuting time colors, a pulse which deliberately
rotates both source polarizations, a nonlinear gate-depleting trajectory, an
approximately synchronized finite sheet, a metric coupling different
daughter or source modes, cancellation using cross-daughter paths, or a
non-normal gate geometry.  No finite-band, collar, pressure, wake, viscosity,
or unforced Navier--Stokes stage is asserted.

## 1. The full complex-polarization blocks

Use the C162--C164 limiting-cone coordinates

\[
 p=A\left(e_r+{e_z\over\sqrt3}\right),\qquad
 g={Ay\over\sqrt3}e_z,\qquad
 D_y=y^2+2y+4,
 \tag{1.1}
\]

where \(A>0\), and use the real orthonormal source and daughter frames

\[
 (e_\sigma,e_t),\qquad
 e_\sigma={e_r-\sqrt3e_z\over2},\qquad
 (e_\perp,e_t),\qquad
 e_\perp=-{1+y\over\sqrt{D_y}}e_r
          +{\sqrt3\over\sqrt{D_y}}e_z.
 \tag{1.2}
\]

Let one positive gate coefficient have unit horizontal polarization

\[
 E=u e_r+v e_t,\qquad u,v\in\mathbb C,
 \tag{1.3}
\]

and let reality supply \(\overline E\) at \(-g\).  Direct use of the
symmetric projected Euler symbol gives the dimensionless forward and
physical reverse blocks

\[
 \boxed{
 F_y(E)=
 \begin{pmatrix}
 -\dfrac{(4-y^2)u}{2\sqrt{D_y}}&0\\[2mm]
 -\dfrac{yv}{2}&u
 \end{pmatrix},\qquad
 R_y(\overline E)=
 \begin{pmatrix}
 -\dfrac{(y+2)\overline u}{\sqrt{D_y}}&0\\[2mm]
 -\dfrac{y\overline v}{\sqrt{D_y}}&\overline u
 \end{pmatrix}.}
 \tag{1.4}
\]

If the Fourier gate coefficient is \(\gamma\), the dimensional blocks are
\(A\gamma F_y\) and \(A\overline\gamma R_y\); their return is therefore
\(A^2|\gamma|^2M_y\), where

\[
 M_y=R_yF_y=
 \begin{pmatrix}
 |u|^2r(y)&0\\ q_y&|u|^2
 \end{pmatrix},
 \qquad
 r(y)={(4-y^2)(y+2)\over2D_y},
 \tag{1.5}
\]

and

\[
 q_y={y(4-y^2)\over2D_y}u\overline v
       -{y\over2}\overline u v.
 \tag{1.6}
\]

For a reality-complete family indexed by its positive heights, write

\[
 W=\sum_{a>0}\left(
   \gamma_aE_a e^{ig_az}
  +\overline{\gamma_aE_a}e^{-ig_az}\right),
 \qquad w_a=|\gamma_a|^2.
 \tag{1.7}
\]

With unit \(E_a\), the natural bare-coefficient scaling corresponding to
the C161 signed-daughter normalization is

\[
 2\sum_{a>0}w_a=1.
 \tag{1.8}
\]

Equal signed gate weights give \(w_a=1/q\) for \(q/2\) positive heights.
The Leray symbol still changes the actual edge norms, so (1.8) is not being
claimed as a physical star normalization.  None of the conclusions below
uses more than \(w_a>0\), and any later common rescaling cannot alter them.

## 2. The common line and common-metric equations

Assume \(u\ne0\), put

\[
 \rho={v\over u}=a+ib,
 \tag{2.1}
\]

and exclude the repeated or singular algebraic heights
\(y\in\{0,-4\}\).  Besides the universal tangential eigenline \(e_t\),
the return (1.5) has the non-tangential eigenline

\[
 e_\sigma+m_y(\rho)e_t,
 \qquad m_y(\rho)={q_y\over |u|^2r(y)-|u|^2}.
 \tag{2.2}
\]

The identities

\[
 r(y)-1=-{y^2(y+4)\over2D_y}
 \tag{2.3}
\]

and (1.6) give the exact slope formula

\[
 \boxed{
 m_y(\rho)
 ={2(y+1)\over y+4}\,\mathop{\rm Re}\rho
  +{2i\over y}\,\mathop{\rm Im}\rho.}
 \tag{2.4}
\]

The same number is forced by a weighted-unitary formulation.  Indeed, let

\[
 Q=\begin{pmatrix}q_{11}&z\\\overline z&q_{22}\end{pmatrix}>0
 \tag{2.5}
\]

be one \(2\)-by-\(2\) source-fiber metric shared edgewise by the daughter
blocks.  If a separate positive daughter metric \(P_y\) makes one edge block
weighted self-adjoint, then

\[
 QR_y=F_y^*P_y
 \quad\Longrightarrow\quad
 QM_y=M_y^*Q.
 \tag{2.6}
\]

For the triangular matrix (1.5), its off-diagonal equation is

\[
 z\bigl(|u|^2-|u|^2r(y)\bigr)=\overline{q_y}\,q_{22},
 \qquad
 \boxed{{z\over q_{22}}=-\overline{m_y(\rho)}.}
 \tag{2.7}
\]

Thus a common edgewise source-fiber metric for several edges requires exactly
the same slope condition as a common non-tangential source line.  This is only
a necessary condition for the independent-edge physical star, which is enough
for that no-go.  It says nothing about a global metric with off-diagonal
couplings between distinct source or daughter modes.

## 3. Reality forces circular polarization on an angular arc

The negative daughter uses \((-y,\overline E)\), hence
\(\rho\mapsto\overline\rho\).  Formula (2.4) gives

\[
 m_y(\rho)-m_{-y}(\overline\rho)
 ={12y\over16-y^2}\,\mathop{\rm Re}\rho.
 \tag{3.1}
\]

Therefore, for \(y\ne0,\pm4\), the two members of one reality pair share a
non-tangential source line or a common source metric only if

\[
 \mathop{\rm Re}{v\over u}=0.
 \tag{3.2}
\]

This condition is already rigid on an open angular arc.  Write one fixed
global horizontal polarization as

\[
 E=(X,Y,0),\qquad
 u=X\cos\phi+Y\sin\phi,qquad
 v=-X\sin\phi+Y\cos\phi.
 \tag{3.3}
\]

Where \(u\ne0\), (3.2) is equivalent to
\(\mathop{\rm Re}(v\overline u)=0\).  Direct expansion gives

\[
 \mathop{\rm Re}(v\overline u)
 ={ |Y|^2-|X|^2\over2}\sin(2\phi)
  +\mathop{\rm Re}(Y\overline X)\cos(2\phi).
 \tag{3.4}
\]

If this vanishes on any open angular interval, both coefficients vanish:

\[
 |X|=|Y|,qquad \mathop{\rm Re}(Y\overline X)=0.
 \tag{3.5}
\]

For nonzero \(E\), necessarily \(X\ne0\); the ratio \(Y/X\) has unit modulus
and zero real part.  Consequently

\[
 \boxed{Y=\sigma iX,qquad \sigma\in\{+1,-1\}.}
 \tag{3.6}
\]

Conversely, (3.6) gives

\[
 u=Xe^{\sigma i\phi},\qquad v=\sigma iu,qquad
 \rho=\sigma i,qquad
 \boxed{m_y={2\sigma i\over y}.}
 \tag{3.7}
\]

Thus circular polarization is not merely one convenient choice: it is the
only fixed horizontal polarization which makes a non-tangential line survive
one nonzero reality pair uniformly on an angular arc.

## 4. Distinct charge heights are incompatible

Let the positive physical gate heights be \(g_a>0\), so that at source
radius \(A\)

\[
 y_a={\sqrt3g_a\over A}.
 \tag{4.1}
\]

After Section 3, a common non-tangential line for two reality pairs would
require

\[
 {2\sigma_ai\over y_a}={2\sigma_bi\over y_b}
 \quad\Longleftrightarrow\quad
 {\sigma_a\over g_a}={\sigma_b\over g_b}.
 \tag{4.2}
\]

For distinct positive magnitudes \(g_a\ne g_b\) and
\(\sigma_a,\sigma_b\in\{\pm1\}\), (4.2) is impossible.  Opposite
helicities at one height are incompatible as well.  Duplicate coefficients
at the same Fourier height simply combine into one coefficient and do not
create distinct charge daughters.

Hence a reality-complete static pure-normal palette with at least two
distinct positive charge magnitudes has no common non-tangential physical
source line and no common edgewise \(2\)-by-\(2\) source metric on an open
angular arc.
The only common scalar line is the tangential line \(e_t\), whose exact full
ladder is the phase multiplier of C164 Section 5 and therefore has constant
pointwise modulus.

This kills the strict C161 independent-edge star in this class.  In that
star each natural daughter line must return to the same source line (and its
dark combinations must not create the second source polarization); merely
arranging cancellation after summing all daughters is not the same
operator.

## 5. A scalar-return two-polarization tight frame also fails

One might retain both source polarizations and ask for a scalar pair-summed
source return.  This proposal is also impossible with positive physical
intensities.

The two diagonal entries of \(M_y+M_{-y}\) are

\[
 |u|^2s(y),\qquad 2|u|^2,
 \tag{5.1}
\]

where the exact even rational function is

\[
 \boxed{
 s(y)=r(y)+r(-y)
 ={8(4-y^2)\over y^4+4y^2+16}.}
 \tag{5.2}
\]

For every nonzero real \(y\),

\[
 \boxed{
 2-s(y)={2y^2(y^2+8)\over y^4+4y^2+16}>0.}
 \tag{5.3}
\]

For an arbitrary finite palette, the source return through all signed first
daughters is lower triangular and has diagonal gap

\[
 \begin{aligned}
 (\mathcal M)_{22}-(\mathcal M)_{11}
 &=\sum_{a>0}w_a|u_a(\phi)|^2\,[2-s(y_a)]\\
 &>0
 \end{aligned}                                      \tag{5.4}
\]

at every source for which at least one radial projection \(u_a\) is
nonzero.  The gate polarizations can change the lower-left entry but cannot
change the sign in (5.4).  If every \(u_a\) vanishes, the block is
rank-degenerate rather than oscillatory.

Consequently no positive-weight finite polarization palette, including an
angular tight frame, makes this pair-summed two-polarization return a scalar
multiple of the identity.  In particular it cannot be the Euclidean
normalized skew star: if the physical reverse is the forward adjoint and the
forward map is an isometry up to one common scale, its return is a scalar
multiple of the identity.  This conclusion does not exclude a non-scalar
aggregate return or a metric coupling different daughter modes.  The
normalization (1.8) only rescales the positive sum in (5.4).

## 6. Exact radial synchronization is impossible for the static palette

In the oscillatory safe band \(0<|y_a|<2\), at an angle for which at least
one \(u_a(\phi)\ne0\), the two nonzero eigenfrequency squares of the
pair-summed static source return are

\[
 \begin{aligned}
 \Lambda_o^2(A,\phi)
 &=A^2\sum_{a>0}w_a|u_a(\phi)|^2
      s\left({c_a\over A}\right),\\
 \Lambda_t^2(A,\phi)
 &=2A^2\sum_{a>0}w_a|u_a(\phi)|^2,
 \qquad c_a=\sqrt3g_a>0.
 \end{aligned}                                      \tag{6.1}
\]

An angular tight frame can remove some \(\phi\)-dependence, but it cannot
remove the radial dependence.  Put \(x=c^2/A^2\in(0,4)\).  Then

\[
 A^2s(c/A)=c^2h(x),
 \qquad
 h(x)={8(4-x)\over x(x^2+4x+16)},                  \tag{6.2}
\]

and

\[
 h'(x)={16(x^3-4x^2-16x-32)
              \over[x(x^2+4x+16)]^2}<0,
 \tag{6.3}
\]

because

\[
 x^3-4x^2-16x-32=x^2(x-4)-16x-32<0
 \quad(0<x<4).                                    \tag{6.4}
\]

As \(A\) increases, \(x\) decreases, so every nonzero term in
\(\Lambda_o^2\) increases strictly.  The tangential value in (6.1) is
already proportional to \(A^2\).  If one fixed pulse time \(T>0\) depleted
the source branch on a radial interval, continuity would require
\(T\Lambda_o(A)\), or \(T\Lambda_t(A)\), to stay in the discrete set
\(\pi/2+\pi\mathbb Z\) on that interval.  It would therefore be constant,
contradicting strict monotonicity.  Thus neither branch has one common
quarter-period on an open radial interval under a fixed static gate palette.

## 7. Consequence for the one-cell target

Height-dependent static polarizations do not close C161 through either its
independent-edge scalar architecture or the scalar-return full-fiber
relaxation tested here.  The exact obstruction is not a missing choice of a
few polarization vectors: reality first forces each viable non-tangential
edge to a circular helicity, distinct heights then demand incompatible
edgewise source lines/metrics, and the scalar-return aggregate retains a
sign-definite rate anisotropy and radial detuning.

The same-geometry route still not ruled out is now narrower: a genuinely
time-ordered, noncommuting sequence may rotate the two source polarizations
between charge edges, or a full aggregate may use cross-daughter couplings or
a non-block-diagonal metric (or a nonlinear leading trajectory), while also
creating physical endpoint coherence with only the C146--C161 action and wake
budgets.  The independent-edge static palette, its scalar-return tight-frame
relaxation, and the universal tangential phase branch cannot do it.
