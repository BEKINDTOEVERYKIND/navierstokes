# C163: a dual-helicity charge gate has a negative compressed product

**Date:** 2026-08-05
**Status:** exact limiting-cone selected-compression sign classification,
signed-height interval, leakage formula, and equal-radius synchronization
obstruction; a physical normalized charge star, finite-band errors, and the
one-cell stage remain open
**Checker:**
[checks/paired_polarization_charge_star_c163.py](../checks/paired_polarization_charge_star_c163.py)

## 0. Claim boundary

C162 found a non-real forward/reverse product for one circular polarization.
That phase is not an intrinsic obstruction of the pure-normal charge geometry.
An equal-magnitude coherent superposition of the two opposite circular
helicities makes the **selected** forward/reverse product strictly negative at
every nonnegative gate/source height ratio, and also for \(-1\le y<0\), away
from one angular rank-zero line.

The selected two-coordinate compression is nevertheless not the C161
Hilbert-space star.  Its two off-diagonal entries have unequal magnitudes.  A
source-dependent diagonal similarity makes that compression skew-Hermitian,
but the full Leray symbol does not preserve the selected source/daughter
plane: its reverse leg generally leaks into the second source polarization.
Consequently the weighted norm below is conserved only by the compressed
matrix, not by the physical Fourier system.  Even at the compressed level, an
equal-radius angular continuum cannot be synchronized by any nonzero set of
shared gate intensities.

All calculations below are for the limiting cone symbol.  No finite-band
packet, repeated-charge collision estimate, localization, pressure estimate,
or nonlinear Navier--Stokes stage is asserted.

## 1. Cone coordinates and the selected source polarization

Write

\[
 p=A\left(e_r+{e_z\over\sqrt3}\right),\qquad
 e_r=(\cos\phi,\sin\phi,0),\qquad
 e_t=(-\sin\phi,\cos\phi,0),                         \tag{1.1}
\]

with \(A>0\).  Put \(c=\cos\phi\), \(s=\sin\phi\), and

\[
 \alpha={c+s\over2\sqrt2},\qquad
 \beta={c-s\over\sqrt2}.                            \tag{1.2}
\]

The C162 selected unit source amplitude is

\[
 a=\alpha e_r+\beta e_t-\sqrt3\alpha e_z,
 \qquad 4\alpha^2+\beta^2=1,
 \qquad a\cdot p=0.                                  \tag{1.3}
\]

Let the pure-normal gate wavevector be

\[
 g={Ay\over\sqrt3}e_z,
 \qquad y={\sqrt3\,g_z\over A}\in\mathbb R.         \tag{1.4}
\]

Thus the tuned C162 height is \(y=2\).  For a horizontal complex gate
polarization \(E=u e_r+v e_t\), the symmetric projected Euler symbol is

\[
 {1\over A}\mathcal S(p,a;g,E)
 =-\alpha u{4-y^2\over\sqrt D}\,e_\perp
   +(\beta u-\alpha yv)e_t,                           \tag{1.5}
\]

where

\[
 D=y^2+2y+4,
 \qquad
 e_\perp=-{1+y\over\sqrt D}e_r+{\sqrt3\over\sqrt D}e_z
 \perp p+g.                                           \tag{1.6}
\]

At \(y=2\), the first component vanishes.  The forward scalar is then

\[
 f=\beta u-2\alpha v={E_x-E_y\over\sqrt2},            \tag{1.7}
\]

which is independent of the source angle.  The reverse coefficient from
\(e_t\) back to \(a\), using the reality companion \(\overline E\) at
\(-g\), is

\[
 r=\beta\overline u.                                  \tag{1.8}
\]

Consequently the tuned selected two-coordinate compression is

\[
 L_E=-iA\begin{pmatrix}0&\beta\overline u\\f&0\end{pmatrix},
 \qquad
 (L_E)_{12}(L_E)_{21}=-A^2\beta\overline u
   (\beta u-2\alpha v).                              \tag{1.9}
\]

This also gives the complete classification of this selected compression.
For
\(\alpha\beta\ne0\), the product is negative real precisely when
\(u\ne0\) and

\[
 {v\over u}\in\mathbb R,
 \qquad
 \beta\left(\beta-2\alpha{v\over u}\right)>0.        \tag{1.10}
\]

The selected compression is skew-Hermitian in its unweighted two-coordinate
metric precisely
when

\[
 \beta\overline u=\overline f
 \quad\Longleftrightarrow\quad \alpha v=0.           \tag{1.11}
\]

Thus the selected compression for a nondegenerate single source is repaired
by choosing the gate radial to it, \(v=0\).  This choice depends on \(\phi\),
and it does not show that the full two-polarization symbol has an invariant
skew block.  At \(\beta=0\), equivalently
\(\phi=\pi/4\pmod\pi\), no nonzero tuned skew edge exists in the selected
compression.

## 2. One shared equal-magnitude dual-helicity polarization

Let

\[
 E_0={e_x-e_y\over\sqrt2}.
 \tag{2.1}
\]

In terms of the circular polarizations

\[
 E_+={e_x+i e_y\over\sqrt2},\qquad
 E_-={e_x-i e_y\over\sqrt2},                         \tag{2.2}
\]

this is the equal-magnitude coherent pair

\[
 E_0={1+i\over2}E_+ + {1-i\over2}E_-.                \tag{2.3}
\]

Put

\[
 h=c-s,\qquad \ell=c+s,
 \qquad h^2+\ell^2=2.                                \tag{2.4}
\]

Then \(u=h/\sqrt2=\beta\), \(v=-\ell/\sqrt2=-2\alpha\).
For \(F=A^{-1}\mathcal S(p,a;g,E_0)\), (1.5) becomes

\[
 F_\perp=-{\ell h(4-y^2)\over4\sqrt D},
 \qquad
 F_t={2h^2+y\ell^2\over4}.                           \tag{2.5}
\]

The reverse row on the same orthonormal daughter frame is

\[
 r_\perp=-{\ell h\over\sqrt D},
 \qquad
 r_t={h^2\over2}.                                    \tag{2.6}
\]

Define the forward squared norm and the directed return by

\[
 N=|F|^2,
 \qquad R=r_\perp F_\perp+r_tF_t.                    \tag{2.7}
\]

Direct simplification gives

\[
 N={\ell^2h^2(4-y^2)^2\over16D}
    +{(2h^2+y\ell^2)^2\over16},                      \tag{2.8}
\]

and

\[
 \boxed{
 R={h^2\over2}\left[
 1+{\ell^2y^2(y-2)\over4D}
 \right].}                                           \tag{2.9}
\]

The Fourier ODE off-diagonal product is therefore

\[
 \boxed{P=-A^2R.}                                    \tag{2.10}
\]

In the notation \(x=g_z/A=y/\sqrt3\) and
\(D_x=1+(x+1/\sqrt3)^2=D/3\), (2.9) is equivalently

\[
 R={h^2\over2}\left[
 1+{\sqrt3x^2\ell^2(x-2/\sqrt3)\over4D_x}
 \right].                                            \tag{2.11}
\]

Equation (2.9) is strictly positive for every \(y\ge0\) when \(h\ne0\).
For \(y\ge2\), the bracket is at least one.  For \(0\le y\le2\), use

\[
 y^2(2-y)\le{32\over27},\qquad D\ge4,
 \qquad \ell^2\le2                                  \tag{2.12}
\]

to obtain

\[
 1+{\ell^2y^2(y-2)\over4D}\ge {23\over27},
 \qquad
 R\ge {23\over54}h^2.                               \tag{2.13}
\]

The reality-symmetric C161 shifts require both signs of \(g_z\).  On the
clean signed interval \(-1\le y\le0\),

\[
 y^2(2-y)\le3,\qquad D=(y+1)^2+3\ge3,
\]

and hence

\[
 1+{\ell^2y^2(y-2)\over4D}\ge {1\over2},
 \qquad R\ge {h^2\over4}.                           \tag{2.14}
\]

Thus the negative selected ODE product is certified on
\(-1\le y<\infty\), in particular for reality-paired shifts with
\(|y|\le1\).  No sign claim is made here for \(y<-1\).

Thus unequal source radii do not restore the C162 complex phase for a
positive gate height: \(y\) varies with \(A\), but the selected product
remains negative.  A signed reality pair has the same conclusion whenever
its negative member obeys \(y\ge-1\).  The angular line \(h=0\) remains an
exact rank defect.

## 3. The exact compressed weighted rotation and its leakage

For \(h\ne0\), let \(b=F/\sqrt N\) be the unit forward-bright daughter and
let \(P_{a,b}\) denote orthogonal projection onto the selected source and
daughter lines.  The two-coordinate compression \(P_{a,b}LP_{a,b}\) is

\[
 L_0=-iA
 \begin{pmatrix}
 0&R/\sqrt N\\
 \sqrt N&0
 \end{pmatrix}.                                      \tag{3.1}
\]

With

\[
 S=\operatorname {diag}\!\left(\sqrt{R/N},1\right),
 \tag{3.2}
\]

one has the exact similarity

\[
 S^{-1}L_0S=-iA\sqrt R
 \begin{pmatrix}0&1\\1&0\end{pmatrix}.              \tag{3.3}
\]

Equivalently, this compressed matrix is skew-Hermitian for the weighted
metric

\[
 W=S^{-*}S^{-1}=\operatorname {diag}(N/R,1),
 \qquad
 L_0^*W+WL_0=0.                                      \tag{3.4}
\]

At the tuned height \(y=2\), the formulas reduce to

\[
 N=1,qquad R={h^2\over2}=\beta^2,qquad
 L_0=-iA\begin{pmatrix}0&\beta^2\\1&0\end{pmatrix}.
 \tag{3.5}
\]

The similarity \(\operatorname {diag}(|\beta|,1)\) converts (3.5) to
\(-iA|\beta|\sigma_x\), but the compressed unweighted norm is replaced by

\[
 {|z|^2\over\beta^2}+|w|^2.                          \tag{3.6}
\]

In the isolated compressed ODE, a unit source quarter-rotates at
\(\pi/(2A|\beta|)\) and produces daughter coordinate \(1/|\beta|\).
This is not a statement about the full Leray ODE.  On any compact set

\[
 |h|\ge h_0>0,\qquad -1\le y\le Y,                  \tag{3.7}
\]

the compressed weighted metric is uniformly equivalent to its unweighted
two-coordinate norm.  No uniform condition number is claimed as
\(h\to0\) or \(y\to\infty\).

The failure of invariance is explicit.  Put

\[
 e_s=-{1\over2}e_r+{\sqrt3\over2}e_z,\qquad
 a_\perp=\beta e_s+2\alpha e_t.                     \tag{3.7a}
\]

Then \(a_\perp\) is the unit source polarization orthogonal to \(a\).  If
\(F=A^{-1}\mathcal S(p,a;g,E_0)\), direct projection gives

\[
 \boxed{
 \left\langle a_\perp,
 \mathcal S(p+g,AF;-g,E_0)\right\rangle
 =A^2\alpha\beta\,{y^2[(2-\beta^2)y+2(1+\beta^2)]\over D}.} \tag{3.7b}
\]

At \(y=2\) this is \(2A^2\alpha\beta\), which is generically nonzero.
Therefore (3.1)--(3.6) are exact identities for a compression, not a closed
physical subsystem or a conserved weighted norm for the full symbol.

For \(q\) charge shifts, the same mismatch already survives in the selected
collective compression.  If the gate weights are \(\gamma_j\), the
compressed forward-bright norm and directed product are proportional to

\[
 \sum_j|\gamma_j|^2N_j,
 \qquad
 \sum_j|\gamma_j|^2R_j,                              \tag{3.8}
\]

respectively.  The compressed collective block is unweighted-skew only if
these two sums agree.  Moreover, \(N_j,R_j\) depend on the source angle and
radius.

There is no additional first-step cycle obstruction in the **compressed
graph**.  In the C161 support chart every daughter \((i,j)\) retains its
source label \(i\), so the source--first-daughter graph is a disjoint union
of stars.  Since every compressed edge product is negative real, one may fix
a weight on each source and choose the phase and positive diagonal weight of
every daughter independently.  This makes each selected
source--forward-bright compression diagonally skew-Hermitian.  If

\[
 P_{ij}=-A_i^2R(\phi_i,y_{ij}),
 \qquad y_{ij}={\sqrt3(g_j)_z\over A_i},              \tag{3.9}
\]

the resulting compressed weighted bright frequency at source \(i\) is

\[
 \boxed{\Lambda_i^2
 =\sum_j|\gamma_j|^2(-P_{ij})
 =A_i^2\sum_j|\gamma_j|^2R(\phi_i,y_{ij}).}           \tag{3.10}
\]

Thus edgewise phase consistency is solved only after compression.  The full
symbol has the additional source-polarization vertices exposed by (3.7b).

There is also an exact synchronization obstruction already inside the
compression.  Consider equal-radius sources \(A_i=A\) whose angles fill an
interval, so \(y_{ij}=y_j\) is independent of \(\phi_i\).  Define

\[
 W_0=\sum_j|\gamma_j|^2>0,\qquad
 C_0=\sum_j|\gamma_j|^2
 {y_j^2(y_j-2)\over4(y_j^2+2y_j+4)}.                \tag{3.11}
\]

Since \(h^2=1-\sin2\phi\) and \(\ell^2=1+\sin2\phi\), (2.9) gives, with
\(x=\sin2\phi\),

\[
 \boxed{{\Lambda(\phi)^2\over A^2}
 ={1\over2}\big[(W_0+C_0)-W_0x-C_0x^2\big].}       \tag{3.12}
\]

The coefficient of \(x\) is the rigid value \(-W_0/2\).  Thus no nonzero
set of shared intensities makes the compressed rate constant on any angular
interval, regardless of the signed gate heights.  This does not rule out a
carefully selected finite set of source angles, approximate synchronization,
time-dependent polarizations, or dynamics through the leaked source line.
Later charge steps also create collisions and cycles, so the compressed
forest argument does not apply to the full ladder.

C163 therefore supplies a uniformly signed, diagonally weighted
**compression** on a restricted cone, not a physical normalized \(q\)-way
star.

## 4. Why two time colors do not remove the selected angular rank defect

At \(y=2\), the two circular colors have

\[
 f_\pm={e^{\pm i\phi}\over\sqrt2}(\beta\mp2i\alpha),
 \qquad
 r_\pm={\beta e^{\mp i\phi}\over\sqrt2}.             \tag{4.1}
\]

If the colors are made orthogonal in time/Floquet label so that cross
terms average away, with nonnegative intensities \(w_\pm\), their directed
product and forward norm are

\[
 R_{\rm pair}={\beta^2\over2}(w_++w_-)
 +i\alpha\beta(w_--w_+),
 \qquad
 N_{\rm pair}={w_++w_-\over2}.                       \tag{4.2}
\]

Equal intensities remove the imaginary part, but

\[
 {R_{\rm pair}\over N_{\rm pair}}=\beta^2.           \tag{4.3}
\]

Thus the balanced two-color selected compression gives the same weighted
defect as (3.5).  It is unweighted-skew only on the exceptional line
\(\beta^2=1\).

At \(\beta=0\) the tuned reverse coefficient into the selected source line
vanishes for **every** horizontal gate polarization at every time.  If the
orthogonal source polarization is discarded, the resulting time-dependent
two-coordinate compression has the triangular form

\[
 \dot z=0,
 \qquad
 \dot w=-iA f(t)z,
 \qquad
 U(T)=\begin{pmatrix}
 1&0\\-iA\int_0^T f(t)\,dt&1
 \end{pmatrix}.                                      \tag{4.4}
\]

If this compressed propagator is unitary, its lower off-diagonal entry is
zero and the propagator is the identity.  Hence no two-color modulation of
this selected compression gives a nontrivial unweighted-skew transfer
uniformly over the full angular circle.  This is not a no-go for the full
two-source-polarization system: the reverse vector can enter the discarded
source line, and that larger block remains open.

## 5. Remaining one-cell obligations

C163 repairs the explicit C162 complex phase of the selected compression on
every cone sector separated from \(c=s\), including signed heights
\(-1\le y<0\).  It does **not** close the terminal splitter.  The exact
remaining obligations are:

1. choose a \(q^2\)-source microlocal sheet separated from \(c=s\) while
   retaining the C154 shear and C147 cardinality requirements;
2. overcome the rigid equal-radius angular synchronization obstruction, or
   prove a factorial-safe approximate transfer on a specially selected
   finite source sheet;
3. close the leaked source polarization and construct a physical
   \(L^2\)-normalized collective star; the compressed weighted identity is
   not a conservation law for the full system;
4. control daughter--gate repeated charges, finite radial/angular bands,
   collars, pressure, retained wakes, C125 relative residuals, and BAFL.

No one-cell Navier--Stokes stage or Millennium conclusion is claimed.
