# C164: the full pure-normal gate fiber has one unitary branch and one signed branch

**Date:** 2026-08-06
**Status:** exact frozen limiting-cone linear symbol, full-polarization
classification, and tangential-fiber no-go; a time-dependent
polarization repair and the nonlinear one-cell stage remain open
**Checker:**
[checks/full_fiber_pure_normal_gate_c164.py](../checks/full_fiber_pure_normal_gate_c164.py)

## 0. Claim boundary

C163 compressed the pure-normal gate interaction to the particular source
line inherited from C149.  That compression hides a simpler exact branch.
This note computes the complete two-polarization source-to-daughter and
daughter-to-source blocks at an arbitrary **signed** normal height.

There are three conclusions.

1. The round trip has the exact eigenvalues

   \[
   \beta^2,
   \qquad
   \beta^2{(4-y^2)(y+2)\over2(y^2+2y+4)}.
   \]

   The first is a tangential branch whose forward and reverse coefficients
   are both exactly \(\beta\), independent of height.  When \(\beta\ne0\),
   the second is elliptic for \(y<2\) except for its zero at \(y=-2\) and
   its collision with the tangential branch at \(y=-4\), which is a genuine
   Jordan defect when \(v\ne0\), and is
   hyperbolic for \(y>2\).  The two full blocks are positively symmetrizable for
   \(|y|<2\), away from the already known angular defect \(\beta=0\).
2. For a generic source angle \(uv\ne0\), the oblique invariant line depends
   injectively on \(y\).  A nonzero reality pair \(\{y,-y\}\) therefore has
   only the tangential source line in common.  This is true throughout the
   safe band \(0<|y|\leq1\).
3. The apparent tangential repair cannot supply C161's missing pointwise
   gain.  On the whole normal-charge fiber, an arbitrary real pure-normal
   gate bundle acts on that polarization by a unit-modulus phase generated
   by a real scalar potential.
   It preserves the physical pointwise modulus exactly.  Starting from one
   charge mode, a large Fourier coefficient spread can be generated, but
   its phases cancel to constant physical amplitude.

The last statement is an exact no-go for realizing C161 through the
**fixed tangential source polarization** using any real pure-normal gate
bundle.  It
does not rule out noncommuting time colors, a pulse that actively rotates
source polarizations, a genuinely two-source-coordinate terminal map, or a
nonlinear leading trajectory.  Nothing here proves a finite-band,
localized, viscous, or nonlinear Navier--Stokes stage.

## 1. The full forward and reverse matrices

Use the C162--C163 cone coordinates

\[
 \begin{aligned}
 e_r&=(\cos\phi,\sin\phi,0),
 &e_t&=(-\sin\phi,\cos\phi,0),\\
 p&=A\left(e_r+{e_z\over\sqrt3}\right),
 &g&={Ay\over\sqrt3}e_z,
 \end{aligned}                                      \tag{1.1}
\]

where \(A>0\) and now \(y\in\mathbb R\) is signed.  Put

\[
 D=y^2+2y+4=(y+1)^2+3.                              \tag{1.2}
\]

Orthonormal bases of the source and daughter divergence-free planes are

\[
 e_\sigma={e_r-\sqrt3e_z\over2},\qquad e_t,
 \qquad
 e_\perp=-{1+y\over\sqrt D}e_r+{\sqrt3\over\sqrt D}e_z,
 \qquad e_t.                                        \tag{1.3}
\]

The equal-magnitude dual-helicity C163 gate is

\[
 E_0={e_x-e_y\over\sqrt2}=u e_r+v e_t,
 \qquad
 u={\cos\phi-\sin\phi\over\sqrt2}=\beta,
 \qquad
 v=-{\cos\phi+\sin\phi\over\sqrt2}=-2\alpha,
 \qquad u^2+v^2=1.                                  \tag{1.4}
\]

For the symmetric projected Euler symbol

\[
 \mathcal S(k,a;q,b)
 =P_{k+q}\{(a\cdot q)b+(b\cdot k)a\},               \tag{1.5}
\]

define the dimensionless operators

\[
 F a=A^{-1}\mathcal S(p,a;g,E_0),
 \qquad
 R b=A^{-1}\mathcal S(p+g,b;-g,E_0).                \tag{1.6}
\]

Rows and columns in the bases (1.3) give the complete matrices

\[
 \boxed{
 F=\begin{pmatrix}
 -\dfrac{u(4-y^2)}{2\sqrt D}&0\\[2mm]
 -\dfrac{yv}{2}&u
 \end{pmatrix},
 \qquad
 R=\begin{pmatrix}
 -\dfrac{u(y+2)}{\sqrt D}&0\\[2mm]
 -\dfrac{yv}{\sqrt D}&u
 \end{pmatrix}.}                                    \tag{1.7}
\]

This uses the physical reality companion \((-g,E_0)\); \(E_0\) is real.
The symbol \(e_\sigma\) is chosen to avoid conflict with the oppositely
oriented auxiliary \(e_s\) in C163 Section 3.  For example, on
\(e_\sigma\) the unprojected forward vector divided by \(A\) is
\(-yE_0/2+ue_\sigma\).  Dotting it with \(e_\perp\) gives the first entry of
\(F\).  Dotting the reverse vector with \(e_\sigma\) gives the first entry of
\(R\).  The remaining entries follow without pressure because \(e_t\) is
orthogonal to every wavevector in this normal-charge fiber.

The full frozen Fourier ODE on source and daughter coordinates is

\[
 {d\over dt}\binom z w
 =-iA
 \underbrace{\begin{pmatrix}0&R\\F&0\end{pmatrix}}_{H_y}
 \binom z w.                                        \tag{1.8}
\]

## 2. Exact spectrum and invariant polarizations

Multiplication of (1.7) gives

\[
 RF=\begin{pmatrix}
 \lambda_o&0\\[1mm]
 -\dfrac{uvy^2(y+1)}D&u^2
 \end{pmatrix},                                     \tag{2.1}
\]

where

\[
 \boxed{
 \lambda_t=u^2,
 \qquad
 \lambda_o
 =u^2{(4-y^2)(y+2)\over2D}
 =u^2{(2-y)(y+2)^2\over2D}.}                        \tag{2.2}
\]

Thus the numerical eigenvalue conjecture is exact.  The target round trip
\(FR\) has the same diagonal and lower entry

\[
 (FR)_{21}={uvy^2\over2\sqrt D}.                    \tag{2.3}
\]

The tangential vectors themselves form an invariant branch:

\[
 Fe_t=ue_t,
 \qquad
 Re_t=ue_t.                                         \tag{2.4}
\]

For \(u\ne0\) and \(y\ne-4\), set

\[
 \xi_y={2v(y+1)\over u(y+4)},
 \qquad
 \eta_y=-{v\sqrt D\over u(y+4)},                   \tag{2.5}
\]

and define

\[
 \sigma_y=e_\sigma+\xi_y e_t,
 \qquad
 d_y=e_\perp+\eta_y e_t.                            \tag{2.6}
\]

Then the second branch is not merely a round-trip eigenline: forward and
reverse maps close separately,

\[
 \boxed{
 F\sigma_y=f_y d_y,
 \qquad
 Rd_y=r_y\sigma_y,}                                 \tag{2.7}
\]

with

\[
 f_y=-{u(4-y^2)\over2\sqrt D},
 \qquad
 r_y=-{u(y+2)\over\sqrt D},
 \qquad f_yr_y=\lambda_o.                           \tag{2.8}
\]

Equivalently, with the source and daughter shear matrices

\[
 S_y=\begin{pmatrix}1&0\\\xi_y&1\end{pmatrix},
 \qquad
 T_y=\begin{pmatrix}1&0\\\eta_y&1\end{pmatrix},   \tag{2.9}
\]

one has the simultaneous bi-diagonalization

\[
 T_y^{-1}FS_y=\operatorname {diag}(f_y,u),
 \qquad
 S_y^{-1}RT_y=\operatorname {diag}(r_y,u).           \tag{2.10}
\]

At \(y=0\), \(RF=u^2I\), so its eigenspace is not unique, but
(2.7) remains a valid branch identity.  The other repeated point is
load-bearing: at \(y=-4\),

\[
 RF=\begin{pmatrix}u^2&0\\4uv&u^2\end{pmatrix}.    \tag{2.11}
\]

For \(uv\ne0\) it is a genuine Jordan block.  The frozen ODE has secular
growth multiplying its oscillation and admits no positive symmetrizer.
When \(v=0\), the lower entry vanishes and this exceptional Jordan defect
is absent.

## 3. Weighted ellipticity, endpoint ranks, and hyperbolicity

For \(u\ne0\) and \(|y|<2\), both \(f_y\) and \(r_y\) are nonzero and have
the same sign, and

\[
 {f_y\over r_y}={2-y\over2}=:w_y>0.                 \tag{3.1}
\]

In the sheared coordinates (2.9), give the oblique source coordinate weight
\(w_y\), and give the other three coordinates weight one.  Then

\[
 \operatorname {diag}(w_y,1)\operatorname {diag}(r_y,u)
 =\operatorname {diag}(f_y,u),                      \tag{3.2}
\]

so \(H_y\) is self-adjoint in this positive metric and the Fourier
generator \(-iAH_y\) is skew-adjoint.  In the original orthonormal frames,
one explicit metric is

\[
 \begin{aligned}
 W_s&=S_y^{-T}\begin{pmatrix}w_y&0\\0&1\end{pmatrix}S_y^{-1}
 =\begin{pmatrix}w_y+\xi_y^2&-\xi_y\\-\xi_y&1\end{pmatrix},\\
 W_d&=T_y^{-T}T_y^{-1}
 =\begin{pmatrix}1+\eta_y^2&-\eta_y\\-\eta_y&1\end{pmatrix}.
 \end{aligned}                                      \tag{3.3}
\]

Their determinants are \(w_y\) and one.  Hence this is an exact positive
symmetrizer for every \(|y|<2\) when \(u\ne0\).  On
\(|y|\leq1\),

\[
 {1\over14}\leq{\lambda_o\over u^2}\leq{9\over2},
 \qquad {1\over2}\leq w_y\leq{3\over2}.            \tag{3.4}
\]

Together with \(|u|\geq u_0>0\), (2.5) and (3.4) make the metric uniformly
equivalent to physical energy.  No uniformity is possible as \(u\to0\).
Indeed, when \(u=0\) and \(yv\ne0\), both round trips vanish but the full
block is a nonzero square-zero operator, not a rotation.

The radial determinants are

\[
 \det F=-{u^2(4-y^2)\over2\sqrt D},
 \qquad
 \det R=-{u^2(y+2)\over\sqrt D}.                    \tag{3.5}
\]

When \(u\ne0\), consequently:

* at \(y=-2\), both maps have rank one and the oblique branch is completely
  dark: \(f_{-2}=r_{-2}=0\);
* at \(y=2\), \(F\) has rank one but \(R\) is invertible.  The oblique
  block is nonzero nilpotent, so this is a one-way secular endpoint rather
  than a rotation;
* for \(y>2\), \(\lambda_o<0\).  The oblique eigenvalues of the Fourier
  generator are real, with exact growth/decay rate

  \[
  A|u|\sqrt{{(y-2)(y+2)^2\over2D}};                 \tag{3.6}
  \]

* for \(y<2\), apart from the zero at \(-2\) and, when \(v\ne0\), the Jordan
  collision at \(-4\), the two branches are elliptic.  In particular every
  reality-safe shift \(|y|\leq1\) is full-rank elliptic when \(u\ne0\).

## 4. A reality pair leaves only the tangential scalar line

For \(uv\ne0\) and \(y\notin\{0,-4\}\), the two source eigenlines of
\(RF\) are \(e_t\) and the line of slope

\[
 m(y)={-uvy^2(y+1)/D\over\lambda_o-u^2}
 ={2v(y+1)\over u(y+4)}
 =-{4\alpha(y+1)\over\beta(y+4)}.                  \tag{4.1}
\]

This slope is injective in \(y\), since

\[
 (y_1+1)(y_2+4)-(y_2+1)(y_1+4)=3(y_1-y_2).         \tag{4.2}
\]

Thus two distinct nonexceptional heights have only \(e_t\) as a common
source invariant line.  In particular, for every \(0<|y|\leq1\), the
reality-required pair \(\{y,-y\}\) has only \(e_t\) in common.  A closed
one-source-coordinate star needs a source line invariant under the return
through every daughter, so an autonomous fixed-\(E_0\) realization is
forced onto this tangential branch.

This argument deliberately concerns a scalar-source star.  A larger map
that retains both source polarizations need not have a common invariant
line and is not ruled out.

## 5. The common branch is an exact phase multiplier

The tangential identity (2.4) extends from one edge to the entire normal
charge lattice.  Fix a horizontal wavevector

\[
 p_h=Ae_r,
 \qquad \kappa=p_h\cdot E_0=Au,                     \tag{5.1}
\]

let \(\theta=e_z\cdot x\), and take an arbitrary real horizontal gate shear,
periodic and smooth in \(\theta\) and time-integrable on the interval of
interest,

\[
 W(x,t)=W_h(\theta,t),
 \qquad W_h(\theta,t)\in\mathbb R^3,
 \qquad W_h\cdot e_z=0.                            \tag{5.2}
\]

This includes an arbitrary bundle
\(W_h=\sum_jE_jF_j(\theta,t)\) of real pure-normal shears, not only the
single fixed polarization \(E_0F\) used to derive (1.7).

Every tangential-fiber perturbation has the form

\[
 V(x,t)=e^{ip_h\cdot x}e_t Z(\theta,t).             \tag{5.3}
\]

Because \(e_t\perp p_h,e_z\), it is divergence free at every normal
charge.  Moreover,

\[
 (V\cdot\nabla)W=0,
 \qquad
 (W\cdot\nabla)V=i(p_h\cdot W_h)V,                  \tag{5.4}
\]

and the right side is already divergence free.  The exact inviscid
linearized equation is therefore the scalar multiplication equation

\[
 \boxed{\partial_tZ=-i\,p_h\cdot W_h(\theta,t)Z.}  \tag{5.5}
\]

For the single-polarization specialization
\(W_h=E_0F\), equivalently, if

\[
 F=\sum_a\widehat F_a e^{ia\theta},
 \qquad
 Z=\sum_mz_m e^{im\theta},
 \qquad
 \widehat F_{-a}=\overline{\widehat F_a},           \tag{5.6}
\]

then

\[
 \dot z_m=-i\kappa\sum_a\widehat F_a z_{m-a}.      \tag{5.7}
\]

The convolution in (5.7) is Hermitian, but its stronger physical-space
form (5.5) gives

\[
 Z(t,\theta)
 =\exp\left[-i\int_0^t p_h\cdot W_h(\theta,s)\,ds\right]Z(0,\theta),
 \qquad
 \boxed{|Z(t,\theta)|=|Z(0,\theta)|.}              \tag{5.8}
\]

Starting with one normal-charge coefficient
\(Z(0,\theta)=ce^{im_0\theta}\), any later charge spread therefore sums
to a function of modulus exactly \(|c|\).  After the \(-p_h\) reality
completion, the real-field supremum remains at most the same fixed factor
\(2|c|\).  The coefficient \(\ell^1\) norm can grow, but it cannot become
C161's missing physical \(\sqrt q\) point-amplitude gain: the charge phases
cancel exactly.

Viscosity cannot turn this particular invariant fiber into a focusing
mechanism.  On the same unlocalized shear, the scalar equation becomes

\[
 \partial_tZ
 =\nu(\partial_\theta^2-|p_h|^2)Z
  -i(p_h\cdot W_h)Z.                               \tag{5.9}
\]

The imaginary potential cancels from the modulus identity,

\[
 \partial_t|Z|^2
 =\nu\partial_\theta^2|Z|^2
  -2\nu|\partial_\theta Z|^2-2\nu|p_h|^2|Z|^2,      \tag{5.10}
\]

so the periodic maximum principle gives

\[
 \|Z(t)\|_{L^\infty_\theta}
 \le e^{-\nu|p_h|^2t}\|Z(0)\|_{L^\infty_\theta}.   \tag{5.11}
\]

This is a statement about the exact scalar fiber only; it does not construct
the time-dependent gate as an unforced background or control localization
collars.

Thus the full-fiber calculation closes the simplest paired-polarization
repair negatively.  The remaining same-geometry target is precise: use an
unfolded or time-dependent gate sequence that genuinely mixes the two
polarization branches while retaining physical endpoint coherence, or
prove that every such sequence inherits an analogous unitary-fiber
obstruction.
