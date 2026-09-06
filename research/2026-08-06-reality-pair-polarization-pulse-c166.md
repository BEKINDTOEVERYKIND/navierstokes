# C166: a reality-paired pulse mixes the full first-neighbour polarization fiber

**Date:** 2026-08-06
**Status:** exact limiting-cone, first-neighbour Galerkin Lie algebra and
two-pulse construction; no full-ladder pulse, uniform charge star, or
nonlinear Navier--Stokes stage
**Checker:**
[checks/reality_pair_polarization_pulse_c166.py](../checks/reality_pair_polarization_pulse_c166.py)

## 0. Claim boundary

C164 leaves open whether time-dependent gate polarizations can escape the
tangential phase-multiplier fiber.  They can in the exact first-neighbour
Galerkin compression, already for one source and a reality-required signed
charge pair.  On the complete two-polarization source and first-daughter
space, the compressed radial and tangential pure-normal gate operators have
a nonzero exact commutator.  For signed height pair \(\{y,-y\}\), its source
block contains

\[
 [H_r,H_t]e_{\sigma}
 ={2y^2(4-y^2)\over y^4+4y^2+16}\,e_t,
 \qquad 0<y<2.                                      \tag{0.1}
\]

This is not merely an infinitesimal effect inside that compression.  At
\(y=1\), a radial quarter-pulse followed by one tangential pulse gives an
explicit finite endpoint with nonzero source-polarization mixing.  After
reality completion, the exhibited point amplitude is exactly
\(\sqrt{3/2}\) times the initial supremum.  The raw endpoint coefficient
energy is \(113/112\) times the initial energy; after equal-energy
renormalization the exhibited amplitude ratio is still
\(\sqrt{168/113}>1\).  The total normalized pulse action is less than
\(5/2\), independent of any cascade parameter.

The positive result is deliberately finite-dimensional.  The space keeps
the source and its first daughters \(p\pm g\), but the same physical gate
also sends those daughters to \(p\pm2g\).  Thus it is not an invariant
subspace of the full charge ladder.  Moreover, the radial gate direction is
adapted to one horizontal source direction, the second pulse is an active
nonnormal stretch rather than an energy-preserving rotation, and the gain
is a fixed constant rather than the \(\sqrt q\) required by C161.  No
unforced pulse, repeated-charge estimate, source sheet, localization,
viscosity, pressure, BAFL estimate, or one-cell stage is claimed.

## 1. The physical signed-pair operators

Fix one C164 cone source

\[
 p=A\left(e_r+{e_z\over\sqrt3}\right),
 \qquad
 g={Ay\over\sqrt3}e_z,
 \qquad 0<y<2.                                      \tag{1.1}
\]

The real gate contains both Fourier coefficients \(+g\) and \(-g\).
Equivalently, its physical spatial factor is \(2\cos(g\cdot x)\), with the
normalization that each signed Fourier coefficient has amplitude one.  Put

\[
 \eta\in\{+y,-y\},\qquad
 D_\eta=\eta^2+2\eta+4.                             \tag{1.2}
\]

Use the source basis \((e_\sigma,e_t)\) from C164 and, at daughter
\(p+\eta Ae_z/\sqrt3\), its basis
\((e_{\perp,\eta},e_t)\).  For a radial horizontal gate polarization
\(E=e_r\), C164's full forward and reverse blocks are

\[
 F^r_\eta=
 \begin{pmatrix}a_\eta&0\\0&1\end{pmatrix},
 \qquad
 R^r_\eta=
 \begin{pmatrix}c_\eta&0\\0&1\end{pmatrix},       \tag{1.3}
\]

where

\[
 a_\eta=-{4-\eta^2\over2\sqrt{D_\eta}},
 \qquad
 c_\eta=-{\eta+2\over\sqrt{D_\eta}}.              \tag{1.4}
\]

For the orthogonal tangential gate polarization \(E=e_t\), they are

\[
 F^t_\eta=b_\eta E_{21},\qquad
 R^t_\eta=d_\eta E_{21},qquad
 b_\eta=-{\eta\over2},\qquad
 d_\eta=-{\eta\over\sqrt{D_\eta}}.                \tag{1.5}
\]

On

\[
 \mathcal X_y=S\oplus D_+\oplus D_-,               \tag{1.6}
\]

define the complete first-neighbour signed-pair matrices

\[
 H_\bullet=
 \begin{pmatrix}
 0&R^\bullet_+&R^\bullet_-\\
 F^\bullet_+&0&0\\
 F^\bullet_-&0&0
 \end{pmatrix},\qquad \bullet\in\{r,t\}.           \tag{1.7}
\]

The frozen Fourier generators are \(-iAH_r\) and \(-iAH_t\).  All entries
in (1.7) include the physical return through the reality companion of the
gate; no one-sided complex gate is being used.  However, (1.7) is precisely
the Galerkin compression \(\Pi_{\mathcal X_y}L\Pi_{\mathcal X_y}\) of the
frozen linearized operator, not its restriction to an invariant physical
subspace.  Calling the following propagator exact always means exact for
this displayed finite matrix.

## 2. Energy metric and the exact Lie bracket

The radial operator has the positive symmetrizer

\[
 M_y=\operatorname {diag}\left(
 1,1,{2\over2-y},1,{2\over2+y},1\right),            \tag{2.1}
\]

in the ordered coordinates

\[
 (s_\sigma,s_t,d_{+,\perp},d_{+,t},
                    d_{-,\perp},d_{-,t}).          \tag{2.2}
\]

Indeed, \(c_\eta/a_\eta=2/(2-\eta)\), so

\[
 M_yH_r=H_r^TM_y.                                   \tag{2.3}
\]

Consequently \(e^{-i\tau H_r}\) is exactly unitary in this positive
metric.  On \(0<y\le1\), the metric is uniformly equivalent to the
physical coefficient \(\ell^2\) energy.

The tangential gate is qualitatively different:

\[
 H_t^2=0,qquad H_t\ne0.                            \tag{2.4}
\]

Thus \(e^{-isH_t}=I-isH_t\).  A nonzero nilpotent matrix cannot be
self-adjoint in any positive inner product, so this pulse is an active
stretching operation, not another weighted unitary rotation.

The source block of the Lie bracket is

\[
 ([H_r,H_t])_{SS}
 =\sum_{\eta=\pm y}
 \left(R^r_\eta F^t_\eta-R^t_\eta F^r_\eta\right).
                                                               \tag{2.5}
\]

Its only nonzero entry is the lower-left entry.  For one signed height,

\[
 b_\eta-d_\eta a_\eta
 =-{\eta(\eta+4)\over D_\eta}.                    \tag{2.6}
\]

After adding the reality pair,

\[
 \begin{aligned}
 \chi(y)
 &:=-{y(y+4)\over y^2+2y+4}
   +{y(4-y)\over y^2-2y+4}\\
 &={2y^2(4-y^2)\over y^4+4y^2+16}>0,
 \qquad 0<y<2.                                    \tag{2.7}
 \end{aligned}
\]

Therefore

\[
 \boxed{[H_r,H_t]e_\sigma=\chi(y)e_t.}             \tag{2.8}
\]

Reality cancels the order-\(y\) part but leaves order-\(y^2\)
source-polarization mixing.  At the convenient safe height \(y=1\),
\(\chi(1)=2/7\).
Hence the oblique and tangential branches are not protected by a common
invariant decomposition once the gate polarization is allowed to change in
time.

## 3. Exact finite two-pulse ordering defect

Set \(A=1\) by normalized time and define

\[
 \lambda_y^2=\sum_{\eta=\pm y}a_\eta c_\eta>0.
                                                               \tag{3.1}
\]

Starting from the oblique source coordinate, the radial pulse is exactly

\[
 U_r(\tau)e_\sigma
 =\cos(\lambda_y\tau)e_\sigma
  -i{\sin(\lambda_y\tau)\over\lambda_y}
       \sum_{\eta=\pm y}a_\eta e_{\eta,\perp}.     \tag{3.2}
\]

The identities

\[
 b_y+b_{-y}=0,qquad
 \sum_{\eta=\pm y}d_\eta a_\eta=-\chi(y)          \tag{3.3}
\]

then give the exact two-pulse composition

\[
 \begin{aligned}
 U_t(s)U_r(\tau)e_\sigma
 ={}&\cos(\lambda_y\tau)e_\sigma
  -i{\sin(\lambda_y\tau)\over\lambda_y}
       \sum_\eta a_\eta e_{\eta,\perp}\\
 &+{s\chi(y)\sin(\lambda_y\tau)\over\lambda_y}e_t
  -is\cos(\lambda_y\tau)
       \sum_\eta b_\eta e_{\eta,t}.              \tag{3.4}
 \end{aligned}
\]

Here \(U_t(s)=I-isH_t\).  Reversing the order leaves
\(\sum b_\eta e_{\eta,t}\) dark under \(H_r\), and hence

\[
 \begin{aligned}
 &(U_t(s)U_r(\tau)-U_r(\tau)U_t(s))e_\sigma\\
 &\quad={s\chi(y)\sin(\lambda_y\tau)\over\lambda_y}e_t
 -is\{\cos(\lambda_y\tau)-1\}
       \sum_\eta b_\eta e_{\eta,t}.               \tag{3.5}
 \end{aligned}
\]

Equation (3.5) is a finite-action, exact ordering defect; it does not rely
on truncating a Baker--Campbell--Hausdorff series.  Its mixed derivative at
the origin recovers (2.8).

## 4. One explicit compressed active focus at \(y=1\)

At \(y=1\),

\[
 \lambda_1^2={8\over7},qquad
 \tau_*={\pi\over2\lambda_1}={\pi\sqrt{14}\over8}. \tag{4.1}
\]

Choose the tangential pulse action \(s=1\).  Formula (3.4) becomes

\[
 \boxed{
 U_t(1)U_r(\tau_*)e_\sigma
 ={\sqrt{14}\over14}e_t
  +i{3\sqrt2\over8}e_{+,\perp}
  +i{\sqrt{42}\over8}e_{-,\perp}.}                 \tag{4.2}
\]

The normalized total action is

\[
 1+\tau_*=1+{\pi\sqrt{14}\over8}<{5\over2}.       \tag{4.3}
\]

The endpoint genuinely leaves the scalar phase fiber: it has nonzero
source \(e_t\) created from an incoming \(e_\sigma\).  Its physical
half-lattice coefficient energy is

\[
 {1\over14}+{18\over64}+{42\over64}
 ={113\over112}.                                   \tag{4.4}
\]

Before the active second pulse, the two radial daughters have physical
coefficient energy \(15/16\), while their \(M_1\)-energy is exactly one.
After the second pulse the same radial part remains and the new tangential
source raises the \(M_1\)-quadratic form to \(15/14\).  Thus (4.2) is not
being mislabeled as an energy-preserving C161 star.  In particular,
\(113/112>1\) records work done on the perturbation by the prescribed
frozen gate.  It is ordinary nonunitary linearized energy transfer from the
base pump, not energy created by an unforced Euler or Navier--Stokes
solution.  A full nonlinear realization would have to retain the
compensating pump deformation, higher charge modes, and wake.

There is nevertheless an exact pointwise concentration.  At \(z=0\), the
two daughter polarizations are

\[
 e_{+,\perp}=-{2\over\sqrt7}e_r+{\sqrt3\over\sqrt7}e_z,
 \qquad e_{-,\perp}=e_z.                           \tag{4.5}
\]

Their coherent complex amplitude in (4.2) is

\[
 i\left(-{3\sqrt{14}\over28}e_r
          +{5\sqrt{42}\over28}e_z\right),         \tag{4.6}
\]

whose squared norm is

\[
 {9\cdot14+25\cdot42\over28^2}={3\over2}.         \tag{4.7}
\]

After adding the \(-p\) reality companion, set \(z=0\) and choose the common
horizontal phase \(e^{ip_h\cdot x}=-i\).  Then (4.6) is real after reality
completion, while the real source-\(e_t\) coefficient in (4.2) contributes
exactly zero at that point.  The initial real source has supremum
\(2|e_\sigma|=2\), while the displayed endpoint point has size
\(2\sqrt{3/2}\).  Therefore the first-neighbour Galerkin circuit has the
exact exhibited raw ratio

\[
 \boxed{G_{\rm exhibited}=\sqrt{3/2}>1.}           \tag{4.8}
\]

This ratio includes the small active energy increase (4.4).  If the whole
endpoint is rescaled by \(\sqrt{112/113}\) back to the initial
half-lattice coefficient energy, the same point instead certifies

\[
 \boxed{G_{\rm exhibited,energy\ normalized}
   =\sqrt{{3/2\over113/112}}=\sqrt{168/113}>1.}    \tag{4.9}
\]

Thus the displayed coherence is not explained solely by the gate-supplied
increase in perturbation energy.  Neither ratio is a claim about the
uncomputed full charge-ladder propagator.

This explicitly shows why the C164 tangential phase-multiplier theorem
cannot be extended to arbitrary polarization pulses.

## 5. What remains load-bearing

The result chooses \(E=e_r\) relative to one source.  Across the C161
two-dimensional source sheet, the radial direction varies with source
angle, whereas a physical gate frequency has only one shared polarization.
In addition, the finite signed-pair circuit omits the repeated charge
outputs \(p\pm2g,p\pm3g,\ldots\) created during a finite pulse.  Those are
not small merely because (4.2) is exact in the first-neighbour matrix.

The same-geometry target is therefore narrower than before:

> construct a shared polarization/time code whose **full charge-ladder**
> propagator approximates the source-adapted radial quarter-pulses
> uniformly over the \(q^2\) source sheet, while producing the normalized
> \(q\)-daughter coherence, controlling the active-stretch energy/wake and
> repeated shifts at the C146 powers, and using an unforced gate trajectory.

C166 removes a broad obstruction inside the first-neighbour algebra:
active time colors can mix the two C164 branches and the compressed circuit
can focus one reality-completed source with order-one action.  It does not
claim that the same pulse focuses the actual full gate flow, and it does
not remove the uniform shared-control, normalization, pump-work, or
full-ladder obstruction that separates the current construction from a
one-cell Navier--Stokes stage.
