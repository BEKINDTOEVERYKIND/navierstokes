# Compact-ring profile flexibility after the C67 mismatch

**Date:** 2026-08-02

**Status:** exact Gavrilov--Constantin--La compatibility identities and
rational-interval obstruction; literature scope audit self-derived

**Scope:** which straight-column profiles can arise from known smooth
compact Euler-ring constructions.  This note is not an existence theorem
for a new annular ring and not a Navier--Stokes singularity claim.

## 1. Outcome

The Gavrilov--Constantin--La hodograph family has a sharp straight-limit
compatibility law.  In the circular straight-tube limit, the ratio of axial
to azimuthal velocity must obey

\[
 \left({W\over V}\right)^2={1\over2}+{C\over s^2}. \tag{1.1}
\]

For the standard **stagnation-core, finite-pitch** seed used in the compact
wake, `W/V` remains finite as both components vanish at the central circle.
That extra property forces `C=0`, hence

\[
                         W=\pm{V\over\sqrt2}.       \tag{1.2}
\]

This is intrinsic to the finite-pitch stagnation-core branch, not merely the
constant-`M` example.  Smooth core regularity alone does **not** force
`C=0`: a nonzero axial core velocity permits `C>0`.  The arbitrary pressure
cutoff still multiplies both components together and cannot change the
underlying ratio.  Thus the free hodograph function does not realize the C65
Batchelor profile inside the standard vanishing-core branch, while broader
core-through-flow branches remain available.

There are two surviving alternatives.

1. The locked-pitch profile can itself be redesigned.  This is done in
   `2026-08-02-locked-pitch-gavrilov-edge-profile.md`, which constructs a
   compact profile in (1.2) with a unique full BAS edge.
2. A `C != 0` seed escapes (1.2).  It may be smooth through the core with a
   nonzero axial velocity; after a pressure cutoff the final compact support
   is still toroidal and hollow.  No current theorem gives the prescribed
   profile and uniform thin-ring estimates required here.

The strongest published prescribed-generator thin-ring theorem also cannot
simply import the C65 Batchelor profile: its monotonicity hypotheses fail at
the target ring by a strict interval-certified gap.

## 2. Exact hodograph constraint

Constantin--La's formulation of Gavrilov's construction has

\[
 p'(\psi)=2M(\psi),\qquad
 A'M-3M^2-2AM'=0,                                  \tag{2.1}
\]

where `A=|u|^2/2` is a flux function.  Equivalently,

\[
                 \left({A\over M^2}\right)'={3\over M}. \tag{2.2}
\]

In a circular straight-tube limit, let `s` be the minor radius, `R` the
major radius, `V(s)` the azimuthal/poloidal speed, and `W(s)` the
axial/toroidal speed.  The leading identities are

\[
 \psi_s=RV,\qquad p_s={V^2\over s},\qquad
 M={V\over2Rs},\qquad A={V^2+W^2\over2}.            \tag{2.3}
\]

Put `h=W/V`.  Substitution of (2.3) into (2.1) cancels every derivative of
the freely selected amplitude `V` and leaves

\[
                       s(h^2)'+2h^2=1.              \tag{2.4}
\]

Integrating gives (1.1).  If `h` remains finite as `s -> 0`, then `C=0`.
This applies to the standard Gavrilov stagnation core, where both leading
components vanish with finite pitch.  It is important not to replace that
hypothesis by smoothness alone.  For every `C>0`, the local profile

\[
 V(s)=s,\qquad W(s)=\sqrt{C+{s^2\over2}}            \tag{2.5}
\]

is smooth as an axisymmetric velocity through `s=0` and satisfies (1.1):
`V e_theta` is linear in Cartesian transverse coordinates and `W` is a
smooth even function of `s`.  Here the axial core speed is nonzero and
`W/V` diverges.  Thus higher hodograph freedom cannot repair an order-one
mismatch within the finite-pitch stagnation branch, but a core-through-flow
branch is a genuine escape.

The primary derivations are in Gavrilov,
<https://arxiv.org/abs/1810.08020>, and Constantin--La,
<https://arxiv.org/abs/1903.11699>.

## 3. Why separate component cutoffs do not help

On a fixed toroidal flux foliation, write the toroidal contribution as
`F(psi)/r`.  If the poloidal and toroidal components were multiplied by
separate functions `a(psi)` and `b(psi)`, the new kinetic energy would be

\[
 \widetilde A
 =a^2A+{(b^2-a^2)F^2\over2r^2}.                    \tag{3.1}
\]

Localizability requires \(\widetilde A\) to depend only on `psi`.  Since `r`
varies along a nondegenerate toroidal flux surface and `F` is nonzero,
(3.1) forces

\[
                            b^2=a^2.                \tag{3.2}
\]

Thus Gavrilov's pressure modulation supplies one common envelope, up to a
component sign.  Independent axial/swirl profile design requires changing
the flux geometry, not merely changing the cutoff.

## 4. The `C != 0` escape

For `C != 0`, equation (1.1) gives the regular local law

\[
 h(s)=\pm\sqrt{{1\over2}+{C\over s^2}}              \tag{4.1}
\]

wherever its radicand is positive.  One scalar amplitude `V(s)` remains free
and `W=hV`.  For `C>0`, (4.1) can be smooth at the core when `V=O(s)` and
`W` has a nonzero even limit, as in (2.5).  For `C<0`, the real profile is
naturally restricted away from the axis.  Either case can match more local
data than the finite-pitch branch and is not excluded by C67, which was
correctly scoped to the standard seed family.

What is missing is global, not formal: one must select and continue a
hodograph solution around a closed toroidal pressure tube, keep the relevant
levels nested through the cutoff, prescribe the needed straight profile,
and establish uniform Gevrey and curvature estimates as the ring becomes
thin.  No theorem in the sources above supplies that package.  The symmetry theorem of
Peralta-Salas--Slobodeanu, arXiv:2606.13462,
<https://arxiv.org/abs/2606.13462>, additionally constrains analytic
localizable regular pressure domains toward the axisymmetric toroidal
setting; it does not construct the missing annular branch.

## 5. Why the prescribed-generator thin-ring theorem does not import C65

Cao--Zhan, *Desingularization of vortex rings in 3 dimensional Euler
flows*, arXiv:2009.13210, <https://arxiv.org/abs/2009.13210>, construct thin
axisymmetric rings for prescribed Grad--Shafranov generators under separate
positivity and monotonicity hypotheses.  Their result gives compact
vorticity and qualitative local profile convergence, but not compact
velocity, a quantitative high-`C^k` convergence rate, or a semigroup
transfer.

There is also an exact incompatibility with C65.  Let

\[
 x=s^2,\qquad E=e^{-x},
\]

and allow the spectrally harmless axial shift

\[
 V={Q(1-E)\over s},\qquad W=E-c.                   \tag{5.1}
\]

For major radius `R`, the required Grad--Shafranov generators are

\[
 HH_\Psi={2RxE(E-c)\over Q(1-E)},                  \tag{5.2}
\]

\[
 R^2(-B')=2RE\left[Q-{x(E-c)\over Q(1-E)}\right]. \tag{5.3}
\]

Nonnegativity of (5.2) at the AO ring forces `c<=E`.  Put

\[
 A={x(E-c)\over1-E},\qquad
 G=E\left(Q-{A\over Q}\right).
\]

Then

\[
 G_x={E\over Q}(A-A'-Q^2),                         \tag{5.4}
\]

and the exact identity

\[
 A-A'-{xE\over1-E}
 ={(E-c)(x-(1-E))\over(1-E)^2}                     \tag{5.5}
\]

is nonnegative at C65.  Its certified parameters satisfy

\[
 {x_*E_*\over1-E_*}-Q_*^2>0.04579,
 \qquad x_*-(1-E_*)>0.1473.                        \tag{5.6}
\]

Thus `G_x>0`, whereas `Psi_x<0`.  Therefore `-B'` decreases as a function
of `Psi`, contrary to the required monotonicity.  If `c>E`, (5.2) is
negative instead.  Every axial shift fails at the target ring itself.

## 6. Literature boundary

No published result located in this audit closes the compact Batchelor or
annular-profile bridge.

- Abe, arXiv:2008.09345, constructs traveling Beltrami rings with compact
  vorticity, but the available regularity/profile constraints do not give a
  smooth compact Batchelor carrier: <https://arxiv.org/abs/2008.09345>.
- Constantin--Drivas--Ginsberg, arXiv:2007.09103, deform Grad--Shafranov
  flows in periodic corrugated pipes, not bounded Euclidean rings:
  <https://arxiv.org/abs/2007.09103>.
- Domínguez-Vázquez--Enciso--Peralta-Salas, arXiv:2005.04380, broaden
  compact support for piecewise-smooth flows with an interface, not the
  smooth carrier required here: <https://arxiv.org/abs/2005.04380>.
- Enciso--Peralta-Salas, arXiv:1210.6271, construct decaying Beltrami vortex
  tubes rather than compact Euler bubbles:
  <https://arxiv.org/abs/1210.6271>.

The actionable order is therefore:

1. use the exact locked-pitch carrier already available in the standard
   Gavrilov family;
2. audit the q-free AO proof transfer and finite-curvature residual; and
3. retain the broader `C != 0` core-through-flow/annular branches as a
   fallback if a different local pitch or axial-shear sign is later required.

The interval calculations in (5.6) and the algebraic ledgers are reproduced
by `checks/compact_ring_profile_gate.py`.
