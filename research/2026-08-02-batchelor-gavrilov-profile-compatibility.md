# The Batchelor--Gavrilov base-profile compatibility gate

**Date:** 2026-08-02

**Status:** exact profile comparison and rational interval certificate;
independently identified in audit

**Scope:** the standard pressure-modulated Gavrilov thin-tube seed and the
full-edge-matched Batchelor column.  This note neither excludes other compact
steady Euler seeds nor constructs a curved unstable ring.

## 1. Why curvature is not yet the first error

The spectral module selects the straight Batchelor column

\[
 V_B(r)=Q_*{1-e^{-r^2}\over r},
 \qquad W_B(r)=e^{-r^2}.                            \tag{1.1}
\]

The geometric module uses the standard pressure-modulated Gavrilov bubble.
Its leading straight-tube profile, recorded in the three-scale rectifier
ledger, is

\[
 V_G(s)=\kappa_0s g(s),
 \qquad W_G(s)={\kappa_0s g(s)\over\sqrt2}.          \tag{1.2}
\]

Consequently

\[
                         W_G'={V_G'\over\sqrt2}.     \tag{1.3}
\]

Uniform spatial rescaling multiplies both derivatives by the same radial
factor, uniform velocity rescaling multiplies both by the same amplitude,
and an axial Galilean shift changes neither derivative.  Axial reflection
can only reverse the sign in (1.3).  Thus the invariant magnitude condition
for this standard family is

\[
                         \left|{W_G'\over V_G'}\right|
                         ={1\over\sqrt2}.            \tag{1.4}
\]

## 2. Exact mismatch at the certified AO ring

Let `x_*=r_*^2` and `Q_*` be the jointly optimized parameters from the
full-edge certificate.  Differentiating (1.1) gives

\[
 V_B'(r_*)
 =Q_*{(2x_*+1)e^{-x_*}-1\over x_*},
 \qquad
 W_B'(r_*)=-2\sqrt{x_*}e^{-x_*}.                    \tag{2.1}
\]

The certified intervals

\[
 0.59671214<x_*<0.59671216,
 \qquad
 0.8278572<Q_*<0.8278581                            \tag{2.2}
\]

give, using only rational Taylor enclosures for the exponential,

\[
 V_B'(r_*)>0>W_B'(r_*),
 \qquad
 2.95<-{W_B'(r_*)\over V_B'(r_*)}<2.96.             \tag{2.3}
\]

In particular, (2.3) disagrees with both signs of (1.4).  An axial
Galilean shift may match the velocity values at one radius, but it cannot
match this first jet.  There is also a global mismatch: the standard
Gavrilov modulation is compactly supported, whereas the Batchelor swirl
has a `Q_*/r` tail.

## 3. Consequence for the quasimode program

The previous architecture informally treated the finite-curvature error as
the first difference between the certified straight mode and the compact
ring.  Equations (1.3)--(2.3) show that this is not justified.  For the
standard seed, the base coefficients already differ by order one at the AO
ring.  Therefore an `O(M^-1)` normalized residual cannot be inferred merely
by bending the Batchelor mode with curvature `M^-1`.

The spectral and geometric modules can be joined only after one of the
following genuinely new results.

1. **Compact Batchelor realization.**  Construct a compact steady Euler
   tube or ring whose rescaled local column converges, in the coefficient
   norm needed by the PDE block, to (1.1) on the complete mode-localization
   region.
2. **Gavrilov-family spectral certificate.**  Select an actually realizable
   compact Gavrilov-type profile, prove the AO high-frequency eigenvalue
   expansion for it, identify its full BAS edge, and establish the same
   uniform sector bound.
3. **Quantified interpolation family.**  Produce compact steady profiles
   approaching the edge-matched Batchelor coefficients with mismatch at most
   the required quasimode residual, uniformly through the curvature limit.

This gate is narrower than an impossibility theorem.  It rules out only the
unproved direct identification of the two existing profiles; it does not
rule out any of the three repairs above.

## 4. Subsequent repair

The second option has now been carried out at principal-symbol level.  The
companion note
`2026-08-02-locked-pitch-gavrilov-edge-profile.md` constructs an explicit
compact pressure modulation satisfying the actual relation
`W_G=V_G/sqrt(2)` whose resonant ring is the unique full BAS edge.  It does
not identify this carrier with the Batchelor profile; it replaces the
Batchelor carrier by a compatible locked-pitch one.

The remaining spectral issue is sharper.  Literal AO Assumption A uses a
quotient that must blow up for every compact locked-pitch bump, although it
cancels from the physical Rayleigh equation.  The companion note records a
self-audited q-free transfer of the AO gluing proof.  Independent audit,
the PDE-sector uniform bound, and the curved residual remain open.

## 5. Reproducibility

`checks/batchelor_gavrilov_profile_compatibility.py` verifies (2.2)--(2.3)
with rational interval arithmetic and checks the invariance of the derivative
ratio under the allowed uniform scalings.
