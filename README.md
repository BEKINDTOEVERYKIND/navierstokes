# Navier--Stokes breakdown research checkpoint

This repository records a theorem-oriented attempt to resolve the
three-dimensional Navier--Stokes Millennium problem.  It is **not** a
claimed solution.

The target is Clay alternative (D): smooth periodic initial data and a
smooth, rapidly time-decaying force for which no global smooth solution
exists.  The force is part of the official problem, but it may not hide the
singularity; its full residual must extend smoothly and be flat at the
terminal cascade time.

## Current outcome

The work has reduced the constructive route to one local analytic theorem.
Several objects that were previously only formal are now exact.

1. The affine Kelvin solution realizes the required amplification law
   \(a\sim\ell^{-\gamma}\), and all scalar constraints meet in the sharp
   window
   \[
   1<\gamma<\frac32.
   \]
2. Scaled, disjoint copies of Gavrilov's compact steady Euler flow give an
   exact discretely self-similar stationary wake.  Packing
   \(N_j\asymp K_j^3\) microbubbles of diameter
   \(\delta_j=\ell_j/K_j\) realizes the complete polynomial-carrier
   energy, derivative, Reynolds, and viscous-loss ledger with **zero**
   leading Euler defect.
3. A sign-changing Gavrilov pressure modulation makes each bubble
   individually zero-angular-momentum.  One improper reflected copy then
   cancels helicity.  The complete finite-dimensional interpolation defect
   has a componentwise compact symmetric anti-divergence with the correct
   inertial scale.
4. The natural equal-shell two-wave transition has an exact rank defect.
   Unequal-radius doublets recover both daughter polarizations.  Three
   explicit integer beat directions synthesize the
   \(\gamma=5/4\) affine strain exactly and give a rank-five rational chart
   of every symmetric trace-free strain; an exact checker verifies a
   \(5\times5\) minor of \(-1\,214\,003\,700\).

There is also a new mandatory correction.  A bare thin Gavrilov endpoint
cannot be dressed by correctors confined to its own disjoint ball.  At the
second velocity jet, flat Navier--Stokes residual forces

\[
 {\mathsf G}_{ij}(W)=\int\partial_kW_i\,\partial_kW_j=c\delta_{ij},
\]

whereas the thin seed has

\[
 G_\parallel-G_\perp
 =\pi^2\varepsilon\int\chi^2+o(\varepsilon)>0.
\]

Its first differentiated pressure therefore has an unavoidable
\(|x|^{-3}\) quadrupole and the second velocity jet an \(|x|^{-4}\) tail.
The final construction must retain that global pressure/centre wake or a
nonzero endpoint oscillatory bath.

## The single missing theorem

The remaining prize-level target is a **three-beat forward
Kelvin--Reynolds parametrix**:

- evolve the complete finite charged material-phase lattice rather than a
  static or finitely closed triad;
- keep the explicit rank-five child-strain chart uniformly invertible;
- route unmatched and nonzero-charge products into the annular wake;
- retain the nonlocal pressure jets and a viscously admissible endpoint
  bath;
- solve through
  \(M_j\asymp j^2/\log j\) at the polynomial carrier \(K_j=j^A\); and
- obtain Gevrey bounds and residual/seam errors \(e^{-c j^2}\).

That last rate beats every physical derivative cost and would make the
terminal force \(C^\infty\)-flat.  No theorem in this repository yet
establishes this parametrix, so no Millennium conclusion is claimed.

## Exact notes

- [Exact Gavrilov wake and packed carrier](research/2026-07-29-gavrilov-dss-wake-construction.md)
- [Localized Gavrilov transition ledger](research/2026-07-29-gavrilov-active-transition-ledger.md)
- [Two-colour rank theorem and three-beat repair](research/2026-07-29-two-colour-endpoint-rank.md)
- [Viscous endpoint-jet obstruction](research/2026-07-29-gavrilov-viscous-endpoint-jet-obstruction.md)
- [Forward one-phase inverse audit](research/2026-07-29-one-phase-cell-inverse-audit.md)
- [Polynomial-carrier Gevrey ledger](research/2026-07-29-polynomial-carrier-ledger.md)
- [Flat-force/Borel attack](research/2026-07-29-flat-force-borel-attack.md)
- [All-support Laurent pole no-go](research/2026-07-29-laurent-null-pole-no-go.md)
- [Audit of the claimed 2026 Navier--Stokes profile](research/2026-07-29-shahmurov-ns-preprint-audit.md)
- [Audit of the claimed axisymmetric reduction](research/2026-07-29-shahmurov-axisymmetric-reduction-audit.md)

## Reproducible algebra

No extra dependencies are needed:

```bash
python checks/forced_recurrence_ledger.py
python checks/two_colour_endpoint_rank.py
```

All `PASS` statements in these programs concern exact algebraic identities
or scaling inequalities.  They do not assert existence of the missing
Navier--Stokes parametrix.
