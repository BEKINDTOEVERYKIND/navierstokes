# Navier--Stokes breakdown research checkpoint

This repository records a theorem-oriented attempt to resolve the
three-dimensional Navier--Stokes Millennium problem.  It is **not** a
claimed solution.

The target is Clay alternative (D): smooth periodic initial data and a
smooth, rapidly time-decaying force for which no global smooth solution
exists.  The force is part of the official problem, but it may not hide the
singularity: the residual must extend smoothly through the terminal cascade
time.  We impose the stronger convenient design that, after subtracting a
fixed smooth terminal background, the shrinking-scale tail is flat there;
it can then be extended with rapid future decay.

## Current outcome

The work has reduced the constructive route to a short linked chain of
local analytic estimates.  Several objects that were previously only formal
are now exact.

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
5. One material phase cannot itself realize a positive-definite covariance,
   but modulo pressure it gives a full stress chart near every tensor with a
   simple top eigenvalue: \(Q=\lambda_{\max}(R)I-R\) has rank two, and two
   kernel rotations plus its transverse block span the trace-free quotient.
   Two exact localized phase controls realize those rotations in one fixed
   torus cohomology class and transport as covectors.  The actual Kelvin
   velocity map is also invertible on each transverse plane, so the full
   chart persists under every finite common affine deformation, with exact
   determinant \(ab/|F^{-T}n|^8\).  For the intended compression this costs
   \(e^{-8t}\).  The loss is physical rather than coordinatic: wave-number
   gain \(s\) and a fixed terminal covariance gap \(q\) require launch
   covariance energy at least \(2qs\).  A bounded terminal pulse restores a
   uniformly conditioned pointwise chart only if a source, reservoir, or
   nonlinear transfer genuinely launches a fresh bath.  Three fixed
   orthogonal phases remain a controlled-splice alternative.
6. The spatially affine viscous carrier/parent/child system has a rigorous
   capture orbit with child/parent ratio \(\asymp r^{-1}\), parent
   depletion \(\asymp r^{-2}\), and a uniform terminal damping gap.
   Its robust version exposes the correct weighted affine-defect and
   exponentially-small-seed norms that a PDE construction must satisfy.
7. The full BAS reduction proves that non-resonant column exponents vanish.
   The Batchelor calculation first exposed and then repaired a transverse
   profile mismatch.  More importantly, an explicit hollow profile in the
   **actual locked-pitch Gavrilov family** is now exactly edge-matched.  Its
   compact pressure modulation has a unique full squared BAS edge
   \((5+\sqrt{22})\Omega_0^2/3\), a strict local fixed-sector maximum, and an
   isolated helical-phase minimum.  Its complete BAS propagator obeys the
   uniform edge bound \(Ce^{\lambda_*|t|}\) with no polynomial prefactor.
   The quotient in literal AO Assumption A necessarily blows up for every
   compact locked-pitch bump, but cancels from the physical Rayleigh
   coefficients.  A source-level audit indicates that AO's eigenmode proof
   transfers to those regularized coefficients; this proof transfer awaits
   independent audit.  The separate Gevrey-2 ledger retains the strict
   cascade choice \(A>4\).
8. The separated second velocity wake is uniform as a Gavrilov torus becomes
   thin when velocity is normalized by its actual size: the
   \(\varepsilon^{-2}\) Laplacian cost cancels the \(\varepsilon^2\)
   tube volume.  The extracted exterior wake remains summable; a newly
   identified same-stage nearest-neighbor term costs \(K^3\), rather than
   the inter-stage \(K^2\), but geometric scale decay still dominates it.
9. The third jet identifies the first genuine thin-aspect loss:
   \(4\nu^2\int\Delta U\otimes\Delta U\) produces a generic
   \(\varepsilon^{-2}\) exterior quadrupole.  It is still perturbative in
   the minor-scale viscous parameter \(\Theta\); polynomial thinning is
   compatible with the cascade.  Its coefficient alone permits
   \(\beta<\gamma-1\), while the stronger minor-scale Reynolds condition
   \(\Theta_j\to0\) requires \(\beta<(\gamma-1)/2\).
10. The fixed-subedge q-free scalar equation has a genuine analytic
    pseudospectral characteristic.  More strongly, after scalar localization
    there is an exact velocity--pressure reconstruction: divergence and the
    azimuthal/axial momentum equations vanish identically, while the only
    residual is \(F_r=\gamma G/(in^2)\).  This removes the unsupported global
    Hodge correction.  The compact Piola comparison is now exact as well:
    it preserves divergence, has no long-ring norm loss, and gives ambient
    residual ratio
    \(C[p^{-2}\|G\|/\|H\|+\varepsilon p+\varepsilon]\).
    On the unchanged analytic ring core, the coalescing scalar operator now
    has an explicit uniform complex-WKB packet with residual
    \(C\eta^2e^{-c\eta^3/s}\).  Conditional on the previously stated
    normalized curved-base and semigroup comparisons, its curvature,
    viscosity, and gain-time parameter balance closes at
    SELF/analytic-proof status.
11. Tracking the pseudomode growth to within \(j^{-g}\) of the full BAS edge
    removes the exponential Duhamel penalty.  The coalescing radial roots have
    effective semiclassical parameter \(s/\eta^3\); both analytic WKB action
    and material-phase coherence close at the single sharp ledger condition
    \(A>3g\).  For the intended \(g=2\) gain, any \(A>6\) gives long-gain
    fidelity and summability.  If the scalar residual itself is retained as
    terminal force, direct \(C^\infty\)-flat domination needs \(A>8\); the
    convenient \(e^{-cj^2}\) rate follows at \(A\ge10\).
12. The full PDE semigroup gate has an energy-level repair.  The symmetric
    strain norm of a locked-pitch column equals the BAS edge exactly at the
    designed ring.  Its only larger log-normal strain maximum lies on the
    spectrally stable inner flank and can be removed through a hollow Gevrey
    ramp inside a strict strain ellipse.  The compact redesigned carrier
    then obeys \(\|e^{tL_U}\|_{2\to2}\le e^{\lambda_*|t|}\) for every
    Fourier sector at once; curvature changes the exponent by only
    \(O(\varepsilon)\).
13. The raw finite-curvature Gavrilov seed is not linear: its axial
    component is radial at the seed circle.  The exact pressure is quadratic
    to leading order, however, and a hollow modulation
    \(g(p/\varepsilon^2)\) avoids the nonsmooth center.  On that fixed annulus,
    tube rescaling and speed normalization make the curved field
    \(O_{C^k}(\varepsilon)\)-close to the chosen locked-pitch column for every
    fixed \(k\).  This supplies the local \(O(\varepsilon p)\) coefficient
    ledger and \(O(\varepsilon)\) strain-edge shift; the exact Piola
    commutator converts it into an aspect-uniform ambient estimate.
14. The exact coalescing dilation produces the fixed limiting symbol
    \(A_*(\Xi^2-1)+ic_*X\) with nonzero bracket and effective parameter
    \(s/\eta^3\).  Uniform holomorphic coefficient bounds follow from the
    unchanged log-normal core; an explicit eikonal phase, transport
    recursion, and optimal truncation give \(e^{-c\eta^3/s}\).  Exact
    pressure/velocity reconstruction has only polynomial loss.  The thin
    axial period requires the corrected winding
    \(m=\operatorname{round}(\beta p/\varepsilon)\); its rounding error is
    lower order, and adjoining the conjugate charge gives a real periodic
    pseudoorbit.
15. The all-order phase combinatorics are also benign for the explicit
    three-phase choice.  Transported coordinate phases obey
    \(\nabla\Phi_k=F^{-T}k\), so every nonzero integer charge has a uniform
    lattice gap under bounded deformation.  There are only \(O(r^3)\)
    charges through order \(r\), and their exact convolution stays within
    \(C^r(r!)^2\).  The zero charge is now isolated as the sole
    non-elliptic channel carrying positive stress and the global wake.
16. On the periodic branch, that zero charge has an explicit global
    symmetric anti-divergence of order minus one.  Zero mean is its only
    compatibility condition, and the inverse preserves the same Gevrey-2
    Fourier norm without a new factorial.  What cannot be discarded is its
    terminal heat state: erasing the wake requires a separate weighted
    time-moment condition for every nonzero spatial Fourier mode.  The
    remaining zero-charge problem is therefore a nonlinear endpoint
    carry-forward/controllability theorem, not an elliptic inverse.
17. The viscous endpoint jet itself now has a canonical target.  Start the
    exact unforced Navier--Stokes evolution from the steady Euler bubble
    instead of holding that bubble fixed.  On a uniform fraction of one
    turnover, the resulting collar differs from the bubble and dissipates
    only by \(O(\mathrm{Re}^{-1})\); its full forward jet is generated by
    an explicit recurrence.  The open task is to match the active
    three-phase transition to that growing-order Gevrey jet.
18. The collar's growing-order factorial cost also fits the existing
    window.  An order-\(n\) term with \(q\) viscous insertions has exactly
    \(n+q\) spatial derivatives, whose extra Gevrey-2 price is at most
    \((4n^2)^q\).  Thus all endpoint jets through \(M_j\) retain a
    \(C^n(n!)^2\) bound when
    \(\mathrm{Re}_j^{-1}M_j^2\to0\), a consequence of the already stronger
    carrier heat gate.  Replacing \(M_j\) by \(2M_j+2\) changes the heat
    parameter by a ratio tending to four (uniformly at most \(16\) for
    \(M_j\ge1\)).  The remaining splice is an operator construction, not a
    scalar factorial obstruction.
19. The one-phase two-harmonic covariance bath is compatible with the
    rank-five endpoint chart at pulse onset.  All new bath--partner
    coefficients are \(O(1)\), not \(O(K)\), and short sequential pulses
    preserve the finite-dimensional Jacobian.  But an exact extreme-charge
    determinant proves that no nonzero finite Fourier corrector can make the
    linearized output purely zero charge, even if a zero-charge corrector is
    admitted.  The viable architecture must
    retain an infinite or growing charged wake; harmonic separation is not a
    finite-cell cure.  The resulting infinite ladder is nevertheless a
    bounded nearest/next-nearest shift on fixed-slow-sector analytic/Gevrey
    charge spaces, uniformly in carrier and cutoff.  Its Dyson tail at depth
    \(M_j\asymp j^2/\log j\) is \(e^{-(1-o(1))cj^2}\).  The correct forward
    endpoint must retain its charged wake: backward viscous erasure is not a
    bounded all-charge two-ended map on any finite-radius Gevrey scale.
20. Literal whole-space affine wake export has a finite-energy barrier.  For
    every material set \(A\),
    \[
      \int_A|X_3(T,a)|\,da
      \le\int_A|a_3|\,da
       +|A|^{1/2}\int_0^T\|u(t)\|_2\,dt.
    \]
    Uniformly bounded energy therefore forbids exponential export of a fixed
    positive volume over arbitrarily long gain windows.  The remaining
    spatial alternatives are slower separation with a better stretching
    ledger, super-small stage volume, prepositioned storage, or Eulerian
    multipole cancellation.

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

The finite central rectifier is not a surviving shortcut.  Reality forces an
infinite sideband ladder, and no nonzero finite ladder closes exactly.  A
second correction is that a lateral carrier envelope destroys exact
pressure-purity of the common high sum; its leading Leray remainder is small
but must be retained.

Nor may the certified Batchelor column simply be bent into the standard
Gavrilov bubble.  The latter has straight-limit relation
\(W_G=V_G/\sqrt2\), while at the Batchelor AO ring
\(V_B'>0>W_B'\), with
\(-W_B'/V_B'\in(2.95,2.96)\).  This order-one first-jet mismatch survives
uniform rescaling, axial reflection, and an axial Galilean shift.  A separate
compatible base-profile certificate therefore had to precede any
\(O(M^{-1})\) curvature estimate.

That mismatch has a principal-symbol repair: redesign the common Gavrilov
pressure modulation rather than importing the Batchelor profile. The
explicit locked-pitch carrier has the required unique full BAS edge, and
the scalar edge-tracking normal form retains the threshold \(A>3g\).
An adversarial audit nevertheless found two errors in the claimed physical
lift. The raw Gavrilov axial component is radial, not linear at the seed
circle, and the thin torus has axial Fourier lattice
\(\varepsilon\mathbb Z\). Both admit explicit corrections: the hollow
multiplier gives the required annular \(C^k\) comparison, and the winding
integer is \(m\asymp p/\varepsilon\). The former global Hodge-projector
comparison on a changing long tube is withdrawn. The spectral module now
has an exact straight velocity--pressure lift and an aspect-uniform compact
Piola residual lemma in the ambient domain.  Its remaining gap is the
uniform coalescing-edge analytic scalar construction and the gain-time
parameter balance; it is not an internally closed carrier chain.

## Current load-bearing theorem: one \(A_2\) cell

The program is now frozen on one geometry. C107–C136 were reconstructed
from the surviving registry and auditor notes, then re-derived with
dependency-free checkers; they are not represented as byte-for-byte recovery
of the lost drafts.

The exact finite-dimensional and scalar-ledger layer now contains:

1. an equal-shell integer \(A_2\) hexagon and a dual-helicity paired gate
   that cancels its first reality sideband;
2. an unavoidable named \(\sqrt6\)-mode at second Picard order once the
   terminal child is active;
3. exact pump-depletion orbits and the weighted ladder polynomial
   \[
   \lambda^6-9\lambda^4+18\lambda^2-9,\qquad
   {633\over250}<\sigma_*<{2533\over1000};
   \]
4. an exact unforced heat-decaying Beltrami pump with ideal gain
   \[
   G(t_*)={d\over\delta}(R-1-\log R);
   \]
5. the factorial schedule
   \[
   \ell_j=(j!)^{-8},\quad a_j=(j!)^{10},\quad
   q_j=(j+1)^8,\quad b_j=(j+1)^{-2},\quad
   F_j=q_j^{3/2},
   \]
   on which the recorded scalar budgets are summable.

The remaining package is **backward-weighted active-focus leakage
(BAFL)**. Put \(n=j+1\). The one-cell trajectory must realize the unforced
depletion/gain blocks, route the dormant seed before activation, actively
concentrate its energy by \(q_j^{3/2}\), and satisfy the two-channel response
bound

\[
 \mathfrak L_j^{\rm act}\le Cn^{-6},\qquad
 \mathfrak L_j^{\rm wake}\le Cn^{-4},\qquad
 L_j\le Cn^2,
\]

so that \(L_j\mathfrak L_j^{\rm act}=O(n^{-4})\).  Constructing the
localized trajectory is part of the BAFL package; the displayed inequalities
are the BAFL estimate on that trajectory.

The split is load-bearing. C140--C141 now close the raw fixed-Fourier
timing calculation.  The C116 gate is exactly 2D3C: its planar pump
heat-decays autonomously and its written \(N=(1,1,1)\) component is a
passive scalar.  A child of size \(b_j=n^{-2}\) creates a nonzero
second-Picard wake of size \(b_j^2=n^{-4}\), but that wake first returns
to the active child after one additional interaction, at
\(b_j^3=n^{-6}\).  Phase tuning within the whole cancelled four-mode
family cannot remove the cubic coefficient.

Thus the homogeneous short gate matches the two raw powers reserved by
BAFL.  It also exposes a necessary conversion problem: every evolution
of smooth periodic data confined to the fixed planar \(A_2\) root lattice
is globally regular 2D3C; within that split the passive scalar cannot
deplete the autonomous planar pump or supply the scheduled
three-dimensional volume focus.  This does not cover the intended
off-plane C118/C119 leaves.  The localized conversion-exposure estimate
(LCE) asks whether the construction can exit the passive-scalar plane and
perform the active focus without directly promoting the \(O(n^{-4})\)
wake into the active chart above \(O(n^{-6})\).  LCE does not replace full
BAFL, which must also control newly generated leakage, nonlinear wake
feedback, and closure of the entrance and retained-wake state classes.

The active-focus audit C142--C148 now rules out the material-envelope
shortcut and identifies the precise failure of the generic stationary
cutoff estimate.  A global
affine Kelvin map can amplify the child by \(q^{3/2}\) while keeping the two
named wake velocities bounded, but a transported three-dimensional
envelope is heat-killed at action \(\asymp q^2\log n\); a stationary
finite-energy cutoff instead has backward-focused parent/collar exposure
\(O(\log h)\) under the direct product bound.  Slowing the strain changes
neither the heat action nor that bound, but the latter is not a universal
lower-bound no-go.  The affine branch survives only conditionally on a
co-moving Lagrangian collar whose integrated wake-to-active block is
\(O(b)\).

There is an exact coherent packet with \(q^3\) conjugate pairs
(\(2q^3\) nonzero Fourier wavevectors) on the existing
\(\mathbb ZN+\mathbb Zk_1+\mathbb Zk_2\) lattice, but it repeats in every
parent period cell and has order-one radial bandwidth.  The scheduled
localized point focus needs \(\asymp(qK)^3\) Fourier degrees, whereas an
exact or narrow C121 shell has only \(O((qK)^2)\).  Thus the next pump cannot
be the exact global six-root orbit; it must be a localized C121-like core
with a non-shell collar.  A direct quadratic collapse of the coherent
packet also needs a macroscopically large gate/companion family.

The surviving gain target is LBRG: a localized broad-band relative-gain
theorem for the exact variable-fiber axial operator C148, including
eikonal rephasing, a positive amplitude cocycle, packet edges, nonlinearity,
and C125 relative to the \(n^{-28}\) growing seed.  The scalar
microseed/reservoir factorization and conditional unfolded-gate exponent
ledger are compatible with the budgets, but they do not prove LBRG or
BAFL.

No localized BAFL trajectory is constructed here. Consequently there is no
unforced Navier–Stokes stage map, cascade theorem, or Millennium conclusion.

## One-cell artifacts

- [Transported-writer reconcile, Clebsch, and Weber (reconstructed C107–C113)](research/2026-08-05-transported-writer-reconcile-c107-c113.md)
- [Terminal A2 paired gate and second-Picard leakage (reconstructed C114–C117)](research/2026-08-05-terminal-triad-hexagon-c114-c117.md)
- [A2 depletion and enclosed ladder eigenvalue (reconstructed C118–C120)](research/2026-08-05-hexagon-depletion-eigenmode-c118-c120.md)
- [Unforced decaying pump and gain gate (reconstructed C121–C125)](research/2026-08-05-unforced-decaying-pump-c121-c125.md)
- [Factorial one-cell schedule (reconstructed C126–C131)](research/2026-08-05-factorial-stage-schedule-c126-c131.md)
- [Just-in-time active focus and BAFL (reconstructed C132–C136)](research/2026-08-05-active-focus-activation-c132-c136.md)
- [One-cell BAFL reduction and explicit A2 leakage](research/2026-08-05-one-cell-stage-map-obstruction.md)
- [Exact 2D3C paired gate and cubic wake return (C140–C141)](research/2026-08-05-paired-gate-passive-return-c140-c141.md)
- [Affine child/wake selector and localization obstruction (C142–C143)](research/2026-08-05-affine-selector-localization-c142-c143.md)
- [Coherent packet, shell endpoint no-go, and relative-gain gate (C144–C148)](research/2026-08-05-coherent-packet-relative-gain-c144-c148.md)
- [Independent audit of reconstructed C107–C136](audit/AUDIT-2026-08-05-c107-c136-reconstruction.md)
- [Exact Gavrilov wake and packed carrier](research/2026-07-29-gavrilov-dss-wake-construction.md)
- [Localized Gavrilov transition ledger](research/2026-07-29-gavrilov-active-transition-ledger.md)
- [Two-colour rank theorem and three-beat repair](research/2026-07-29-two-colour-endpoint-rank.md)
- [Viscous endpoint-jet obstruction](research/2026-07-29-gavrilov-viscous-endpoint-jet-obstruction.md)
- [Forward one-phase inverse audit](research/2026-07-29-one-phase-cell-inverse-audit.md)
- [Polynomial-carrier Gevrey ledger](research/2026-07-29-polynomial-carrier-ledger.md)
- [Flat-force/Borel attack](research/2026-07-29-flat-force-borel-attack.md)
- [All-support Laurent pole no-go](research/2026-07-29-laurent-null-pole-no-go.md)
- [Embedded-quasimode nonlinear-instability gate](research/2026-08-01-embedded-quasimode-instability-gate.md)
- [AO long-gain quasimode ledger](research/2026-08-02-ao-long-gain-quasimode-ledger.md)
- [Explicit AO Batchelor global-b certificate](research/2026-08-02-ao-batchelor-global-bas-certificate.md)
- [Full Batchelor BAS cocycle and transverse counterexample](research/2026-08-02-ao-batchelor-full-bas-cocycle.md)
- [Full-edge-matched Batchelor AO profile](research/2026-08-02-ao-batchelor-full-edge-matched-profile.md)
- [Uniform finite-time BAS propagator](research/2026-08-02-bas-uniform-finite-time-propagator.md)
- [Batchelor--Gavrilov profile compatibility gate](research/2026-08-02-batchelor-gavrilov-profile-compatibility.md)
- [Locked-pitch Gavrilov full-edge carrier](research/2026-08-02-locked-pitch-gavrilov-edge-profile.md)
- [Compact-ring profile flexibility audit](research/2026-08-02-compact-ring-profile-flexibility-audit.md)
- [Corrected positive-stress chart](research/2026-08-01-three-phase-positive-stress-chart.md)
- [Affine-core viscous capture](research/2026-08-01-affine-core-viscous-capture.md)
- [Reality-sideband/Bloch audit](research/2026-08-01-reality-sideband-bloch-audit.md)
- [Corrected multicolour core-buffer audit](research/2026-08-01-multicolour-core-buffer-audit.md)
- [Characteristic-envelope pressure ledger](research/2026-08-02-characteristic-envelope-pressure-ledger.md)
- [Thin-torus exterior-wake uniformity](research/2026-08-02-thin-torus-exterior-wake-uniformity.md)
- [Thin-torus third-jet pressure](research/2026-08-02-thin-torus-third-jet-pressure.md)
- [Integrated transition checkpoint](research/2026-08-02-integrated-transition-checkpoint.md)
- [Recovered forced-subedge and curved pressure block](research/2026-08-03-recovered-subedge-curved-block.md)
- [Edge-tracking coalescing pseudomode](research/2026-08-03-edge-tracking-coalescing-pseudomode.md)
- [Strain-capped locked-pitch semigroup](research/2026-08-03-strain-capped-locked-pitch-semigroup.md)
- [Corrected hollow Gavrilov conormal comparison](research/2026-08-03-gavrilov-normalized-taylor-comparison.md)
- [Three material phases and the all-order lattice majorant](research/2026-08-03-three-material-phase-lattice-majorant.md)
- [Periodic zero-charge routing and terminal wake](research/2026-08-03-periodic-zero-charge-wake-routing.md)
- [Exact Navier--Stokes endpoint collar](research/2026-08-03-exact-navier-stokes-endpoint-collar.md)
- [Endpoint-collar Gevrey jet majorant](research/2026-08-03-endpoint-collar-gevrey-jet.md)
- [One-material-phase pressure-gauge boundary chart](research/2026-08-03-one-material-phase-pressure-gauge-chart.md)
- [Localized fixed-cohomology material-phase core](research/2026-08-03-local-material-phase-core-chart.md)
- [Exact q-free velocity--pressure lift](research/2026-08-03-qfree-velocity-pressure-lift.md)
- [Kelvin covariance-chart persistence and conditioning](research/2026-08-03-kelvin-covariance-chart-persistence.md)
- [Aspect-uniform compact Piola packet transport](research/2026-08-03-aspect-uniform-piola-packet-transport.md)
- [Finite endpoint-jet interpolation and order doubling](research/2026-08-03-endpoint-jet-interpolation.md)
- [Nonlinear periodic wake carry and scale-shift gate](research/2026-08-03-nonlinear-periodic-wake-carry.md)
- [Affine-cone wake export and finite-energy localization gap](research/2026-08-03-affine-cone-wake-export.md)
- [One-phase two-harmonic sideband coupling and finite-support no-go](research/2026-08-03-one-phase-two-harmonic-sideband-coupling.md)
- [Uniform analytic coalescing-edge scalar packet](research/2026-08-03-uniform-analytic-coalescing-scalar.md)
- [Long-window Kelvin passive no-go and terminal relaunch gate](research/2026-08-03-kelvin-terminal-relaunch-gate.md)
- [Infinite two-harmonic charge-ladder semigroup](research/2026-08-03-two-harmonic-charge-ladder-semigroup.md)
- [Audit of the claimed 2026 Navier--Stokes profile](research/2026-07-29-shahmurov-ns-preprint-audit.md)
- [Audit of the claimed axisymmetric reduction](research/2026-07-29-shahmurov-axisymmetric-reduction-audit.md)

## Reproducible algebra

No extra dependencies are needed:

```bash
for f in checks/*.py; do python "$f"; done
```

The current dependency-free checkers all pass.  Every `PASS` statement
concerns an exact algebraic identity or scaling inequality; none asserts
existence of the missing Navier--Stokes parametrix.
