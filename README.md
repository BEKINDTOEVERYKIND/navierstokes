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

## 2026-08-30 frontier correction

The repository no longer describes any implementation as “the one realistic
path left.”  C195 certifies strictly mapped forward and sign-reflected inverse
coefficient cones on an explicit finite-horizon \(C^0\) off-ray tube, with
gain above \(3000\) in both directions.  It does not construct an invariant
or canonical stable/unstable bundle.  C196 then identifies the actual
endpoint discriminator:

* exact real periodic divergence-free profiles with fixed-energy
  \(q^{3/2}\) concentration exist kinematically on a fixed aperture;
* C180's retained band supports at most
  \(8\delta^{3/2}q^{3/2}/J^2\) concentration one-sided; and
* C194 controls one beam, so a \(q^2\)-carrier endpoint needs a new uniform
  Fourier-integral almost-orthogonality and band-retention theorem.

The active portfolio now tests four distinct mechanisms: the \(A_2\)
same-witness cascade, a direct full-wake nonautonomous Navier--Stokes
invariant graph, a nonaxisymmetric Euler-dominant Type-II inner-profile
discovery that still requires certified finite-energy outer matching, and
an all-order terminally-flat forced construction.  Their fixed failure
criteria are in [FRONTIER.md](FRONTIER.md), and the counterfactual success
logic is in
[research/2026-08-30-counterfactual-success-portfolio.md](research/2026-08-30-counterfactual-success-portfolio.md).
None is a claimed singular solution.

## Historical modular-route chronology (C1--C194)

The numbered body below is retained as a provenance record of the earlier
module-by-module program.  Its phrases such as “remaining,” “next,” and
“the route” describe the local frontier at the checkpoint where each
paragraph was written; they are not the current architecture verdict and do
not assert that any module was globally last or unique.  The authoritative
current state is the 2026-08-30 correction above together with
[FRONTIER.md](FRONTIER.md).  Within that historical program, several objects
that were previously only formal became exact.

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

At that checkpoint, the modular program grouped its unresolved package
under **backward-weighted active-focus leakage (BAFL)**. Put \(n=j+1\). The
one-cell trajectory was required to realize the unforced
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

C171 separates the two parts of that condition.  Transporting a compact
curl potential as a one-form gives the exact divergence-free Piola field
\(U\circ X=Fw_0\), whose parent cross-residual is
\(2(U\cdot\nabla)V\), not zero.  Conditional on its full Piola/F-jet factor
obeying \(M_F\log n=o(n^2)\), the normalized backward parent-cross scale is
\(M_Fq^{-1}\log h=o(n^{-6})\), and one \(O(n^2)\) chart charge is
\(o(n^{-4})\).  A rigid center-following collar already has the same
parent-cross arithmetic, but still owes boundary-crossing control.  Leray
projection is \(L^2\)-contractive while its pressure tail remains global.
None of these facts gives the retained wake an automatic extra factor
\(b\): C171 exhibits a bright real zero-charge \(A_2\) wake-to-child triad
and a bright allowed Piola ordering.  The actual C140 fixed-projector wake
does gain one interaction and returns through C141's nonzero cubic term, but
a localized moving block still needs a quantitative kernel or export
theorem.

C172 closes exact affine pressure darkness as a local substitute for that
theorem on the C142 selector.  For a constant matrix, the dark compression
has determinant \(k^T\operatorname{adj}(A)k/|k|^2\).  Its rank-two
characteristic set is the union of the left- and right-kernel planes, so no
nonzero \(L^2(\mathbb R^3)\) localized field can satisfy
\({\mathbb P}(AU)=0\) when \(\operatorname{rank}A\ge2\).  The fixed
trace-free rank-one curl-kernel exception has no Piola amplitude gain.  In
the displayed C142 launch frame, the C140 child and named wakes are bright;
making all three \(N\)-fibers dark for one constant selector forces
\(AN=0\) and removes their Kelvin velocity gain.  Universal zero-order
affine multipliers reduce to scalar Piola, while a universal
translation-invariant finite-order velocity rule cannot replace the
rational Hodge symbol on an open frequency cone.  This last statement does
not cover variable-coefficient auxiliary systems or mode-tailored
antidivergences.

C173 removes a finite-dimensional portion of the collar problem without
claiming the broad-band theorem.  Two concentric symmetric affine curl
collars with weights \(128/127\) and \(-1/127\) retain the full inner affine
strain while canceling the degree-three moment.  Their transform vanishes
to at least fifth order at low frequency, giving a \(q^{-2}\) improvement
over one collar on fixed chart-tracked parent-scale coordinates at fixed
core normalization.  An order-12 material-label polynomial notch, with
coefficient \(\ell^1\)-norm \(343/64\), kills the named squared radii \(2\)
and \(6\) through two Fourier derivatives and hence kills the formal first
curvature variation there.  The same robust construction on \(q\) distinct
axial labels needs degree at least \(3q\) (radial differential order at
least \(6q\)), outside the current Gevrey budget.  The full \(L^2\)
residual, actual C140 profile normalization, moving broad bundle, curvature
remainder, and periodic pressure tails remain in MCKC/BAFL.

C174 gives a complementary prescribed time-slicing ledger.  With
\(J=O(\log h)\) bounded-variation rigid-collar slices, its switches,
parent-cross term, self term, factorially small viscosity, and
trajectory-separated instantaneous pressure tails fit the scheduled
powers under the note's stated hypotheses.  Exact recharting telescopes to
the same physical child cocycle; it preserves Floquet multipliers only for
matching endpoint charts and does not erase C154's physical shear.  At the
central first-jet level, an isotropic \(J^{-1}\) residence shrink is
sufficient at return times, but costs \(J^3\) in reciprocal volume,
\(J^{3/2}\) in fixed-energy amplitude, and \(J\) in bandwidth; finite-tube
second jets and between-return residence are open.  A physical reset with
relative mismatch \(\delta_m\) costs \(\delta_m b\) in each individual
terminal-weighted norm, so a separated triangle-inequality proof needs
\(\sum_m\delta_m=O(b^2)\).  A signed Duhamel cancellation is not excluded.
This is the rechart-or-rewrite obstruction (RRR), not an unforced stage.

C176 closes the exact-background residence question left in that historical
C174 ledger.  Regular action--angle--axial coordinates give first and second
flow jets \(O(J)\) and \(O(J^2)\) at every intermediate time.  A correlated
packet fits in an \(O(Jr)\) tube, while a fixed-aperture
\(q\times q\times(q/J)\) carrier slab contains
\(\Omega(q^3/J)\) reality pairs and stays at \(O(q)\) frequency under the
central-fiber C154 shear.  An enclosing spherical collar costs only
\(q^{-1}\operatorname{polylog}q\) in the backward ledger.  This closes
background residence and arithmetic capacity, not localized finite-frequency
propagation or the unforced collar correction.

C177 and C178 attack that unforced correction from complementary sides.
A same-curl Beltrami reservoir can be preloaded into the pump as an exact
unforced heat-decaying background: it contributes a homogeneous tangent and
no additive writer.  Different curl values instead create the exact cross
source \((\kappa_1-\kappa_2){\mathbb P}(U\times R)\), while a common shell's
single heat envelope supplies no polarization timing.  For the collar, a
fixed compact monotone ramp has no exact \(L^2\) backward-heat preparation,
but a forward analytic buffer has an explicit contractive entrance datum
with exact pure-heat terminal cancellation.  Local Laplacian-polynomial
preparations also retain compact entrance support and converge terminally.
These are genuine autonomous/scalar repairs; neither proves their required
full \(A_2\), C125, pressure-tail, or nonlinear propagation estimates.

C179--C184 isolate the terminal-reservoir branch without promoting an
instantaneous Fourier edge to a stage map. The exact 2D3C background
\(v_{A_2}+\Theta N\) supports arbitrary planar passive frequencies and its
off-plane one-edge Leray symbol is uniformly two-polarization invertible
away from the exact normal-charge and equal-radius defects. The factorial
C121 shell is richer still: it contains a reality-complete octahedral
palette of \(q+O(1)\) modes whose normalized forward rows form an almost
isometric full-polarization star on a narrowed \(q^2\)-source C176 slab,
with only constant-cost first-collision pruning. Its leading term,
however, is scalar transport and its normalized strain is only
\(O(q^{-1/2})\). C181 closes the corresponding affine principal test: the
actual Euler Kelvin multiplier of a terminal vertical shear is uniformly
bounded on the fixed C159/C176 input cone for arbitrarily long shear action.
The missing \(\sqrt q\) gain appears only when that input cone degenerates
at least as fast as \(q^{-1/2}\). The remaining finite-frequency velocity
question is narrowed further by C182. On a fixed \(O(1)\) normalized
interval, any complete linearized C180 propagator whose entrance packet
obeys the explicit unnormalized \(L^2\), \(L^\infty\), and scale-\(q\)
\(W^{1/4,16}\) tail bounds has endpoint at most \(Cbq^{9/8}\). This is
short of the raw coherent target by \(q^{-3/8}\), and remains short after
the C180 \(J^2\) normalization tax. The result is conditional: C176/C180
do not construct those entrance/tail bounds, and the large-torus Hölder and
generic \(B^0_{\infty,1}\) shortcuts are invalid.

C183 makes the long 2D3C evolution exact in Lagrangian coordinates without
turning it into an endpoint deformation formula. For each vertical Fourier
charge, a unit-modulus gauge removes scalar transport and gives
\(D^Tp=p_0-mt g_0\) together with the exact Piola constraint. The remaining
Kelvin velocity is a time-ordered two-polarization connection. At a planar
return with \(L=D_T^T\), a common periodic coefficient requires
\(Lg_0=g_0\) and \(Lp_0=p_0-mTg_0\); unipotent defects generate explicit
linear or quadratic covector drift. Stationary C152 scalars \(H(f)\) do
admit returning covectors after the exact retuning
\(u_h\cdot p_0=-mTc\). If a common periodic frame is constructed and its
connection differs from C159 by \(\varepsilon_q=O(q^{-1/2})\) per return,
the extra factor over \(O(\log q)\) returns is only \(1+o(1)\). C180's
spatial coefficient ledger does not prove this common-frame hypothesis,
and a broad nonstationary PPRG scalar need not satisfy the return resonance.

C184 tests the most natural return-coherent PPRG repair, a stationary
polynomial \(\Theta=H(f)\). Its diagonal Fourier coefficients are genuine
triangular coordinates and can produce a normalized full-polarization
tangent row, so there is no selected-row rank obstruction. The complete
polynomial is fatal instead: an outer diagonal coefficient \(\tau_q\)
forces an off-line binomial face of coefficient-\(\ell^2\) size at least
\(|\tau_q|\{(9/4)^q-1\}/\sqrt q\). At the required star scale this is
\(\gtrsim b(9/4)^q/q^2\), exponentially above every polynomial reservoir
or wake budget, and the high polynomial parameters have zero first jet on
the C159 zero loop. This **stationary-polynomial collateral explosion
(SPCE)** closes only the degree-\(q\) polynomial line-palette shortcut.
Nonpolynomial profiles, deliberate recycling of the collateral face, and
non-returning time-ordered PPRG remain open.

C185 inverts the proof architecture and lands the first explicit positive
infinite-dimensional PDE-operator estimate in this branch.  Combining
C159's independently certified periodic cone with the abstract-level
Shvydkoy spectral inclusion gives the robust operator-norm statement

\[
 \|G_{nT}\|_{2\to2}\ge e^{n/5}>\left({6\over5}\right)^n.
\]

C189's independent cross-audit holds the stronger
\(r_{\rm ess}(G_T)\ge e^{1/5}\) citation pending a paper-body check of
Shvydkoy Theorem 4.1.  It is not used downstream until that check lands.

This is genuine infinite-dimensional inviscid operator growth, not a
Galerkin or frozen-symbol statement.  It does not give a viscous scale
return.  The accompanying fixed-point audit maps the program to
Albritton--Bru\'e--Colombo, Elgindi, and Chen--Hou and isolates the one
component with no completed analogue: an **unforced viscous
physical-velocity scale return (UVSR)** for the complete active-plus-wake
state.  The repository now treats BAFL, C125, and RIGM as projections of one
renormalization residual/trapping estimate rather than successive proof
architectures.  The root `FRONTIER.md` is the mandatory session boot state
and records the pre-registered PPRG kill trigger.

C186 gives a conditional abstract finite-dimensional PPRG witness, with a
sharp realization boundary.  C183's planar return is unipotent, but its
Kelvin polarization monodromy is only known to have determinant one.  If
two candidate polarization generators are square-zero,
the scalar \(\tau=\operatorname{tr}(N_1N_2)\) vanishes exactly when their
unipotent episode maps share a flag.  The explicit transverse pair
\(I+E_{12},I+E_{21}\) has per-episode exponent greater than \(12/25\), and
entrywise \(1/100\) perturbations in \(SL(2)\) retain exponent greater than
\(9/20\).  A smooth bounded exact passive 2D3C Euler example also rotates a
scalar gradient through transverse directions, so passivity alone cannot
force gradient alignment; it says nothing by itself about invariant lines
of the Kelvin connection and does not realize the two polarization blocks.
What remains is the load-bearing PDE test: realize two
such robust polarization blocks along one unforced passive orbit in a
common C183 return fiber and validate the finite-frequency evolution.

C187 reconciles the two outstanding audit items.  Ordinary \(H^3\)
linearized Navier--Stokes Duhamel continuity holds with the explicit
constant

\[
 44\sqrt{2T/\nu}\exp(7744V^2T/\nu),
\]

but this explodes on the factorial-viscosity schedule and cannot prove the
old abstract (5.1); that display is now folded into the full trapping
premise until its structured state space is defined.  The session-only
infinite-ladder values \(0.66855\ldots\) and \(2.63707\ldots\) are withdrawn:
no operator, tail enclosure, proof, or checker exists for them.  C120's
finite \(6\times6\) enclosure is unaffected.

C188 completes the demand-side UVSR scalar-corridor audit.  Under the
equal-normalized-shape convention, same-energy transfer into a \(q^{-3}\)
child fixes the ledger focus multiplier at \(q^{3/2}\); unequal physical
shapes contribute the explicit square root of their \(L^2\)-constant
ratio.  The optimizable quantity is the net gain \(g=bF\); in the
equal-shape power ledger, \(g=bq^{3/2}=q^\gamma\).  The bare
high-Re/finite-energy corridor has
infimum \(1\).  If the bounded-profile C176 worst-case upper majorant
\(C(1+\Lambda)q^{-1}J^{7/2}\) must remain separately below the \(b^3\)
active tolerance, its sharp polynomial infimum is \(7/6\), with a
sufficient logarithmic correction at the boundary.  This is a result for
the declared envelope, not a lower bound on the physical collar.  The
equal-shape all-sequence forms are \(g=q\rho\),
\(1<\rho=bq^{1/2}<q^{1/2}\), with
\(\log\rho/\log q\to0\) at the direct floor, in
the direct class and
\(g=K_*^{1/3}q^{7/6}(\log q)^{7/6}\omega\), with
\(\omega\to\infty\), \(\log\omega/\log q\to0\), and \(0<b<1\), for
the declared envelope, where
\(K_*=C_{\rm col}C_J^{7/2}\) and \(C_{\rm col}\) absorbs the uniform
bound on \(1+\Lambda_j\).  Arbitrary sequences must additionally satisfy
C188's exact global energy-product condition.  The scalar schedule
\(q=n^4,g=2q,b=2n^{-2}\) attains direct polynomial order one, doubles
Reynolds each return, and has total normalized energy at most \(16/3\),
but it does not reuse C180's shell theorem.  The current
\(5/4\) ledger is therefore not envelope-optimal: keeping C180's
proved \(q=n^8\) shell and taking \(b=n^{-5/2}\) gives the exact lower
scalar point \(g=q^{19/16}\), after respecifying C161 with
\(H=n^{51/2}\) and \(J_{\rm split}=\lceil n^{5/2}\rceil\).

The constraint-theory audit also corrects the proposed use of partial
regularity.  CKN epsilon regularity constrains cylinders centered at an
actual singular point; at every smooth finite-stage center the standard
scaled quantities instead have explicit \(r^2,r^3,r^4\) upper bounds and
tend to zero.  This does not control stage radii in cylinders centered at
the eventual singular point.  The local energy inequality contains
indefinite inward pressure/advective flux, and Tao's quantitative \(L^3\)
lower rate is global, subsequential, and has an unspecified exponent.  No
numerical per-stage power follows from the landed scalar ledger without
additional center, occupancy, pressure, and wake hypotheses.  Finally,
the exact law \(\mu_{j+1}=(q_j/g_j)\mu_j\) shows that a fixed supercritical
scaling has no positive-viscosity autonomous fixed point.  The actual
factorial construction is the nonautonomous augmented map
\((X,\mu,n)\mapsto(\mathcal R_{n,\mu}X,(q_n/g_n)\mu,n+1)\); an autonomous
inviscid limiting map remains open.

C190 completes the pre-registered two-episode test on C186's
rotating-gradient orbit.  Consecutive determinant-one episodes force
\(m=0\); in the common co-rotating frame both quarter maps are the same
parabolic shear, their product trace is exactly two, and full-return powers
grow at most linearly.  This is a theorem for that orbit, not a no-go for
all passive 2D3C dynamics.

C191 tests the proposed direct reuse of C185 against C182's raw
\(q^{3/8}\) deficit.  The certified exponent is scalar-power capable:
\[
R_\Delta=\left\lceil{15\over8}\log q\right\rceil
\quad\Longrightarrow\quad
\|G_{R_\Delta T}\|_{2\to2}\ge q^{3/8},
\]
with required pump action
\[
{45\over8}\log q<TR_\Delta
<{57\over10}\log q+{76\over25}.
\]
The declared scalar collar and heat budgets remain feasible under an
explicit logarithmic allocation, but no landed theorem supplies that
lower stage coverage.  More decisively, C182 is an \(L^\infty\) upper
bound while C185 is an unrestricted supremal \(L^2\) operator-norm lower
bound with no common finite band.  A common \(q^{3/8}\) multiplier overruns
child energy by \(q^{3/4}\); fixed-energy rescaling cancels it but does not
resolve C125.  In C147's coherent-packet template, uniform rescaling leaves
the divergent ratio at \((57/2)n^8\log n\).  Thus PPRG is
not obsolete: the missing object is same-witness, retained-band,
fixed-final-energy concentration, not bare norm growth.

C191 also catches a class-scope inconsistency.  The C152/C159/C185
background is itself passive 2D3C and supplies an \(m\ne0\) exact return
with growth, so universal secular no-return is false on C179/C183's broad
class.  It is periodic/resonant, however, and therefore does not satisfy
the accepted non-resonant PPRG witness form.  No third PPRG form is
introduced; its unresolved accepted target is one genuinely incommensurate
non-fixed orbit with transported return or growth and the same-witness
physical endpoint estimate.  This is not the only program route.  Failure
of a search is not promoted to a class no-go.

The same-curl reservoir is preloaded and has no autonomous heat switch, so
its possible action through a logarithmic C159 window is not covered by
C182.  C191 forbids treating the inherited \(O(\log q)\) upper charge as
actual lower residence and forbids subtracting C185's \(L^2\) exponent from
the concentration demand.  Finite-band propagation, viscosity, depletion,
C125, RIGM, BAFL, and the stage remain open.

C192 sharply strengthens the positive asset on the same certified orbit.
Instead of integrating the unstable amplitude column, it bounds the raw
cooperative generator on 2048 cells and propagates a directed positive
lower product. For \(w=(1,3/20)^T\),
\[
                 Mw>3000w,
 \qquad \|G_{rT}\|_{2\to2}\ge3000^r>e^{8r},
\]
using only C189's approved operator-norm source scope.  The raw deficit now
uses
\[
R_\Delta=\left\lceil{3\over64}\log q\right\rceil,
\qquad
TR_\Delta<{57\over400}\log q+{76\over25};
\]
the logarithmic action coefficient is exactly forty times smaller than in
C191. Consequently a first-order \(q^{-1/2}e^{\Gamma t}\) remainder has
the nonempty power window \(\Gamma<350/57\), and an order-\(-1\) remainder
has \(\Gamma<550/57\).  C194 now proves the first local pressure-resolved
remainder with \(\Gamma=6\), but not its periodic/off-ray completion.  A
narrow-cone Gaussian at this scale also has width \(q^{-1/4}\), not the
required \(q^{-1}\) child width, and scalar growth still cancels under
fixed-energy normalization.  C192 therefore makes a direct finite-frequency
bridge plausible without claiming the physical endpoint.

C193 removes that last scalar-normalization collision at the exact
principal level.  It certifies the C159 contracting line and its physical
conditioning (projector norm \(<5/4\)), then places an explicit localized
profile in the expanding line and a broad profile in the contracting line.
After
\[
 R_{\rm filt}=\left\lceil\frac38\log n\right\rceil+1,
 \qquad q=n^8,
\]
the entrance and endpoint \(L^2\) norms agree, no discrete return exceeds
that energy, and the concentration quotient improves by more than
\((1750/251)n^3>q^{3/8}\).  This is a complex principal two-fiber lemma,
not yet a real finite-frequency solution.

C194 supplies the complementary upper-error theorem.  Its exact curl
ansatz is divergence free, its pressure cancels the leading Kelvin
equation, and the full local linearized-Euler error is bounded by
\[
\frac{\hbar e^{6t}|b_*|}{|k_0|}
\left[
4{,}199{,}040\,\varepsilon^{-1}(1+t)^3\|\nabla\chi\|_2+
2{,}898{,}006{,}000{,}000{,}000(1+t)^7\|\chi\|_2
\right].
\]
At \(\hbar=q^{-1}\), the \(q^{-1/4}\) concentration width has conditional
relative margin \(q^{-27/100}\) against a hypothetical same-packet signal;
the \(q^{-1/2}\) stress test retains \(q^{-1/50}\).  At the C194 checkpoint,
the local same-witness composition problem was to certify off-ray
forward/backward cone transport, periodic reality and retained band, and
combine the results on one solution.  C195 and C196 later separated this
into a proved finite-horizon \(C^0\) dominated-cone field, exact kinematic
endpoints, and an unproved uniform
multi-beam/band-retention theorem; these are distinct obligations rather
than one bridge.
C193's extra return pays the fixed factor \(e^{912/25}\).

There is an exact coherent packet with \(q^3\) conjugate pairs
(\(2q^3\) nonzero Fourier wavevectors) on the existing
\(\mathbb ZN+\mathbb Zk_1+\mathbb Zk_2\) lattice, but it repeats in every
parent period cell and has order-one radial bandwidth.  The scheduled
localized point focus needs \(\asymp(qK)^3\) Fourier degrees, whereas an
exact or narrow C121 shell has only \(O((qK)^2)\).  Thus the next pump cannot
be the exact global six-root orbit; it must be a localized C121-like core
with a non-shell collar.  A direct quadratic collapse of the coherent
packet also needs a macroscopically large gate/companion family.

C149 supplies a local positive amplitude cocycle on the same pump: a
nonsymmetric three-root Beltrami field has a vertical elliptic ray with
bounded periodic covector and dimensionless inertial-time Floquet exponent
\((9/16)\sqrt3(2+\delta)\epsilon+O(\epsilon^2)\).  C152 supplies a
same-geometry alternative to that ray's axial-export defect.  For
\(\delta=4/5\), the regular level \(f_\delta=0\) is one contractible closed
streamline with zero axial drift, an exact rank-one linearized return, and
a two-plane of periodic
covectors containing an explicit off-plane direction.  The corresponding
Kelvin amplitude calculation has a very large numerical multiplier
(C153), and C159 replaces that floating orientation with an
outward-rounded cooperative-cone certificate.  C192 strengthens the same
certificate to \(Mw>3000w\), hence \(\rho(M)>3000\), without integrating
an amplitude column or using the numerical trace as a premise.  C154 also
shows that a fiber band with
nondegenerate width in the shearing direction widens from \(q\) to
\(q\log q\) over the gain window.  Correlation with the periodic-covector
plane removes that shear but leaves only a two-dimensional fiber; a
three-dimensional candidate must instead absorb a logarithmic narrowing
or supply a different correlation mechanism.

C150--C151 expose the nonlinear gate.  Normal charge makes the retained
quadratic projection exactly dark and delays support-level return to degree
three, but an explicit pair of limiting unstable rays has instantaneous
wake \(QA^2\) and growing-line return \(Q^2A^3\).  C155 now closes the
zero-incoming-wake, first-period averaging loophole for the quarter pair:
after one complete limiting Kelvin period the exact selected coefficient is
\(-3.097745\ldots<-21/8\), and it retains its sign for sufficiently small
positive \(\epsilon\).  The endpoint wake is nonzero, so the parent-only
period map is not Markov and must be enlarged by a wake variable or
invariant graph.

C156 still proves that on one fixed limiting resonant ring, cancelling the
entire quadratic complement forces a single Fourier ray.  Its stronger
persistent-wake coercivity claim failed adversarial audit.  The reverse
parent does not see a separately reinitialized wake, so its coefficient is
not \(\mathfrak K(-t)\): with the one shared neutral difference wake it is
exactly \(-\mathfrak K(t)\).  At quarter separation the two coefficients are
\(-9\pi/8\) and \(+9\pi/8\), and their retained-energy secular sum vanishes
as energy conservation requires.  Thus the claimed fixed-ring scalar-phase
no-go is withdrawn.  The corrected antisymmetric normal form instead has
an exact positive balanced three-ray state, so fixed-ring nonlinear balance
is a genuine repair candidate, not merely an absence of obstruction.

C157 gives the first exact thickening test.  Two reflected decompositions
on distinct radii share one quadratic output, but their projected wake
vectors are linearly independent; their two scalar pair products cannot
cancel that output.  The determinant loses one power of the radial gap, so
adjacent \(q^{-1}\)-spaced layers can cost \(O(q)\) in the inverse.  This
does not classify outputs with more than two decompositions.  C160 shows
that unequal-radius periodic sum and difference returns do cancel in their
Cesàro long-time mean, but adjacent layers need time \(O(q)\) to detune.
That mechanism is therefore too slow for the \(O(\log q)\) gain window.

The actual \(q^3\)-mode packet therefore remains undecided.  C158 fixes the
normalization in one-representative-per-reality-pair coordinates: a bounded
diagonal cubic kernel is measured by \(Q^2A_{\ell^1}^2/M\), not by applying
the two-ray coefficient directly to C147's total point seed.  The selected
cone makes physical point amplitude comparable to this coefficient scale
within fixed factors, but not equal to it.  For \(M=q^3\) the upper
parameter is \(q^{-5}\) at the seed and \(q^{3/2}\) at the target; this is
not a coercive threshold.

C161 gives a sharper just-in-time architecture.  Carry only \(q^2\) modes
through the long gain, matching the real two-dimensional capacity suggested
by C154's periodic-covector plane, and add the third Fourier coordinate at
the endpoint with a normalized \(q\)-way pure-normal charge star.  This
dimension count is not a theorem that the exact plane contains a rank-two
torus lattice; C170 below records that arithmetic obstruction.  On the
cofinal even-\(n\) schedule, \(q\) nonzero reality-symmetric shifts give
exactly \(q^3\) first daughters.  The abstract star preserves \(\ell^2\)
and meets the conditional \(b,b^2,b^3\) target/wake/return powers.  These are
half-lattice coefficient scales; physical point coherence is not yet a
theorem.  C162 rules out the simplest tuned circular gate because its
forward/reverse product is \(-1/4+i/4\), not negative real.  C163 repairs
that phase in the selected compression with the fixed equal-magnitude
dual-helicity gate \(E_0=(e_x-e_y)/\sqrt2\), uniformly for a signed reality
pair of charge heights.  The repair is not closed: the reverse leg forces
the orthogonal source polarization.  Moreover every nonzero static
shared-intensity mixture has
a rigid \(-W\sin(2\phi)/2\) rate term, so it cannot synchronize an angular
interval.  C164 computes the complete two-polarization fiber.  A reality
pair generically shares only its tangential line, and that entire charge
ladder is an exact phase multiplier whose pointwise modulus is conserved
inviscidly and decreases in the linearized viscous scalar equation.  C165
shows that a prescribed common chirped two-level control can synchronize
populations only with \(O(J)\) action and with a generally rate-dependent
phase, outside the C161 \(O(1)\)-action ledger.

Time ordering is nevertheless algebraically active.  C166 gives an audited
nonzero radial/tangential commutator and a fixed first-neighbour Galerkin
point ratio.  C168 proves that a fixed **point** gain survives on the
complete integer charge ladder: a unit tangential pulse gives \(\sqrt2\),
and a nonzero radial-first pulse retains gain above \(13/10\) after all
repeated charge walks are included; its half-lattice coefficient energy is
bounded but is not conserved.  These are prescribed-pump linear gains, not
unforced or energy-preserving splitter theorems.  C167 rules out the static
independent-edge scalar star and its scalar-return full-fiber relaxation on
an open angular/radial sheet; cross-daughter aggregates and time ordering
remain open.  C169 supplies the full nonlinear boundary: in the
fixed 2D3C plane the total transverse field is passive.  A \(q\)-signed-mode
gate with coefficient norm \(G_t\) can produce transverse supremum at most
\(\sqrt qG_t\); the proposed \(bq^{3/2}\) endpoint requires
\(G_t\gtrsim bq\), not the scheduled \(O(b)\).  Here \(G_t\) is a
coefficient-\(\ell^2\) norm, while its energy is \(G_t^2\).  A surviving
converter of the C161 charge-star type must leave the fixed plane through a
genuinely three-dimensional noncommuting evolution.  A wholly different
in-plane focusing theorem remains logically open, but would replace rather
than realize that pure-normal source-to-bright block.  C170 also closes the
tempting aligned-packet shortcut:
the C159 periodic plane projects isomorphically to the horizontal plane,
and a \(q^{-1}\)-angular sector about the C159 ray at frequency \(O(q)\)
has only \(O(q)\) compatible integer modes.  C154 allows at most one normal
lift per projected mode through \(O(\log q)\) returns.  Moreover, the exact
plane's integer rank depends on the unresolved rationality of \(\sigma\); the rank-two
ambient projected-lattice count is only an upper bound for exact-plane
occupancy.  The required \(q^2\) packet must therefore be a genuinely
two-directional finite-frequency lattice/microlocal construction, with
C154 shear and near-plane arithmetic controlled, rather than a single
source-adapted ray.  The surviving terminal target is no longer a
pure-normal source-to-bright block: it must supply physical endpoint
coherence and depletion through an unforced, genuinely three-dimensional,
time-ordered full-polarization converter while satisfying C125/BAFL.
Its collar subproblem is now the two-clause **material-collar kernel closure
(MCKC)**: control the full Piola chart or rigid-collar boundary crossing,
and prove an \(O(b)\) integrated return for the actual retained C140 wake or
export it with uniform pressure/heat tails.  MCKC is necessary inside
LCE/BAFL; it is not a stage theorem.

C172--C176 sharpen that target without closing it.  Exact local pressure
darkness is unavailable for the rank-two selector.  C173 pays the named
finite material-label fibers but not the moving (q^2)-to-(q^3) active
bundle.  C174 replaces the rigid-collar part of MCKC(i) by the RRR fork:
either prove finite-tube residence and coherence under the accumulated
C154 shear, including second jets and between-return control, or prove a
signed physical-rewrite cancellation stronger than the separated condition
\(\sum_m\delta_m=O(b^2)\).  In either branch, the actual retained wake still
needs the broad-band Duhamel kernel/export estimate in MCKC(ii), and the
prescribed collar sources must be absorbed into one unforced trajectory.
C176 resolves the exact-background residence, second-jet, and carrier-count
portion of the first branch; base-dependent finite-frequency propagation
remains open.  C175 identifies the remaining kernel requirement as RIGM:
common material phases preserve the dangerous \(A_2\) resonance, and an
\(O(b)\)-tilted invariant graph over a specified admissible wake class already
requires an \(O(b)\) first-order Melnikov transfer.  The surviving choices are
therefore signed full-packet cancellation, \(O(b)\)-action just-in-time
exposure, or controlled spatial export--not an automatic eikonal, heat, or
graph-transform gain.

C177 removes the fatal additive-writer ratio on an exact same-curl manifold,
conditional on the actual propagated tangent satisfying C125.  Its remaining
same-curl reservoir realization gate (SCRG) is to recover the full physical
polarization and charge action from one static shell, or cancel every
noncommon-curl cross source.  C178 separately shows that prescribed collar
ramps can be absorbed exactly at the pure-heat endpoint after analytic
buffering, and approximately by compact local entrance data.  The tuned C142
rate makes that scalar preparation smaller than the \(n^{-28}\) seed, but the
full \(A_2\) evolution and RIGM/MCKC response are still the load-bearing
estimates.

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
- [Elliptic A2 gain cocycle and derivative-amplified cubic return (C149–C151)](research/2026-08-05-a2-elliptic-gain-c149-c151.md)
- [Zero-drift A2 orbit and covector shear (C152–C154)](research/2026-08-05-zero-drift-a2-orbit-c152-c154.md)
- [Floquet-averaged DACR and fixed-ring cancellation boundary (C155–C156)](research/2026-08-05-floquet-averaged-dacr-c155-c156.md)
- [Two-radius additive quartet and packet normalization (C157–C158)](research/2026-08-05-two-radius-quartet-normalization-c157-c158.md)
- [Certified zero-drift cooperative Kelvin cone (C159)](research/2026-08-05-zero-drift-cooperative-cone-c159.md)
- [Unequal-radius detuned DACR cancellation (C160)](research/2026-08-05-two-radius-detuned-dacr-c160.md)
- [Normalized terminal charge-star ledger (C161)](research/2026-08-05-normalized-terminal-charge-star-c161.md)
- [Tuned circular normal-gate obstruction (C162)](research/2026-08-05-circular-normal-gate-star-no-go-c162.md)
- [Equal-magnitude dual-helicity compressed gate and synchronization obstruction (C163)](research/2026-08-05-paired-polarization-charge-star-c163.md)
- [Full-fiber pure-normal gate classification (C164)](research/2026-08-06-full-fiber-pure-normal-gate-c164.md)
- [Charge-star synchronization and prescribed chirped control (C165)](research/2026-08-06-charge-star-synchronization-control-c165.md)
- [Reality-paired prescribed Galerkin polarization pulse (C166)](research/2026-08-06-reality-pair-polarization-pulse-c166.md)
- [Static reality-paired gate-bundle obstruction (C167)](research/2026-08-06-reality-paired-gate-bundle-rank-no-go-c167.md)
- [Complete prescribed-pump charge-ladder pulse (C168)](research/2026-08-06-complete-charge-ladder-pulse-c168.md)
- [Fixed-plane passive-star resource obstruction (C169)](research/2026-08-06-fixed-plane-passive-star-no-go-c169.md)
- [Periodic-plane aligned-sheet capacity obstruction (C170)](research/2026-08-06-periodic-plane-aligned-sheet-no-go-c170.md)
- [Co-moving Piola collar and MCKC boundary (C171)](research/2026-08-06-comoving-piola-collar-c171.md)
- [Affine pressure-dark material-transport obstruction (C172)](research/2026-08-06-affine-pressure-dark-material-transport-c172.md)
- [Paired multipole Piola collar and finite-shell notch (C173)](research/2026-08-06-paired-multipole-piola-collar-c173.md)
- [Time-sliced rigid collar and rechart-or-rewrite obstruction (C174)](research/2026-08-06-time-sliced-rigid-collar-c174.md)
- [Lagrangian resonance and invariant-graph Melnikov gate (C175)](research/2026-08-13-lagrangian-resonant-graph-obstruction-c175.md)
- [Correlated A2 tube residence and carrier capacity (C176)](research/2026-08-13-a2-correlated-tube-residence-c176.md)
- [Same-curl unforced reservoir and cross-curl gate (C177)](research/2026-08-13-same-curl-unforced-reservoir-c177.md)
- [Compact and analytic-buffer heat preparation (C178)](research/2026-08-13-compact-heat-preparation-c178.md)
- [Planar passive reservoir and full-rank one-edge chart (C179)](research/2026-08-13-planar-passive-reservoir-c179.md)
- [Factorial-shell full-polarization transport star (C180)](research/2026-08-14-factorial-shell-transport-star-c180.md)
- [Static vertical-shear Euler Kelvin propagator (C181)](research/2026-08-14-static-vertical-shear-kelvin-c181.md)
- [Conditional finite-p terminal no-focus theorem (C182)](research/2026-08-14-holder-terminal-no-focus-c182.md)
- [Exact 2D3C Lagrangian gauge and common-Floquet return obstruction (C183)](research/2026-08-14-2d3c-lagrangian-gauge-c183.md)
- [Stationary polynomial star and exponential collateral obstruction (C184)](research/2026-08-14-stationary-polynomial-star-c184.md)
- [Fixed-point literature map and operator growth (C185)](research/2026-08-23-fixed-point-literature-map-c185.md)
- [Conditional unipotent PPRG dichotomy (C186)](research/2026-08-23-unipotent-pprg-dichotomy-c186.md)
- [Explicit Duhamel constant and ladder audit (C187)](research/2026-08-23-duhamel-and-ladder-audit-c187.md)
- [UVSR schedule corridor and epsilon-regularity audit (C188)](research/2026-08-25-uvsr-schedule-corridor-c188.md)
- [Adversarial audit of C188](audit/AUDIT-2026-08-25-c188.md)
- [Independent cross-audit of C185--C187 (C189)](audit/AUDIT-2026-08-25-c185-c187.md)
- [Rotating-gradient return obstruction (C190)](research/2026-08-25-rotating-gradient-return-obstruction-c190.md)
- [Adversarial audit of C190](audit/AUDIT-2026-08-25-c190.md)
- [C185 deficit accounting and passive-class correction (C191)](research/2026-08-27-c185-deficit-accounting-c191.md)
- [Adversarial audit of C191](audit/AUDIT-2026-08-27-c191.md)
- [Strong zero-drift gain and short finite-frequency clock (C192)](research/2026-08-28-strong-zero-drift-gain-c192.md)
- [Adversarial audit of C192](audit/AUDIT-2026-08-28-c192.md)
- [Hyperbolic fixed-energy polarization filter (C193)](research/2026-08-28-hyperbolic-polarization-filter-c193.md)
- [Pressure-resolved local curl/WKB bridge (C194)](research/2026-08-28-a2-curl-wkb-bridge-c194.md)
- [Adversarial audit of C193--C194](audit/AUDIT-2026-08-28-c193-c194.md)
- [Quantitative off-ray \(A_2\) finite-horizon dominated cones (C195)](research/2026-08-30-a2-offray-hyperbolic-bundle-c195.md)
- [Periodic endpoint construction and phase-space ceiling (C196)](research/2026-08-30-periodic-endpoint-phase-space-c196.md)
- [Adversarial audit of C195--C196 and the success portfolio](audit/AUDIT-2026-08-30-c195-c196.md)
- [Counterfactual success portfolio](research/2026-08-30-counterfactual-success-portfolio.md)
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
