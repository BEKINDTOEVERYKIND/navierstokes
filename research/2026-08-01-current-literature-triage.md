# Current-literature triage for the surviving Navier--Stokes routes

**Date:** 2026-08-01
**Status:** primary-source triage, not a literature review and not a claim of
priority or correctness for unreviewed preprints.

The purpose of this note is strategic: identify machinery that can change a
load-bearing step in the repository's two surviving branches, and distinguish
it from results that concern forward self-similar singular data, weak
solutions, models, or excluded Type-I scenarios.

## 1. Machinery worth importing now

### 1.1 Coercive plus compact plus certified finite rank

Hou, Wang and Yang give a computer-assisted construction of an exact forward
self-similar Navier--Stokes profile and unstable eigenpair by decomposing the
linearized operator into a coercive part and a compact perturbation, then
approximating the compact part by a finite-rank operator with a controlled
tail:

* T. Hou, Y. Wang and C. Yang,
  [*Nonuniqueness of Leray--Hopf solutions to the unforced incompressible 3D
  Navier--Stokes Equation*](https://arxiv.org/abs/2509.25116).

This is directly relevant methodologically to two open tasks here:

1. certify the full finite-frequency AO residual/semigroup data if the
   quasimode gate cannot be closed analytically; and
2. validate the finite-dimensional complement of a numerically constructed
   three-phase endpoint map after the infinite-dimensional part has a
   coercive estimate.

It is not itself a Clay blow-up result.  The paper explicitly starts from a
forward scale-invariant \(|x|^{-1}\) singular datum; its localized datum is
compactly supported but belongs to
\(C^\infty(\mathbb R^3\setminus\{0\})\cap L^q\), \(q<3\), rather than being
smooth at the origin.  The authors also distinguish this setting from
backward self-similarity arising from smooth initial data.

### 1.2 Analytic finite-rank correction of exact local constraints

Chen and Hou isolate a reusable analytic correction mechanism: numerics
determine approximate profiles and low defect modes, while explicit
finite-rank Taylor corrections impose the exact vanishing conditions required
by singularly weighted estimates:

* J. Chen and T. Hou,
  [*Analytic finite-rank corrections for singularly weighted estimates in a
  computer-assisted proof of 3D Euler singularity*](https://arxiv.org/abs/2607.15256).

The plausible import is narrow but valuable.  Apply the correction principle
to the finite list of endpoint moments, pressure-centre modes, and phase
gauges **after** the three-phase evolution supplies a tame approximate right
inverse.  It does not create the missing infinite-dimensional coercivity or
Gevrey-2 estimate.

### 1.3 Embedded quasimodes instead of isolated spectrum

The older sources become newly useful at the recovered AO endpoint:

* S. Friedlander, W. Strauss and M. Vishik,
  [*Nonlinear instability in an ideal fluid*](https://www.numdam.org/article/AIHPC_1997__14_2_187_0.pdf);
* R. Shvydkoy,
  [*The essential spectrum of advective equations*](https://arxiv.org/abs/math-ph/0412019).

The repository's C50--C51 continuation extracts a quantitative quasimode
gate from these ideas.  If its normalized hypotheses hold, the exact
finite-frequency isolated eigenvalue is unnecessary.

## 2. New results that constrain rather than close the route

### 2.1 Type-I rotated self-similarity is further restricted

Pineau and Vicol prove Liouville theorems for rotated backward self-similar
and rotated discretely self-similar solutions under Type-I bounds in several
parameter regimes:

* B. Pineau and V. Vicol,
  [*On rotated backwards self-similar solutions of the incompressible 3D
  Navier--Stokes equations*](https://arxiv.org/abs/2607.09619).

This makes a rotated Type-I shortcut less attractive.  It does not exclude
the repository's Type-II/high-Reynolds window, whose leading scaling is
Eulerian and non-tight.

### 2.2 Homothetic forward profiles do not provide a singular-limit shortcut

Binz and Coiculescu rule out nontrivial sufficiently regular homothetic
forward self-similar solutions in 3D and identify the homothetic class as the
one compatible with a particular singular-limit route to Leray--Hopf
nonuniqueness:

* T. Binz and M. Coiculescu,
  [*Homothetic Self-Similar Solutions to the Incompressible Navier--Stokes
  Equations*](https://arxiv.org/abs/2607.12159).

This is another reason not to base the blow-up program on amplitude
homothety of one forward profile.  Their approximate unstable-eigenvalue
example is interesting for quasimode methodology, not a backward blow-up.

### 2.3 Smooth 1D limiting profiles lift only to low-regularity 3D Euler

Chen constructs smooth one-dimensional limiting profiles which, in companion
work, lift to 3D Euler velocities of class \(C^{1,\alpha}\),
\(0<\alpha<1/3\):

* J. Chen,
  [*Asymptotically Self-Similar Blowup for 3D Incompressible Euler with
  \(C^{1,1/3-}\) Velocity I*](https://arxiv.org/abs/2605.15149).

This is a major Euler development but not a direct Navier--Stokes seed: the
lifted velocity does not have the classical second derivatives needed for a
smooth viscous solution, and the earlier repository audit finds the relevant
Euler-II scaling super-parabolic for a direct viscosity perturbation.

### 2.4 Forward unstable self-similar Navier--Stokes profiles concern
singular initial data

Ionescu, Jia and Palasek numerically construct global axisymmetric
swirl-free forward self-similar profiles with pointwise residual about
\(10^{-10}\) and unstable modes:

* A. Ionescu, H. Jia and S. Palasek,
  [*On the non-uniqueness of solutions of the axi-symmetric swirl-free
  Navier--Stokes equations, I*](https://arxiv.org/abs/2606.07501).

The numerical and spectral tools may help a certified finite-dimensional
block.  Forward profiles homogeneous of degree \(-1\) at infinity address
weak nonuniqueness from critical/singular data, not finite-time breakdown
from smooth Clay data.

## 3. Model and finite-scale results

Palasek's new shell model has genuine viscous Type-II blow-up and is designed
to make the embedding step more tractable:

* S. Palasek,
  [*Finite-time blow-up in an elementary model of the 3D Navier--Stokes
  equations*](https://arxiv.org/abs/2605.13827).

It strengthens Step 1 of a shell-model program.  The paper itself emphasizes
that embedding Step 2 into the full Navier--Stokes nonlinearity remains the
hard part.  The repository's sideband, recurrence, and finite-cell audits are
precisely tests of that gap; no new full-PDE embedding theorem is supplied.

Yu's supply--tax reduction gives a rigorous finite-scale accounting
alternative for suitable weak solutions:

* R. Yu,
  [*Critical Ledgers and Scale-Defect Cascades for Navier--Stokes*](https://arxiv.org/abs/2606.13887).

It is useful language for falsifying a candidate cascade, but the paper
explicitly remains finite-scale/conditional.  It neither supplies the
recurrent transition nor rules out every profitable critical mechanism.

## 4. Strategic conclusion

No source above closes smooth-data global regularity or finite-time breakdown.
The literature changes the implementation ranking as follows.

1. **Constructive Type-II route:** three positive material phases + the
   exact Gavrilov storage wake + the now-summable global second-jet pressure
   wake.  Use analytic finite-rank corrections only for the finite moment
   block and coercive/compact certification for a finite residual block.
2. **Spectral trigger route:** recompute the AO \(L^p\) growth ratio and
   normalized residual.  Apply the quasimode gate before attempting an exact
   embedded eigenvalue/global Weber theorem.
3. **Computer-assisted alternative:** search directly for one local
   three-phase return stage with interval-certified residual and a coercive
   infinite-dimensional complement.  This is meaningful only after an
   analytic tail estimate prevents finite truncation from assuming the
   desired answer.
4. **De-prioritize:** Type-I/RDSS profile searches, homothetic forward
   limits, direct low-regularity Euler-to-Navier--Stokes perturbation, and
   shell-model optimization without a new embedding mechanism.

The single load-bearing theorem remains an all-order, localized,
three-phase material-transition inverse with a global pressure-wake state.
The new sources offer credible tools for its finite-dimensional and
certification layers, but not its principal analytic estimate.
