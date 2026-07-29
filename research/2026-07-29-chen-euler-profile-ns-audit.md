# Chen's 2026 Euler profile versus full Navier--Stokes: a decisive scaling audit

Date: 2026-07-29

## Verdict

Jiajie Chen's 2026 construction is a major exact Euler result, but its
singular profile cannot be imported directly, perturbatively, or by a change
of clock into the Clay Navier--Stokes problem.

The decisive obstruction is not merely that the theorem treats Euler rather
than Navier--Stokes.  Its core focuses **faster than the parabolic scale**.
If \(s=T-t\) and \(c=c_{x,\alpha}\), Chen's profile has

\[
 \omega=s^{-1}\Omega(x/s^c),\qquad
 \ell(s)=s^c,\qquad
 a(s)=s^{c-1}.
\]

Near \(\alpha=1/3\), the paper proves

\[
 c=\frac{8}{9}(1/3-\alpha)^{-1}
   +o\big((1/3-\alpha)^{-1}\big),
\]

so in particular \(c\gg1\).  At that scale the ratio of viscosity to the
Euler terms is

\[
 \mu(s)=\nu s^{1-2c}\longrightarrow\infty,
\]

the local Reynolds number is

\[
 \operatorname{Re}(s)=\frac{a(s)\ell(s)}{\nu}
 =\frac{s^{2c-1}}{\nu}\longrightarrow0,
\]

and the core lies far below the heat length:

\[
 \frac{\ell(s)}{\sqrt{\nu s}}
 =\nu^{-1/2}s^{c-1/2}\longrightarrow0.
\]

Thus heat acts on a time much shorter than one nonlinear turnover.  There is
no small-viscosity expansion about this Euler profile as \(t\uparrow T\);
the nominal perturbation parameter diverges.

There are two independent regularity mismatches as well.  Chen's initial
velocity is \(C^{1,\alpha}\), not smooth, and its axisymmetric no-swirl cusp
is the singular mechanism.  A smooth finite-energy Navier--Stokes
realization in the same symmetry class is globally regular.  Moreover, for
the near-critical profiles \(c>1\), the core velocity amplitude tends to
zero even while its gradient blows up.  A localized smooth Navier--Stokes
singularity cannot be produced by adding only a bounded smooth exterior to
such a core.

The genuinely reusable contribution is proof technology: an exact
nonlinear fixed point lifting a one-dimensional profile into a nonlocal
three-dimensional fluid equation; anisotropic weighted Biot--Savart
estimates; outgoing-flow stability; trajectory integration by parts; and a
finite-codimension treatment of unstable modes.  These tools may be useful
for a different, high-Reynolds, non-precompact cascade.  They do not rescue
Chen's profile itself.

## 1. What the paper actually proves

Chen constructs exact self-similar vorticity profiles for three-dimensional
axisymmetric Euler without swirl.  For any regularity label
\(\alpha\in(0,1/3)\), the theorem supplies a nontrivial \(C^\alpha\)
profile.  The quantitatively constructed family is near the endpoint
\(\alpha=1/3\), where

\[
 \omega^\theta_\alpha(t,x)
 =\frac1{1-t}
 \Omega^\theta_{*,\alpha}
 \left(\frac{x}{(1-t)^{c_{x,\alpha}}}\right)
\]

and

\[
 c_{x,\alpha}\asymp(1/3-\alpha)^{-1}.
\]

The statement for all lower Hölder exponents uses the inclusion of the
near-critical profile in every lower \(C^\gamma\) class; it does not produce
a Clay-compatible smooth profile by sending \(\gamma\) upward.

The exact profile has slow spatial decay and infinite kinetic energy.
Chen then proves finite-codimension stability and obtains asymptotically
self-similar blow-up from compactly supported \(C^\gamma\) vorticity and
\(C^{1,\alpha}\cap L^2\) velocity.  This repairs the energy at spatial
infinity, but it does not repair the cusp at the symmetry axis.  The paper
explicitly identifies that low regularity as coming from the vanishing order
near \(r=0\).

These distinctions matter for the Clay statement.  The official problem
requires a \(C^\infty\) divergence-free initial velocity and a \(C^\infty\)
force.  To prove breakdown by prescribing a classical branch, that branch
must be smooth on every compact subinterval before \(T\).  Chen's branch is
only \(C^{1,\alpha}\) at every such time.

Primary sources:

* J. Chen, [*Asymptotically Self-Similar Blowup for 3D Incompressible Euler
  with \(C^{1,1/3-}\) Velocity II*](https://arxiv.org/abs/2605.15130),
  especially Theorems 1.1--1.2.
* C. Fefferman, [official Clay problem
  statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf).

## 2. Full-Laplacian scaling

Set \(s=T-t\).  The Biot--Savart law is of order \(-1\), so Chen's vorticity
ansatz implies the velocity scaling

\[
 u_E(x,t)=s^{c-1}U(y),\qquad y=x/s^c.
\]

Consequently,

\[
 \begin{array}{c|c}
 \text{quantity} & \text{size}\\ \hline
 u_E & s^{c-1}\\
 \nabla u_E,\ \omega_E & s^{-1}\\
 \partial_tu_E & s^{c-2}\\
 u_E\cdot\nabla u_E & s^{c-2}\\
 \nu\Delta u_E & \nu s^{-c-1}.
 \end{array}
\]

The velocity equation therefore gives

\[
 \frac{\|\nu\Delta u_E\|}
 {\|\partial_tu_E\|+\|u_E\cdot\nabla u_E\|}
 \asymp \nu s^{1-2c}.
\]

The same computation in vorticity form is even more transparent:

\[
 \partial_t\omega,\quad
 u\cdot\nabla\omega,\quad
 \omega\cdot\nabla u
 \asymp s^{-2},
\qquad
 \nu\Delta\omega\asymp\nu s^{-1-2c}.
\]

Again the ratio is \(\nu s^{1-2c}\).  For every fixed positive viscosity and
every \(c>1/2\), diffusion eventually dominates.  The crossover satisfies

\[
 s_\nu^{\,2c-1}\asymp\nu.
\]

Making \(\nu\) numerically small merely postpones the crossover; it cannot
remove it from a cascade reaching \(s=0\).

Equivalently, the nonlinear turnover and heat times are

\[
 t_{\rm nl}\asymp\frac{\ell}{a}=s,\qquad
 t_{\rm heat}\asymp\frac{\ell^2}{\nu}
 =\frac{s^{2c}}{\nu},
\]

so

\[
 \frac{t_{\rm heat}}{t_{\rm nl}}
 =\operatorname{Re}(s)\longrightarrow0.
\]

The same conclusion appears directly in the heat semigroup.  A packet at
frequency \(k(s)\asymp s^{-c}\), left to evolve for the remaining time
\(s\), acquires the damping factor

\[
 \exp(-\nu k(s)^2s)
 =\exp(-\nu s^{1-2c})\longrightarrow0.
\]

Maintaining that packet would require continuous nonlinear replenishment on
a time shorter than \(t_{\rm heat}\), while its own turnover time is
\(t_{\rm nl}\gg t_{\rm heat}\).

This is a structural obstruction, not a question of insufficient numerical
resolution.

### 2.1 The prescribed-force residual is singular

Because \(u_E\) already solves Euler, inserting it into Navier--Stokes leaves

\[
 f_E=-\nu\Delta u_E
\]

after applying the Leray projection.  The Laplacian of a divergence-free
field is already divergence-free, so pressure cannot hide it.  For every
spatial derivative order \(m\), a nondegenerate profile has

\[
 \nabla^m f_E
 \asymp \nu s^{-c-1-mc}.
\]

Time derivatives add further negative powers of \(s\).  Thus the required
force is not merely non-flat at \(T\): it fails to extend continuously.
Also, since Chen's velocity is only \(C^{1,\alpha}\), this residual is not a
smooth classical force even at a fixed preterminal time.

Trying to cancel \(-\nu\Delta u_E\) with a small velocity corrector reverses
the perturbative hierarchy.  The corrector equation is driven by a term
larger than the Euler operator by the factor
\(\mu(s)\to\infty\).  A formal expansion in viscosity therefore has
increasing, rather than decreasing, coefficients.

There are only two elementary same-scale alternatives, and neither keeps
the Chen core as the leading object.  Linear cancellation of Laplacians at
frequency \(k\) cancels the total velocity at that frequency as well,
because \(\Delta\) has the scalar multiplier \(-|k|^2\).  Nonlinear
cancellation requires an added amplitude \(A\) with

\[
 A^2k\asymp\nu Ak^2,\qquad\text{hence}\qquad A\asymp\nu k.
\]

At \(k=s^{-c}\), this is larger than Chen's amplitude
\(s^{c-1}\) by the factor \(\nu s^{1-2c}=\mu(s)\to\infty\).
Such a component is a new viscosity-scale singular mechanism, not a
corrector to Chen's Euler solution.  For \(c>1/2\), its leading
nonlinearity--viscosity balance is quasistationary; a localized
finite-energy fixed profile is then killed by the standard energy identity
for stationary unforced Navier--Stokes.  Escaping that conclusion again
requires a non-precompact wake or a genuinely time-dependent multi-scale
balance.

## 3. General similarity threshold

For comparison, take a fixed-profile ansatz

\[
 u=s^{-A}U(x/s^B).
\]

Time, nonlinear, and viscous terms have sizes

\[
 s^{-A-1},\qquad s^{-2A-B},\qquad
 \nu s^{-A-2B}.
\]

An Euler time--nonlinearity balance requires

\[
 A+B=1.
\]

On this Euler line the viscous-to-Euler ratio is

\[
 \nu s^{1-2B}.
\]

Hence:

* \(B<1/2\): Euler-dominated and potentially perturbative;
* \(B=1/2\): full parabolic, time--nonlinearity--viscosity balance;
* \(B>1/2\): diffusion-dominated.

Chen has \(B=c\gg1\), deep in the third region.

This also explains why simply seeking another localized fixed Euler profile
does not give an easy escape.  The Euler profile equation on the line
\(A+B=1\) is

\[
 (1-B)U+B\,y\cdot\nabla U+
 \mathbb P\nabla\cdot(U\otimes U)=0.
\]

Assuming enough decay to justify integration by parts, its \(L^2\) identity
is

\[
 \left(1-\frac{5B}{2}\right)\|U\|_2^2=0.
\]

A nonzero finite-energy fixed profile therefore forces \(B=2/5\).
The energy-scale Type-I/DSS concentration is itself subject to the
Chae--Wolf exclusion.  A serious high-Reynolds alternative consequently
has to be non-precompact, non-fixed-profile, or retain an active wake rather
than merely replacing Chen's exponent:

* D. Chae and J. Wolf, [*Energy concentrations and Type I blow-up for the
  3D Euler equations*](https://arxiv.org/abs/1706.02020).

## 4. Why the obvious modifications do not evade the audit

### 4.1 Smoothing the axis cusp

Chen's finite-energy initial vorticity has the form

\[
 \omega^\theta_0=r^\gamma g,\qquad
 0<\gamma<1/3,\qquad g\in C_c^\infty.
\]

The factor \(r^\gamma\) is precisely the nonsmooth part.  Mollifying it
inside a fixed radius \(\delta\) changes the solution once the active length
\(s^c\) reaches \(\delta\).  For every fixed \(\delta>0\), that occurs before
the proposed singular time.  A sequence \(\delta_n\downarrow0\) gives a
sequence of smooth problems, not one smooth datum that blows up.

A time-dependent smoothing radius does not fix the scale conflict.  Let
\(\rho(s)\) be its physical radius.  If
\(\rho(s)\lesssim\ell(s)=s^c\), the Chen core is retained, but its largest
curvature is at least the already sub-parabolic curvature
\(\ell(s)^{-2}\), and is worse when \(\rho\ll\ell\).  If instead
\(\rho(s)\gtrsim\sqrt{\nu s}\), then, because
\(\ell(s)/\sqrt{\nu s}\to0\), the smoothing removes the entire similarity
core.  The intermediate cutoff annulus also creates time, Laplacian, and
Biot--Savart commutator errors.  Making all of those errors smooth through
\(T\) would be a new all-order inverse construction, not a mollified
instance of Chen's solution.

There is also a theorem-level barrier: smooth finite-energy axisymmetric
no-swirl Euler and Navier--Stokes data are globally regular.  The classical
whole-space result is:

* M. R. Ukhovskii and V. I. Yudovich, [*Axially symmetric flows of ideal and
  viscous fluids filling the whole
  space*](https://doi.org/10.1016/0021-8928(68)90147-0).

Therefore a smooth version must break the symmetry or add swirl.  Once it
does so, Chen's profile and stability theorem no longer supply the needed
solution.  Swirl by itself also does not alter the factor
\(\nu s^{1-2c}\); it would have to create a new exponent, a new leading
Laplacian cancellation, or a non-self-similar cascade.

An infinite nested smoothing that is \(C^\infty\) at the axis but resembles
\(r^\gamma\) on selected annuli is not ruled out by this elementary
argument.  In the no-swirl class it remains globally regular; outside that
class it becomes a new multiscale cascade problem rather than a smoothing
of the proved profile.

### 4.2 Changing the clock

The similarity exponent is not a freely adjustable clock.  Let

\[
 u(x,t)=a(t)U(x/\ell(t))
\]

use the same nondegenerate Euler profile with exponent \(c\).  Matching the
coefficients of its profile equation requires

\[
 a'=(1-c)\frac{a^2}{\ell},
\qquad
 \ell'=-ca.
\]

Eliminating \(t\) gives

\[
 a\propto\ell^{(c-1)/c},\qquad
 \ell\propto(T-t)^c.
\]

Thus a reparameterization changes only multiplicative constants.  It cannot
move \(c\) from the diffusion-dominated side to \(c<1/2\).

More directly, substituting the same profile at an exponent \(B\ne c\)
leaves the leading Euler residual

\[
 (c-B)(U-y\cdot\nabla U),
\]

times the singular physical prefactor.  It vanishes only for a degree-one
homogeneous/affine degeneracy, which is neither the Chen profile nor a
localized finite-energy field.  A smooth terminal force cannot absorb this
leading residual.

A time-dependent spatial dilation can be implemented only by changing the
flow, for example by adding an affine strain.  That is a potentially useful
Kelvin-wave mechanism, but the affine field is nonperiodic and has infinite
energy, and cutting it off creates the unresolved Reynolds-stress/wake
problem.  It is therefore a bridge to the separate localized return-cell
program, not a reparameterization covered by Chen's theorem.

### 4.3 Using the paper's anisotropy

Chen's profile is anisotropic in shape, decay, and characteristic flow, but
the similarity coordinate in the theorem is

\[
 y=x/s^c
\]

with one spatial exponent.  The paper's anisotropic weighted estimates do
not weaken the full Laplacian scaling.

For a genuinely anisotropic ansatz with lengths
\(\ell_j=s^{c_j}\), an active \(j\)-directional Laplacian divided by a
turnover term carries the factor

\[
 \nu s^{1-2c_j}.
\]

Every direction with \(c_j>1/2\) and a nonzero leading fast-variable
Laplacian is therefore diffusion-dominated.  Cancellation among directions
with the same largest exponent would require the profile to be harmonic in
those fast variables.  In a localized finite-energy class that harmonic
nullspace is zero; in a periodic class it is constant.  Exact independence
or affinity can also remove a directional second derivative, but:

1. independence reduces the active dynamics toward a lower-dimensional
   globally regular class;
2. a nonzero affine dependence is incompatible with finite energy and
   periodicity;
3. localizing either construction reintroduces curvature and heat.

A new anisotropic profile for which every genuinely curved direction has
\(c_j\le1/2\) remains logically possible.  It is not a modification
controlled by Chen's theorem and must solve a different full profile or
return-cell problem.

### 4.4 Hiding the singularity in a bounded velocity

For the constructed near-critical family \(c>1\), on any fixed bounded
similarity region (or a localized version of the core),

\[
 \|u_{\rm core}\|_\infty\asymp s^{c-1}\to0,
\qquad
 \|u_{\rm core}\|_3\asymp s^{2c-1}\to0,
\]

while \(\|\nabla u_{\rm core}\|_\infty\asymp s^{-1}\).
This is a valid low-regularity Euler singularity mechanism.  It is the wrong
shape for a localized smooth Navier--Stokes singularity: bounded velocity on
a finite interval lies in a Serrin regularity class, and bounded
\(L^\infty_tL^3_x\) is also an endpoint regularity criterion.  Adding a
bounded smooth exterior does not change that conclusion.  A second
component would have to become critical or unbounded; then that second
component, not the Chen core, carries the Navier--Stokes singularity.

For an accessible primary proof of the \(L^\infty_tL^3_x\) criterion, see:

* I. Gallagher, G. Koch, and F. Planchon, [*A profile decomposition approach
  to the \(L^\infty_t(L^3_x)\) Navier--Stokes regularity
  criterion*](https://arxiv.org/abs/1012.0145).

## 5. What is genuinely reusable

The following parts of Chen's work are worth importing into the surviving
cascade program.

### 5.1 Lifting a reduced model into the full nonlocal equation

The construction begins with a smooth one-dimensional profile, regards it
as an approximate three-dimensional profile with poor radial decay, and
uses a Schauder fixed point to obtain an exact 3D Euler profile.  This is a
real precedent for turning a reduced singular mechanism into an exact
nonlocal fluid solution rather than merely comparing formal equations.

### 5.2 Anisotropic weighted Biot--Savart estimates

The proof uses different radial and axial weights, angular averaging, and
sharp control of the Biot--Savart operator when the approximate profile has
little radial decay.  A packet-plus-annular-wake construction has the same
kind of mismatch: good decay in one geometric variable and a long,
structured tail in another.

### 5.3 Outgoing flow as a stability mechanism

The self-similar characteristic field satisfies a quantitative outgoing
condition.  Perturbations are carried from the core toward infinity, which
turns the local transport operator into a damped part of the linearization.
For a return cascade, the analogous goal is not to delete the outgoing
remainder but to place it in a controlled wake space.

This part is not plug-and-play.  Much of Chen's strong radial damping is
powered by the same large parameter \(c\asymp(1/3-\alpha)^{-1}\) that makes
the physical core super-parabolic.  A high-Reynolds cascade would need a
different source of coercivity, such as carrier separation or a
stage-dependent annular escape estimate.  In particular, the affine Kelvin
strain is a saddle rather than a globally outgoing field, so an analogue of
Chen's outgoing inequality has to be proved, not assumed.

### 5.4 Integration by parts along trajectories

Chen controls an apparently unbounded radial derivative by integrating
along characteristics and using the equation a second time, effectively
trading it for a better-decaying axial derivative plus lower-order terms.
This is directly suggestive for off-shell and pressure-tail terms in a
localized Kelvin packet.

### 5.5 Local transport plus compact nonlocal perturbation

The stability operator is split into a local transport/multiplication part
and a compact nonlocal Biot--Savart perturbation.  Semigroup estimates then
reduce instability to finitely many modes.  This is close in spirit to the
projected endpoint submersion needed by a packet-plus-wake return cell:
prescribe finitely many child coordinates and leave the infinite-dimensional
wake free.

## 6. Prize-relevant surviving variants

Only variants that abandon the diffusion-dominated Chen scaling survive.

1. **High-Reynolds annular cascade.**  Use Chen's weighted fixed-point and
   outgoing-flow methods to control the global wake of stages satisfying
   \(\operatorname{Re}_j\to\infty\), rather than using the Chen profile as
   the active core.
2. **Anisotropic heat-scale return cell.**  Seek a genuinely new profile
   with every curved direction at or above the heat length and solve the
   full Navier--Stokes, rather than Euler, endpoint map.  Any direction
   below the heat length must be exactly dark to the Laplacian, a severe
   algebraic gate.
3. **Symmetry-breaking multiscale cusp replacement.**  Replace the single
   \(r^\alpha\) cusp by an infinite smooth, non-axisymmetric hierarchy whose
   active stages reset before heat wins.  Chen's finite-codimension
   machinery might control modulation, but the construction would be a new
   cascade and must retain all discarded wake.

These are long shots, not consequences of the 2026 theorem.  No GPU
experiment on the Chen profile can answer the missing question: the
full-Laplacian scaling already decides that profile's fate.  Computation
becomes relevant only after a concrete high-Reynolds endpoint operator or
anisotropic full-NS profile equation has been derived.
