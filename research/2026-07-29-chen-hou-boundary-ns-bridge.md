# Chen--Hou smooth boundary Euler blowup versus periodic Navier--Stokes

Date: 2026-07-29

## Decisive verdict

The smooth-data Chen--Hou boundary blowup is **not** a viable leading core for
the Clay periodic Navier--Stokes problem.  Three independent obstructions
close the direct bridge.

1. The core contracts in both active meridional directions like
   \[
   \ell(t)\asymp (T-t)^\beta,\qquad
   \beta=-\frac{\bar c_l}{\bar c_\omega}
   \approx 2.920561.
   \]
   This is far below the heat scale.  In Chen--Hou's dynamic variables, a
   fixed positive viscosity has coefficient
   \[
   \mu_2(t)=\nu\,\frac{C_\omega}{C_l^2}
   \asymp \nu (T-t)^{1-2\beta}
   =\nu (T-t)^{-4.841122}.
   \]
   Thus viscosity grows without bound relative to every Euler term.  It
   cannot be included as a small correction, even formally to all orders.

2. Smooth reflection already fails at the first boundary trace.  A geometric
   reflection across a wall requires the normal velocity to be odd and the
   tangential velocities to be even.  It consequently requires the
   tangential vorticity, in particular \(\omega^\theta\), to vanish at the
   wall.  Chen--Hou explicitly use
   \(\omega(x,0)\ne0\), where \(\omega=\omega^\theta/r\), and observe that its
   odd elliptic extension across the wall is not Hölder in the normal
   direction.  That odd extension is a Green-function device, not a smooth
   Euler or Navier--Stokes reflection.

3. Reflection across the cylindrical surface \(r=1\) is not an exact
   Euclidean symmetry.  The map \(r\mapsto2-r\) changes the coefficients
   \(1/r\) and \(1/r^2\).  Doubling the cylinder instead produces a curved,
   and under the naive double only Lipschitz-metric, manifold rather than the
   flat three-torus.  The boundary therefore cannot be removed by parity
   without an equation error.

The remaining escape routes all abandon the Chen--Hou core as the leading
object: a new parabolic-scale viscous profile, a non-precompact multiscale
cascade whose second derivatives are not close to the Euler profile, or a
force singular at \(T\).  The last option is inadmissible for Clay
alternative (D).

What survives is a useful conditional no-go theorem, stated in Section 7,
and a plausible fractional-dissipation stability target below order
\(\sigma<0.1712\).  Neither resolves ordinary Navier--Stokes.

## 1. Primary-source facts being audited

Chen and Hou prove stable, nearly self-similar finite-time blowup for smooth
finite-energy data in the axisymmetric Euler equations on

\[
0\le r\le1,\qquad z\in\mathbb T,
\]

with a solid boundary at \(r=1\).  Their variables are

\[
\Gamma=ru^\theta,\qquad
\theta=\Gamma^2,\qquad
\omega=\frac{\omega^\theta}{r}.
\]

The wall condition is impermeability, expressed through the angular stream
function as

\[
\widetilde\phi(1,z)=0.
\]

It is not a no-slip or stress-free Navier--Stokes boundary condition.
The angular velocity and angular vorticity are odd and periodic in \(z\),
and the blowup occurs at the boundary circle \(r=1,z=0\).

The rescaled coordinates in Part I, equations (6.8)--(6.13), are

\[
X=C_l^{-1}z,\qquad Y=C_l^{-1}(1-r).
\]

There is one common length \(C_l\).  The profile is strongly anisotropic in
its *values and derivatives*, but there are not different axial and normal
similarity exponents.  In three-dimensional geometry the anisotropy is:

| direction | physical core length |
|---|---:|
| axial \(z\) | \(C_l\asymp s^\beta\) |
| wall-normal \(1-r\) | \(C_l\asymp s^\beta\) |
| azimuthal \(r\vartheta\) | \(O(1)\), because the solution is axisymmetric |

Here and below \(s=T-t\).

Part I gives the rigorously enclosed approximate-profile parameters

\[
\bar c_l\approx3.00649898,\qquad
\bar c_\omega\approx-1.02942516,\qquad
\frac{\bar c_l}{\bar c_\omega}\approx-2.9205600.
\]

Part II verifies the stability inequalities and error bounds needed by the
Part I analysis.  It does not introduce a viscous estimate or a
Navier--Stokes profile.

The dynamic scales satisfy

\[
C_\omega(\tau)\asymp e^{\bar c_\omega\tau},\qquad
C_l(\tau)\asymp e^{-\bar c_l\tau},\qquad
s\asymp C_\omega(\tau).
\]

Consequently

\[
C_l\asymp s^\beta,\qquad
\beta=-\frac{\bar c_l}{\bar c_\omega}
=2.920561005\ldots.
\]

Primary sources:

- J. Chen and T. Y. Hou,
  [*Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler
  equations with smooth data I: Analysis*](https://arxiv.org/abs/2210.07191),
  especially equations (2.6)--(2.12), (2.23), (6.2)--(6.13), and Theorem 4.
- J. Chen and T. Y. Hou,
  [*Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler
  equations with smooth data II: Rigorous Numerics*](https://arxiv.org/abs/2305.05660),
  especially Section 2 and the rigorous-numerics appendices.
- C. Fefferman,
  [official Clay problem statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf),
  equations (8)--(11) and alternative (D).

## 2. Component amplitudes and all active similarity exponents

The vorticity rescaling is

\[
\omega_{\rm ss}(X,Y,\tau)
=C_\omega\,\omega_{\rm phy}(C_lX,C_lY,t),
\]

so the angular-vorticity variable and meridional velocity have sizes

\[
\omega_{\rm phy}\asymp s^{-1},\qquad
u_{\rm mer}\asymp\frac{C_l}{C_\omega}
\asymp s^{\beta-1}=s^{1.920561\ldots}.
\]

Although the meridional velocity tends to zero in the shrinking core, its
gradient grows like \(s^{-1}\).

The swirl has a different amplitude.  Chen--Hou's relation

\[
C_\theta=\frac{C_\omega^2}{C_l}
\]

gives

\[
\theta_{\rm phy}=\Gamma^2
\asymp \frac{C_l}{C_\omega^2}
\asymp s^{\beta-2},
\]

and hence

\[
\Gamma,\ u^\theta
\asymp s^{\beta/2-1}
=s^{0.4602805\ldots}.
\]

The swirl also tends to zero, but its meridional derivatives have size

\[
\nabla_{r,z}u^\theta
\asymp s^{-\beta/2-1}
=s^{-2.4602805\ldots}.
\]

Thus it is important not to assign one amplitude exponent to all three
velocity components.  The relevant Euler and viscous sizes are:

| quantity | Euler/time scale | full-Laplacian scale |
|---|---:|---:|
| \(\omega^\theta/r\) | \(s^{-2}\) | \(\nu s^{-1-2\beta}\) |
| \(u_{\rm mer}\) | \(s^{\beta-2}\) | \(\nu s^{-\beta-1}\) |
| \(\Gamma\) or \(u^\theta\) | \(s^{\beta/2-2}\) | \(\nu s^{-3\beta/2-1}\) |

Every row has the same ratio

\[
\frac{\text{viscosity}}{\text{Euler}}
\asymp \nu s^{1-2\beta}
=\nu s^{-4.841122\ldots}.
\]

In particular, inserting the Euler velocity into Navier--Stokes requires

\[
f_E=-\nu\Delta u_E
\]

after Leray projection.  The Laplacian of a divergence-free field is already
divergence-free, so this term cannot be hidden in the pressure.  For a
nonharmonic Chen--Hou profile the meridional part of \(f_E\) grows like
\(\nu s^{-\beta-1}\), and the swirl part grows like
\(\nu s^{-3\beta/2-1}\).  Numerically these are

\[
\nu s^{-3.920561\ldots},
\qquad
\nu s^{-5.380842\ldots}.
\]

The required force is not merely non-flat at \(T\); it does not extend
continuously there.

## 3. Exact viscous audit in cylindrical variables

For ordinary three-dimensional axisymmetric Navier--Stokes, the circulation
and angular-vorticity variable obey the viscous operators

\[
\begin{aligned}
\Gamma_t+u^r\Gamma_r+u^z\Gamma_z
&=\nu\left(\Gamma_{rr}-\frac1r\Gamma_r+\Gamma_{zz}\right),\\
\omega_t+u^r\omega_r+u^z\omega_z
&=\frac1{r^4}\partial_z(\Gamma^2)
 +\nu\left(\omega_{rr}+\frac3r\omega_r+\omega_{zz}\right).
\end{aligned}
\]

These follow as the \(n=3\) case of equation (1.1) in Hou's generalized
axisymmetric formulation, which explicitly states that \(n=3\) recovers
ordinary Navier--Stokes.

Under

\[
z=C_lX,\qquad r=1-C_lY,\qquad
\frac{dt}{d\tau}=C_\omega,
\]

the viscous coefficients sort into three levels:

| cylindrical term | coefficient in dynamic variables | behavior |
|---|---:|---:|
| \(\partial_{XX}+\partial_{YY}\) | \(\displaystyle\mu_2=\nu C_\omega C_l^{-2}\) | \(\nu s^{-4.841122}\to\infty\) |
| \(r^{-1}\partial_r\) | \(\displaystyle\mu_1=\nu C_\omega C_l^{-1}r^{-1}\) | \(\nu s^{-1.920561}\to\infty\) |
| vector-basis \(r^{-2}\) terms | \(\displaystyle\mu_0=\nu C_\omega r^{-2}\) | \(\nu s\to0\) |
| azimuthal derivatives | zero | axisymmetry |

The factors \(3/r\) in the \(\omega\) equation and \(-1/r\) in the
\(\Gamma\) equation only change signs and constants in the second row.
The two genuine heat directions, \(z\) and \(1-r\), have the same fatal
coefficient \(\mu_2\).  The first-order curvature term also diverges, though
it is smaller than the principal heat term by \(C_l\).  Only the algebraic
vector-basis curvature is asymptotically negligible.

Equivalently,

\[
\frac{\ell(s)}{\sqrt{\nu s}}
=\frac{s^{\beta-1/2}}{\sqrt{\nu}}\longrightarrow0,
\]

and

\[
\frac{t_{\rm heat}}{t_{\rm Euler}}
=\frac{\ell^2/\nu}{s}
=\frac{s^{2\beta-1}}{\nu}\longrightarrow0.
\]

For any fixed \(\nu>0\), the crossover is

\[
s_\nu\asymp\nu^{1/(2\beta-1)}
=\nu^{0.2065637\ldots}.
\]

Choosing a very small viscosity only postpones the crossover.  It cannot
remove it from a construction that reaches \(s=0\).

For reference, the closest current numerical work by Hou does not contradict
this audit.  It obtains a nearly parabolic exponent \(c_l\approx0.5233\) only
for a **generalized** axisymmetric system with dimension about \(3.188\) and
solution-dependent viscosity; its constant-viscosity experiment uses a
generalized Boussinesq system with unequal viscosities and dimension about
\(4.73\).  It is not the \(n=3\), equal, constant-viscosity Clay equation.
See T. Y. Hou,
[*Nearly self-similar blowup of generalized axisymmetric Navier--Stokes and
Boussinesq equations*](https://arxiv.org/abs/2405.10916), Sections 1.1--1.3.

## 4. Why an all-order viscous correction does not repair the core

The parameter of a putative viscous perturbation is not \(\nu\) by itself.
It is

\[
\mu_2(\tau)=\nu C_\omega C_l^{-2}
\asymp\nu e^{(2\bar c_l+\bar c_\omega)\tau}.
\]

Since

\[
2\bar c_l+\bar c_\omega
\approx4.9835728>0,
\]

successive rescaled times make the effective perturbation larger.  A formal
\(\nu\)-series acquires powers

\[
\left(\nu s^{1-2\beta}\right)^m,
\]

which grow with \(m\).  It is not an asymptotic expansion near the proposed
singularity.

There is also no profile-level cancellation available.  If the leading
rescaled vorticity \(\Omega\) remained nontrivial and \(C^2\)-close to the
Chen--Hou profile, dividing the rescaled equation by
\(\mu_2\to\infty\) would give

\[
\Delta_{X,Y}\Omega=0
\]

at leading order.  The Chen--Hou profile is nontrivial, localized in the
physical construction, and decays in the similarity far field; it is not a
decaying harmonic function.  A same-frequency linear corrector that cancels
its Laplacian must cancel the profile itself at leading amplitude.  A
high-frequency, small-amplitude corrector might be small in \(C^0\), but not
in the rescaled \(C^2\) topology tested by the force.  Such a construction
would be a new non-precompact multiscale mechanism, not an all-order
perturbation of Chen--Hou.

The only honest inner reorganization is to move to a scale at least
\(\sqrt{\nu s}\) and solve a new viscous profile problem.  The exponent
\(\beta=2.92056\) and the Chen--Hou stability theorem then cease to be the
leading description.

## 5. The smooth-reflection obstruction

First flatten the wall locally and let \(Y\) be its normal coordinate.
A geometric reflection of a vector field is

\[
u_Y(X,-Y)=-u_Y(X,Y),\qquad
u_{\rm tan}(X,-Y)=u_{\rm tan}(X,Y).
\]

For a \(C^\infty\) reflected field this requires the complete jet conditions

\[
\partial_Y^{2m}u_Y|_{Y=0}=0,\qquad
\partial_Y^{2m+1}u_{\rm tan}|_{Y=0}=0,
\qquad m=0,1,2,\ldots.
\]

Impermeability supplies only the first condition with \(m=0\).  It gives
none of the tangential stress or higher-jet conditions.

Already at first order, geometric parity forces tangential vorticity to be
odd across the wall and therefore zero on it.  In the axisymmetric
coordinates this includes

\[
\omega^\theta|_{r=1}=0.
\]

Indeed, \(u^r(1,z)=0\) for every \(z\), so
\(\partial_z u^r(1,z)=0\), while

\[
\omega^\theta=\partial_z u^r-\partial_r u^z.
\]

Thus a nonzero wall value of \(\omega^\theta\) is exactly a nonzero normal
derivative of the tangential velocity \(u^z\).  Even reflection of \(u^z\)
then reverses that derivative and creates a first-derivative jump.

Chen--Hou's analysis explicitly says, in the discussion following equation
(2.33) in Part I,

\[
\omega(x,0)\ne0,\qquad \eta(x,0)\ne0,
\]

and \(\omega=\omega^\theta/r\).  In Section 3.4 they further state that the
odd extension \(W\) used for the half-plane Biot--Savart law is not Hölder
in the normal direction across the boundary.

At the single symmetry point \(z=0\), oddness in \(z\) makes
\(\omega^\theta(1,0)=0\).  This does not help: the boundary profile is
nonzero for neighboring rescaled \(X\), so there is no open wall
neighborhood on which the reflected velocity is smooth.  Therefore:

- the Green-function odd extension is discontinuous in vorticity;
- a reflected velocity cannot be \(C^1\), much less \(C^\infty\);
- its Laplacian contains a wall-supported distribution;
- removing that distribution would require a surface force, which is not a
  smooth periodic Clay force.

The standard viscous alternative does not help.  A no-slip wall requires
all velocity components to vanish, while Chen--Hou impose only no flow.  A
stress-free reflection wall requires the tangential normal derivatives, and
hence the relevant tangential vorticity traces, to vanish.  Chen--Hou satisfy
neither Navier--Stokes boundary class.

## 6. Curvature and periodicization

Even if the trace mismatch were redesigned, the cylindrical reflection
would not solve the same flat equation.  In the Euclidean metric

\[
ds^2=dr^2+r^2d\vartheta^2+dz^2,
\]

the map \(r\mapsto2-r\) is not an isometry.  It changes \(r^2\), \(1/r\), and
\(1/r^2\).  This is visible directly in Chen--Hou's exact rescaled elliptic
operator:

\[
-\Delta_{X,Y}\phi
+\frac{C_l}{r}\phi_Y
+\frac{2C_l^2}{r^2}\phi
=r\omega.
\]

The \(C_l/r\) and \(C_l^2/r^2\) terms are small in the Euler blowup analysis,
which is why the half-plane Boussinesq profile is relevant asymptotically.
They are not identically symmetric across \(Y=0\).

There are only three elementary ways to try to remove the wall:

1. **Reflect into \(r>1\) in Euclidean space.**  The reflected field fails
   the exact cylindrical equations because the coefficients do not reflect.
2. **Double the cylinder as a manifold.**  The doubled metric is not the
   flat metric of \(\mathbb T^3\); the naive even double of
   \((1-Y)^2d\vartheta^2\) is not smooth at \(Y=0\).
3. **Cut off and embed a wall patch in \(\mathbb T^3\).**  The cutoff,
   divergence correction, and pressure matching create a residual.  Exterior
   residuals might be made smooth, but the wall trace and the divergent
   principal viscous residual occur in the shrinking core and cannot.

The original \(z\)-periodicity and physical \(\vartheta\)-periodicity do not
solve the missing radial periodicity.  The axisymmetric singular set is a
circle, not a local point that can simply be copied into a Cartesian torus
without breaking the proved symmetry class.

## 7. A surviving conditional theorem

The following proposition is not a solution of Navier--Stokes, but it is a
useful route-closing theorem that can be made fully rigorous.

**Proposition (super-parabolic reflected-profile obstruction).**
Let \(s=T-t\), \(\beta>1/2\), and let a smooth forced periodic
Navier--Stokes solution with fixed \(\nu>0\) have, around a shrinking
axisymmetric ring or a locally reflected wall core, the vorticity form

\[
\omega_{\rm phy}=s^{-1}
\left(\Omega\!\left(\frac{z-z_*}{s^\beta},
\frac{r-r_*}{s^\beta}\right)+o_{C^2_{\rm loc}}(1)\right),
\]

where \(\Omega\) is bounded, nonzero, and tends to zero at spatial infinity.
Assume the force and its curl remain bounded near \(T\), as they must for a
smooth Clay force.  Then no such one-profile asymptotic can hold.

**Proof sketch.**  The time, transport, and stretching terms are
\(O(s^{-2})\).  The two principal meridional diffusion terms are

\[
\nu s^{-1-2\beta}
\left(\Delta_{X,Y}\Omega+o(1)\right).
\]

The curvature first- and zeroth-order terms are smaller than the principal
diffusion by \(s^\beta\) and \(s^{2\beta}\), respectively.  The smooth force
is negligible after vorticity rescaling.  Divide the equation by
\(\nu s^{-1-2\beta}\).  Since \(s^{2\beta-1}\to0\), the limit is

\[
\Delta_{X,Y}\Omega=0.
\]

After a genuinely smooth reflection the profile is defined on the full
cross-sectional plane.  A bounded harmonic function that tends to zero is
zero, contradicting the hypothesis.  For Chen--Hou, the reflection premise
already fails because \(\omega^\theta\) has nonzero wall trace.  \(\square\)

This proposition deliberately assumes \(C^2\) profile compactness.  Its
escape hatch is exactly a non-precompact cascade or a new inner heat scale.
That is why it closes the Chen--Hou *profile bridge* without claiming global
regularity of Navier--Stokes.

## 8. What, if anything, is reusable

### 8.1 A non-Clay fractional-dissipation target

For dissipation \(\nu(-\Delta)^\sigma\), the rescaled coefficient is

\[
\mu_\sigma
=\nu\frac{C_\omega}{C_l^{2\sigma}}
\asymp\nu s^{1-2\sigma\beta}.
\]

It decays, rather than grows, precisely when

\[
\sigma<\frac1{2\beta}
=0.17119998\ldots.
\]

Subject to choosing the initial effective coefficient small, Chen--Hou's
weighted \(L^\infty\)-\(C^{1/2}\) stability framework may plausibly prove
persistence under such weak hypodissipation.  This would be a legitimate
new theorem, but it is far from the ordinary Laplacian \(\sigma=1\).

### 8.2 Proof technology

The reusable ingredients are the computer-assisted construction of an
approximate steady state, rigorous residual enclosure, finite-rank treatment
of nonlocal terms, and weighted \(L^\infty\)-Hölder stability on an
unbounded similarity domain.  They could be applied to a genuinely viscous,
near-parabolic profile if one is found.

### 8.3 Strategic implication

The Chen--Hou mechanism focuses at a fixed super-parabolic exponent.  It
does not benefit from replacing a geometric carrier sequence by a slowly
varying schedule such as \(K_j=j^A\): the physical effective viscosity of
the proposed profile is already forced to diverge as \(s^{-4.841122}\).
A slow bounded-ratio Kelvin handoff is therefore a distinct cascade
architecture, not a modulation of this boundary blowup.

## Bottom line

The Chen--Hou papers supply a genuine smooth Euler singularity with boundary
and powerful stability machinery.  They do **not** supply a wall-removable
or viscosity-perturbative seed for periodic Navier--Stokes.  The two numbers

\[
\boxed{\beta=2.920561\ldots},\qquad
\boxed{\mu_2=\nu s^{-4.841122\ldots}}
\]

close the viscous route, while

\[
\boxed{\omega^\theta|_{r=1}\ne0}
\]

closes smooth parity reflection at its first necessary trace.
