# Glass boundary control is not an unforced exterior wake

Date: 2026-07-29

## Verdict

Glass's exact controllability theorem is a powerful bounded-domain Euler
theorem, but it does **not** provide an analytic transition and its auxiliary
exterior domain is not an Euler wake.  The proof extends fields outside the
physical domain for transport bookkeeping, permits nonzero exterior
divergence, inserts exterior vorticity, and deletes vorticity packets after
they have left the physical domain.

There are three sharp consequences for the proposed transmutation.

1. Smooth gluing across the boundary of the ball requires the exterior
   solution to reproduce the **complete velocity trace, pressure trace, and
   all interface jets**.  Glass chooses only the normal velocity and incoming
   tangential vorticity as controls.  Its theorem is not a right inverse for
   this overdetermined exterior Cauchy problem.
2. A true unforced Euler wake cannot create the vorticity that Glass ejects
   and injects.  A transition from zero to a vortical target requires
   preloaded exterior vorticity, and a finite-energy wake must initially
   contain at least the final interior kinetic energy.  Thus a wake initially
   at rest, or one added stage by stage at zero cost, is impossible.
3. A generic Glass Euler path cannot be embedded exactly into unforced
   Navier--Stokes even locally: equality of the two equations would force
   \(\Delta\omega=0\) in the ball.  Positive-time spatial analyticity of
   regular Navier--Stokes solutions gives an independent endpoint
   obstruction.

Iterating Glass on annuli does not repair these defects.  An annulus has two
boundary components, while the theorem requires the control set to meet
both.  Controlling only the new outer sphere is outside the theorem and is
generically forbidden by Kelvin circulation.  Allowing control on the inner
sphere leaves the original actuator in place rather than pushing it
outward.

There is a conditional Euler route left open: prescribe a nontrivial,
preloaded exterior vortex reservoir and solve a new full-trace exterior
Euler matching problem with scale-uniform estimates.  That is the actual
``control-transmutation lemma'' needed below.  It is not a consequence of
Glass and is essentially the missing global transition theorem itself.
Nothing here proves a singularity or resolves the Navier--Stokes problem.

---

## 1. What Glass actually proves

The primary source is O. Glass,
[Exact boundary controllability of 3-D Euler equation, ESAIM COCV 5
(2000), 1--44](https://www.numdam.org/item/COCV_2000__5__1_0.pdf).

Let \(\Omega\subset\mathbb R^3\) be bounded, connected, and
\(C^\infty\), and let

\[
 \Gamma_0\subset\partial\Omega
\]

be a nonempty relatively open set meeting every connected component of
\(\partial\Omega\).  For \(\alpha\in(0,1)\), Glass assumes

\[
 y_0,y_1\in C^{2,\alpha}(\overline\Omega;\mathbb R^3),\qquad
 \operatorname{div}y_i=0,
\tag{1.1}
\]

and only the endpoint boundary condition

\[
 y_i\cdot n=0
 \quad\hbox{on }\partial\Omega\setminus\Gamma_0.
\tag{1.2}
\]

For every \(T>0\), Theorem 1.1 produces

\[
 y\in
 C([0,T];C^{1,\alpha}(\overline\Omega))
 \cap
 L^\infty([0,T];C^{2,\alpha}(\overline\Omega))
\tag{1.3}
\]

which solves Euler distributionally in time, has
\(y(0)=y_0,\ y(T)=y_1\), and satisfies

\[
 y\cdot n=0
 \quad\hbox{on }
 (\partial\Omega\setminus\Gamma_0)\times[0,T].
\tag{1.4}
\]

The pressure is asserted in the distribution space
\(\mathcal D'(\Omega\times(0,T))\).  The sign convention in the paper is
\(\partial_t y+(y\cdot\nabla)y=\nabla p\).

Three details matter.

### 1.1 This is not an analytic controllability theorem

The theorem is at \(C^{2,\alpha}\) endpoint regularity, with the trajectory
only continuous into \(C^{1,\alpha}\) and bounded into \(C^{2,\alpha}\).
It does not assert \(C^\infty\), Gevrey, or analytic dependence in space or
time.

Indeed, Proposition 1.3 states

\[
 y=p=0\qquad\hbox{for }t\in[1/2,1].
\tag{1.5}
\]

There is a minor internal discrepancy worth recording: the fixed-point
operator in equation (3.44) is explicitly set to zero on \([3/4,1]\), not
on \([1/2,1]\).  Either displayed interval gives what matters here: the
nontrivial return trajectory vanishes on a terminal time interval.  A
nonzero time-analytic trajectory cannot have this property.

There is also a spatial boundary obstruction to localized analytic control.
Suppose \(\Omega=B\), the complement of \(\Gamma_0\) contains a relatively
open subset of the connected sphere, and \(y(t,\cdot)\) is real analytic in
a neighborhood of \(\overline B\).  Then

\[
 g(t,\cdot)=y(t,\cdot)\cdot n
\]

is analytic on \(\partial B\).  Since \(g=0\) on an open subset, analytic
unique continuation on the sphere gives \(g=0\) on all of \(\partial B\).
Thus the localized through-boundary flushing used by the return method
cannot itself be analytic up to the boundary.  Taking
\(\Gamma_0=\partial B\) removes this particular observation, but does not
upgrade Glass's regularity theorem.

### 1.2 The controls are not the full velocity trace

Glass notes after Proposition 1.3 that one may regard the controls as

\[
 y\cdot n\quad\hbox{on }\Gamma_0,
\qquad
 \omega\wedge n\quad\hbox{where }y\cdot n<0,
\tag{1.6}
\]

where \(\omega=\operatorname{curl}y\).  The second datum specifies incoming
vorticity.  Neither the tangential velocity trace nor the boundary pressure
is a prescribed output of the theorem.

This distinction is harmless for bounded-domain controllability.  It is
decisive for gluing to an exterior classical solution.

### 1.3 The exterior domain in the proof is bookkeeping

The return trajectory is a potential flow

\[
 \overline y=\nabla\theta,\qquad
 \Delta\theta=0\ \hbox{in }\Omega,\qquad
 \partial_n\theta=0\ \hbox{on }\partial\Omega\setminus\Gamma_0,
\tag{1.7}
\]

whose flow flushes every chosen packet out through \(\Gamma_0\).  The proof
introduces a larger domain \(\widetilde\Omega\) and a bounded linear
extension \(\pi\), with compact support in a still larger bounded set.  For
an iterate \(u\), the advecting field outside \(\Omega\) is

\[
 \widetilde u=\overline y+\pi(u-\overline y).
\tag{1.8}
\]

Only its restriction to \(\Omega\) is divergence free.  Consequently the
auxiliary exterior vorticity is transported by

\[
 \partial_t\omega+(\widetilde u\cdot\nabla)\omega
 =
 (\omega\cdot\nabla)\widetilde u
 -
 \omega\,\operatorname{div}\widetilde u
\quad\hbox{in }\widetilde\Omega,
\tag{1.9}
\]

not by the vorticity equation of a global incompressible Euler solution.

More decisively, in equations (3.19)--(3.24) Glass decomposes the initial
vorticity into packets.  When packet \(w_i\) has been transported so that

\[
 \operatorname{supp}w_i(t_i^-)\cap\Omega=\varnothing,
\tag{1.10}
\]

the auxiliary field is redefined at \(t_i\) with \(w_i\) omitted.  The
restriction to \(\Omega\) stays continuous because the omitted packet is
already outside.  The field on \(\widetilde\Omega\) is not a continuous
unforced Euler evolution.

For multiply connected \(\Omega\), Lemma 2.3 also starts with vorticity
\(\operatorname{curl}\aleph_i\) supported in
\(\widetilde\Omega\setminus\Omega\), passes it through \(\Omega\), and
expels it again in order to adjust the harmonic circulation coordinates.
This is exactly the sort of preloaded exterior vorticity a physical wake
would need, but the paper uses it only as an auxiliary transport device.

Finally, arbitrary-time controllability is obtained from small-data
null-controllability by Euler scaling and time reversal.  If a small
trajectory on \([0,1]\) starts from \(\varepsilon y_0\), then

\[
 \widetilde y(x,t)
 =
 \varepsilon^{-1}y(x,t/\varepsilon),
\qquad 0\le t\le\varepsilon.
\tag{1.11}
\]

The endpoint is restored to \(y_0\), but the return-flow amplitude is
multiplied by \(\varepsilon^{-1}\).  The theorem therefore contains no
uniform small-control estimate as \(T\downarrow0\).

---

## 2. The exact interface problem

Let \(\Sigma=\partial\Omega\), with normal \(n\) pointing from \(\Omega\)
to its exterior.  For piecewise fields write

\[
 [f]=f^+-f^-,
\]

where minus denotes the interior trace and plus the exterior trace.

### Proposition 2.1 (Euler interface conditions)

Let \((u^\pm,q^\pm)\) be \(C^1\) Euler solutions up to the fixed interface
\(\Sigma\).  Define the piecewise field

\[
 U=u^-\mathbf1_\Omega+u^+\mathbf1_{\Omega^c},
\qquad
 Q=q^-\mathbf1_\Omega+q^+\mathbf1_{\Omega^c}.
\tag{2.1}
\]

Then \((U,Q)\) is a distributional Euler solution across \(\Sigma\) if and
only if

\[
 [u\cdot n]=0,
\qquad
 [(u\cdot n)u+qn]=0.
\tag{2.2}
\]

If the common normal velocity is \(g\), these conditions become

\[
 [q]=0,\qquad g[u_\tau]=0.
\tag{2.3}
\]

In particular, wherever the boundary control has \(g\ne0\), even a weak
vortex-sheet gluing must match the full velocity trace.  Where \(g=0\),
(2.2) permits a tangential jump, but that is a vortex sheet and not a
classical solution.

#### Proof

In conservation form Euler is

\[
 \partial_tU+\operatorname{div}(U\otimes U+QI)=0.
\tag{2.4}
\]

Distributional differentiation of the two characteristic functions
gives the surface terms

\[
 [u\cdot n]\,\delta_\Sigma
\]

in incompressibility and

\[
 [(u\cdot n)u+qn]\,\delta_\Sigma
\]

in momentum.  Their vanishing is exactly (2.2).  With equal normal traces,
the tangential and normal projections give (2.3). \(\square\)

A global \(C^1\) field additionally needs \([u]=[\nabla u]=0\).  A global
\(C^\infty\) field needs equality of every spatial and time-compatible
jet.  Pressure may be shifted by a function of time on either connected
piece, but this only adjusts one scalar constant; it cannot match a
nonconstant pressure trace.

### Proposition 2.2 (Navier--Stokes interface conditions)

For viscosity \(\nu>0\), a piecewise classical field cannot have
\([u]\ne0\), since its second derivatives then contain derivatives of
surface measures.  After imposing \([u]=0\), the distributional momentum
condition is continuity of total traction:

\[
 \left[
 \left(
 u\otimes u+qI
 -\nu(\nabla u+\nabla u^{T})
 \right)n
 \right]=0.
\tag{2.5}
\]

Since the velocity trace already matches, this reduces to

\[
 \left[
 \left(
 qI-\nu(\nabla u+\nabla u^{T})
 \right)n
 \right]=0.
\tag{2.6}
\]

A smooth gluing again requires all jets, not merely (2.6).

### The missing control-transmutation lemma

The exact lemma needed by the proposed strategy is therefore the following.

> **Exterior-wake right inverse.**  Given a controlled interior Euler path
> \((u^-,q^-)\) on \(B_\ell\times[0,\tau]\), construct a finite-energy
> exterior Euler path \((u^+,q^+)\) such that:
>
> 1. \(u^+\), \(q^+\) solve unforced Euler in
>    \((\mathbb R^3\setminus\overline B_\ell)\times[0,\tau]\);
> 2. their complete compatible jets equal those of
>    \((u^-,q^-)\) on the sphere;
> 3. the global pressure induced by \(U\otimes U\) gives exactly \(q^-\)
>    in the ball, up to a function of time;
> 4. the exterior initial vorticity, energy, helicity, and outgoing wake
>    are included in the state rather than discarded; and
> 5. the construction has scale-uniform tame bounds strong enough to
>    iterate.

If this lemma holds, Proposition 2.1 makes the gluing immediate.  But the
lemma is a nonlinear exterior Euler Cauchy problem with both velocity and
pressure data on a timelike boundary.  Glass's endpoint theorem supplies
none of items 2--5.

There is also a nonlocal pressure condition hidden in item 3.  In whole
space,

\[
 -\Delta q=\partial_i\partial_j(U_iU_j).
\tag{2.7}
\]

Changing \(U\) in a remote annulus changes the pressure inside the ball by
a harmonic function.  Preserving the prescribed interior acceleration
requires that harmonic correction to be spatially constant, equivalently
that all of its nonconstant interior multipoles vanish.  Separation of
supports does not provide exact dynamical decoupling.

---

## 3. A concrete transmutation no-go theorem

The following analytic endpoint lies within Glass's hypotheses and makes
the obstruction explicit.  On the unit ball \(B\), fix
\(a\in\mathbb R^3\setminus\{0\}\) and set

\[
 v(x)=(1+|x|^2)(a\times x).
\tag{3.1}
\]

It is polynomial, divergence free, and tangent to every sphere.  Its
vorticity is

\[
 \omega_v
 =
 2a+4|x|^2a-2(a\cdot x)x,
\qquad
 \Delta\omega_v=20a.
\tag{3.2}
\]

Thus Glass's theorem supplies, for every \(T>0\), a controlled Euler path
from \(0\) to \(v\) on \(B\), for any admissible control patch.

### Theorem 3.1 (rest-wake Euler and all-wake Navier--Stokes obstructions)

For the endpoint pair \(0\to v\):

1. there is no finite-energy smooth whole-space Euler embedding whose
   exterior initial vorticity vanishes;
2. any finite-energy whole-space Euler embedding must have initial
   exterior energy at least

   \[
   E_B(v)
   =
   \frac12\int_B|v|^2\,dx
   =
   \frac{752\pi}{945}|a|^2;
   \tag{3.3}
   \]

3. there is no smooth whole-space Navier--Stokes embedding with viscosity
   \(\nu>0\) in the standard finite-energy, rapidly decaying mild class.

#### Proof of 1

Let \(X(t,\cdot)\) be the volume-preserving flow of a smooth whole-space
Euler solution.  Cauchy's vorticity formula is

\[
 \omega(t,X(t,\alpha))
 =
 D_\alpha X(t,\alpha)\,\omega_0(\alpha).
\tag{3.4}
\]

If the interior datum is zero and the exterior initial vorticity also
vanishes, then \(\omega_0=0\) globally.  Equation (3.4) gives
\(\omega(T)=0\), contrary to (3.2).  Finite energy also excludes a
nonzero global curl-free, divergence-free constant remainder. \(\square\)

#### Proof of 2

For a smooth decaying Euler solution, total kinetic energy is conserved.
At \(t=0\), the interior energy is zero, whereas at \(t=T\) it is
\(E_B(v)\).  Hence

\[
 E_{\rm ext}(0)
 =
 E_{\rm total}(0)
 =
 E_{\rm total}(T)
 =
 E_B(v)+E_{\rm ext}(T)
 \ge E_B(v).
\tag{3.5}
\]

The value in (3.3) follows by spherical integration. \(\square\)

#### Proof of 3

A smooth Navier--Stokes solution in the standard whole-space
finite-energy/mild class is spatially real analytic at every positive
regular time; see
Grujić--Kukavica,
[Space Analyticity for the Navier--Stokes and Related Equations with
Initial Data in \(L^p\)](https://doi.org/10.1006/jfan.1997.3167).
If a whole-space solution agreed with the controlled path on \(B\), then
at time \(T>0\) its analytic velocity would equal the polynomial (3.1) on
an open set.  The identity theorem forces

\[
 U(T,x)=(1+|x|^2)(a\times x)
\quad\hbox{for every }x\in\mathbb R^3,
\]

which is not in \(L^2(\mathbb R^3)\).  This contradicts finite energy.
\(\square\)

The first conclusion is not a no-go for every imaginable Euler embedding.
It identifies the exact price of escape: the exterior initial datum must
already contain the vorticity that will later enter the ball.  Boundary
control has become global initial-data engineering.

Because \(D_\alpha X\) is invertible, (3.4) also gives the exact support
identity

\[
 \operatorname{supp}\omega(t)
 =
 X(t,\operatorname{supp}\omega_0).
\tag{3.6}
\]

Thus the preload must have the correct transported support and vortex-line
geometry, not merely enough scalar energy.  Glass's incoming-vorticity
control suppresses this global Lagrangian compatibility condition.

### Proposition 3.2 (local viscous incompatibility)

Let a smooth velocity \(u\) solve Euler on an open cylinder
\(\mathcal O\times I\).  If a smooth unforced Navier--Stokes solution with
\(\nu>0\) agrees with \(u\) there, then

\[
 \Delta\operatorname{curl}u=0
\quad\hbox{on }\mathcal O\times I.
\tag{3.7}
\]

#### Proof

Allowing different pressures, subtract Euler from Navier--Stokes:

\[
 \nu\Delta u=\nabla\phi
\tag{3.8}
\]

for some scalar distribution \(\phi\).  Taking curl gives (3.7).
\(\square\)

For (3.1), equation (3.2) violates (3.7).  Thus even without using spatial
analyticity, no smooth Euler transition which reaches this endpoint
through an open time interval can be the restriction of a
Navier--Stokes solution.  An exterior correction cannot repair a residual
\(\nu\Delta\omega\) located strictly inside the ball.

This proposition leaves the irrotational return flow
\(\overline y=\nabla\theta\) itself untouched: harmonic potential flows
have \(\omega=0\).  What fails is the generic vortical correction needed
to connect arbitrary endpoints.

---

## 4. Why annular iteration does not move the actuator

Let

\[
 A_{r,R}=\{x:r<|x|<R\}.
\]

Its boundary has two connected components.  Glass's hypothesis requires
the control set \(\Gamma_0\) to meet both of them.  This is not a technical
artifact: Remark 1.2 derives the necessity from Kelvin's circulation law.

Consequently:

* choosing \(\Gamma_0\) only on \(\{|x|=R\}\) does not satisfy the theorem;
* including \(\{|x|=r\}\) in \(\Gamma_0\) leaves a freely chosen actuator
  on the old interface;
* prescribing the full inner trace from a previous ball is much stronger
  than allowing Glass to choose controls there; and
* controlling endpoints on a larger ball does not preserve a previously
  selected trajectory on the smaller ball.

This gives a simple induction failure.  Suppose a first stage produces
\(u_0\) on \(B_{R_0}\).  To push its control to \(R_1>R_0\), the next stage
would need to solve Euler on \(A_{R_0,R_1}\) with the complete inner trace
of \(u_0\) fixed and with only the outer trace free.  That is the
exterior-wake right inverse of Section 2, not Glass's theorem.  Applying
Glass to the entire larger ball instead chooses a new interior trajectory
and generally changes \(u_0\) instantly through pressure.

An infinite-radius limiting argument has the same issue.  A compatible
sequence

\[
 U_{N+1}|_{B_{R_N}\times[0,T]}
 =
 U_N
\tag{4.1}
\]

would already define one global unforced solution on the union of the
balls.  Producing (4.1) is exactly the missing theorem.  Approximate
agreement on compact sets is insufficient for an exact construction
unless one proves convergence through the nonlinear equation, pressure,
all interface jets, and the global invariant budgets.

For Navier--Stokes the situation is stricter.  Parabolic smoothing forbids
nonzero solutions with exact open annular gaps at positive times, and
Proposition 3.2 says that moving the cutoff residual outward does not
remove the viscous residual still present in the Euler core.

### 4.1 A cutoff converts boundary control into body force

One can always choose a smooth solenoidal extension \(V\) which agrees
with the interior velocity in a collar of the ball and is cut off farther
out.  This is a useful kinematic operation, but dynamically it gives

\[
 F
 =
 \partial_tV
 +
 \mathbb P\operatorname{div}(V\otimes V)
\tag{4.2}
\]

for Euler, or

\[
 F_\nu
 =
 \partial_tV
 +
 \mathbb P\operatorname{div}(V\otimes V)
 -
 \nu\Delta V
\tag{4.3}
\]

for Navier--Stokes.  Here \(\mathbb P\) is the whole-space Leray
projection.  The unprojected mismatch can be arranged in the cutoff
collar, but \(\mathbb P\) is nonlocal.  Unless (4.2) or (4.3) vanishes
identically, the construction is a forced solution.  Replacing \(V\) by
\(V+w\) and cancelling \(F\) is precisely the nonlinear exterior matching
problem of Section 2.

There is a logically valid limiting reduction.  Suppose one could build
forced fields \(V_N\) and pressures \(Q_N\) such that

\[
 \operatorname{supp}F_N
 \subset\mathbb R^3\setminus B_{R_N},
\qquad R_N\longrightarrow\infty,
\tag{4.4}
\]

and, on every compact spacetime set,

\[
 (V_N,Q_N)\longrightarrow(U,Q)
\tag{4.5}
\]

strongly enough to pass the quadratic equation and the required
derivatives.  Then \((U,Q)\) is an unforced whole-space solution.  This is
the precise sense in which forcing may be pushed to infinity.

But (4.4)--(4.5) are the conclusion that the proposed annular iteration
must prove.  Glass gives neither a correction which preserves the old
interior trajectory nor uniform compactness as the outer radius changes.
If one additionally requires exact consistency
\(V_{N+1}=V_N\) on \(B_{R_N}\), the union is already the desired global
solution, so the formulation has not simplified the construction.

---

## 5. Quantitative cost of the Glass return flow

Although the exterior matching theorem is absent, the scaling bill of the
return method can be computed exactly.

Let \((\bar y,\bar q)\) be a fixed nonzero unit-ball, unit-time return flow.
On a ball of radius \(\ell\) and a time interval of length \(\tau\), its
Euler scaling is

\[
 y_{\ell,\tau}(x,t)
 =
 \frac{\ell}{\tau}
 \bar y\!\left(\frac{x}{\ell},\frac{t}{\tau}\right),
\qquad
 q_{\ell,\tau}(x,t)
 =
 \frac{\ell^2}{\tau^2}
 \bar q\!\left(\frac{x}{\ell},\frac{t}{\tau}\right).
\tag{5.1}
\]

At a time when \(\bar y\ne0\),

\[
\begin{aligned}
 \|y_{\ell,\tau}\|_\infty
 &\asymp \frac{\ell}{\tau},\\
 \|\operatorname{curl}y_{\ell,\tau}\|_\infty
 &\asymp \frac1{\tau},\\
 \|\nabla^m y_{\ell,\tau}\|_\infty
 &\asymp \frac{\ell^{1-m}}{\tau},\\
 \|y_{\ell,\tau}\|_2^2
 &\asymp \frac{\ell^5}{\tau^2},\\
 \|q_{\ell,\tau}\|_\infty
 &\asymp \frac{\ell^2}{\tau^2}.
\end{aligned}
\tag{5.2}
\]

The boundary power obeys the exact scaling

\[
 {\cal P}_{\ell,\tau}(t)
 :=
 \int_{\partial B_\ell}
 \left(\frac{|y|^2}{2}+q\right)y\cdot n\,dS
 =
 \frac{\ell^5}{\tau^3}
 \overline{\cal P}(t/\tau).
\tag{5.3}
\]

For a nontrivial zero-to-zero loop, the interior energy is not constant,
so \(\overline{\cal P}\) is not identically zero.  Hence the supremum
power is of order \(\ell^5/\tau^3\), and the time-integrated absolute work
is of order \(\ell^5/\tau^2\), the same as its peak kinetic energy.  The
net signed work can cancel, but a physical wake must still accommodate the
intermediate energy.

The viscous dissipation scale of the same profile is

\[
 \nu\int_0^\tau
 \|\nabla y_{\ell,\tau}(t)\|_2^2\,dt
 \asymp
 \nu\frac{\ell^3}{\tau}.
\tag{5.4}
\]

Write the active velocity amplitude as

\[
 a_\ell=\frac{\ell}{\tau}=\ell^{-\gamma},
\qquad
 \tau=\ell^{1+\gamma}.
\tag{5.5}
\]

Then

\[
 E_\ell\asymp\ell^{3-2\gamma},
\qquad
 D_\ell\asymp\nu\ell^{2-\gamma},
\qquad
 \|\omega_\ell\|_\infty\asymp\ell^{-1-\gamma}.
\tag{5.6}
\]

For geometrically shrinking, disjoint retained wakes, the scalar energy
series is summable only when

\[
 \gamma<\frac32.
\tag{5.7}
\]

At the parabolic critical scaling \(\gamma=1\),

\[
 E_\ell\asymp\ell,\qquad
 D_\ell\asymp\nu\ell,
\tag{5.8}
\]

so energy and dissipation are summable.  The scalar ledger therefore does
not exclude a critical shrinking wake.  It says that such a wake cannot be
discarded and that its \(C^m\) and vorticity norms grow rapidly.

At fixed \(\ell>0\), by contrast, compressing infinitely many Glass return
flows into \(\tau_j\downarrow0\) requires peak energy
\(\gtrsim\tau_j^{-2}\).  This is incompatible with finite conserved Euler
energy and with the nonincreasing energy of unforced Navier--Stokes.
Shrinking physical support is essential even before interface matching is
considered.

Glass gives no scale-uniform constants for the exterior-wake problem.  Its
smallness threshold depends on the domain and control geometry, while the
time rescaling (1.11) amplifies the return flow.  Thus (5.7)--(5.8) are
necessary dimensional budgets, not an iteration theorem.

---

## 6. The surviving theorem target

The Glass route becomes prize-relevant only after proving a statement of
the following form.

> **Preloaded Euler wake theorem.**  Construct one normalized, smooth
> interior transition together with an exterior Euler reservoir whose
> incoming vorticity is already present at the initial time, whose full
> interface jets and pressure match, and whose outgoing state is retained
> in a weighted annular phase space.  Prove a tame right inverse with
> constants uniform under the intended scaling, and include the wake in
> the next endpoint map.

If such a theorem existed with the bounds in (5.2), geometrically shrinking
copies with \(\gamma<3/2\) would pass the elementary energy test.  The
critical choice \(\gamma=1\) is especially compatible with a summable
finite-energy, non-tight wake.  This is a legitimate conditional
reduction; it agrees with the scalar escape ledger in
[`2026-07-29-critical-return-cell-escape-boundary.md`](2026-07-29-critical-return-cell-escape-boundary.md).
The exact stationary bubbles in
[`2026-07-29-gavrilov-dss-wake-construction.md`](2026-07-29-gavrilov-dss-wake-construction.md)
can provide storage pieces, but their zero boundary flux does not supply
the active full-trace transition required here.

What cannot be used as that theorem is:

* the auxiliary extension \(\pi\) in Glass;
* deletion of packets after they leave the ball;
* a wake initialized at zero;
* Glass controllability on shells with only outer-boundary control;
* cutoff and divergence correction followed by moving the residual to a
  farther annulus; or
* an exact Euler core inside an unforced Navier--Stokes solution, except in
  the special harmonic-vorticity class (3.7).

The conceptual gain from Glass is therefore narrower but still useful.  It
shows that a bounded Euler core can flush and replace vorticity if the
boundary is treated as an actuator and reservoir.  The global problem is
to realize that reservoir as one dynamically compatible, finite-energy
Euler state.  For Navier--Stokes, the Euler transition cannot simply be
reused: a genuinely viscous full-space transition must be constructed
from the start.

## Sources

1. O. Glass, [Exact boundary controllability of 3-D Euler
   equation](https://www.numdam.org/item/COCV_2000__5__1_0.pdf), ESAIM:
   COCV **5** (2000), 1--44.
2. Z. Grujić and I. Kukavica,
   [Space Analyticity for the Navier--Stokes and Related Equations with
   Initial Data in \(L^p\)](https://doi.org/10.1006/jfan.1997.3167),
   Journal of Functional Analysis **152** (1998), 447--466.
