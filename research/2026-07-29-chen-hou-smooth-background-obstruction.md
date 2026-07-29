# A smooth terminal background cannot viscously dress the Chen--Hou core

Date: 2026-07-29

## Result and claim boundary

This note tests a specific possible bridge from the Chen--Hou smooth
boundary Euler singularity to Clay alternative (D):

1. retain the Chen--Hou field \(v\) as the leading shrinking core;
2. add a divergence-free background \(b\) which is smooth through the
   terminal time;
3. use the cross interactions with \(b\), a pressure adjustment, and a
   smooth periodic force to cancel \(-\nu\Delta v\); and
4. smoothly remove the wall and embed the result in a flat three-torus.

There is a sharp obstruction.  The exact Euler/Navier--Stokes subtraction
identity is

\[
\boxed{
 -\nu\Delta v+(b\cdot\nabla)v+(v\cdot\nabla)b+\nabla q
 =
 f-b_t-(b\cdot\nabla)b+\nu\Delta b .
}
\tag{0.1}
\]

At a tubular core of transverse length \(\ell\) and velocity amplitude
\(a\), divide (0.1) by the viscous size \(a/\ell^2\).  A terminally smooth
background contributes only at relative orders \(\ell\) and \(\ell^2\);
a bounded force contributes at relative order \(\ell^2/a\).  Therefore,
whenever

\[
 \ell\longrightarrow0,\qquad \frac{\ell^2}{a}\longrightarrow0,
\tag{0.2}
\]

every strong one-profile limit must solve the homogeneous transverse
Stokes equation.  Equivalently, its transverse vorticity must be harmonic.
For a smooth full-plane periodic embedding, a bounded decaying harmonic
vorticity profile is zero.

Chen--Hou's meridional scales obey

\[
 \ell=C_l\asymp s^\beta,\qquad
 a_{\rm mer}=\frac{C_l}{C_\omega}\asymp s^{\beta-1},\qquad
 \beta=2.920561\ldots,\qquad s=T-t,
\tag{0.3}
\]

so

\[
 \frac{\ell^2}{a_{\rm mer}}
 =C_lC_\omega
 \asymp s^{\beta+1}\longrightarrow0.
\tag{0.4}
\]

Thus a smooth terminal background cannot cancel viscosity while preserving
a nonzero localized full-plane Chen--Hou profile.  More quantitatively, a
non-Stokes profile forces

\[
 \|f(t)\|_{L^\infty}
 \gtrsim
 \nu\frac{a_{\rm mer}}{\ell^2}
 =
 \frac{\nu}{C_lC_\omega}
 \asymp
 \nu s^{-\beta-1}
 =
 \nu s^{-3.920561\ldots}
\tag{0.5}
\]

along a sequence approaching blowup.  Such a force is not terminally
smooth.

The only leading-order loophole is a transverse harmonic/Stokes profile.
It is nontrivial on a half-plane only because boundary data are available.
That loophole does not produce a torus construction: the actual Chen--Hou
wall trace has no smooth geometric reflection, while any hypothetical
smooth decaying full-plane extension is killed by the harmonic Liouville
theorem.

The theorem below assumes a strong tubular profile limit.  Chen--Hou prove
a nearly self-similar boundary singularity in weighted
\(L^\infty\)-Hölder spaces, not a periodic Navier--Stokes solution.  The
claim here is conditional and route-closing:

> any proposed periodic conversion which retains their core in the stated
> \(H^1_{\rm loc}\) profile sense cannot be completed by a background and
> force smooth through \(T\).

A non-precompact boundary layer, a second singular scale, or a background
with unbounded derivatives lies outside the theorem.  Each is a new
singular mechanism, not a smooth dressing of the proved Euler core.  No
Navier--Stokes singularity or Millennium resolution is claimed.

---

## 1. Primary-source scales and geometry

Chen and Hou consider smooth finite-energy axisymmetric Euler flow in

\[
 0\le r\le1,\qquad z\in\mathbb T,
\tag{1.1}
\]

with an impermeable solid boundary at \(r=1\).  Their variables are

\[
 \Gamma=ru^\theta,\qquad
 \theta=\Gamma^2,\qquad
 \omega=\frac{\omega^\theta}{r}.
\tag{1.2}
\]

Part I, equations (6.8)--(6.11), uses

\[
 X=C_l(\tau)^{-1}z,\qquad
 Y=C_l(\tau)^{-1}(1-r),
\tag{1.3}
\]

and

\[
 \begin{aligned}
 \omega_{\rm ss}(X,Y,\tau)
 &=
 C_\omega(\tau)\,
 \omega_{\rm phy}(1-C_lY,C_lX,t(\tau)),\\
 C_\omega(\tau)
 &=
 C_\omega(0)\exp\left(\int_0^\tau c_\omega(q)\,dq\right),\\
 C_l(\tau)
 &=
 C_l(0)\exp\left(-\int_0^\tau c_l(q)\,dq\right),\\
 \frac{dt}{d\tau}&=C_\omega(\tau).
 \end{aligned}
\tag{1.4}
\]

The rigorously enclosed approximate-profile parameters are

\[
 \bar c_l\approx3.00649898,\qquad
 \bar c_\omega\approx-1.02942516.
\tag{1.5}
\]

Consequently, with \(s=T-t\),

\[
 C_\omega\asymp s,\qquad
 C_l\asymp s^\beta,\qquad
 \beta=-\frac{\bar c_l}{\bar c_\omega}
 =2.920561\ldots.
\tag{1.6}
\]

The numerical exponent records the enclosed nearly steady profile.  The
theorem below does not require an exact power-law limit: it uses only
\(C_l\to0\), \(C_\omega\to0\), and the exact dynamic-scale identities in
(1.7).

The physical component amplitudes are

\[
 \begin{aligned}
 \omega^\theta/r&\asymp s^{-1},\\
 u_{\rm mer}&\asymp \frac{C_l}{C_\omega}
 =s^{\beta-1},\\
 \Gamma,\ u^\theta&\asymp
 \left(\frac{C_l}{C_\omega^2}\right)^{1/2}
 =s^{\beta/2-1}.
 \end{aligned}
\tag{1.7}
\]

Both active meridional directions, \(z\) and \(1-r\), have the common
length \(C_l\).  The azimuthal direction has length \(O(1)\), because the
solution is axisymmetric.  The blowup set is the boundary circle
\(\{r=1,z=0\}\).

The exact rescaled elliptic equation in Part I, equation (6.13), is

\[
 -(\partial_{XX}+\partial_{YY})\phi
 +\frac{C_l}{r}\partial_Y\phi
 +\frac{C_l^2}{r^2}\phi
 =r\omega,
 \qquad r=1-C_lY,
\tag{1.8}
\]

with \(\phi(X,0)=0\).  This displays both facts needed below:

* the frozen principal cross-sectional geometry is the flat half-plane;
* cylindrical curvature enters one and two powers of \(C_l\) later.

Part II rigorously verifies the stability constants used in Part I.  It
does not introduce a viscous profile or a Navier--Stokes stability
estimate.

## 2. Exact subtraction from a smooth background

Let \(v\) be a divergence-free Euler path,

\[
 v_t+(v\cdot\nabla)v+\nabla p_E=0.
\tag{2.1}
\]

Let \(b\) be another divergence-free path and set

\[
 u=v+b.
\tag{2.2}
\]

If \(u\) solves forced Navier--Stokes,

\[
 u_t+(u\cdot\nabla)u+\nabla p-\nu\Delta u=f,
\tag{2.3}
\]

subtracting (2.1) gives (0.1), with \(q=p-p_E\).  This identity has no
asymptotic or projection error.

More generally, if a cutoff or periodicization makes \(v\) solve

\[
 v_t+(v\cdot\nabla)v+\nabla p_E=e,
\tag{2.4}
\]

then the right side of (0.1) gains \(-e\).  A defect \(e\) which is
bounded through \(T\) has exactly the same lower-order status as a smooth
force.  A defect of size \(a/\ell^2\) can cancel viscosity, but is then
itself the forbidden singular force.

Suppose

\[
 b\in C^1([0,T];C^2),\qquad
 e\in C([0,T];L^\infty),\qquad
 f\in C([0,T);L^\infty).
\tag{2.5}
\]

Then

\[
 G:=-e-b_t-(b\cdot\nabla)b+\nu\Delta b
\tag{2.6}
\]

is bounded through \(T\), and (0.1) becomes a stationary Stokes equation
at every fixed time:

\[
 -\nu\Delta v+(b\cdot\nabla)v+(v\cdot\nabla)b+\nabla q=f+G.
\tag{2.7}
\]

The key point is the order of this elliptic equation.  The term to be
cancelled has two derivatives of \(v\).  A smooth background supplies at
most one derivative of \(v\).

In flat rescaled coordinates \(x_\perp=\ell y\), \(v=aV\), equation (2.7)
has the exact schematic normalization

\[
 -\nu\Delta_yV
 +\ell\,b\cdot\nabla_yV
 +\ell^2(V\cdot\nabla)b
 +\nabla_yQ
 =
 \frac{\ell^2}{a}(f+G).
\tag{2.8}
\]

Smooth metric and cylindrical connection terms are \(O(\ell)\) relative
to the principal Laplacian.  Formula (2.8) already shows why neither a
smooth background nor pressure can remove a non-Stokes principal profile.

## 3. A tubular Stokes-limit theorem

The straight periodic tube states the argument without geometric
notation.  Curved rings are treated immediately afterward.

Let

\[
 \Sigma=\{x_1=x_2=0\}\subset\mathbb T^3
\]

be a periodic straight curve, and write \(x_\perp=(x_1,x_2)\) and
\(\sigma=x_3\).

> **Theorem 3.1 (smooth-background Stokes gate).**
> Let \(t_n\uparrow T\), \(\ell_n\downarrow0\), and \(a_n>0\) satisfy
> \[
>  \frac{\ell_n^2}{a_n}\longrightarrow0.
> \tag{3.1}
> \]
> Let \(v,b,u,f,e\) obey (2.2)--(2.5) on \(\mathbb T^3\).  Assume that,
> for every \(R<\infty\),
> \[
>  V_n(y,\sigma)
>  :=
>  a_n^{-1}v(\ell_n y,\sigma,t_n)
>  \longrightarrow V(y)
> \tag{3.2}
> \]
> in \(H^1(B_R\times\mathbb T_\sigma)\), where \(V\) is independent of
> \(\sigma\).
>
> If \(f\) is bounded along \(t_n\), then
> \[
>  \int_{\mathbb R^2}\nabla_yV:\nabla_y\Phi\,dy=0
> \tag{3.3}
> \]
> for every compactly supported \(\sigma\)-independent
> divergence-free test field \(\Phi\).
>
> Equivalently, \(V\) solves the transverse homogeneous Stokes system
> \[
>  -\Delta_yV+\nabla_yQ=0,\qquad
>  \nabla_y\cdot V_\perp=0
> \tag{3.4}
> \]
> in distributions.  In particular,
> \[
>  \Delta_y(\partial_{y_1}V_2-\partial_{y_2}V_1)=0,
>  \qquad
>  \Delta_yV_3=0.
> \tag{3.5}
> \]

### Proof

Fix a compactly supported smooth \(\Phi(y)\) satisfying

\[
 \partial_{y_1}\Phi_1+\partial_{y_2}\Phi_2=0.
\]

The field

\[
 \Phi_n(x)=\Phi(x_\perp/\ell_n)
\tag{3.6}
\]

is an exactly divergence-free periodic test field for all sufficiently
large \(n\).  Test (2.7) at time \(t_n\) against \(\Phi_n\), integrate the
viscous term by parts, and divide by \(a_n\).  The pressure disappears.

Changing variables \(x_\perp=\ell_ny\) gives

\[
 \frac{\nu}{a_n}\int\nabla v:\nabla\Phi_n\,dx
 \longrightarrow
 \nu|\mathbb T_\sigma|
 \int_{\mathbb R^2}\nabla_yV:\nabla_y\Phi\,dy.
\tag{3.7}
\]

The background transport term can be integrated by parts because
\(\operatorname{div}b=0\):

\[
 \begin{aligned}
 \frac1{a_n}
 \left|\int(b\cdot\nabla)v\cdot\Phi_n\,dx\right|
 &=
 \frac1{a_n}
 \left|\int v\cdot(b\cdot\nabla)\Phi_n\,dx\right|\\
 &\le C_\Phi\ell_n.
 \end{aligned}
\tag{3.8}
\]

Local \(H^1\) convergence in (3.2) supplies the uniform local \(L^1\)
bound used here.  Similarly,

\[
 \frac1{a_n}
 \left|\int(v\cdot\nabla)b\cdot\Phi_n\,dx\right|
 \le C_\Phi\ell_n^2,
\tag{3.9}
\]

while

\[
 \frac1{a_n}
 \left|\int G\cdot\Phi_n\,dx\right|
 \le
 C_\Phi\frac{\ell_n^2}{a_n}
 \longrightarrow0.
\tag{3.10}
\]

If \(f\) is bounded along the sequence, it obeys the same estimate.
Passing to the limit proves (3.3).  The distributional de Rham theorem
gives the momentum equation in (3.4).  Rescaling
\(\operatorname{div}v=0\) gives
\(\nabla_y\cdot(V_n)_\perp+\ell_n\partial_\sigma(V_n)_3=0\), so the
assumed \(H^1\) convergence gives the divergence equation in (3.4).
Taking the transverse curl gives the first equation in
(3.5); testing with \((0,0,\phi)\) gives the second. \(\square\)

### Corollary 3.2 (localized full-plane no-go)

Under the hypotheses of Theorem 3.1, suppose that the transverse
vorticity

\[
 \Omega=\partial_{y_1}V_2-\partial_{y_2}V_1
\tag{3.11}
\]

is bounded, tends to zero as \(|y|\to\infty\), and is nonzero.  Then no
bounded terminal force exists.

Indeed, (3.5) makes \(\Omega\) a bounded harmonic function on
\(\mathbb R^2\).  Liouville's theorem and its decay give
\(\Omega=0\), a contradiction.

The same proof applies to a smooth closed curve in a flat torus.  In
Fermi coordinates the volume form and inverse metric are their frozen
values plus \(O(\ell_n)\).  The divergence-free test (3.6) has a local
correction of relative size \(O(\ell_n)\).  All curvature and frame
terms therefore join the vanishing right sides of (3.8)--(3.9).  The
limit remains (3.3).

## 4. Quantitative force lower bound

The proof gives more than a compactness contradiction.  Suppose \(V\)
does not satisfy (3.3).  Choose a compact divergence-free \(\Phi\) for
which

\[
 J_\Phi
 :=
 \int_{\mathbb R^2}\nabla_yV:\nabla_y\Phi\,dy
 \ne0.
\tag{4.1}
\]

Equations (3.7)--(3.10) give

\[
 \frac1{a_n}\int_{\mathbb T^3}f(t_n,x)\cdot\Phi_n(x)\,dx
 =
 \nu|\mathbb T_\sigma|J_\Phi+o(1).
\tag{4.2}
\]

Since

\[
 \|\Phi_n\|_{L^1(\mathbb T^3)}
 =
 \ell_n^2|\mathbb T_\sigma|
 \|\Phi\|_{L^1(\mathbb R^2)},
\]

one obtains

\[
\boxed{
 \liminf_{n\to\infty}
 \frac{\ell_n^2}{a_n}
 \|f(t_n)\|_{L^\infty}
 \ge
 \nu\frac{|J_\Phi|}{\|\Phi\|_{L^1}}>0.
}
\tag{4.3}
\]

Thus the force needed to preserve a non-Stokes core has the sharp
principal size

\[
 \|f(t_n)\|_\infty\gtrsim\nu a_n\ell_n^{-2}.
\tag{4.4}
\]

The estimate also quantifies what a purported background would have to
do.  Pointwise balance in (2.7) requires at least one of

\[
 \|b\|_\infty\gtrsim\frac{\nu}{\ell},\qquad
 \|\nabla b\|_\infty\gtrsim\frac{\nu}{\ell^2},\qquad
 \|b_t\|_\infty+\nu\|\Delta b\|_\infty
 \gtrsim\frac{\nu a}{\ell^2}.
\tag{4.5}
\]

Every alternative violates terminal smoothness as \(\ell\to0\), unless
the leading transverse Stokes defect is already zero.

## 5. Application to the Chen--Hou exponents

For the meridional velocity, take the exact dynamic scales

\[
 \ell=C_l,\qquad
 a=a_{\rm mer}=\frac{C_l}{C_\omega}.
\tag{5.1}
\]

Thus

\[
 \frac{\ell^2}{a_{\rm mer}}
 =C_lC_\omega
 \asymp s^{\beta+1}\longrightarrow0.
\tag{5.1a}
\]

The terms in (2.7) then have the sizes

\[
\begin{array}{c|c|c}
\text{term}&\text{physical size}&
\text{ratio to }\nu a_{\rm mer}/\ell^2\\ \hline
\nu\Delta v_{\rm mer}
 &\nu s^{-\beta-1}&1\\
(b\cdot\nabla)v_{\rm mer}
 &O(s^{-1})&O(\nu^{-1}s^\beta)\\
(v_{\rm mer}\cdot\nabla)b
 &O(s^{\beta-1})&O(\nu^{-1}s^{2\beta})\\
f+G\ \text{bounded}
 &O(1)&O(\nu^{-1}s^{\beta+1}).
\end{array}
\tag{5.2}
\]

All background and force terms are asymptotically lower order.  In
particular, (4.4) becomes

\[
 \|f(t)\|_\infty
 \gtrsim\frac{\nu}{C_lC_\omega}
 \asymp\nu s^{-\beta-1}
 =
 \nu s^{-3.920561\ldots}.
\tag{5.3}
\]

The swirl amplitude is larger than the meridional amplitude:

\[
 a_\theta
 =\frac{C_l^{1/2}}{C_\omega}
 \asymp s^{\beta/2-1},\qquad
 \frac{a_\theta}{a_{\rm mer}}
 =C_l^{-1/2}
 \asymp s^{-\beta/2}.
\tag{5.4}
\]

This does not spoil the meridional test.  In exact cylindrical geometry,
the axisymmetric vector Laplacian preserves the meridional and azimuthal
component blocks.  Under a general smooth tubular frame, the worst
first-order frame coupling has relative size

\[
 \ell\frac{a_\theta}{a_{\rm mer}}
 =
 C_l^{1/2}
 \asymp s^{\beta/2}\longrightarrow0,
\tag{5.5}
\]

and a zeroth-order background coupling has relative size

\[
 \ell^2\frac{a_\theta}{a_{\rm mer}}
 =
 C_l^{3/2}
 \asymp s^{3\beta/2}\longrightarrow0.
\tag{5.6}
\]

Thus the different swirl amplitude cannot cancel the principal
meridional Laplacian while the added background remains smooth.

The same ledger in the azimuthal component is even more singular:

\[
 \nu\frac{a_\theta}{\ell^2}
 =
 \frac{\nu}{C_\omega C_l^{3/2}}
 \asymp
 \nu s^{-3\beta/2-1}
 =
 \nu s^{-5.380841\ldots}.
\tag{5.7}
\]

Only one non-Stokes component is needed for the contradiction, so (5.3)
already closes the smooth-background route.

## 6. The harmonic and lower-order loopholes

Theorem 3.1 isolates the only possible leading cancellation:

\[
 -\Delta_yV+\nabla_yQ=0.
\tag{6.1}
\]

This criterion is exact.  It includes any leading cancellation that might
be called "the Laplacian is harmonic" or "viscosity is only pressure."
Because \(\Delta v\) is divergence-free whenever \(v\) is, a decaying
full-plane instance of (6.1) has no hidden pressure freedom: its vorticity
is harmonic and zero.

There are three distinct settings.

### 6.1 Full-plane localized profile

A smooth periodic embedding produces a full transverse plane after
blowup.  If its profile is bounded and decays in the similarity far field,
(6.1) is trivial by Corollary 3.2.  This is the setting needed for the
proposed torus conversion.

### 6.2 Half-plane profile with boundary data

On \(\mathbb R\times\mathbb R_+\), nonzero decaying harmonic functions can
be generated by boundary values.  Therefore the harmonic Liouville step
alone does not exclude a wall profile.  But a half-plane solution is not a
Clay periodic solution; the wall must still be removed.

Chen--Hou impose only impermeability through their stream-function
condition.  Their Part I analysis explicitly uses a vorticity profile
whose odd elliptic extension across the wall is not Hölder in the normal
direction.  A geometric reflection of velocity would require normal
velocity odd and tangential velocity even, hence tangential vorticity zero
at the wall.  The Chen--Hou wall vorticity is nonzero for neighboring
rescaled axial points.  Thus their Green-function extension is not a
smooth Euler or Navier--Stokes extension.

If one replaces the jump by a boundary layer, then either:

* the layer has thickness comparable to \(\ell\) and changes the leading
  profile; or
* its thickness is \(o(\ell)\), in which case the rescaled \(H^1\)
  compactness in Theorem 3.1 fails and its second derivatives are even
  larger.

Both alternatives are new singular mechanisms.

### 6.3 Curvature after a hypothetical harmonic principal profile

The exact ordinary axisymmetric Navier--Stokes diffusion operators are

\[
\begin{aligned}
 D_\Gamma
 &=
 \partial_{rr}-\frac1r\partial_r+\partial_{zz},\\
 D_\omega
 &=
 \partial_{rr}+\frac3r\partial_r+\partial_{zz}.
\end{aligned}
\tag{6.2}
\]

With

\[
 r=1-\ell Y,\qquad z=\ell X,
\]

they become

\[
\begin{aligned}
 D_\Gamma
 &=
 \ell^{-2}\Delta_{X,Y}
 +\frac1{r\ell}\partial_Y,\\
 D_\omega
 &=
 \ell^{-2}\Delta_{X,Y}
 -\frac3{r\ell}\partial_Y.
\end{aligned}
\tag{6.3}
\]

Thus curvature is smaller than the principal heat operator by exactly one
power of \(\ell\).  A smooth background drift also first appears at that
order.  If one formally assumes
\(\Delta_{X,Y}\Gamma=\Delta_{X,Y}\omega=0\), and attempts to cancel the
remaining curvature solely by a smooth axisymmetric meridional drift with
no other same-order profile coupling, coefficientwise cancellation would
require

\[
 b^r=-\frac{\nu}{r}\quad\text{in the \(\Gamma\) equation},
\qquad
 b^r=\frac{3\nu}{r}\quad\text{in the \(\omega\) equation}.
\tag{6.4}
\]

These requirements are incompatible wherever the corresponding normal
derivatives are both active.  They are also individually incompatible
with impermeability at \(r=1\), which requires \(b^r=0\) if the background
is to preserve the wall.

Equation (6.4) is only a conditional lower-order test; this note does not
assert that the Chen--Hou profiles are transversely harmonic.  Its purpose
is to show that even the half-plane harmonic loophole does not turn the
curvature terms into one freely tunable scalar correction.

## 7. Why a same-scale corrector removes the core

One might drop terminal smoothness of \(b\) but keep its amplitude small.
At the same transverse scale, write

\[
 b_{\rm in}=aB(x_\perp/\ell).
\tag{7.1}
\]

Linear cancellation of the leading operator requires the Stokes system

\[
 -\Delta_y(V+B)+\nabla_yQ=0,\qquad
 \nabla_y\cdot(V+B)_\perp=0.
\tag{7.2}
\]

If \(V+B\) is bounded and decays on the full transverse plane, then

\[
 B=-V.
\tag{7.3}
\]

Indeed, the transverse curl and the tangential component are decaying
harmonic functions.  They vanish; divergence, curl, and decay then force
the meridional component to vanish as well.

Thus a same-frequency linear viscous corrector cancels the leading
Chen--Hou core itself.  It is not a dressing which preserves the core.

A higher-frequency corrector can have small velocity amplitude and large
second derivatives.  Such a family has no strong one-profile
\(H^1_{\rm loc}\) limit and falls outside Theorem 3.1.  It is precisely a
new non-precompact cascade or oscillatory bath.  Its nonlinear
interactions, energy, and terminal force must be constructed independently;
none is supplied by the Chen--Hou theorem.

## 8. Constructive criterion and remaining escape routes

The exact criterion furnished by this note is:

> **Smooth-background criterion.**  
> A super-parabolic shrinking Euler core with
> \(\ell^2/a\to0\) can be converted into a forced Navier--Stokes path by a
> terminally smooth background and bounded force only if every strong
> transverse profile lies in the kernel of the frozen Stokes operator.
> For a localized smooth full-plane profile, that kernel has zero
> decaying vorticity.

For the Chen--Hou exponents the criterion fails for every nonzero localized
periodic profile.  The surviving possibilities all change the leading
problem:

1. a background with velocity \(O(\nu/\ell)\) or gradient
   \(O(\nu/\ell^2)\), hence not smooth through \(T\);
2. a force of size at least (5.3), hence not a Clay force;
3. a same-scale corrector cancelling the Chen--Hou profile itself;
4. an \(o(\ell)\) boundary layer or higher-frequency bath with no strong
   profile compactness;
5. a new parabolic-scale viscous core rather than the
   \(s^{2.920561\ldots}\) Euler core; or
6. retention of a physical boundary, which is outside periodic
   alternative (D).

This is a no-go theorem for a smooth terminal dressing, not for
Navier--Stokes singularity formation in general.

## Primary sources

* J. Chen and T. Y. Hou,
  [*Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler
  equations with smooth data I: Analysis*](https://arxiv.org/abs/2210.07191),
  especially equations (2.23), (6.6)--(6.13), and Theorem 4.
* J. Chen and T. Y. Hou,
  [*Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler
  equations with smooth data II: Rigorous
  Numerics*](https://arxiv.org/abs/2305.05660).
* C. Fefferman,
  [*Existence and smoothness of the Navier--Stokes
  equation*](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf),
  especially periodic alternative (D) and its smooth-force condition.
