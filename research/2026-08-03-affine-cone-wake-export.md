# Affine-cone wake export on \(\mathbb R^3\): exact linear decay and a finite-energy displacement barrier

**Date:** 2026-08-03
**Status:** exact global-affine linear propagator, inviscid core-decoupling
estimate, and finite-energy material-displacement no-go; viscous Gaussian
leakage estimate self-derived; no nonlinear finite-energy export theorem.
**Scope:** a possible one-way-spatial-export repair of C90.  Nothing here
constructs a Navier--Stokes singularity or replaces the missing localized
transition theorem.

## 1. Outcome

There is a clean mechanism that the periodic wake problem cannot use.  On
\(\mathbb R^3\), linearize about the global trace-free affine strain

\[
 U_*(x)=Sx,\qquad
 S=\operatorname{diag}(-\alpha,-\beta,\gamma),\qquad
 \gamma=\alpha+\beta,\qquad \alpha,\beta>0.                 \tag{1.1}
\]

An inviscid vorticity packet initially contained in the positive outgoing
cone \(x_3\ge R\) is carried to \(x_3\ge e^{\gamma t}R\).  Its \(L^1\)
vorticity can grow like \(e^{\gamma t}\), but the Biot--Savart kernel gains
two powers of distance.  Consequently, for every \(m\ge0\),

\[
 \|\nabla^m u(t)\|_{L^\infty(B_r)}
 \le C_m R^{-2-m}e^{-(1+m)\gamma t}\|\omega _0\|_{L^1},
 \qquad 0<r\le R/2.                                        \tag{1.2}
\]

This is a genuine gain-window-uniform wake-to-core bound: its constant is
independent of the terminal time, even if the affine amplification interval
grows like \(j^2\) or longer.  It is stronger than a bounded wake
propagator because the old wake becomes exponentially invisible to the
active core.

The global \(L^1\) hypothesis is not essential.  A non-tight wake may instead
be measured shellwise by

\[
 \|\omega\|_{{\cal W}_m(R)}
 :=\int_{x_3\ge R}x_3^{-2-m}|\omega(x)|\,dx.                \tag{1.3}
\]

The right side of (1.2) can be replaced by
\(C_me^{-(1+m)\gamma t}\|\omega_0\|_{{\cal W}_m(R)}\).
This is the more relevant version for a critical noncompact wake whose total
vorticity need not lie in \(L^1\).

Constant viscosity does not destroy the mechanism, but it destroys exact
support.  The exact solution is an anisotropic Gaussian plume.  A standard
Gaussian split gives, schematically,

\[
 \|\nabla^m u(t)\|_{L^\infty(B_r)}
 \le C_{m,S}\|\omega _0\|_{L^1}
 \left[
 R^{-2-m}e^{-(1+m)\gamma t}
 +\ell_\nu^{-2-m}e^{-c_S(R-r)^2/\ell_\nu^2}
 \right],                                                   \tag{1.4}
\]

where

\[
 \ell_\nu=max\left\{
 \sqrt{\nu/\alpha},\sqrt{\nu/\beta},\sqrt{\nu/\gamma}
 \right\}.                                                   \tag{1.5}
\]

The second term is the time-uniform Gaussian leakage floor.  It is
exponentially small at high P\'eclet number \(R^2/\ell_\nu^2\).

The decisive limitation is geometric, not algebraic: the global affine
field \(Sx\) has infinite energy.  Section 7 proves more than a missing
construction.  For every fixed positive-volume material set in a uniformly
finite-energy flow, its average distance can grow at most linearly in time;
literal exponential export for arbitrarily long times is impossible.  If
the affine field is merely cut off to a finite active core, the wake also
stops receiving exponential outward transport as soon as it exits that
core, while the cutoff creates a collar error and a global pressure field.
Thus (1.2)--(1.4) identify a viable **model estimate**, but a finite-energy
repair must use weaker separation, exponentially small transported volume,
or a non-material cancellation/storage mechanism.

## 2. Exact inviscid propagator

Let \(u\) be a divergence-free perturbation and
\(\omega=\nabla\times u\).  Linearized Navier--Stokes about (1.1) has
vorticity equation

\[
 \partial_t\omega+(Sx\cdot\nabla)\omega-S\omega=\nu\Delta\omega.
                                                                    \tag{2.1}
\]

For \(\nu=0\), characteristics are

\[
 X(t,a)=e^{tS}a
 =\big(e^{-\alpha t}a_1,e^{-\beta t}a_2,e^{\gamma t}a_3\big),
                                                                    \tag{2.2}
\]

and the exact solution is

\[
 \boxed{\omega(t,x)=e^{tS}\omega _0(e^{-tS}x).}             \tag{2.3}
\]

Because \(\operatorname{tr}S=0\), \(\det e^{tS}=1\).  Thus each component
obeys the exact identity

\[
 \|\omega_i(t)\|_{L^1}=e^{s_i t}\|\omega_{0,i}\|_{L^1},
 \qquad (s_1,s_2,s_3)=(-\alpha,-\beta,\gamma),              \tag{2.4}
\]

and, for the Euclidean vector norm,

\[
 \|\omega(t)\|_{L^1}\le e^{\gamma t}\|\omega_0\|_{L^1}. \tag{2.5}
\]

The symmetry of \(S\) makes (2.1) preserve
\(\nabla\cdot\omega=0\).  Therefore a divergence-free initial vorticity
still recovers a divergence-free velocity through the whole-space
Biot--Savart law.

## 3. Exact cone export

For \(R>0\) and \(\kappa>0\), define the positive cone

\[
 \Gamma_{R,\kappa}^{+}
 =\{a:a_3\ge R,\ |a_1|\le\kappa a_3,
                    \ |a_2|\le\kappa a_3\}.                \tag{3.1}
\]

If \(\operatorname{supp}\omega_0\subset\Gamma_{R,\kappa}^{+}\), then
(2.2) gives

\[
 \begin{aligned}
 x_3&\ge e^{\gamma t}R,\\
 {|x_1|\over x_3}&\le
       \kappa e^{-(\alpha+\gamma)t},\\
 {|x_2|\over x_3}&\le
       \kappa e^{-(\beta+\gamma)t}.
 \end{aligned}                                               \tag{3.2}
\]

The packet not only moves outward; its physical cone narrows.  The same
statement holds in the negative cone \(a_3\le-R\).  A two-sided outgoing
wake may therefore be placed in the union of the positive and negative
cones without changing the estimate below.

## 4. Biot--Savart beats vortex stretching

On \(\mathbb R^3\),

\[
 u(x)={1\over4\pi}\int_{\mathbb R^3}
             {\omega(y)\times(x-y)\over|x-y|^3}\,dy.        \tag{4.1}
\]

Its differentiated kernel satisfies

\[
 |\nabla^mK(x-y)|\le C_m|x-y|^{-2-m}.                       \tag{4.2}
\]

Take \(|x|\le r\le R/2\).  Equations (3.2) and (4.2) imply

\[
 \begin{aligned}
 |\nabla^m u(t,x)|
 &\le C_m(e^{\gamma t}R-r)^{-2-m}\|\omega(t)\|_{L^1}\\
 &\le C_mR^{-2-m}e^{-(1+m)\gamma t}\|\omega_0\|_{L^1}.
 \end{aligned}                                               \tag{4.3}
\]

This proves (1.2).  The exponent has a transparent ledger:

\[
 \underbrace{e^{\gamma t}}_{L^1\text{ stretching}}
 \underbrace{e^{-(2+m)\gamma t}}_{\text{distance and kernel}}
 =e^{-(1+m)\gamma t}.                                       \tag{4.4}
\]

No spectral gap and no viscosity are used.  In particular, if
\({\cal E}_T\) denotes the terminal map from incoming cone vorticity to the
velocity jet in \(B_r\), then

\[
 \|{\cal E}_T\|_{L^1_\omega(\Gamma_{R,\kappa})
             \to C^m(B_r)}
 \le C_mR^{-2-m}e^{-(1+m)\gamma T},                         \tag{4.5}
\]

and

\[
 \sup_{T\ge0}\|{\cal E}_T\|\le C_mR^{-2-m}.               \tag{4.6}
\]

This is the precise linear one-way-export analogue of the uniform
wake-restricted constant requested in C90.

There is also a non-tight version.  Before taking the crude minimum distance
in (4.3), retain the initial outgoing coordinate inside the integral.  This
gives

\[
 \|\nabla^m u(t)\|_{L^\infty(B_r)}
 \le C_me^{-(1+m)\gamma t}
       \int_{a_3\ge R}a_3^{-2-m}|\omega_0(a)|\,da.          \tag{4.7}
\]

For example, a dyadic shell at radius \(D\) with velocity size
\(D^{-\delta}\) has vorticity \(L^1\) size of order
\(D^{2-\delta}\), so its contribution to the \(m=0\) integral is
\(O(D^{-\delta})\).  Geometrically separated outer shells are summable for
every \(\delta>0\).  Thus the estimate does not require the full retained
wake to have finite total \(L^1\) vorticity.  It does require every relevant
shell to lie in the same outgoing cone; this is an additional geometric
constraint not supplied by the present Gavrilov storage construction.

## 5. Constant viscosity: exact anisotropic Gaussian

Put \(y=e^{-tS}x\) and write

\[
 \omega(t,x)=e^{tS}W(t,y).                                  \tag{5.1}
\]

Then (2.1) becomes

\[
 \partial_tW
 =\nu\sum_{i=1}^3 e^{-2s_it}\partial_{y_i}^2W.             \tag{5.2}
\]

Define

\[
 A_i(t)=\int_0^te^{-2s_i\tau}\,d\tau,\qquad
 B_i(t)=e^{2s_it}A_i(t).                                    \tag{5.3}
\]

Explicitly,

\[
 \begin{array}{c|c|c}
 i&A_i(t)&B_i(t)\\ \hline
 1&\dfrac{e^{2\alpha t}-1}{2\alpha}
   &\dfrac{1-e^{-2\alpha t}}{2\alpha}\\[2mm]
 2&\dfrac{e^{2\beta t}-1}{2\beta}
   &\dfrac{1-e^{-2\beta t}}{2\beta}\\[2mm]
 3&\dfrac{1-e^{-2\gamma t}}{2\gamma}
   &\dfrac{e^{2\gamma t}-1}{2\gamma}.
 \end{array}                                                 \tag{5.4}
\]

Let

\[
 H_t(x,z)=\prod_{i=1}^3{1\over\sqrt{4\pi\nu B_i(t)}}
 \exp\left[-{(x_i-e^{s_it}z_i)^2\over4\nu B_i(t)}\right]. \tag{5.5}
\]

Since \(\prod_i e^{s_it}=1\), this kernel has unit integral in \(x\), and
the exact viscous solution is

\[
 \boxed{
 \omega_i(t,x)=e^{s_it}\int_{\mathbb R^3}
                         H_t(x,z)\omega_{0,i}(z)\,dz.}      \tag{5.6}
\]

Formula (5.6) shows the geometry directly.  The physical transverse
variances saturate:

\[
 2\nu B_1(t)\le\nu/\alpha,qquad
 2\nu B_2(t)\le\nu/\beta,                                  \tag{5.7}
\]

whereas both the center and standard deviation in the outgoing direction
grow like \(e^{\gamma t}\):

\[
 e^{\gamma t}z_3,qquad
 \sqrt{2\nu B_3(t)}
 =e^{\gamma t}\sqrt{{\nu(1-e^{-2\gamma t})\over\gamma}}.  \tag{5.8}
\]

Their ratio is time-uniform.  If \(z_3\ge R\), the probability of reaching
the material half-space \(x_3/e^{\gamma t}<R/2\) is bounded by

\[
 C\exp\left(-c{\gamma R^2\over\nu}\right).                 \tag{5.9}
\]

To obtain (1.4), split the Gaussian in the material outgoing coordinate
\(\zeta_3=e^{-\gamma t}x_3\).

* On \(\zeta_3\ge R/2\), apply the far-field kernel bound exactly as in
  (4.2).  This gives the first term of (1.4).
* On \(\zeta_3<R/2\), use (5.9).  The physical longitudinal density has the
  additional Jacobian \(e^{-\gamma t}\), which cancels the worst stretching
  multiplier \(e^{\gamma t}\).  Convolving the locally integrable
  \(|x|^{-2}\) kernel with the transverse Gaussians costs
  \(\ell_\nu^{-2}\).  For \(m>0\), move derivatives from the singular
  Biot--Savart kernel onto the Gaussian; each derivative costs at most one
  further \(\ell_\nu^{-1}\).  Small-time powers are absorbed by the stronger
  Gaussian \(e^{-cR^2/(\nu t)}\).

Assuming the high-P\'eclet separation \(R-r\ge4\ell_\nu\), this proves the
schematic estimate (1.4), with constants depending on the fixed strain
ratios and on \(r/R<1\).  Without that separation the two terms can be
replaced by a non-small bound at the larger of the geometric and viscous
core scales.  The argument is a direct heat-kernel estimate, but it has not
received an independent analytic cross-audit.  For a non-tight wake, the
same split is applied shell by shell; the leakage term is weighted by the
Gaussian moment
\(\int e^{-c a_3^2/\ell_\nu^2}|\omega_0(a)|\,da\), rather than by the
possibly infinite global \(L^1\) norm.

For an incoming wake, (1.4) is uniform for every gain length.  For a
time-distributed vorticity source \(f(s)\), Duhamel gives the same kernel
with \(t\) replaced by \(t-s\).  The inviscid far term is integrable in
\(t-s\), but the viscous leakage floor is multiplied by
\(\int_0^t\|f(s)\|_1ds\).  Thus a persistent source still needs either a
bounded total source ledger or a P\'eclet number large enough that the
exponential floor beats the gain-window length.

## 6. Nonlinear stability condition, not a theorem

For the full perturbation about \(Sx\), inviscid vorticity satisfies

\[
 \partial_t\omega+(Sx+u)\cdot\nabla\omega
 -(S+\nabla u)\omega=0.                                    \tag{6.1}
\]

Along a nonlinear characteristic,

\[
 {d\over dt}\big(e^{-\gamma t}X_3(t)\big)
 =e^{-\gamma t}u_3(t,X(t)).                                \tag{6.2}
\]

Hence the one-sided separation survives if, for example,

\[
 \int_0^\infty e^{-\gamma t}\|u_3(t)\|_{L^\infty}\,dt
 <R/2.                                                       \tag{6.3}
\]

Similarly, the linear stretching ledger survives with a bounded loss if

\[
 \int_0^\infty\|\nabla u(t)\|_{L^\infty}\,dt<\infty.       \tag{6.4}
\]

Under (6.3)--(6.4), a perturbed version of (4.3) follows with a constant
depending exponentially on the integral in (6.4).  These are useful
bootstrap hypotheses, not consequences of the current construction.  An
order-one work-carrying wake need not have small or time-integrable
self-strain.  It can translate, tilt out of the cone, or feed vorticity back
toward the core.  This is the first nonlinear obstruction.

## 7. Finite-energy material-displacement no-go

There is a sharper reason that the literal exponential export law (2.2)
cannot survive a finite-energy globalization on a fixed positive-volume
material set.

Let \(U\) be a smooth divergence-free velocity on \([0,T]\times\mathbb R^3\)
with flow map

\[
 \partial_tX(t,a)=U(t,X(t,a)),\qquad X(0,a)=a.               \tag{7.1}
\]

Let \(A\) be a bounded measurable set of volume \(V>0\).  Incompressibility
gives \(\det\nabla_aX=1\), and hence

\[
 \begin{aligned}
 {d\over dt}\int_AX_3(t,a)\,da
 &=\int_{X(t,A)}U_3(t,x)\,dx,\\
 \left|{d\over dt}\int_AX_3(t,a)\,da\right|
 &\le V^{1/2}\|U(t)\|_{L^2(\mathbb R^3)}.                  \tag{7.2}
 \end{aligned}
\]

After integration,

\[
 \left|\int_AX_3(T,a)\,da-\int_Aa_3\,da\right|
 \le V^{1/2}\int_0^T\|U(t)\|_2\,dt.                       \tag{7.3}
\]

There is no cancellation loophole from using both outgoing cones.  Apply
the same calculation to smooth approximations of \(|X_3|\) and pass to the
limit:

\[
 \int_A|X_3(T,a)|\,da
 \le\int_A|a_3|\,da
 +V^{1/2}\int_0^T\|U(t)\|_2\,dt.                           \tag{7.4}
\]

The radial version with \(|X|\) is identical.  Thus the average distance of
a finite-volume material set grows at most like the time-integrated kinetic
energy norm.

If every label in \(A\) satisfies \(X_3(T,a)\ge D(T)\), then (7.3) gives the
exact finite-window gate

\[
 VD(T)\le\int_Aa_3\,da
 +V^{1/2}\int_0^T\|U(t)\|_2\,dt.                            \tag{7.5}
\]

If only a fixed fraction \(\theta>0\) satisfies
\(|X_3(T,a)|\ge D(T)\), (7.4) gives the same inequality with
\(VD(T)\) replaced by \(\theta VD(T)\).  Consequently no fixed
positive-volume fraction can obey \(|X_3(t,a)|\ge Re^{\gamma t}\) for
arbitrarily long times when \(\|U(t)\|_2\) is uniformly bounded.

For a smooth unforced Navier--Stokes solution, the energy inequality gives

\[
 \|U(t)\|_2\le E_0:=\|U(0)\|_2,\qquad
 \int_0^T\|U(t)\|_2dt\le E_0T.                             \tag{7.6}
\]

The same conclusion holds for smooth Euler during its lifespan by energy
conservation.  If \(D(T)=Re^{\gamma T}\) and the initial material
barycenter is negligible relative to \(D(T)\), (7.5) requires

\[
 V^{1/2}\lesssim {E_0T\over R}e^{-\gamma T},\qquad
 V\lesssim {E_0^2T^2\over R^2}e^{-2\gamma T}.              \tag{7.7}
\]

This is also a quantitative finite-window obstruction.  Under a normalized
stage model in which \(E_0/R\) is uniformly controlled and
\(G_j\asymp c_Gj^2\), a material set exported with the literal affine law
would need

\[
 V_j\lesssim \operatorname{poly}(j)
              e^{-2\gamma c_Gj^2}.                          \tag{7.8}
\]

A merely geometric-in-\(j\) stage volume, even with polynomial carrier
factors, is much larger than the right side.  This cascade observation is
**conditional**: the physical normalization of \(E_{0,j}/R_j\), the actual
gain time, and the volume of the material wake set must all be checked
before applying (7.8).

The no-go has important limits.

* It concerns a material set transported by a smooth flow.  Viscous
  vorticity has no material support, so (7.5) does not directly constrain
  the Gaussian vorticity plume in Section 5.  It still constrains any
  proposed material carrier or passive-label export mechanism in a smooth
  Navier--Stokes velocity.
* With external forcing, replace \(E_0T\) by
  \(\int_0^T\|U(t)\|_2dt\).  Exponential material export requires that
  integral, or the reciprocal square root of the transported volume, to
  grow exponentially.
* A stage-dependent material volume may be as small as (7.7), and a finite
  interval alone gives no contradiction.  The result forbids a fixed
  positive-volume exponential channel for arbitrarily long time.
* It does not forbid linear or polynomial separation whose Biot--Savart
  gain beats a weaker vorticity-growth ledger.  Nor does it forbid
  prepositioned storage shells, multipole cancellation, or another
  Eulerian mechanism that makes the old wake invisible without transporting
  a positive-volume set exponentially far away.

This no-go explains exactly how the global affine model evades finite
energy: maintaining \(U_3\simeq\gamma X_3\) on a fixed-volume exponentially
outgoing packet forces an exponentially large \(L^2\) velocity budget, and
the full field \(Sx\) has infinite budget from the outset.

## 8. Why this does not yet repair C90

### 8.1 The affine background is not Clay-admissible

The kinetic energy of \(Sx\) on \(\mathbb R^3\) is infinite:

\[
 \int_{\mathbb R^3}|Sx|^2dx=\infty.                         \tag{8.1}
\]

It is an exact Euler and Navier--Stokes solution because \(\Delta(Sx)=0\)
and \((Sx\cdot\nabla)Sx=S^2x\) is a gradient, but it is only a local model.

### 8.2 A localized affine core does not keep exporting

Suppose the designed velocity agrees with \(Sx\) only in \(B_L\).  The
packet receives the exponential law (2.2) only until it exits \(B_L\).
Outside, a compactly supported or decaying background supplies no such
uniform strain.  The distance then need not continue to grow, so (4.5) is
replaced at best by a fixed \(L^{-2-m}\) far-field bound.  Over a growing
gain window this is bounded but generally not contracting under the desired
stage return.

Cutting off \(Sx\) also violates the equation in the collar.  Correcting
that defect creates vorticity, work, and a nonlocal pressure/velocity tail.
No finite-energy unforced flow is currently known that provides the same
one-way cone for all relevant times while retaining the affine capture in
the core.

### 8.3 Pressure is absent only after the hard geometry is granted

Vorticity removes pressure from (2.1), and whole-space Biot--Savart already
accounts for the instantaneous velocity tail of a separated packet.  This
is why (4.3) is legitimate despite the lack of finite propagation.
However, the pressure that enforces a **localized** affine core is part of
the missing global construction.  It can change the exterior transport and
cannot be prescribed independently of the cutoff correction.

### 8.4 Viscosity gives infinite propagation

Equation (5.6) has no compact support at positive time.  The leakage is
exponentially small, not zero.  A single incoming packet has the uniform
floor (1.4), but repeated or continuously generated wakes require a summable
source ledger.  Derivative constants also contain powers of
\(\ell_\nu^{-1}\), so the all-order Gevrey estimate must check that the
P\'eclet exponential dominates those powers through the growing order
\(M_j\).

### 8.5 The torus has no export direction

The vector field \(Sx\) is not periodic.  More fundamentally, a compact
torus has no spatial infinity: a transported wake wraps around and remains
visible to the lowest periodic modes.  Thus this mechanism cannot repair
the periodic heat no-go in C90.  It requires the \(\mathbb R^3\) Clay
alternative and a genuinely noncompact retained wake.

## 9. Sharpened open target

The calculation suggests replacing a generic periodic endpoint estimate by
the following geometric theorem.

> **Finite-energy cone-export theorem (open).** Construct a smooth,
> divergence-free, finite-energy unforced Navier--Stokes transition on
> \(\mathbb R^3\) with an active affine-capture core and an outgoing wake
> channel such that: (i) the wake acquires a distance or cancellation gain
> \(D(t)\) whose Biot--Savart decay beats its vorticity stretching without
> demanding forbidden fixed-volume exponential material displacement;
> (ii) nonlinear self-advection preserves the channel; (iii) the localized
> pressure/collar correction belongs to the same channel; and (iv) the
> resulting wake-to-core endpoint map is uniform through the full gain
> window in the required Gevrey orders.

The global-affine model proves that the local exponent ledger is favorable:
distance supplies two Biot--Savart powers while stretching costs only one.
Section 7 proves that its literal fixed-volume exponential transport cannot
be realized by a uniformly finite-energy flow for arbitrarily long times.
What remains open is a weaker finite-energy decoupling mechanism: slower
separation with still slower vorticity growth, exponentially small
stage-dependent transported volume with a verified physical normalization,
or an Eulerian cancellation/storage construction that does not rely on
material export.
