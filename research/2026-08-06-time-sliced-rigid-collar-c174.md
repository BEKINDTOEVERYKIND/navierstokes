# C174: time-sliced rigid collars have a rechart-or-rewrite dichotomy

**Date:** 2026-08-06
**Status:** exact prescribed-ramp, scale, pressure-moment, recharting, and
central-return ledgers; unforced realization, finite-tube residence, and
the integrated wake kernel remain conditional
**Checker:**
[checks/time_sliced_rigid_collar_c174.py](../checks/time_sliced_rigid_collar_c174.py)

## 0. Claim boundary

C171 shows that a center-following rigid curl collar has the favorable
backward parent-cross scale

\[
                         q^{-1}\log h.              \tag{0.1}
\]

Unlike a full material Piola tube, it has no accumulated \(F\)-jet
polynomial.  This note tests whether one can exploit that fact through the
whole C142/C159 gain window by dividing it into

\[
             J=O(\log h)                            \tag{0.2}
\]

subintervals of \(O(1)\) strain action and replacing/reorienting the compact
collar completion at every boundary.

The collar ledger is favorable as a residual budget. If adjacent
completions agree with the same
prescribed affine converter on a common child core, a smooth switch is
exactly divergence free and its time-integrated \(L^2\) cost is the
difference of two collar profiles.  One switch has backward scale
\(O(\Lambda q^{-1})\), where

\[
                     \Lambda={\lambda\ell\over a},  \tag{0.3}
\]

and all \(J\) switches cost

\[
             O(\Lambda q^{-1}\log h)
             =O(\Lambda n^{-8}\log n)=o(n^{-6})    \tag{0.4}
\]

whenever \(\Lambda\log n=o(n^2)\).  Parent-cross and self terms have the
same or smaller power; viscosity is absorbed by the factorial Reynolds
number.  Thus an \(O(1)\) collar-profile change does **not** recreate
C143's fatal backward \(O(\log h)\) scale.  Its amplitude is
\(\lambda r\), not the parent/child amplitude.

The ramp is prescribed and creates the residual displayed below. It is
not, by itself, an unforced Navier--Stokes operation. Turning this
residual budget into an exact unforced trajectory still requires the
stage correction/closure theorem.

The child ledger has the opposite conclusion.  A coordinate rechart of the
same physical child introduces no equation residual and preserves the exact
product cocycle, but it does not reset physical deformation.  C154's
central return remains

\[
                  F^m=I+mR,\qquad R^2=0,           \tag{0.5}
\]

and the covector bandwidth still shears by \(O(J)\).  If instead one
physically rewrites the child to an undeformed reference profile, an
\(O(1)\) relative mismatch at any slice has the same individual backward
terminal size, independently of when it is made. With active seed \(b\), a
triangle-inequality proof charges \(O(b)\) per rewrite and \(O(bJ)\) over
the window, whereas BAFL permits \(O(b^3)\). Such a robust sum-of-norms
proof therefore needs average relative rewrite mismatch

\[
                       O(b^2/J),                  \tag{0.6}
\]

not \(O(1)\).

This is the exact **rechart-or-rewrite obstruction (RRR)**:

> A pure rechart preserves the intended C159 cocycle but leaves the
> accumulated child shear and finite-tube residence problem unchanged.
> A physical reset removes that deformation only by rewriting the child.
> Each order-one rewrite has an individual \(b\)-scale weighted response,
> so separated sum-of-norms control loses the two extra powers required by
> BAFL. A specially phased signed cancellation is not ruled out.

At the exact central first-jet level, the remaining deformation is only
linear in \(J\), so an initial core smaller by \(O(J^{-1})\) would suffice.
Whether a finite-radius packet stays in that core requires a uniform
second-flow-jet bound; C154 proves only the central first jet.  Time slicing
therefore removes the selector's Piola-chart growth but does not close
MCKC(i), child residence, coherence, the pressure/wake evolution, LCE,
BAFL, or the full unforced stage.

## 1. A sliced rigid-curl construction

Use the C143 scales

\[
 r={\ell\over q},\qquad q=n^8,\qquad
 h=q^{3/2}=n^{12},\qquad b=n^{-2}.                 \tag{1.1}
\]

Let \(c(t)\) follow the central parent trajectory,

\[
                         c'(t)=V(t,c(t)).          \tag{1.2}
\]

Split the gain interval into \(I_m=[t_m,t_{m+1}]\),
\(0\le m<J\), so that

\[
 \lambda |I_m|\le C_0,\qquad J\le C_1(1+\log h).  \tag{1.3}
\]

For each slice choose a trace-free affine core \(S_m(t)y\),
\(\|S_m\|\le C\lambda\), and a fixed-scale cutoff/tube
\(\chi_m(t,y/r)\). In addition to uniform spatial seminorms, assume the
profile path has \(O(1)\) dimensionless total variation on a slice:

\[
 \int_{I_m}\left(\|S_m'(t)\|
 +\lambda\|\partial_t\chi_m(t,\cdot)\|_{C^2}\right)dt
 \le C\lambda.                                    \tag{1.3a}
\]

Put

\[
 {\cal A}_m(t,y)=
 \chi_m(t,y/r){(S_m(t)y)\times y\over3},
 \qquad
 U_m(t,x)=\nabla_x\times{\cal A}_m(t,x-c(t)).      \tag{1.4}
\]

Every \(U_m\) is exactly divergence free and compactly supported in a ball
or tube of diameter \(O(r)\).  On its inner core,

\[
                         U_m(t,c+y)=S_m(t)y.       \tag{1.5}
\]

Uniform cutoff and orientation bounds give

\[
 \begin{aligned}
 \|U_m\|_2&\le C\lambda r^{5/2},&
 \|\nabla U_m\|_2&\le C\lambda r^{3/2},\\
 \|\Delta U_m\|_2&\le C\lambda r^{1/2},&
 \|U_m\|_\infty&\le C\lambda r.                   \tag{1.6}
 \end{aligned}
\]

There is no long-time \(F\)-jet in (1.6): each field is built directly in
the current rigid frame.  Rotation of an anisotropic tube is harmless in
these norms as long as its aspect ratio is uniformly bounded.  A spherical
cutoff is orientation independent.

Without (1.3a), an arbitrarily oscillatory prescribed orientation or
cutoff would have an uncontrolled time-derivative residual; the spatial
bounds (1.6) alone do not control that term.

## 2. Exact ramps and the collar switching cost

Near \(t_{m+1}\), let \(\theta_m\) increase smoothly from zero to one and
set

\[
 U=(1-\theta_m)U_m+\theta_m U_{m+1}.               \tag{2.1}
\]

Linearity of divergence and curl gives

\[
                         \nabla\cdot U=0            \tag{2.2}
\]

exactly.  The switching part of the time derivative is

\[
 G_m^{\rm sw}=\theta_m'(U_{m+1}-U_m).              \tag{2.3}
\]

For a monotone ramp,

\[
 \boxed{
 \int\|G_m^{\rm sw}(t)\|_2\,dt
 \le\sup_{\rm ramp}\|U_{m+1}-U_m\|_2
 \le C\lambda r^{5/2}.}                           \tag{2.4}
\]

The estimate is independent of ramp duration.  Making a ramp faster does
not reduce its integrated cost; making it slower merely moves the same
variation through time.

The same scale covers smooth within-slice profile motion satisfying
(1.3a). It is a bound for the residual of a prescribed approximate
trajectory. No external forcing is admitted in the target problem, so
this source must ultimately be canceled by the unforced stage
construction; C174 only proves that its norm fits the scalar budget.

The useful case is when the two completions realize the same prescribed
affine field on a common child core during the overlap.  Then

\[
             U_{m+1}-U_m=0
             \quad\hbox{on the common core},       \tag{2.5}
\]

and (2.3) is a collar source.  If the prescribed matrix itself rotates,
choose a smooth symmetric trace-free interpolation \(S(t)\).  The global
affine field \(S(t)y\) has

\[
 \partial_t(Sy)+(Sy\cdot\nabla)(Sy)
       =\{S'+S^2\}y
       =\nabla\left({1\over2}y\cdot(S'+S^2)y\right), \tag{2.6}
\]

so its core ramp is pressure-exact.  Only the compact completion produces
the collar residual.  A ramp with \(|S'|\lesssim\lambda^2\) over
time \(O(\lambda^{-1})\) again has the integrated scale (2.4).
This sentence applies directly to C142's symmetric affine selector.  On
the nonsymmetric C159 orbit one instead uses the actual continuous parent
jet; its curvature/viscous Taylor remainder is not declared pressure-exact
and remains in the parent-residual ledger.

Normalize (2.4) by \(a\ell^{3/2}\), then apply the worst full backward
factor \(h\).  One switch costs

\[
 h{\lambda r^{5/2}\over a\ell^{3/2}}
 =\Lambda hq^{-5/2}
 =\boxed{\Lambda q^{-1}}.                         \tag{2.7}
\]

Here \(\lambda r^{5/2}\) is the physical \(L^2\) size of the selector
collar difference, \(a\ell^{3/2}\) is the parent normalization, and the
C159/C142 gain \(h\) is applied exactly once.  This is not a reset of the
physical child, whose different normalization is treated in Section 6.

All switches therefore obey

\[
 \varepsilon_{\rm sw}^{\rm back}
 \le C\Lambda Jq^{-1}
 \le C\Lambda q^{-1}(1+\log h).                   \tag{2.8}
\]

On (1.1), its logarithmic part is

\[
             12C\Lambda n^{-8}\log n=o(n^{-6})
 \quad\Longleftrightarrow\quad
             \Lambda\log n=o(n^2).                \tag{2.9}
\]

The bound uses the full backward factor at every boundary.  Weighting each
switch only by its remaining gain can improve it, but is unnecessary for
the exponent conclusion.

## 3. Parent-cross, self, and viscous ledgers

For a frozen rigid profile, C171's exact identity is the last two terms
below. For the time-dependent profiles in Section 1, the complete identity
is

\[
 D_t^VU+(U\cdot\nabla)V
 =\partial_t^{\rm prof}U
  +[V(t,x)-V(t,c)]\cdot\nabla U+(\nabla V)U.       \tag{3.1}
\]

If \(\|\nabla V\|_\infty\le C a/\ell\), (1.6) gives

\[
 \|[V-V(c)]\cdot\nabla U+(\nabla V)U\|_2
 \le C{a\lambda\over\ell}r^{5/2}.                 \tag{3.2}
\]

On one slice, \(\lambda|I_m|\le C_0\), hence

\[
 {1\over a\ell^{3/2}}
 \int_{I_m}\|G^{\rm par}\|_2dt
 \le Cq^{-5/2}.                                   \tag{3.3}
\]

The self and viscous terms satisfy

\[
 \begin{aligned}
 {1\over a\ell^{3/2}}
 \int_{I_m}\|(U\cdot\nabla)U\|_2dt
 &\le C\Lambda q^{-5/2},\\
 {1\over a\ell^{3/2}}
 \int_{I_m}\|\nu\Delta U\|_2dt
 &\le C{\rm Re}^{-1}q^{-1/2},
 \qquad {\rm Re}={a\ell\over\nu}.                 \tag{3.4}
 \end{aligned}
\]

After summing \(J\) slices and applying \(h=q^{3/2}\),

\[
 \boxed{
 \begin{aligned}
 \varepsilon_{\rm par}^{\rm back}
     &\le CJq^{-1},\\
 \varepsilon_{\rm self}^{\rm back}
     &\le C\Lambda Jq^{-1},\\
 \varepsilon_{\nu}^{\rm back}
     &\le C{\rm Re}^{-1}Jq.
 \end{aligned}}                                   \tag{3.5}
\]

The first two have the same favorable power as (2.8).  On the factorial
schedule \({\rm Re}^{-1}=\nu(j!)^{-2}\), the last term is
\(O(\nu(j!)^{-2}n^8\log n)\) and is \(o(n^{-6})\).  These are raw
backward \(L^2\) ledgers.  They do not prove contraction of the localized
active evolution family, nor do they close the retained wake.

The first term in (3.1) is charged in Section 2 under (1.3a); without that
hypothesis it would be uncontrolled.

## 4. Instantaneous pressure-tail ledger

At each fixed time, let \(G_m(t,\cdot)\) denote any compact residual above,
including a switch. Every such instantaneous residual has zero spatial
mean:

* compact curls have zero mean, so the time derivative in (2.3) does too;
* each advective term integrates to zero for divergence-free fields; and
* the Laplacian integrates to zero.

On \(\mathbb R^3\), outside the instantaneous support, the Leray kernel is
homogeneous of degree \(-3\). Subtracting its value at the instantaneous
collar center and using \(\int G_m(t)=0\) gives, for
\(d(t)=|x(t)-c(t)|\ge4r\),

\[
 |{\mathbb P}G_m(t,x(t))|
 \le {C\over d(t)^4}
        \int |y-c(t)|\,|G_m(t,y)|\,dy.             \tag{4.1}
\]

Integrating (4.1) in time is legitimate along any observation trajectory
which remains separated by \(d(t)\asymp\ell\). It avoids a false fixed
support assertion: the time-integrated Eulerian source can occupy the
whole tube swept out by \(c(t)\). The integrated first-moment ledgers per
slice are

\[
 \begin{array}{c|c|c}
 \text{source}&
 \displaystyle\int_{I_m}\!\int |y-c(t)|\,|G_m(t,y)|\,dy\,dt&
 \displaystyle {1\over a}\int_{I_m}|{\mathbb P}G_m(t,x(t))|dt\\ \hline
 \text{parent cross}&Ca r^5/\ell&Cq^{-5}\\
 \text{self}&C\lambda r^5&C\Lambda q^{-5}\\
 \text{viscosity}&C\nu r^3&C{\rm Re}^{-1}q^{-3}.
 \end{array}                                      \tag{4.2}
\]

The pure switching source (2.3) is itself divergence free at each time,
so Leray leaves it unchanged and it has **zero** exterior pressure tail
outside its instantaneous support. Its compact \(L^2\) cost is already
charged in (2.7)--(2.8). Parent-cross and self-interactions generated
during a ramp are covered by the corresponding rows of (4.2).

After \(J\) slices and the backward factor \(h\), the separated pointwise
scales are

\[
 \boxed{
 CJq^{-7/2},\qquad
 C\Lambda Jq^{-7/2},\qquad
 C{\rm Re}^{-1}Jq^{-3/2}.}                        \tag{4.3}
\]

These powers are far below (3.5). The degree-\(-4\) image sum is
absolutely summable over a three-dimensional periodic lattice only when
zero mean holds for each translated source (or an equivalent cellwise
multipole cancellation) and the observation path remains uniformly
separated from every copy on the stated scale. If cancellation holds only
after summing several sources, their individual degree-\(-3\) tails cannot
be discarded before the periodic sum. The estimate also does not cover a
growing number of nearby collars without charging that multiplicity.

Equations (4.1)--(4.3) are an **instantaneous,
trajectory-separated pressure estimate integrated in time**, not a BAFL
pressure theorem. During the actual gain, pressure is inserted into a
non-normal active propagator; supports can move, heat spreads them, and
the retained wake is not just a sum of separated tails. Those dynamical
obligations remain open.

## 5. Exact preservation of the child cocycle under recharting

Let \(P_m\) be the exact physical child propagator on slice \(I_m\), and
let \(Q_m\) be any invertible reorientation used to express the same
physical child in slice coordinates.  The coordinate propagator is

\[
                     \widetilde P_m=Q_{m+1}^{-1}P_mQ_m.      \tag{5.1}
\]

The product telescopes:

\[
 \boxed{
 \widetilde P_{J-1}\cdots\widetilde P_0
 =Q_J^{-1}(P_{J-1}\cdots P_0)Q_0.}               \tag{5.2}
\]

Thus exact reorientation/recharting does not alter the **physical**
propagator and produces no PDE source because the physical field has not
changed. If \(Q_J=Q_0\), (5.2) is a similarity and preserves Floquet
multipliers. If the endpoint charts are uniformly conditioned (in
particular, orthogonal reorientations), physical and coordinate norm gains
are uniformly comparable.

For arbitrary \(Q_0,Q_J\), however, the two-sided coordinate matrix
\(Q_J^{-1}PQ_0\) need not have the same eigenvalues, singular values, or
phase convention as \(P\). Equation (5.2) alone therefore does not prove a
uniform coordinate gain/coherence chart; that endpoint conditioning is an
additional requirement.

But (5.2) is only a change of coordinates.  It cannot change a physical
support, physical gradient, or Fourier bandwidth.  In particular C154's
exact central return derivative is

\[
 F_*=I+R,\qquad R=u\otimes g_0,\qquad
 g_0\cdot u=0,\qquad R^2=0.                       \tag{5.3}
\]

Therefore

\[
 \boxed{
 F_*^m=I+mR,\qquad
 F_*^{-Tm}=I-mR^T.}                               \tag{5.4}
\]

The central physical first jet and generic covector width retain their
linear \(m\)-growth after every coordinate reset.  The C154 requirement of
a correlated packet or an initial shear-direction width
\(O(q/J)\) is unchanged.

At the purely linearized residence level, if

\[
 D_J=\max_{0\le m\le J}\|I+mR\|,
 \qquad \rho_0\le {r\over2D_J},                   \tag{5.5}
\]

then the initial ball \(B_{\rho_0}\) remains in the rigid core \(B_r\) at
all return times. Since \(D_J=O(1+J)\), this sufficient isotropic choice
costs only a logarithmic shrink in radius. It is not free: the physical
core volume is smaller by \(O(J^{-3})\), its reciprocal volume
normalization is \(O(J^3)\), fixed-energy \(L^2\) amplitude normalization
would be \(O(J^{3/2})\), and resolving the smaller radius costs an
additional \(O(J)\) maximum Fourier width in each isotropically narrowed
direction.

This is not a lower bound. A shear-aligned anisotropic tube may narrow
only the active shear direction and pay less volume. Compatibility of
either choice with the C159/C170 \(q^2\) lattice capacity, coherence, and
energy normalization is not proved.

For a finite tube one also needs a Taylor remainder.  Write schematically

\[
 \Phi^m(c+y)=c+F_*^my+{\cal R}_m(y),\qquad
 |{\cal R}_m(y)|\le H_J{|y|^2\over\ell}.           \tag{5.6}
\]

The choice \(\rho_0=r/(4D_J)\) remains inside \(B_r\) if, for example,

\[
                         H_J\le4D_J^2q.            \tag{5.7}
\]

Every fixed polylogarithmic \(H_J\) satisfies (5.7), but C154 does not
prove such a second-flow-jet bound on the finite tube.  Equation (5.7) is
only a sufficient finite-radius bound at the discrete return times.
Residence between returns additionally requires a uniform within-period
first/second-flow-jet bound.  These together form the child-residence gate
left by the pure rechart.

## 6. A physical child reset is too expensive

Suppose instead that at slice \(m\) one replaces the actual deformed child
by an undeformed reference child.  Let its current gain be \(g_m\), with
\(1\le g_m\le h\), and let \(\delta_m\) be the relative mismatch in the
prescribed growing child coordinate (with its scheduled \(L^2\)
normalization).  At the scheduled child volume \(r^3\), its normalized
current ledger size is

\[
              \delta_m\,b\,g_m\,q^{-3/2}.          \tag{6.1}
\]

The scalar remaining backward weight is \(h/g_m\).  Hence the
terminal-weighted growing-coordinate ledger is exactly

\[
 {h\over g_m}
   \delta_m b g_m q^{-3/2}
 =\boxed{\delta_m b},                              \tag{6.2}
\]

independent of the slice time. If \(\delta_m\ge0\) is the norm of the
mismatch in the normalized growing coordinate, the sum of the individual
scalar weighted-response norms is exactly

\[
              \sum_m\|\text{weighted rewrite}_m\|
                   = b\sum_m\delta_m.               \tag{6.3}
\]

The norm of the signed total Duhamel response is at most (6.3) and can be
smaller through phase or vector cancellation; (6.3) is not a lower bound
for that signed total. The active pre-chart
BAFL allowance is \(O(b^3)\).  Therefore a robust separated-source proof
using the triangle inequality requires

\[
                      \sum_m\delta_m=O(b^2).       \tag{6.4}
\]

For \(J\asymp\log h\), a separated sum-of-norms proof needs the average
mismatch in (0.6). If that proof is used and \(\delta_m\asymp1\), then

\[
 \varepsilon_{\rm rewrite}^{\rm back}
 \asymp bJ=O(n^{-2}\log n),
 \qquad
 {\varepsilon_{\rm rewrite}^{\rm back}\over b^3}
 \asymp {J\over b^2}=O(n^4\log n).                \tag{6.5}
\]

Relative to the terminal child itself, the accumulated loss is
\(O(J)=O(\log h)\).  Slowing a ramp cannot repair (6.2), because its
integrated time derivative is the profile variation.

Equation (6.3) is the sum-of-norms ledger for a robust separated-source
argument, not a universal signed lower bound on a specially engineered
sequence of rewrites.  Exact telescoping could cancel profile jumps, but
then the final physical deformation has not been reset; a cancellation
that also preserves an undeformed final child is part of the missing
nonperturbative stage construction.

## 7. Exact surviving obstruction

Time slicing gives a useful positive narrowing:

* rigid collar completions can be changed \(O(\log h)\) times without an
  \(F\)-jet polynomial or a fatal exponent loss;
* their switching, parent-cross, self, viscosity, and
  trajectory-separated pressure powers all fit the factorial schedule;
  and
* exact coordinate reorientation telescopes to the same physical C159
  child cocycle, subject to a separate endpoint chart-conditioning gate.

It does not physically reset the child.  The single remaining gate for this
variant is

> **RRR, the rechart-or-rewrite obstruction.** Either keep the physical
> child and prove finite-tube residence/coherence under the accumulated
> C154 shear, including (5.7), or rewrite it with
> \(\sum_m\delta_m=O(b^2)\) in a separated sum-of-norms proof, or else
> prove the required cancellation directly in the signed Duhamel response,
> while preserving the target cocycle and retained-wake bounds.

The first branch is a narrowed form of MCKC(i); the second requires a
nonperturbative cancellation far stronger than ordinary slice ramps.
Neither branch controls the wake-to-active kernel in MCKC(ii).  The full
unforced stage, LCE, BAFL, and any blow-up conclusion remain conditional.
