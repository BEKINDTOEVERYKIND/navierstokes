# C176: the exact \(A_2\) normal form closes background tube residence

**Date:** 2026-08-13
**Status:** exact background-flow normal form, jet, lattice-capacity,
nilpotent-correlation, enlarged-collar, and fixed-cocycle rewrite ledgers;
finite-frequency localization, nonlinear wake closure, C125, BAFL, and the
unforced stage remain open
**Checker:**
[checks/a2_correlated_tube_residence_c176.py](../checks/a2_correlated_tube_residence_c176.py)

## 0. Claim boundary

C174 left the rechart-or-rewrite obstruction (RRR).  Its first branch asked
for a second-flow-jet and between-return residence theorem for a packet
carried through \(J=O(\log h)\) returns of the C152--C159 zero-drift orbit.
For the **exact heat-decaying \(A_2\) background**, that part is not an
unknown Gronwall problem.  The regular level foliation gives an exact local
action--angle--axial normal form

\[
 E'=0,\qquad \theta'=\omega(E),\qquad z'=-\sqrt2E.       \tag{0.1}
\]

Consequently the full between-return flow has first jet \(O(J)\), second
jet \(O(J^2)\), and, more generally, every fixed-order jet grows only
polynomially in \(J\).  An ordinary child ball of radius
\(r=\ell/q\) remains inside a ball of radius \(O(Jr)\).  For the
phase-space-correlated construction, one may instead choose from the start
an envelope with two widths \(O(r)\) and one width \(O(Jr)\).  Its volume
is \(O(Jr^3)\), whereas the incompressible image of an ordinary \(r\)-ball
still has volume \(O(r^3)\); the larger tube is a convenient cover, not a
claim of physical volume creation.  This proves the missing exact-background
residence estimate at every intermediate time, not only at return times.

The polylogarithmic enlargement does not spoil C174's scalar collar
budget.  A deliberately crude spherical collar of radius \(R=CJr\) has
worst switching/parent/self backward ledger

\[
                       O\!\left(q^{-1}J^{7/2}\right),       \tag{0.2}
\]

up to the displayed strain ratio.  With

\[
 q=n^8,\qquad b=n^{-2},\qquad h=q^{3/2}=n^{12},
 \qquad J=O(\log h),                                      \tag{0.3}
\]

this is \(n^{-8}\operatorname{polylog}n=o(n^{-6})\), and it
remains \(o(n^{-4})\) after one \(n^2\) chart charge.  The factorial
Reynolds factor absorbs the enlarged viscous ledger.

There is also enough torus arithmetic for the **carrier centres**.  A
fixed-aperture box about the C159 ray, of widths

\[
                  q\ \times\ q\ \times\ {q\over J},       \tag{0.4}
\]

where the short frequency direction is normal to the periodic-covector
plane, contains \(\Omega(q^3/J)\) integer points.  It therefore contains
\(q^2\) reality-paired carrier centres for all large stages.  The exact
C154 shear keeps the whole central-fiber box at \(O(q)\) width for all
\(m\le J\).
This does not contradict C170: C170's \(O(q)\) count concerned aperture
\(O(q^{-1})\) and normal tolerance \(O(J^{-1})\); (0.4) uses a fixed
two-directional aperture and normal tolerance \(O(q/J)\).

Strictness of the C159 cooperative inequalities gives, by compact
continuity, a qualitative fixed position--projective-covector tube in
which the **principal** Kelvin cone still expands uniformly.  What is not
supplied here is the
finite-frequency localized parametrix, physical endpoint coherence, the
unforced full-polarization terminal converter, or the C125/BAFL response.
In particular, compact localization convolves the carrier set with
sidebands; (0.4) is a carrier-centre and phase-space compatibility theorem,
not a compact band-limited Navier--Stokes solution.

Finally, signed physical rewrites have no free telescoping mechanism when
the intervening linear cocycle is held fixed.  Their signed terminal
Duhamel sum is exactly the difference between the rewritten endpoint and
the unrewritten cocycle endpoint.  If that sum vanishes, the final physical
child state equals the unrewritten endpoint; a coordinate rechart cannot
then reset its physical C154 deformation.  A final mismatch of order one
in the specified normalized growing coordinate therefore has backward
size \(\Theta(b)\), not \(O(b^3)\).  This is a fixed-cocycle statement; a
nonlinear reservoir which changes the cocycle is not excluded, but then
its dynamics and wake are precisely part of the remaining unforced BAFL
theorem.

Thus C176 closes the **background residence and lattice-capacity part** of
RRR on the existing geometry.  It does not close MCKC(ii), LCE, C125,
BAFL, the one-cell map, or the Millennium problem.

## 1. Exact action--angle--axial normal form

Retain C152's dimensionless field

\[
 U=N\times\nabla f-\sqrt2 fN,
 \qquad
 f(a,b)=\cos a+\cos b+{4\over5}\cos(a+b),             \tag{1.1}
\]

with \(a=r_1\cdot x\), \(b=r_2\cdot x\), and
\(N=(1,1,1)\).  Put

\[
                         z={N\cdot x\over3}.              \tag{1.2}
\]

The planar phase flow is Hamiltonian with first integral \(E=f\), while

\[
                         \dot z=-\sqrt2E.                 \tag{1.3}
\]

C152 proves that \(E=0\) is a compact regular closed level.  Hence there
is \(E_*>0\) such that \(|E|<E_*\) is a regular annulus.  Choose a smooth
transversal and let \(\theta\in\mathbb R/2\pi\mathbb Z\) be normalized
time around each level.  If \(T(E)\) is its period, then

\[
                         \omega(E)={2\pi\over T(E)}        \tag{1.4}
\]

is smooth, and the flow is exactly (0.1).  On the lifted angle coordinate,

\[
 \boxed{
 \Psi_\tau(E,\theta,z)
 =\bigl(E,\theta+\omega(E)\tau,z-\sqrt2E\tau\bigr).}      \tag{1.5}
\]

No linearization has been used in (1.5).

For the unforced heat-decaying Beltrami pump

\[
 V(t,x)=A(t)U(Kx),\qquad A(t)=A_0e^{-2\nu K^2t},           \tag{1.6}
\]

the same formula holds after the inertial-time change

\[
                    \tau(t)=K\int_0^t A(s)\,ds.           \tag{1.7}
\]

The sign or rate of the heat clock changes how much physical time is used,
but not the spatial orbit or any of the following \(\tau\)-flow identities
on an interval where the clock is monotone.

### 1.1 Polynomial flow jets

In the normal coordinates,

\[
 D\Psi_\tau=
 \begin{pmatrix}
 1&0&0\\
 \tau\omega'(E)&1&0\\
 -\sqrt2\tau&0&1
 \end{pmatrix},
 \qquad
 \partial_E^2(\Psi_\tau)_\theta=\tau\omega''(E),          \tag{1.8}
\]

and all other second derivatives vanish.  Let \(\Xi\) be the smooth
coordinate map from \((E,\theta,z)\) back to physical space.  On a compact
subannulus, \(\Xi,\Xi^{-1}\), and their fixed-order derivatives are
bounded.  The chain rule applied to

\[
                  \Phi_\tau=\Xi\Psi_\tau\Xi^{-1}          \tag{1.9}
\]

therefore gives

\[
 \boxed{
 \begin{aligned}
  \|D\Phi_\tau\|+\|D\Phi_\tau^{-1}\|&\le C(1+|\tau|),\\
  \|D^2\Phi_\tau\|+\|D^2\Phi_\tau^{-1}\|
      &\le C(1+|\tau|)^2.
 \end{aligned}}                                           \tag{1.10}
\]

More generally, every fixed \(s\)-jet is \(O_s((1+|\tau|)^s)\).
The quadratic power in the second line comes from the two copies of
\(D\Psi_\tau\) entering \(D^2\Xi[D\Psi_\tau,D\Psi_\tau]\); it is a
uniform upper bound, not a claim that every entry has quadratic growth.

At \(\tau=mT(0)\), (1.10) is consistent with C152's exact derivative

\[
 D\Phi_{mT(0)}=I+mR,
 \qquad R=u\otimes g_0,qquad R^2=0.                       \tag{1.11}
\]

Unlike C174's schematic hypothesis \(H_J\), (1.10) also controls every
time between returns.

## 2. Exact finite-tube residence

Let \(\ell=K^{-1}\), \(r=\ell/q\), and consider
\(0\le\tau\le C_0J\).  In the \((\theta,z)\)-plane the first-order slip
caused by an energy displacement is along

\[
                         v_0=(\omega'(0),-\sqrt2).          \tag{2.1}
\]

Take an initial curvilinear parallelepiped with

* energy width \(O(q^{-1})\);
* width \(O(q^{-1})\) in the direction perpendicular to \(v_0\); and
* width \(O(Jq^{-1})\) in the \(v_0\) direction.

Formula (1.5) gives

\[
 (\omega(E)-\omega(0),-\sqrt2E)\tau
   =E\tau v_0+O(E^2\tau)(1,0).                            \tag{2.2}
\]

Since \(|E|=O(q^{-1})\), \(J/q\to0\), and \(J/q^2\ll q^{-1}\), the
image of this deliberately correlated envelope remains in a comparable
parallelepiped for **every**
\(0\le\tau\le C_0J\).  Returning to physical units, it has two widths
\(O(r)\), one width \(O(Jr)\), and volume

\[
                              O(Jr^3).                     \tag{2.3}
\]

Relative to an \(r^3\) child core, the chosen bounding envelope is a factor
\(J\) larger.  If one fills that full envelope roughly uniformly at fixed
\(L^2\) energy, its point amplitude loses \(J^{1/2}\); conversely, keeping
that filled-envelope point amplitude would cost a factor \(J\) in energy.
These are normalization options, not consequences for the
volume-preserving image of an ordinary ball.  They are polylogarithmic
rather than new powers of \(n\), but are not declared free: the terminal
coherence/conversion theorem must carry the chosen normalization.  The
long spatial width \(Jr\) is reciprocal to the short carrier-centre
frequency width \(q/J\) used in Section 3.

The same conclusion follows more crudely from (1.10): an initial ball
\(B_r\) has linear image diameter \(O(Jr)\), while its second-order
remainder is

\[
                 O\!\left({J^2r^2\over\ell}\right)
       =O\!\left(Jr{J\over q}\right)=o(Jr).               \tag{2.4}
\]

Thus a ball of radius \(R=CJr\) contains the full tube.  In C174's
notation one may take

\[
             D_J\le CJ,\qquad H_J\le CJ^2,                \tag{2.5}
\]

so its sufficient inequality \(H_J\le4D_J^2q\) holds for all large
\(q\).  This proves the missing second-jet and between-return residence
gate for the exact background pump.

It does **not** prove stability of this tube after adding the localized
collar, the growing child, its nonlinear wake, or the unforced correction.
Those perturbations must be controlled in the C125/BAFL norm rather than
silently inserted into (1.5).

## 3. The correlated \(q^2\) carrier packet is arithmetically available

Let

\[
 e_3={u\over|u|},\qquad e_1={g_0\over|g_0|},
 \qquad e_2=e_3\times e_1.                                \tag{3.1}
\]

Because \(u\cdot g_0=0\), this is an orthonormal frame.  Let
\(\bar k\in u^\perp\) be the C159 periodic ray, normalized to fixed
nonzero length.  For a sufficiently small fixed \(\delta>0\), define the
oriented box

\[
 \begin{aligned}
 {\cal B}_{q,J}=\{k\in\mathbb R^3:\;&
 |e_1\cdot(k-q\bar k)|\le\delta q,\\
 &|e_2\cdot(k-q\bar k)|\le\delta q,\\
 &|e_3\cdot(k-q\bar k)|\le\delta q/J\}.
 \end{aligned}                                            \tag{3.2}
\]

### 3.1 A rotation-independent lattice lower bound

Nearest-integer unit cubes partition \(\mathbb R^3\).  Every point in the
Euclidean \(\sqrt3/2\)-erosion of (3.2) has its nearest lattice point in
(3.2).  The eroded box has volume comparable to \(q^3/J\) once
\(q/J\to\infty\).  Consequently

\[
 \boxed{
   \#\bigl(\mathbb Z^3\cap{\cal B}_{q,J}\bigr)
       \ge c_\delta{q^3\over J}}                           \tag{3.3}
\]

for all sufficiently large \(q/J\), independently of the unresolved
rationality of C170's exact plane.  For every integer
\(k\in{\cal B}_{q,J}\), its reality partner lies in
\(-{\cal B}_{q,J}\).  Taking \(\delta\) small enough that the boxes avoid
the origin, and identifying a pair at most twice, (3.3) supplies at least
\((c_\delta/2)q^3/J\) distinct reality pairs.  Since \(q/J\to\infty\), one
may select \(q^2\) pairs for every sufficiently large stage.

This is genuinely different from the shortcut excluded by C170.  In
C170's narrow sector, the projected aperture is \(O(q^{-1})\), so there
are only \(O(q)\) projected points, and the alignment tolerance is
\(O(J^{-1})\), allowing at most one normal lift.  Here the projected
aperture is fixed and the normal thickness is \(q/J\).  The projected
count is \(O(q^2)\), with \(O(q/J)\) possible normal layers.  The lower
bound (3.3) is consistent with both of those upper-count mechanisms.

### 3.2 Exact central-fiber C154 propagation of the slab

C154 gives

\[
 K^m=I-mg_0\otimes u.                                     \tag{3.4}
\]

Because \(\bar k\in u^\perp\), every \(k\in{\cal B}_{q,J}\) satisfies

\[
 |u\cdot k|\le C_\delta{q\over J}.                        \tag{3.5}
\]

Therefore, for \(0\le m\le J\),

\[
 |K^mk-k|
       =m|g_0|\,|u\cdot k|
       \le C_\delta q.                                    \tag{3.6}
\]

For the fiber over the exact central orbit, the entire carrier box remains
at frequency \(O(q)\), inside a fixed slightly enlarged angular
neighborhood when \(\delta\) is small.  Floquet decomposition supplies the
same central-fiber bounded-width conclusion between returns: the
within-period covector propagator is uniformly bounded and multiplies the
return shear (3.4).  Equations (3.4)--(3.6) are not an exact propagation
formula at neighboring base points of a localized packet.

The Kelvin generator is homogeneous of degree zero in the covector.  To
compare (3.2) with C159, normalize every nonzero carrier by a smooth
projective section (for example, unit length at the initial point) and use
the correspondingly rescaled transverse frame.  One must not insert the
integer covector \(k\sim q\bar k\) directly into C159's unnormalised
\((E_1,E_2)\) frame while keeping its fixed coordinate vector \(w\): that
would introduce a spurious \(q\)-dependent basis scaling.  After
projective normalization, (3.6) puts all carrier paths in one fixed compact
neighborhood of the normalized central path.

The strict C159 bounds

\[
 B_{12}>32,\quad B_{21}>{9\over10},\quad
 Bw>{1\over5}w                                      \tag{3.7}
\]

hold on a compact periodic position--projective-covector path and have
positive margins.  Smooth dependence of the Kelvin generator on both
variables therefore gives a fixed open phase-space tube in which a
slightly weaker cooperative inequality still holds.  The action--angle
residence estimate keeps the exact-background base points close, while
(3.6) keeps the central-fiber carrier paths in the projective tube when
\(\delta\) is chosen sufficiently small.  This gives qualitative
principal-cone compatibility for the slowly sheared, generally
nonperiodic covectors over the stated finite \(J\)-window.  It is not a
Floquet multiplier theorem for each off-ray covector, and promoting it to
a uniform localized packet evolution still requires control of the
base-dependent covector map, polarization frame, and finite-frequency
connection terms.  No explicit value of \(\delta\) or such error bound is
claimed.

Frequency thickness \(q/J\) is dual to the physical length \(Jr\) in
(2.3).  Thus the lattice correlation and the material residence repair are
the same phase-space object, not two unrelated assumptions.  C161's
terminal \(q\)-way charge step could kinematically restore an \(O(q)\)
third frequency width and an \(O(r)\) third spatial width.  Its physical,
unforced, full-polarization realization remains open by C163--C170.

## 4. The enlarged rigid collar still fits the factorial budget

To avoid hiding aspect-ratio constants, enclose the correlated tube in a
spherical rigid core of radius

\[
                           R=CJr=CJ{\ell\over q}.           \tag{4.1}
\]

Repeat C174's curl completion with \(r\) replaced by \(R\).  One switch,
parent-cross term, or self term acquires the \(L^2\) volume/amplitude
factor \(J^{5/2}\).  Summing \(J\) slices and applying
\(h=q^{3/2}\) gives

\[
 \boxed{
 \begin{aligned}
  \varepsilon_{\rm sw}^{\rm back}
      &\le C\Lambda J^{7/2}q^{-1},\\
  \varepsilon_{\rm par}^{\rm back}
      &\le CJ^{7/2}q^{-1},\\
  \varepsilon_{\rm self}^{\rm back}
      &\le C\Lambda J^{7/2}q^{-1},\\
  \varepsilon_{\nu}^{\rm back}
      &\le C{\rm Re}^{-1}J^{3/2}q.
 \end{aligned}}                                           \tag{4.2}
\]

Here \(\Lambda=\lambda\ell/a\), exactly as in C174.  The corresponding
trajectory-separated pressure-moment powers are

\[
 C(1+\Lambda)J^6q^{-7/2},
 \qquad C{\rm Re}^{-1}J^4q^{-3/2}.                        \tag{4.3}
\]

On (0.3), with \(J\le C\log h\),

\[
 J^{7/2}q^{-1}
    =O\!\left(n^{-8}(\log n)^{7/2}\right)=o(n^{-6}),      \tag{4.4}
\]

and after one \(n^2\) chart charge it is

\[
 O\!\left(n^{-6}(\log n)^{7/2}\right)=o(n^{-4}).          \tag{4.5}
\]

The switch/self conclusion assumes the explicit, slightly stronger
polylog profile condition

\[
                         \Lambda J^{7/2}=o(n^2).           \tag{4.6}
\]

Every fixed or fixed-polylogarithmic \(\Lambda\) satisfies (4.6), but
C174's weaker condition \(\Lambda\log n=o(n^2)\) does not imply it for an
arbitrary \(n\)-dependent \(\Lambda\).  The factorial schedule
\({\rm Re}^{-1}=\nu(j!)^{-2}\) absorbs the polynomial \(qJ^{3/2}\) in
the viscous line.

Equations (4.2)--(4.6) are raw prescribed-collar response ledgers.  Leray
pressure is global, the backward active propagator need not be
contractive, and the collar source has not been absorbed into an exact
unforced solution.  Those are BAFL obligations, just as in C174.

## 5. Exact no-free-rewrite identity for a fixed cocycle

Let \({\cal U}(t,s)\) be any fixed linear evolution family, and let a
sequence of physical rewrites enter as sources \(f_m\) supported in
disjoint time windows.  Variation of constants gives the exact endpoint

\[
 z(T)={\cal U}(T,0)z_0
       +\sum_m\int {\cal U}(T,s)f_m(s)\,ds.                \tag{5.1}
\]

For ideal impulses \(\xi_m\) at \(t_m\), the second term is simply
\(\sum_m{\cal U}(T,t_m)\xi_m\).  Thus

\[
 \boxed{
 \text{signed rewrite response}
       =\text{rewritten endpoint}-\text{unrewritten endpoint}.} \tag{5.2}
\]

If the signed sum cancels exactly, the final physical state is the
unrewritten cocycle state.  Coordinate recharts may describe it
differently but cannot reset its physical C154 deformation.

Suppose instead that the desired undeformed reference endpoint differs
from the unrewritten endpoint by relative amount \(\delta\) in the
normalized growing coordinate.  With scheduled active size \(b\), (5.2)
forces terminal weighted response

\[
                              \delta b.                    \tag{5.3}
\]

The active BAFL allowance is \(Cb^3\), so a fixed-cocycle rewrite must
satisfy

\[
                         \delta\le Cb^2.                   \tag{5.4}
\]

For \(b=n^{-2}\), an order-one reset exceeds the allowance by
\(b^{-2}=n^4\), independently of how its intermediate signed pieces are
phased.  This strengthens C174's triangle-inequality observation only in
the explicitly fixed-cocycle endpoint setting.

It is not a no-go for a nonlinear reservoir which changes the generator,
exchanges energy with the child, and later restores the target amplitude.
Such a mechanism is no longer a telescoping rewrite of the fixed C159
cocycle.  In unforced Navier--Stokes it must be included as an actual
velocity component, and its complement, pressure, depletion, and retained
wake must satisfy C125/BAFL.

## 6. Exact surviving theorem target

C176 removes two uncertainties from the RRR residence branch:

1. on the exact decaying \(A_2\) background, the full second-flow-jet and
   between-return bounds are polynomial and the child stays in an
   \(O(Jr)\) tube; and
2. a C154-correlated, genuinely two-directional integer carrier set has
   more than the required \(q^2\) capacity at only polylogarithmic spatial
   elongation.

It also shows that signed rewriting does not provide an independent free
repair while the target cocycle is held fixed.  The surviving load-bearing
statement is therefore narrower:

> Construct a spatially localized divergence-free finite-frequency packet
> realizing the correlated carrier-centre bundle (3.2), with a uniform
> right/left C159 growing pair and the C125 relative estimate; realize the
> terminal
> genuinely three-dimensional full-polarization converter; and prove that
> the localized collar, pressure, heat, nonlinear child, and retained C140
> wake obey the two BAFL response channels in one unforced trajectory.

The exact normal form does not control perturbations of the trajectory by
the growing child.  The qualitative principal cone does not control the
finite-frequency C148 connection and boundary terms.  The lattice count
does not prove spatial localization with controlled Fourier sidebands or
endpoint phase coherence.
The raw collar powers do not prove the non-normal Duhamel estimate.
Therefore C125, MCKC(ii), LCE, BAFL, the one-cell theorem, and any blow-up
claim remain open.

## 7. Verification boundary

`checks/a2_correlated_tube_residence_c176.py` verifies the exact normal-form
flow identities, nilpotent physical/covector shear, correlated box
propagation, an axis-aligned representative of the lattice-capacity count,
all enlarged-collar exponents, the factorial schedule comparisons, and the
fixed-cocycle signed rewrite identity.  It does not certify the standard
action-angle coordinate-existence theorem, choose the actual C159
continuity radius, construct a localized packet, or prove C125, MCKC,
LCE, BAFL, an unforced stage, or a singular solution.
