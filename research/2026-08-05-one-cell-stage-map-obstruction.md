# One-cell stage map: reduction to backward-weighted active-focus leakage

Date: 2026-08-05  
Status: **CONDITIONAL REDUCTION / SELF-DERIVED.  Not a Navier--Stokes
stage theorem and not a blow-up result.**

## 1. Purpose

The program no longer needs another cascade geometry.  It needs one
unforced Navier--Stokes map taking the structured state at stage \(j\) to
the structured state at stage \(j+1\).  The finite-dimensional ladder,
depletion, decaying-pump, factorial, focusing, pressure, viscosity, and
wake ledgers determine a candidate trajectory.  This note identifies the
one estimate that is not supplied by those scalar ledgers.

The missing package is called **backward-weighted active-focus leakage**
(BAFL): construction of the localized candidate trajectory together with
the two-channel response estimate stated below.  The phrase *BAFL estimate*
will refer only to that response estimate; the *BAFL package* also includes
construction of the trajectory on which it is measured.  For the selected
geometry, its first
unavoidable source is the explicit reality-symmetric \(A_2\) difference
shell.
It measures the full Leray-projected response to every quadratic mode that
leaves the intended ladder, with the response weighted by how much of the
remaining amplifier it can still traverse.

## 2. Exact error equation

Let \(I_j=[t_j,t_j^+]\), and let \(v_j\) denote the smooth, divergence-free
candidate one-cell trajectory obtained by embedding the intended ladder,
connector, focus, and carried wake.  It is an approximate solution of the
unforced equation:

\[
 \partial_t v_j+\mathbb P\nabla\!\cdot(v_j\otimes v_j)
       -\nu\Delta v_j=F_j.                                      \tag{2.1}
\]

For an exact solution \(u_j=v_j+e_j\), subtraction gives the identity

\[
 \partial_t e_j-L_j(t)e_j
   =-F_j-\mathbb P\nabla\!\cdot(e_j\otimes e_j),                 \tag{2.2}
\]

where

\[
 L_j(t)e=\nu\Delta e-
 \mathbb P\nabla\!\cdot(v_j\otimes e+e\otimes v_j).             \tag{2.3}
\]

If \({\cal U}_j(t,s)\) is the evolution family for \(L_j\), Duhamel's
formula is exactly

\[
 e_j(t)={\cal U}_j(t,t_j)e_j(t_j)
 -\int_{t_j}^t {\cal U}_j(t,s)F_j(s)\,ds
 -\int_{t_j}^t {\cal U}_j(t,s)
       \mathbb P\nabla\!\cdot(e_j\otimes e_j)(s)\,ds.           \tag{2.4}
\]

No mode truncation has been made in (2.2)--(2.4).

## 3. Intended and leaked residuals

Let \(\Pi_j(t)\) denote the moving projection onto the intended \(A_2\)
ladder and its pump/child coordinates, and put \(Q_j=I-\Pi_j\).  The
\(Q_j\)-state includes the explicitly retained wake.  This is bookkeeping,
not an invariant Fourier projection.  Split

\[
 F_j=F_j^{\rm paid}+F_j^{\rm leak},\qquad
 F_j^{\rm leak}=Q_j
   \mathbb P\nabla\!\cdot(v_j\otimes v_j),                      \tag{3.1}
\]

with every curvature, cutoff, heat, pressure, chart, and scheduled-wake
term placed in \(F_j^{\rm paid}\).  Put \(n=j+1\), let \(L_j\) be the
actual terminal active-chart loss, and assume \(1\le L_j\le Cn^2\).
For

\[
 w_F(t)=\int_{t_j}^t {\cal U}_j(t,s)F(s)\,ds,
\]

define the pre-chart active and retained-wake response seminorms

\[
 \begin{aligned}
  \|F\|_{{\cal R}^{\rm act}_j}
   &:=\sup_{t\in I_j}\|\Pi_j(t)w_F(t)\|_{{\cal A}_j},\\
  \|F\|_{{\cal R}^{\rm wake}_j}
   &:=\sup_{t\in I_j}\|Q_j(t)w_F(t)\|_{{\cal W}_j},\\
  \|F\|_{{\cal R}_j}
   &:=n^2\|F\|_{{\cal R}^{\rm act}_j}
      +\|F\|_{{\cal R}^{\rm wake}_j}.
 \end{aligned}                                                \tag{3.2}
\]

The active and wake norms include the remaining dynamical gain between
insertion time \(s\) and the terminal time.  The factor \(n^2\) reserves
the worst allowed finite-dimensional chart loss only for active return;
the fixed constant in \(L_j\le Cn^2\) is absorbed in later constants.
Define the analogous state norm \(\|e\|_{{\cal Y}_j}\) with the same
weighted split.  Thus
(3.2), rather than an unweighted \(L^1_tX\) norm, is the correct measure
for just-in-time activation.  The reconstructed scalar ledgers are designed
to give the still-conditional paid-residual estimate

\[
 \|F_j^{\rm paid}\|_{{\cal R}_j}=o(n^{-4}).                    \tag{3.3}
\]

### 3.1 The first leaked shell is forced by reality

The leading leak can be identified before localization.  Use the exact
integer \(A_2\) roots

\[
 k_1=(1,-1,0),\qquad k_2=(0,1,-1),\qquad
 k_3=(-1,0,1),\qquad k_1+k_2+k_3=0,
\]

put \(N=(1,1,1)\), and set \(t_i=N\times k_i\).  Every complex
polarization transverse to \(k_i\) has the form

\[
 a_i=x_i N+y_i t_i.
\]

For two Fourier coefficients define the symmetrized Euler interaction

\[
 \mathcal B_{p,q}(a,b)=
 P_{p+q}\big((a\mathbin\cdot q)b+(b\mathbin\cdot p)a\big).
\]

Direct integer arithmetic gives, cyclically in \((i,j)=(1,2),(2,3),(3,1)\),

\[
 \begin{aligned}
 \mathcal B_{k_i,k_j}(a_i,a_j)
   &=3\,(y_i x_j-y_jx_i)N,\\
 \mathcal B_{k_i,-k_j}(a_i,\overline{a_j})
   &=-3\,(y_i\overline{x_j}+\overline{y_j}x_i)N.
 \end{aligned}                                                   \tag{3.4}
\]

The first line is the intended hexagon interaction at \(-k_m\).  Reality
forces the second line at the difference wave number
\(q_{ij}=k_i-k_j\), whose length is \(\sqrt3\,|k_i|\).

There is a sharp algebraic distinction between a paired gate and a
three-leaf orbit.  For one active pair, the difference coefficient can
vanish while the sum coefficient is nonzero: one may take
\(x_i/y_i=-\overline{x_j/y_j}\) with nonzero real part.  This is the
algebraic opening for an unfolded, just-in-time paired gate.  But if all
three modes are nonzero and all three difference coefficients in (3.4)
vanish, then all three intended coefficients vanish.  Indeed, when
\(y_iy_j\ne0\), writing \(r_i=x_i/y_i\) gives

\[
 r_i=-\overline{r_j}\quad\hbox{cyclically},
\]

so \(r_1=r_2=r_3\) is purely imaginary and
\(y_ix_j-y_jx_i=y_iy_j(r_j-r_i)=0\).  If some \(y_i=0\), nonzero of all
three modes and the same difference equations force all \(y_i=0\), again
annihilating the intended interactions.

The obstruction is quantitative.  Put

\[
 D_{ij}=y_i\overline{x_j}+\overline{y_j}x_i,
 \qquad S_{ij}=y_ix_j-y_jx_i.
\]

If \(y_1y_2y_3\ne0\), then the exact identity

\[
 r_2-r_1=
 \frac{D_{23}}{y_2\overline{y_3}}
 -\overline{\frac{D_{31}}{y_3\overline{y_1}}}
\]

implies

\[
 |S_{12}|\le
 \frac{|y_1|}{|y_3|}|D_{23}|
 +\frac{|y_2|}{|y_3|}|D_{31}|.                         \tag{3.5}
\]

Thus, when the three in-plane amplitudes are comparable within a factor
\(\kappa\), an intended interaction of size \(s\) forces at least one
difference interaction of size \(s/(2\kappa)\).  A simultaneous balanced
three-leaf stage has order-one raw leakage; only temporal unfolding,
spatial export, or a cancellation in the *propagated* Duhamel response can
meet the active-return \(n^{-6}\) budget.

Thus an active real three-leaf hexagon is not an invariant Fourier cell.
At least one \(q_{ij}\) mode is produced.  The one-cell theorem must either
route paired gates before the third mode becomes exposed, or prove that the
backward-weighted response of these explicit difference modes is small.
This is the first concrete component of (BAFL), not an unspecified
``mode-leakage error.''

The reconstructed terminal calculation C116--C117 makes the timing gain
precise.  Two helicity channels can cancel the first difference output of a
single parent pair while retaining its child.  That child necessarily
creates a named length-\(\sqrt6\) parent--child output on the next Picard
iterate.  For the fixed nondegenerate paired gate in C116, if the two
parents are normalized to order one and the gate is stopped while the
short-time expansion is valid, when the child reaches the scheduled
dormant size

\[
 b_j=(j+1)^{-2},
\]

then at least one forced second-Picard mode has a nonzero leading coefficient
and natural size

\[
 \Theta(b_j^2)=\Theta((j+1)^{-4}).                            \tag{3.6}
\]

This exactly saturates the retained-wake allowance.  Under the earlier
conservative pre-propagation reading, it would be two powers too large if
passed directly through a chart with condition number \(O(j^2)\) as an
active-coordinate error.  Hence the C107 fork and terminal leakage are the
same one-cell issue: prove that (3.6) remains in a transported dark/exported
wake channel, or gain another factor \(j^{-2}\) in its propagated return to
the active chart.  Raw support and size counting alone did not determine
whether the fixed-projector return gains that factor; the desired child
and the second-Picard output scale as the first and second powers of the
same gate time.

C140--C141 resolve this raw fixed-projector alternative exactly.  The
second-Picard mode does not enter the six-root dynamics at order \(b_j^2\):
relative to the retained six-root equation, its first active child return
is the nonzero cubic term

\[
                  \Theta(b_j^3)=\Theta(n^{-6}).               \tag{3.7}
\]

The reason is structural, not an accidental cancellation.  The C116 datum
and every descendant in the fixed \(A_2\) root lattice are 2D3C.  The
planar pump evolves autonomously, while the child and wake are passive
\(N=(1,1,1)\)-directed scalars.  This closes the homogeneous short-gate
power count.  Within this split the passive scalar cannot deplete the
autonomous planar pump or supply the active three-dimensional volume
focus.  This does not refute the abstract C118/C119 normal form or its
off-plane leaves \(q_m\pm Kr_i\), \(m\ne0\), which lie outside the fixed
root lattice.

Consequently the bare fixed-projector timing objection is resolved, and
one necessary estimate is isolated at the exit from that invariant plane.
For a moving active projector \(P_j(t)\), with
\(Q_j=I-P_j\), the direct conversion exposure contains

\[
 \left(\dot P_jQ_j+P_jA_jQ_j\right)w_j,                       \tag{3.8}
\]

where \(A_j\) is the exact Leray-projected linearized conversion/focus
generator and
\(w_j=Q_ju_j\) is the \(O(n^{-4})\) wake.  The localized
conversion-exposure estimate (LCE) is the requirement that the
backward-weighted Duhamel response of (3.8) be \(O(n^{-6})\).  For fixed
Fourier \(P\), the analogous \(P\leftarrow Q\) response first occurs
through the extra parent--wake interaction just proved.  Localization,
off-plane routing, pressure, and focus remain open.  LCE is a necessary
subestimate, not a replacement for the full BAFL estimate below: newly
generated leakage, nonlinear wake feedback, and closure of the
entrance/wake state classes also remain to be controlled.

### 3.2 Active-focus audit: affine failure and the broad-band target

C142--C151 test the active-focus step without changing the \(A_2\)
geometry. The strongest global affine selector is explicit: with
\(h=q^{3/2}\), it amplifies the intended child velocity by
\(\asymp h\) while leaving the two named C117 wake velocities bounded.
This does not survive generic localization. A materially transported
three-dimensional envelope incurs the exact heat action

\[
 {4Lq^2\log n\,(h^2-1)\over3(h^2-h^{-2})}
 \sim {4L\over3}q^2\log n,                                  \tag{3.9}
\]

whereas the direct product estimate for a stationary compact affine core
has the strain-speed-independent backward-focused \(L^2\) scale

\[
                         O(\log h).                         \tag{3.10}
\]

Thus slowing the selector does not improve this estimate. Equation (3.10)
is not a universal lower bound: an algebraic or Lagrangian cancellation
could make the actual interaction smaller. The surviving form is
conditional on a co-moving Lagrangian collar which cancels the leading
parent transport and gives an integrated wake-to-active block \(O(b)\).
The direct \(O(\log h)\) estimate does not place the retained
\(b^2=n^{-4}\) wake below the \(b^3=n^{-6}\) active allowance; the ratio
of its natural scale to that allowance is \(12n^2\log n\).

The [C171 co-moving collar audit](2026-08-06-comoving-piola-collar-c171.md)
and its [checker](../checks/comoving_piola_collar_c171.py) separate the
parent-cross loss from the wake-return loss.  If a compact curl potential is
transported as a one-form by the parent flow \(X\), the resulting Piola
velocity is exactly divergence free and satisfies

\[
 U\circ X=Fw_0,\qquad
 D_t^VU=(U\mathbin\cdot\nabla)V,
 \qquad G_{\rm parent}=2(U\mathbin\cdot\nabla)V.       \tag{3.10a}
\]

Thus material transport removes the large stationary term
\(V\cdot\nabla\chi\), but the velocity-linearized strain terms add rather
than cancel.  With the full Piola/F-jet polynomial denoted by \(M_F\), its
normalized time-integrated and backward-focused scales are

\[
 M_Fq^{-5/2}\log h
 \quad\xrightarrow{\ \times h\ }\quad
 M_Fq^{-1}\log h=12M_Fn^{-8}\log n.                \tag{3.10b}
\]

Conditional on

\[
 M_F\log n=o(n^2),                                  \tag{3.10c}
\]

the last scale is \(o(n^{-6})\); one prescribed-ray chart charge
\(\kappa_{\rm ch}\le Cn^2=Cq^{1/4}\) gives
\(M_Fq^{-3/4}\log h=o(n^{-4})\).  A rigid collar following only the central
trajectory proves the same parent-cross arithmetic without \(M_F\), but it
is not material and therefore still owes boundary-crossing and retained-wake
control.  Leray projection and the pressure complement are \(L^2\)
contractions of the raw residual, but the pressure tail is global and the
backward active propagator need not be contractive.

This gain does not supply the missing factor \(b\) in the wake-to-active
block.  C171 gives an explicit real zero-charge \(A_2\) wake-to-child triad
and a separate allowed polarization for which the Piola ordering is bright;
reality partners add rather than cancel.  For the actual C140 wake, C141's
fixed-projector calculation remains the sharper statement: its first return
is the nonzero cubic term, so its factor \(b\) comes from one additional
interaction, not from charge or reality by itself.  The exact collar
obligation is therefore

> **Material-collar kernel closure (MCKC).** (i) Prove (3.10c), including
> the required first two label derivatives, or control boundary crossing
> for the rigid center-following substitute; and (ii) prove for the actual
> retained C140 wake \(w\) that its integrated wake-to-active block obeys
> \(\|{\cal L}_{WA}w\|\le Cb\|w\|\), or export it with stage-uniform
> pressure/heat tails.

MCKC settles neither the inherited self/viscous response nor LCE/BAFL.  It
is the now-explicit collar/kernel subproblem inside those estimates.

The [C172 pressure-dark audit](2026-08-06-affine-pressure-dark-material-transport-c172.md)
and its [checker](../checks/affine_pressure_dark_material_transport_c172.py)
rule out exact pressure darkness as a local replacement for MCKC on the
existing selector.  For a constant matrix \(A\), the determinant of the
pressure-dark compression on \(k^\perp\) is
\(k^T\operatorname{adj}(A)k/|k|^2\).  A rank-two matrix therefore has
dark frequencies only on the union of the planes normal to its right and
left kernels, and no nonzero \(L^2(\mathbb R^3)\) localized field can obey
\({\mathbb P}(AU)=0\) when \(\operatorname{rank}A\ge2\).  The fixed
trace-free rank-one curl-kernel exception has no material, \(L^p\), or
pointwise Piola amplitude gain.  In the displayed C142 launch frame,
\(A=\sigma(d\otimes N+N\otimes d)\); per displayed \(N\)-coefficient, the
full \(2AU\) residual has projected squared sizes \(216\sigma^2\) at the
child and \(162\sigma^2\) at each named wake.  Making all three start-frame
fibers dark for one constant selector forces \(AN=0\) and removes their
Kelvin velocity gain.  Scalar Piola is the only universal
datum/label-independent zero-order affine multiplier, and no universal
translation-invariant finite-order velocity operator reproduces the
rational Hodge generator on an open frequency cone.  Variable-coefficient
auxiliary systems and mode-tailored antidivergences are outside that last
no-go.

The [C173 paired-multipole audit](2026-08-06-paired-multipole-piola-collar-c173.md)
and its [checker](../checks/paired_multipole_piola_collar_c173.py) supply an
exact finite-coordinate repair, not the broad-band response theorem.  The
paired symmetric collars

\[
             w_{\rm pair}={128\over127}w_r-{1\over127}w_{2r}
                                                               \tag{3.10d}
\]

retain \(Sx\) on the inner core while canceling the full degree-three
moment.  Their Fourier transform vanishes to at least order five at the
origin, two orders beyond one symmetric collar.  At fixed local-core
normalization this gives a \(q^{-2}\) improvement on fixed,
chart-tracked parent-scale coordinates.  A material-label triple-shell
polynomial of spatial order 12 and coefficient \(\ell^1\)-norm \(343/64\)
preserves the core and kills the named \(|\zeta|^2=2,6\) fibers through two
Fourier derivatives, including the affine coefficient and its formal first
curvature variation.  This removes the finite named-fiber part only.  A
robust underlying-jet-independent two-jet notch on \(q\) distinct axial
labels has scalar degree at least \(3q\), or radial differential order at
least \(6q\), outside the current Gevrey budget.  The actual C140 profile
normalization, moving active bundle, curvature remainder, full local
\(L^2\)/non-normal response, and periodic pressure tails remain open.

The [C174 time-sliced rigid-collar audit](2026-08-06-time-sliced-rigid-collar-c174.md)
and its [checker](../checks/time_sliced_rigid_collar_c174.py) test the rigid
branch by slicing the gain window into \(J=O(\log h)\) prescribed pieces.
Under its bounded profile-variation hypothesis,
\(\Lambda\log n=o(n^2)\), factorial Reynolds decay, and the stated
trajectory-separation assumptions, the switch, parent-cross, self,
viscous, and instantaneous pressure-tail ledgers integrated along a
separated trajectory all fit the schedule.  These are residual bounds for
a prescribed approximate trajectory, not an unforced operation.  An exact
coordinate rechart satisfies

\[
 \widetilde P_{J-1}\cdots\widetilde P_0
 =Q_J^{-1}(P_{J-1}\cdots P_0)Q_0.                \tag{3.10e}
\]

It represents the same physical cocycle and produces no PDE source, but it
is a Floquet similarity only when \(Q_J=Q_0\).  It does not reset physical
support, bandwidth, or C154's accumulated shear.  At the central first-jet
return times an isotropic radius shrink \(O(J^{-1})\) is sufficient for
residence, with reciprocal-volume \(O(J^3)\), fixed-energy amplitude
\(O(J^{3/2})\), and bandwidth \(O(J)\) costs.  A uniform second-flow-jet
bound and residence between returns remain open.

A physical rewrite has a different ledger.  If \(\delta_m\) is its relative
mismatch in the normalized growing child coordinate, then its individual
terminal-weighted norm is exactly \(\delta_m b\).  Hence a robust
triangle-inequality proof of the \(O(b^3)\) active allowance requires

\[
                         \sum_m\delta_m=O(b^2).    \tag{3.10f}
\]

The sum of individual norms is not a lower bound for a specially phased
signed Duhamel response.  This leaves the exact **rechart-or-rewrite
obstruction (RRR)**: either prove finite-tube residence/coherence under the
C154 shear, including second jets and between-return control, or prove the
signed physical-rewrite cancellation directly while preserving the target
cocycle.  Neither branch supplies MCKC(ii)'s broad-band wake kernel/export
bound or absorbs the prescribed collar residual into an unforced stage.

C176 resolves the exact-background part of the first branch.  In the
regular action--angle--axial annulus, the full flow jets through the gain
window satisfy \(D\Phi=O(J)\), \(D^2\Phi=O(J^2)\), and a correlated packet
has an \(O(Jr)\) bounding tube at every intermediate time.  A fixed-aperture
\(q\times q\times(q/J)\) carrier slab contains
\(\Omega(q^3/J)\) unordered reality pairs; on the central fiber, C154 keeps
its frequency at \(O(q)\) through \(J\) returns.  Enlarging the rigid collar
to the tube costs only
\(O((1+\Lambda)J^{7/2}q^{-1})=n^{-8}\operatorname{polylog}n\).
This does not propagate the localized base-dependent finite-frequency
packet, and it does not turn the prescribed collar into an unforced field.

C175 sharpens the second clause of MCKC.  A common material eikonal preserves
the exact relation \(e_1-k_1=k_c\), and the allowed-class Leray--heat Born
response on that edge has one sign and an order-one inviscid limit.  For a
specified active/wake splitting, the graph equation

\[
 K'=AK+B-KD-KCK
\]

implies the exact Melnikov identity

\[
 [\ell Kr]_{t_0}^{t_1}
 =\int_{t_0}^{t_1}(\ell Br-\ell KCKr)\,dt.        \tag{3.10g}
\]

Therefore an \(O(b)\) graph on the stated admissible wake class requires
the corresponding first-order transfer to be \(O(b)\) already.  This is
**RIGM**, the resonant invariant-graph Melnikov gate.  It can still be met
by signed packet cancellation, an \(O(b)\)-action just-in-time exposure, or
spatial export with controlled global tails; no splitting-independent
no-go is claimed.

C177 supplies an exact autonomous half-repair to that gate.  If the old
pump and a preloaded reservoir share one curl eigenvalue, their heat-decaying
sum is an exact unforced Beltrami background.  The reservoir then appears in
the perturbation equation only through a homogeneous tangent, so there is no
absolute writer to divide by the \(n^{-28}\) seed.  Under C177's explicit
propagated-tangent hypothesis the relative charge is only
\(O(b\log n)=o(1)\).  The repair is not yet the converter: distinct curl
components create the exact additive block

\[
 \varepsilon(\kappa_1-\kappa_2)
      {\mathbb P}(U_0\times R_0),                 \tag{3.10h}
\]

and one common shell has only a scalar heat envelope, so it cannot supply
the prescribed noncommuting polarization clock.  Recovering full physical
polarization on a same-curl multi-direction shell, or canceling (3.10h), is
the **same-curl reservoir realization gate (SCRG)**.

C178 addresses the separate question of absorbing the prescribed collar
ramps into entrance data.  A fixed nonzero compact profile times a monotone
scalar ramp is not in the exact backward heat range of any \(L^2\) datum.
There are nevertheless two constructive repairs.  Forward buffering

\[
 g=e^{\sigma\Delta}h,
 \qquad
 z_0=\int_0^T\vartheta'(s)e^{(\sigma-\nu s)\Delta}h\,ds,
 \qquad \sigma\ge\nu T,                            \tag{3.10i}
\]

gives exact pure-heat terminal cancellation with an \(L^2\)-contractive
entrance datum, at the price of analytic global tails.  Alternatively,
finite polynomials in \(-\Delta\) give compact divergence-free entrance
data whose terminal errors tend to zero.  On the tuned C142 rate the scalar
preparation is factorially smaller than the C125 seed.  None of these heat
identities propagates through the actual localized \(A_2\) linearization or
controls its nonlinear, pressure, and retained-wake response; those remain
inside RIGM/MCKC/BAFL.

The coherent alternative has an exact kinematic target but an exact
endpoint obstruction. The existing lattice contains a divergence-free
packet with \(q^3\) conjugate pairs (\(2q^3\) nonzero Fourier wavevectors)
and point gain \(cq^{3/2}\) in every parent period cell. However, at the
scheduled energy and point amplitude, one isolated child requires at
least

\[
                         (qK)^3                             \tag{3.11}
\]

coherent Fourier degrees. An exact or relatively narrow C121 shell at
radius \(O(qK)\) has only \(O((qK)^2)\). Consequently the focused endpoint
cannot be the exact global six-root pump; it must be a localized
C121-like core with a non-shell collar/wake. A one-step quadratic
collapse of \(M\) coherent sources to six pump targets needs at least
\(M/6\) oriented gate frequencies and has at least \(M/6\)
support-allowed oriented companions before polarization cancellation
(\(M/12\) if unoriented reality pairs are counted).

Temporal unfolding remains arithmetically possible. If the per-gate
target contributions are aligned with size \(\asymp b\theta\) and the
wake/active returns obey \(O(b\theta^2)\), \(O(b\theta^3)\), then
\(J=b^{-1}\) microgates of strength \(\theta=b\) give

\[
 \text{target}\asymp Jb\theta=b,\qquad
 \text{wake}\lesssim Jb\theta^2=b^2,\qquad
 \text{return}\lesssim Jb\theta^3=b^3.                       \tag{3.12}
\]

The unresolved issue is uniformity over the many simultaneous coherent
channels.

The reservoir variant moves that uniformity into C125. A point seed
\(s=n^{-16}\) on child volume \(q^{-3}\) has growing-coordinate size
\(\varepsilon=n^{-28}\); gain \(H=n^{26}\) reaches the scheduled point,
\(L^2\), and energy targets. But a same-support writer term of size
\(26n^{-20}\log n\) is larger than the seed by
\(26n^8\log n\). Absolute summability therefore gives no retained gain.

The exact axial-layer linearization is the variable-fiber operator C148,
not the finite C120 matrix and not a Toeplitz block. C149 settles one local
ray of its principal high-frequency problem: a nonsymmetric three-root C121
pump has a vertical elliptic ray with bounded periodic covector and a
simple expanding Floquet multiplier. Its dimensionless inertial-time
exponent is

\[
 {9\over16}\sqrt3(2+\delta)\epsilon+O(\epsilon^2),
 \qquad \epsilon={1-\delta\over2+\delta}>0.                  \tag{3.13}
\]

For an inviscid or instantaneously frozen physical background
\(P_0U_\delta(Kx)\), this exponent is multiplied by \(P_0K\) per unit
physical time; the unforced heat-decaying pump instead uses the C149
inertial-time reparametrization.

This removes the hyperbolic frequency-stretch objection at the local
principal-cocycle level, and C127 makes its \(O(\log n)\) heat ledger
harmless. It does not yet provide a localized packet. For fixed small
\(\epsilon\), the selected critical-point ray has nonzero axial speed and
travels \(O(\log n)\) parent lengths during the gain interval, requiring a
co-moving or long core.

C152 supplies an exact zero-drift alternative on this same pump. For
\(\delta=4/5\), the regular level \(f_\delta=0\) is one contractible closed
streamline in the phase torus. Its full three-dimensional linearized return
is the rank-one shear \(F=I+u\otimes g_0\), and its periodic-covector plane
\(u^\perp\) contains an explicit off-plane direction. The Kelvin amplitude
monodromy on that direction has numerical trace
\(16716.8837799\ldots\); C153 retains this number only as orientation.
C159 now supplies the theorem-grade sign. It reduces the periodic covector
to one scalar quadrature and, in a periodic transverse frame, certifies by
outward arithmetic that a fixed positive vector is a strict cooperative
subsolution. Hence \(\rho(M)>1\); the exact determinant-one identity gives
\(\operatorname{tr}M>2\). No floating trace or precomputed amplitude column
is a premise. C154 also gives the exact fiberwise shear
\((F^{-T})^\ell=I-\ell g_0\otimes u\). An \(O(q)\) band with nondegenerate
width in the shearing direction therefore grows to \(\Theta(q\log q)\)
over \(\Theta(\log q)\) returns. Correlation of every displacement with
\(u^\perp\) removes the shear but leaves a two-dimensional fiber; a
three-dimensional candidate must instead pay a \(q/\log q\) initial width
in the shearing direction or provide another correlation mechanism.

Normal charge supplies an exact partial nonlinear protection. If the
retained real charge band is \(S\cup(-S)\), with
\(S=[m_0-\Delta,m_0+\Delta]\) and \(m_0>3\Delta\), then

\[
 P_HB(v_H,v_H)=0,                                           \tag{3.14}
\]

and degree three is the first support-allowed active return. However, C151
shows that this grading does not remove derivatives. Two explicit limiting
unstable rays have the exact principal scales

\[
 \text{wake}\asymp QA^2,\qquad
 \text{growing-line return}\asymp Q^2A^3,                   \tag{3.15}
\]

with nonzero return projection \(-9/40\); the displayed same-charge \(A_2\)
lattice sequences converge to this limiting symbol. For fixed small
\(\epsilon>0\), nearby retuned integer directions retain it by continuity,
but their microlocal PDE realization remains open.

We call (3.15) the **derivative-amplified cubic return (DACR)**. For one
quarter-period pair in the limiting cocycle, C155 closes the
zero-incoming-wake, first-period averaging question. Its complete selected
coefficient is

\[
 \mathfrak C_0={3\over64}\left[-24\pi-3
 +3\cos\!\left({4\sqrt{10}\pi\over5}\right)
 +4\sqrt{10}\sin\!\left({4\sqrt{10}\pi\over5}\right)\right]
 <-{21\over8}. \tag{3.16}
\]

The sign persists for continuously normalized right/dual-left branches at
all sufficiently small positive \(\epsilon\). This rules out Floquet
averaging as a cancellation for that pair, but it does not produce an
iterable scalar map: the sum wake is nonzero after the period and must be
retained in the state. Moreover the resulting \(Q^{-1}\) balance is a
single-Fourier-coefficient scale. C147's \(q^{-2}\) quantity is a total
coherent point amplitude over \(q^3\) modes, so no packet gain threshold
follows without the aggregate normalization and all cross-channel sums.

C156 gives one exact fixed-ring boundary: cancellation of the whole
quadratic complement on one limiting resonant ring forces support onto one
Fourier ray. Its proposed stronger secular-coercivity boundary does not
survive a shared-wake audit. The coefficient \(\mathfrak K(-t)\) comes from
reinitializing the neutral difference wake after exchanging the parents;
it is not the reverse-parent coefficient in the same two-mode evolution.
With the one causal difference wake shared by both equations, exact
reevaluation gives

\[
 \mathfrak K_{0\leftarrow\phi}(t)=\mathfrak K(t),\qquad
 \mathfrak K_{\phi\leftarrow0}(t)=-\mathfrak K(t).           \tag{3.17}
\]

At quarter separation these are \(-9\pi/8\) and \(+9\pi/8\). Their
retained-energy secular sum is zero, as the Euler energy identity requires.
The former negative pair-symmetric formula and its scalar-phase no-go are
therefore withdrawn.  In fact the corrected antisymmetric fixed-ring
normal form has the explicit strictly positive three-ray balance (4.17) of
C156. This is a genuine candidate invariant distribution, although it does
not cancel the quadratic wake or survive finite \(\epsilon\) by itself.

C157 begins the required radial thickening. For two distinct radii, two
reflected pair decompositions have the same quadratic output but linearly
independent projected wake vectors. Their two scalar pair products cannot
cancel that output. The determinant is proportional to \(A-B\), so its
inverse loses one power of radial separation and can cost \(O(q)\) on
adjacent \(q^{-1}\)-spaced layers. A full slab may have more decompositions,
which this rank-two test does not classify.

C158 fixes the reduced-half-lattice diagonal normalization. For \(M\)
equal-modulus positive-frequency representatives and an order-one bounded
diagonal kernel,

\[
 \Lambda_{\rm up}={Q^2A_{\ell^1}^2\over M}.                 \tag{3.18}
\]

At \(M=q^3\), \(Q=q\), this is \(q^{-5}\) at the C147 seed, one at the
formal coefficient scale \(q^{1/2}\), and \(q^{3/2}\) at the target.
Reality doubles the Fourier energy, while the selected cone supplies only
fixed-factor comparability between \(A_{\ell^1}\) and the physical coherent
point component. This remains an upper-bound ledger, not a lower bound or a
recovered cubic no-go.

C160 tests whether radial detuning supplies the missing cancellation. For
two unequal radii and fixed interior nonzero angular separation, the unique
periodic sum- and difference-wake return means cancel exactly, and the causal
directed return has zero Cesàro long-time mean. However,

\[
 \omega_-\asymp {|A-B|\over A+B}.                           \tag{3.19}
\]

Thus adjacent normalized layers separated by \(q^{-1}\) detune only after
time \(O(q)\), far longer than the \(O(\log q)\) gain window. This is not a
useful short-stage cancellation estimate.

C161 avoids carrying all \(q^3\) modes through that nonlinear window. It
retains a \(q^2\)-mode source packet during the long gain, matching the real
two-dimensional capacity suggested by C154's periodic-covector plane, and
creates the third Fourier coordinate only at the endpoint. This dimension
count does not supply \(q^2\) exact torus frequencies in that plane; C170
below makes its arithmetic boundary explicit. On the cofinal even-\(n\)
schedule, a reality-symmetric set of \(q\) nonzero pure-normal shifts gives
\(q^3\) distinct first daughters, while the
normalized abstract bright state

\[
 B_q=q^{-1/2}\sum_{a=1}^q e_a                              \tag{3.20}
\]

preserves \(\ell^2\) and changes the half-lattice coefficient-\(\ell^1\)
scale from \(bq=n^6\) to \(bq^{3/2}=n^{10}\). The pure-normal gate bundle
has zero gate--gate Euler self-interaction. Splitting the abstract quarter
rotation into \(J=b^{-1}\) factors meets the conditional C146 powers
\(b,b^2,b^3\). These facts are support and Hilbert-space arithmetic, not a
physical Leray propagator or a point-coherence lower bound.

C162 rules out the most direct circular realization of (3.20). At its
tuned height the forward edge is nonzero, but the exact forward/reverse
product is

\[
 \kappa_{\rm f}\kappa_{\rm r}=-{1\over4}+{i\over4},          \tag{3.21}
\]

and the reverse vector also leaves the selected source line. No phase or
positive diagonal rescaling turns this specified block into a skew star.
Other polarizations and time-dependent/full gate bundles remain open.

C163 repairs the complex phase, but only after compressing to the selected
source and forward-bright lines.  For the fixed equal-magnitude
dual-helicity gate \(E_0=(e_x-e_y)/\sqrt2\), signed height
\(-1\le y<\infty\),
and \(h=\cos\phi-\sin\phi\), the compressed product is

\[
 -A^2R,\qquad
 R={h^2\over2}\left[1+
 { (\cos\phi+\sin\phi)^2y^2(y-2)
    \over4(y^2+2y+4)}\right].                       \tag{3.22}
\]

It is strictly negative away from \(h=0\), and on the symmetric
reality-safe range \(|y|\le1\) one has \(R\ge h^2/4\).  The selected
compression is diagonally similar to a skew block.  The full reverse
interaction, however, has an explicit same-wavevector component in the
orthogonal source polarization, so this weighted rotation is not a closed
physical subsystem.  Static synchronization also fails on an angular
interval: for any nonzero shared intensity sum \(W\), the squared rate has
the rigid linear coefficient \(-W\sin(2\phi)/2\), independent of the gate
heights.  Full two-polarization evolution or a genuinely time-dependent
unfolding is therefore load-bearing.

C164 resolves the complete frozen two-polarization block.  Its source
round-trip eigenvalues are

\[
 u^2,\qquad
 u^2{(2-y)(y+2)^2\over2(y^2+2y+4)}.                \tag{3.23}
\]

A nonzero reality pair generically shares only the tangential source line.
On the whole charge lattice that line is not a hidden focusing channel: for
an arbitrary real pure-normal shear bundle it satisfies

\[
 Z_t=-i(p_h\mathbin\cdot W_h)Z,\qquad
 \|Z(t)\|_\infty\le e^{-\nu|p_h|^2t}\|Z(0)\|_\infty. \tag{3.24}
\]

Thus coefficient spreading on the common line is phase cancellation in
the exact linearized scalar fiber. C165 also closes the most direct
synchronization controls negatively. A prescribed common chirped two-level
pulse has an explicit \(O(J^{-1})\) population error on a compact rate
interval, but its normalized envelope area is
\(2J/3\), not the \(O(1)\) action in C161, and its final phase is generally
rate dependent.

Time ordering does escape the common-line algebra at finite dimension.
C166 proves

\[
 [H_r,H_t]e_\sigma={2\over7}e_t                 \tag{3.25}
\]

for the reality pair \(y=\pm1\), and its explicit first-neighbour two-pulse
Galerkin endpoint has point ratio \(\sqrt{3/2}\). The raw
coefficient-energy factor \(113/112\) is prescribed-pump work, not a
conservative splitter or an unforced pulse. C167 then rules out the static
independent-edge scalar star and its scalar-return full-fiber relaxation on
an open angular/radial sheet: reality forces the non-tangential common slope
\(2\sigma i/y\), incompatible across distinct heights, while the
pair-summed full-fiber return retains the strict gap

\[
 {2y^2(y^2+8)\over y^4+4y^2+16}>0.                \tag{3.26}
\]

C168 removes C166's first-neighbour truncation for one source.  The radial
gate is a bounded full-integer nearest-neighbour operator and the
tangential gate is globally square-zero. A unit tangential pulse has exact
perturbation **point** gain \(\sqrt2\) and half-lattice coefficient energy
\(3/2\); after a nonzero radial action \(1/100\), the all-walk tail is below
\(1/32\), the point gain remains above \(13/10\), and
the half-lattice endpoint energy is below three.  This is an exact
prescribed-pump linear statement, not a nonlinear or unforced stage.

C169 gives the exact nonlinear interpretation.  Every real fixed-plane
2D3C solution splits into autonomous two-dimensional Navier--Stokes and a
passive transverse scalar.  If the latter initially has \(q_{\rm s}\)
signed modes and coefficient norm \(G_t\), then

\[
 \|\Theta(t)\|_\infty\le\sqrt{q_{\rm s}}G_t.       \tag{3.27}
\]

For C161's literal \(q_{\rm s}=q\) convention, a coherent transverse target
\(cbq^{3/2}\) requires \(G_t\ge cbq\), or half that for the difference of
two full solutions. Here \(G_t\) is a coefficient-\(\ell^2\) norm and its
energy is \(G_t^2\). Thus the fixed-plane completion contradicts the
additional scheduled gate **norm** \(G=O(b)\); the C168 linear gain
measures sensitivity of a pre-existing scalar reservoir.  A viable terminal
converter of the C161 charge-star type must therefore leave the fixed plane
genuinely three-dimensionally. A wholly different theorem focusing an
in-plane component remains logically open, but it would replace rather than
realize the pure-normal source-to-bright block.

C170 closes the aligned-source-sheet shortcut at the C159 amplifier.  The
C152 periodic plane has the exact form

\[
 \mathcal P=\operatorname{span}_{\mathbb R}\{d,r+\sigma N\}, \tag{3.28}
\]

and horizontal projection is an isomorphism on \(\mathcal P\). If
\(\sigma\notin\mathbb Q\), the exact integer intersection is only
\(\mathbb Z d\); if \(\sigma\in\mathbb Q\), it has rank two, and C170 does
not decide which case holds. Independently of that arithmetic, the ambient
projected integer lattice has capacity \(O(\delta q^2+q)\) in the tested
angular aperture \(\delta\) about the C159 horizontal ray at size \(O(q)\).
Thus its \(C/q\)-sector contains at most \(O(q)\) projected modes, and its
exact-plane subset cannot contain more. The central C159 ray has exactly
scale-independent monodromy. Nearby real
exact-plane directions inherit hyperbolicity only as a qualitative smooth
principal-cocycle consequence; along their individual expanding lines the
relative gain varies by \(1+O((\log q)/q)\) through \(O(\log q)\) returns.
C154 forces an off-plane lift that remains aligned through
\(\Theta(\log q)\) returns to satisfy

\[
 |u\mathbin\cdot\Delta k|=O((\log q)^{-1}),         \tag{3.29}
\]

leaving at most one integer normal lift per projected point.  Hence the
required \(q^2\) growing packet must occupy a genuinely two-directional
sheet and cannot share one source-adapted pure-normal polarization.

The actual focused packet has the multi-radius and multi-charge thickness
required by C145, but C161 postpones that thickness until a terminal
charge-star interval. The surviving stage theorem must first construct a
finite-frequency \(q^2\) lattice/microlocal packet in a fixed hyperbolic
neighborhood of the certified C159 ray, including near-plane arithmetic,
C154 shear, and a genuinely two-directional polarization chart. It must
then replace the failed pure-normal source-to-bright block by an unforced,
genuinely three-dimensional, time-ordered full-polarization converter with
physical endpoint coherence, source depletion, and repeated-charge
collision control. A different in-plane focusing theorem would be a
replacement architecture inside the same stage, not a completion of C161's
block. Finally the construction must prove C125 relative to the \(n^{-28}\)
seed, land in the localized next-pump core, and satisfy LCE. This remains a
component of BAFL, not a replacement for it.

## 4. The unresolved BAFL estimate

For a candidate trajectory with the advertised pump, routing, activation,
focus, and wake modules, the **BAFL estimate** is equivalently the pair

\[
 \boxed{
 \begin{aligned}
  \|F_j^{\rm leak}\|_{{\cal R}^{\rm act}_j}
   +\|{\cal U}_j(\cdot,t_j)e_j(t_j)\|_{{\cal Y}^{\rm act}_j}
      &\le c_0n^{-6},\\
  \|F_j^{\rm leak}\|_{{\cal R}^{\rm wake}_j}
   +\|{\cal U}_j(\cdot,t_j)e_j(t_j)\|_{{\cal Y}^{\rm wake}_j}
      &\le c_0n^{-4}.
 \end{aligned}}                                                \tag{BAFL}
\]

Up to a fixed change of \(c_0\), this is the combined bound
\(\|F_j^{\rm leak}\|_{{\cal R}_j}+
\|{\cal U}_j(\cdot,t_j)e_j(t_j)\|_{{\cal Y}_j}\lesssim n^{-4}\).

It must hold for the full spatially localized trajectory, not merely for
the finite ladder ODE.  It includes:

1. leakage created while the weak seed is routed;
2. leakage created during frequency conversion;
3. leakage exposed to the \(R_j^{3/2}\) active focus;
4. pressure-mediated return from spatially separated gates; and
5. feedback of the retained charged and zero-charge wakes.

An exact finite Fourier circuit cannot be substituted for (BAFL); the
off-ladder response is part of the state.  Conversely, (BAFL) does not
require those modes to vanish.  It permits an \(O(n^{-4})\) retained wake,
but requires every component that returns to the active ladder to be
\(O(n^{-6})\) before the chart.  The latter gives an \(O(n^{-4})\)
chart-weighted error uniformly over all \(L_j\le Cn^2\).

## 5. Conditional one-cell theorem

Assume the following uniform statements for all sufficiently large \(j\).

1. The candidate \(v_j\) starts in the stage-\(j\) chart and its intended
   coordinates end at the stage-\((j+1)\) pump, with the scheduled retained
   wake.
2. The terminal chart has condition number at most \(C(j+1)^2\).
3. Equation (3.3) holds.
4. The space measured by \({\cal Y}_j\) is a complete solution space for
   (2.4), and the bilinear Duhamel map obeys

\[
 \left\|\int_{t_j}^{\cdot}{\cal U}_j(\cdot,s)
 \mathbb P\nabla\!\cdot(f\otimes g)(s)\,ds\right\|_{{\cal Y}_j}
 \le K\|f\|_{{\cal Y}_j}\|g\|_{{\cal Y}_j}                  \tag{5.1}
\]

Let

\[
 Z_j:=c_0n^{-4}+\|F_j^{\rm paid}\|_{{\cal R}_j},
 \qquad KZ_j\le\frac18.                                     \tag{5.2}
\]

5. (BAFL) holds.

Then the exact unforced Navier--Stokes solution through that cell exists on
\(I_j\), remains within \(2c_0n^{-4}+o(n^{-4})\) of \(v_j\) in
\({\cal Y}_j\), and both its terminal active coordinates and its retained
wake differ from the stage-\((j+1)\) target by \(O(n^{-4})\).

Indeed, (2.4), (5.1), and the two BAFL channels give

\[
 E_j\le Z_j+K E_j^2,\qquad Z_j=c_0n^{-4}+o(n^{-4}).             \tag{5.3}
\]

On the ball \(E_j\le2Z_j\), the right side is at most \(2Z_j\) when
\(4KZ_j\le1\), and the same estimate makes the Duhamel map a contraction.
The active-chart loss is already included in \({\cal Y}_j\); it is not
charged a second time.  Since \(\sum_j(j+1)^{-4}<\infty\), both the
structured endpoint error and retained wake are summable.

This proves a **conditional one-cell closure theorem**, not (BAFL).

## 6. Why this is the necessary gate for a robust perturbative closure

Introduce a parameter \(\epsilon\) multiplying the leaked source in the
error equation, with all other data fixed.  Differentiating this
parametrized Duhamel problem at \(\epsilon=0\) gives

\[
 \partial_\epsilon e_j(t)|_{\epsilon=0}
 =-\int_{t_j}^t{\cal U}_j(t,s)F_j^{\rm leak}(s)\,ds.            \tag{6.1}
\]

Therefore any stage estimate that is uniformly Lipschitz with respect to
an independently injected leak and has terminal tolerance \(O(n^{-4})\)
requires the corresponding linear response to obey the BAFL split.  The
active component requires the uniform pre-chart bound \(O(n^{-6})\), while
the retained-wake component may be \(O(n^{-4})\).  This is a
derivative-level obstruction for a robust perturbative proof.  It is not a
logical necessity for a single isolated exact trajectory: such a trajectory
could use a nonperturbative cancellation between sources that the separated
BAFL norm intentionally does not credit.

## 7. Exact target for the next proof

The target is now:

> **Prove (BAFL) for one spatially localized \(A_2\) stage, or exhibit one
> explicit leaked interaction whose retained-wake response is larger than
> \(cn^{-4}\), or whose active pre-chart response is larger than
> \(cn^{-6}\), for every allowed
> routing/focusing choice.**

The first outcome closes the one-cell map after the recorded modules are
cross-audited.  The second kills this stage geometry by a named mode and a
quantitative lower bound.  Neither outcome requires opening another
cascade architecture.

Within this target, the active-focus/gain construction is no longer
open-ended: the generic transported high-frequency envelope is excluded by
(3.9), the direct stationary-cutoff estimate fails at (3.10), and the
exact-shell endpoint is excluded by (3.11). The surviving obligation is LBRG producing
the C159-near \(q^2\) finite-frequency packet with C154/C170
arithmetic-and-shear control, followed by an unfolded localized endpoint
conversion that leaves the fixed plane through an unforced genuinely
three-dimensional noncommuting full-polarization evolution. It must include
normalized depletion, physical point coherence, and repeated charge
collisions. Its localization must also prove MCKC, including the global
pressure/heat export alternative, and separately prove the propagated
self/viscous bounds; all modules are measured in BAFL's two response
channels. A specially canceled stationary cutoff is included in that
surviving obligation, not excluded by (3.10).

C172--C176 make that collar obligation finite and explicit.  Exact affine
pressure darkness cannot close the rank-two selector.  C173 removes the
named finite material-label fibers but not the moving broad bundle.  On the
rigid-sliced branch, C174 requires one of the two RRR outcomes: prove
finite-tube residence and coherence under C154's physical shear, with the
second-jet and between-return bounds, or prove a signed rewrite cancellation
which is not available from triangle inequalities unless
\(\sum_m\delta_m=O(b^2)\).  C176 now closes the exact-background residence
part of the first branch: the action--angle--axial flow has first and second
jets \(O(J)\) and \(O(J^2)\), a correlated \(q\times q\times(q/J)\) slab
contains the required \(q^2\) carrier centres, and the enlarged collar costs
only \(q^{-1}\operatorname{polylog}q\).  This is a central-fiber/principal
compatibility theorem, not the finite-frequency localized propagator.

C175 also identifies the remaining broad-band kernel precisely.  Common
material phases preserve the named wake resonance, and a uniformly
\(O(b)\)-tilted invariant graph over a specified admissible wake class can
exist only when the corresponding first-order Melnikov transfer is already
\(O(b)\).  Call this residual operator estimate **RIGM**.  Thus the surviving
MCKC(ii) alternatives are a signed full-packet RIGM cancellation, an
\(O(b)\)-action just-in-time exposure, or genuine spatial export with global
pressure/heat control.  They must still be combined with propagated
self/viscous control, the localized finite-frequency C159 packet, endpoint
coherence, the genuinely three-dimensional converter, and absorption of
every prescribed collar or rewrite source into one unforced Navier--Stokes
trajectory.  The finite-notch, time-sliced, residence, and graph identities
alone do not constitute a stage map.

## 8. Verification boundary

`checks/one_cell_stage_map_obstruction.py` verifies the exponent,
bootstrap, and summability arithmetic used above.  It cannot verify the
PDE evolution-family estimate (BAFL), which remains the umbrella analytic
obstruction in this reduction. LBRG is its now-explicit active-focus/gain
component.  `checks/comoving_piola_collar_c171.py` verifies C171's exact
curl/Piola, Helmholtz, Fourier-triad, and scale arithmetic, but not MCKC,
LCE, BAFL, or the unforced stage.
`checks/affine_pressure_dark_material_transport_c172.py` verifies the C172
fiber classification, selector brightness, Piola uniqueness, and constant
symbol obstruction.  `checks/paired_multipole_piola_collar_c173.py`
verifies the C173 moment/notch and fixed-coordinate scale algebra.
`checks/time_sliced_rigid_collar_c174.py` verifies the C174 prescribed-ramp,
rechart, return-time residence, and rewrite ledgers.  They do not verify the
broad-band Duhamel response, second-flow-jet/between-return residence,
signed rewrite cancellation, RRR, MCKC, LCE, BAFL, or an unforced stage.
`checks/lagrangian_resonant_graph_obstruction_c175.py` verifies C175's
material resonance, allowed-class Leray--heat response, and finite-block
Riccati/Melnikov identities, but not the actual broad-packet RIGM estimate.
`checks/a2_correlated_tube_residence_c176.py` verifies C176's exact
background flow, central-fiber slab, lattice-capacity, enlarged-collar, and
fixed-cocycle ledgers.  It does not verify the localized finite-frequency
packet, converter, C125, MCKC, LCE, BAFL, or an unforced stage.
`checks/same_curl_unforced_reservoir_c177.py` verifies C177's exact
same-curl darkness, homogeneous perturbation structure, cross-curl source,
heat-flatness, and scale ledgers.  Its C125 conclusion remains conditional
on the propagated-tangent bound, and it does not verify SCRG or a stage.
`checks/compact_heat_preparation_c178.py` verifies C178's pure-heat Duhamel,
analytic-buffer, compact polynomial-preparation, and schedule arithmetic.
It does not verify the full \(A_2\) evolution, pressure tails, nonlinear
closure, C125, RIGM, MCKC, LCE, BAFL, or an unforced stage.
