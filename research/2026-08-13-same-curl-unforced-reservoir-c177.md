# C177: a same-curl reservoir is an exact unforced homogeneous microgate, but heat cannot time-order it

**Date:** 2026-08-13
**Status:** exact same-curl unforced-background, homogeneous-tangent,
cross-curl, heat-profile, and scale identities; conditional C125 relative
implication under an explicit propagated tangent bound; full polarization,
the charge band, depletion, RIGM, BAFL, and the stage remain open
**Checker:**
[checks/same_curl_unforced_reservoir_c177.py](../checks/same_curl_unforced_reservoir_c177.py)

## 0. Claim boundary

C165, C166, C168, and C174 use prescribed polarization or collar pulses.
C175 then shows that every uncanceled resonant exposure must have only
\(O(b)\) action, or satisfy an equally strong signed Melnikov/export
estimate.  This note tests the narrowest autonomous repair on the unchanged
A2 geometry: preload a small plane-wave/Beltrami reservoir and let its exact
unforced Navier--Stokes evolution provide the microgate.

There is a genuine exact positive result.  If the reservoir and the old
pump have the same curl eigenvalue, their sum is again Beltrami.  The
complete sum heat-decays as an exact unforced Navier--Stokes background.
After perturbing about that background, the reservoir contributes only a
homogeneous linear tangent \(\varepsilon L_Gv\).  There is no additive
writer independent of the microseed.  Consequently the C147
same-support-writer ratio \(n^8\log n\) is absent on this exact manifold;
under the explicit propagated-tangent hypothesis in Section 2, the C125
relative charge is only \(O(b\log n)=o(1)\).

The repair is narrow.

* A same-curl heat orbit has a fixed spatial shape.  Heat multiplies every
  coefficient by one positive scalar, so it cannot switch, rotate
  polarization, chirp phase, or implement noncommuting time colors.
* On the factorial schedule, every polynomial-frequency heat clock is flat
  to superpolynomial accuracy through polynomial or logarithmic action
  windows.  Order-one heat timing requires a factorial frequency ratio or
  an independently justified cancellation of large reservoirs.
* If two Beltrami components have curl eigenvalues \(\kappa_1\ne\kappa_2\),
  their exact projected cross interaction is

  \[
   {\mathbb P}\{(U\cdot\nabla)R+(R\cdot\nabla)U\}
   =(\kappa_1-\kappa_2){\mathbb P}(U\times R).       \tag{0.1}
  \]

  Thus an \(O(b)\) opposite-helicity or different-shell reservoir
  generically restores an additive source of size \(b\) times the
  normalized curl gap and the displayed Leray block.  It is \(O(b)\) only
  for an order-one normalized gap and can be larger on a \(qK\)-shell.
  “Generically” is only a boundary description: (0.1) can vanish for a
  specially engineered Leray kernel or signed multi-channel cancellation.

This leaves the **same-curl reservoir realization gate (SCRG)**:

> Construct the full physical polarization and terminal charge action from
> a same-curl multi-direction shell, with its actual retained wake in the
> C175 admissible class; or prove exact cancellation/slaving of every
> noncommon-curl cross source.  In either case include reservoir depletion,
> repeated charges, localization, and C125/RIGM/BAFL.

The note does not identify the existing prescribed full-polarization
matrices with the same-curl tangent.  It does not prove a \(q\)-way star,
time ordering, nonlinear depletion orbit, localized pump, C125, RIGM,
BAFL, an unforced one-cell stage, or a singularity.

## 1. Exact same-curl unforced embedding

Work on the torus, or on a periodic chart before localization.  For a
nonzero real number \(\kappa\), let

\[
 {\cal E}_\kappa=\{U:\nabla\cdot U=0,\qquad
                         \nabla\times U=\kappa U\}. \tag{1.1}
\]

Every smooth field in \({\cal E}_\kappa\) is supported on the Fourier shell
\(|k|=|\kappa|\) and one helical sign.  The space is linear.  Therefore,
for arbitrary real \(P_0,G_0\in{\cal E}_\kappa\) and scalar
\(\varepsilon\),

\[
 B_\varepsilon(t)
 =e^{-\nu\kappa^2t}(P_0+\varepsilon G_0)             \tag{1.2}
\]

satisfies

\[
 \nabla\times B_\varepsilon=\kappa B_\varepsilon,
 \qquad
 \Delta B_\varepsilon=-\kappa^2B_\varepsilon.       \tag{1.3}
\]

The vector identity

\[
 (B\cdot\nabla)B
 =\nabla{|B|^2\over2}-B\times(\nabla\times B)         \tag{1.4}
\]

then gives

\[
 {\mathbb P}(B_\varepsilon\cdot\nabla B_\varepsilon)=0.
                                                               \tag{1.5}
\]

Hence (1.2), with pressure \(-|B_\varepsilon|^2/2\) in the convention
\(\partial_tB+B\cdot\nabla B+\nabla p=\nu\Delta B\), up to a function of
time, is an exact unforced Navier--Stokes solution.  This is not a Galerkin
assertion: every cross interaction between \(P_0\) and \(G_0\) is absorbed
into the gradient in (1.4).

Let \(u=B_\varepsilon+v\) be any exact solution and subtract (1.2).  With

\[
 L_Hv=-{\mathbb P}\{(H\cdot\nabla)v+(v\cdot\nabla)H\},
 \qquad
 {\cal N}(v,v)=-{\mathbb P}(v\cdot\nabla v),        \tag{1.6}
\]

the exact perturbation equation is

\[
 \boxed{
 \partial_tv=\nu\Delta v
 +e^{-\nu\kappa^2t}L_{P_0}v
 +\varepsilon e^{-\nu\kappa^2t}L_{G_0}v
 +{\cal N}(v,v).}                                  \tag{1.7}
\]

In particular, by smooth-solution uniqueness on the interval in question,

\[
                         v(0)=0\quad\Longrightarrow\quad v(t)=0. \tag{1.8}
\]

There is no term of size \(\varepsilon\) independent of \(v\).  This is
the exact sense in which the preloaded reservoir is an autonomous
homogeneous gate rather than a prescribed writer.

In normalized inertial time \(\tau=a_jK_jt\), write

\[
 \mu_j={\nu K_j\over a_j},
 \qquad |\kappa|=sK_j.                              \tag{1.9}
\]

If the normalized tangent \(L_{G_0}\) has order-one size, the reservoir
action through a window of length \(T\), for \(\varepsilon\geq0\), is

\[
 \Theta(\varepsilon,\mu_js^2,T)
 =\varepsilon\int_0^T e^{-\mu_js^2\tau}\,d\tau
 =\varepsilon T\,{1-e^{-x}\over x},
 \qquad x=\mu_js^2T,                               \tag{1.10}
\]

with the continuous value \(\Theta=\varepsilon T\) at \(x=0\).  For
\(0\le x\le1\),

\[
                 {\varepsilon T\over2}
                 \le\Theta\le\varepsilon T.        \tag{1.11}
\]

Thus \(\varepsilon=b\) gives one \(O(b)\)-action microgate on an
\(O(1)\) window.  If the same static generator acts for
\(J=b^{-1}\) units while heat remains flat, its total action is order one,
not order \(b\).  This can represent repeated factors of one generator,
but not the switched noncommuting polarizations used in C166/C168.

## 2. What improves in the C125 relative ledger

The C147 microseed has growing-coordinate size

\[
                         z_0=\varepsilon_{\rm seed}=n^{-28}. \tag{2.1}
\]

A prescribed same-support writer of absolute size
\(n^{-20}\log n\) produces the fatal relative ratio
\(n^8\log n\).  Equation (1.7) has no such source.

Time in this section is the normalized inertial time of (1.9).  To state
the exact conditional improvement, let \(\Phi(t)\) be the C125
left trajectory for the intended growing block and let \(G(t)\) be its
gain exponent.  The reservoir contribution to the exact C125 identity is

\[
 I_G={1\over|z_0|}
 \int_0^T e^{-G(t)}
 \left|\Phi(t)\!\left[
   \varepsilon e^{-\mu_js^2t}L_{G_0}v(t)\right]\right|dt. \tag{2.2}
\]

Assume, in the actual C125 packet norms, the propagated homogeneous tangent
bound

\[
 e^{-G(t)}\left|\Phi(t)[L_{G_0}v(t)]\right|
 \le C_G|z_0|,
 \qquad 0\le t\le T.                              \tag{2.3}
\]

Then (2.2) gives exactly

\[
                         \boxed{I_G\le C_G\varepsilon T.}  \tag{2.4}
\]

For \(\varepsilon=b=n^{-2}\) and \(T=O(\log n)\),

\[
                         I_G=O(n^{-2}\log n)=o(1). \tag{2.5}
\]

The seed size cancels from (2.2)--(2.4); there is no \(n^{28}\) inverse
penalty.  This is an exact implication from (2.3), not a proof of (2.3).
Same-order off-ladder descendants, moving projectors, localization, and
the C175 Melnikov exposure can violate (2.3).  The reservoir may instead
be included in the intended leading generator, in which case its retained
effect is not a defect, but its complement response and endpoint
decomposition still require C125/RIGM.

## 3. Energy and depletion are not free

If \(P_0\) and \(G_0\) have disjoint Fourier support, or merely satisfy
\(\langle P_0,G_0\rangle_{L^2}=0\), then the background energy is exactly

\[
 \|B_\varepsilon(t)\|_2^2
 =e^{-2\nu\kappa^2t}
   \left(\|P_0\|_2^2+\varepsilon^2\|G_0\|_2^2\right). \tag{3.1}
\]

There is no inviscid pump--reservoir exchange on the orbit \(v=0\).  When
the supports overlap, only the energy of the sum is intrinsic; the
coefficients still keep their fixed ratio under (1.2).

Once a perturbation is present, its exact energy identity from (1.7) is

\[
 {1\over2}{d\over dt}\|v\|_2^2+\nu\|\nabla v\|_2^2
 =-\int v\cdot(v\cdot\nabla)B_\varepsilon\,dx.     \tag{3.2}
\]

Thus a frozen-background linear gain represents work by the pump/reservoir.
The compensating nonlinear deformation is contained in \(v\), including
its components on the original pump and gate modes; it cannot be omitted
from an exact depletion theorem.

The C161 scales make the issue borderline rather than impossible.  If
\(\|v\|\asymp b\), \(\varepsilon=b\), and normalized derivatives are
order one, the gate portion of the work rate has scale

\[
                              b\|v\|^2\asymp b^3.   \tag{3.3}
\]

One \(O(1)\) micro-window therefore charges \(O(b^3)\), matching the
per-step active-return power.  Over \(J=b^{-1}\) windows the absolute work
can reach

\[
                              Jb^3=b^2,             \tag{3.4}
\]

the entire source energy and the energy of one order-\(b\) reservoir.
Consequently the same-curl embedding makes the gate autonomous, but it does
not justify freezing its amplitude through the full order-one action.
Either the effective source--daughter block must be energy-skew, or the
reservoir/pump depletion and all generated modes must enter the leading
trajectory, as in the C118/C124 philosophy.  Equations (3.3)--(3.4) are
scale upper ledgers, not lower bounds on work and not a depletion orbit.

Using \(J\) mutually \(L^2\)-orthogonal reservoirs of norm \(b\) would
instead preload energy \(Jb^2=b\), before cross interactions and
localization.  Without orthogonality the energy of the sum, not the sum of
the named reservoir energies, is intrinsic.  One shared order-\(b\)
reservoir has the correct \(b^2\) energy scale but supplies only its static
generator.  This is why autonomous time ordering, not raw energy, is the
immediate obstruction.

## 4. Exact cross-curl source

Let \(U\in{\cal E}_{\kappa_1}\) and
\(R\in{\cal E}_{\kappa_2}\).  The exact vector identity

\[
 \nabla(U\cdot R)
 =(U\cdot\nabla)R+(R\cdot\nabla)U
   +U\times(\nabla\times R)+R\times(\nabla\times U) \tag{4.1}
\]

gives

\[
 \boxed{
 {\mathbb P}\{(U\cdot\nabla)R+(R\cdot\nabla)U\}
 =(\kappa_1-\kappa_2){\mathbb P}(U\times R).}      \tag{4.2}
\]

For the separately heat-decaying trial background

\[
 \widetilde B(t)
 =e^{-\nu\kappa_1^2t}U_0
 +\varepsilon e^{-\nu\kappa_2^2t}R_0,             \tag{4.3}
\]

the exact projected Navier--Stokes residual is

\[
 \boxed{
 {\cal S}_{\rm cross}(t)
 =\varepsilon(\kappa_1-\kappa_2)
 e^{-\nu(\kappa_1^2+\kappa_2^2)t}
 {\mathbb P}(U_0\times R_0).}                     \tag{4.4}
\]

The correction on the right side of the evolution equation has the
opposite sign.  When \(\kappa_1=\kappa_2\), (4.4) vanishes and recovers
Section 1.  When \(\kappa_1\ne\kappa_2\), its exact size is
\(\varepsilon\) times the curl gap and the Leray-projected product.  For an
order-one normalized gap this is generically order \(\varepsilon\); for a
\(qK\)-scale gap it may carry an additional \(q\).  It vanishes when the
displayed Leray projection does, or may cancel with other channels.
Otherwise it is an additive source even at \(v=0\), so a microseed cannot
remove it by being made smaller.

At one wavevector the two physical helical polarizations belong to the
opposite curl eigenspaces \({\cal E}_{|k|}\) and
\({\cal E}_{-|k|}\), so their curl gap is \(2|k|\).  Whether their actual
cross source is nonzero still depends on
\({\mathbb P}(U\times R)\).  The existing full-polarization prescribed
controls therefore do not automatically lie on the same-curl dark
manifold.  They may still be realized by several same-helicity directions,
or by an exact \({\mathbb P}(U\times R)\) cancellation; neither statement
is proved here.

## 5. Heat cannot supply the clock

For a fixed divergence-free Fourier coefficient, pure heat evolution is

\[
                  \widehat U(t,k)=e^{-\nu|k|^2t}\widehat U(0,k). \tag{5.1}
\]

The multiplier is real and positive.  It changes neither phase nor
polarization, and both helicities at the same \(k\) have the same rate.
For the entire same-curl shell, \(|k|=|\kappa|\), so (1.2) has exactly one
common scalar envelope.  No heat-generated rotation, chirp, or switching
is possible on that exact unforced manifold.

On C127, with \(n=j+1\),

\[
                         \mu_j=\nu((n-1)!)^{-2}.    \tag{5.2}
\]

Let a reservoir have normalized frequency ratio \(s_j\le n^C\) and act
for \(T_j\le n^D\), where \(C,D\) are fixed.  Its complete heat exponent
obeys

\[
 x_j=\mu_js_j^2T_j
 \le {\nu n^{2C+D}\over((n-1)!)^2}.                \tag{5.3}
\]

For every fixed \(M\),

\[
                         n^Mx_j\longrightarrow0.   \tag{5.4}
\]

Thus every polynomial-band reservoir is temporally flat to
superpolynomial accuracy on polynomial windows; logarithmic windows are
included.  To obtain \(x_j\asymp1\) by heat alone requires

\[
                         s_j\asymp(\mu_jT_j)^{-1/2}, \tag{5.5}
\]

which is factorially larger than the current \(q=n^8\) bands when \(T_j\)
is polynomial.  A signed sum of different heat rates might synthesize a
time profile through cancellation, but no bounded-conditioning theorem is
known; it also leaves the single-eigenspace exact-dark setting and incurs
the cross sources in Section 4.

## 6. Consequences for the existing charge geometry

The literal C161 pure-normal shifts are

\[
                         g_a=aN.                   \tag{6.1}
\]

Their radii are \(|a||N|\).  At one nonzero shell radius there are only the
two collinear choices \(\{a,-a\}\), not \(q\) distinct shifts.  Hence the
literal \(q\)-shift bundle cannot be one same-curl Beltrami reservoir.

This does not rule out replacing those shifts by \(q\) distinct integer
vectors on one sphere.  A same-radius, same-helicity palette is exactly the
strongest continuation left by C177.  It must prove that its
Leray-projected tangent spans the two physical source polarizations
uniformly over the C176 \(q^2\) correlated slab, with the collective
\(q^{-1/2}\) normalization and endpoint coherence.  It must separately
control its interaction with the old A2 pump: exact same-curl darkness
requires the palette to share the old pump shell and helicity, while a
\(qK\)-scale shell has a different curl eigenvalue and incurs (4.4).

Moreover, a bank of flat preloaded reservoirs is present simultaneously.
Its sum is one static generator if all components share \(\kappa\); heat
does not produce the noncommuting time ordering which C166 used to escape
the fixed polarization fiber.  Whether spatial separation, Lagrangian
encounters, or a static multi-direction symbol can replace that timing is
an open same-geometry question.

The exact surviving alternatives are:

1. construct a same-curl, multi-direction, static generator whose full
   charge-ladder propagator provides the required physical polarization
   and coherence;
2. generate time ordering through the unforced nonlinear/depletion orbit,
   not through heat;
3. prove a signed kernel/slaving cancellation of (4.4) for a
   noncommon-curl palette; or
4. export the complete cross source with C175/C178-compatible global tails.

Every branch still owes localization, reservoir depletion, repeated
charges, C125, RIGM, BAFL, and the exact endpoint stage.

## 7. Verification boundary

The dependency-free checker verifies:

* exact Fourier Beltrami and projected same-curl darkness on a nontrivial
  reality-complete multi-direction example;
* the cross-curl identity and a nonzero opposite-helicity source;
* orthogonal-background energy arithmetic and the homogeneous/no-writer
  support structure;
* the action, C125-relative, depletion-power, and preload-energy ledgers;
* fixed-mode scalar heat evolution and the factorial-flatness comparison;
  and
* the pure-normal fixed-radius cardinality boundary.

It cannot verify a same-shell full-polarization symbol, nonlinear
depletion, localization, the physical C176 slab realization, C125, RIGM,
BAFL, or an unforced one-cell stage.
