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

## 8. Verification boundary

`checks/one_cell_stage_map_obstruction.py` verifies the exponent,
bootstrap, and summability arithmetic used above.  It cannot verify the
PDE evolution-family estimate (BAFL), which remains the sole named analytic
obstruction in this reduction.
