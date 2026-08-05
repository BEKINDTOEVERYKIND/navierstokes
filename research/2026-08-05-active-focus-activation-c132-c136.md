# Just-in-time active focus: recovered claims C132--C136

**Date:** 2026-08-05

**Status:** exact material-phase, concentration, and schedule algebra;
conditional activation ordering under the displayed raw-feedback model;
exact statement of the one-cell closure target; **OPEN** localized
Navier--Stokes realization.  No singular solution is claimed.

This note stays inside the existing one-cell A2-ladder program.  It does not
introduce a new carrier geometry.  Its purpose is to say exactly what
``activate just in time and focus'' must mean on the factorial schedule of
`2026-08-05-factorial-stage-schedule-c126-c131.md`, and to reduce the missing
PDE input to one named estimate.

Throughout, \(j\ge1\), \(n=j+1\), and

\[
 q_j=n^8,\qquad b_j=n^{-2},\qquad
 F_j=q_j^{3/2}=n^{12},\qquad
 g_j=b_jF_j=n^{10},\qquad M_j=n^{7/2}.                     \tag{0.1}
\]

Here \(q_j\) is the child frequency ratio, \(b_j\) is the dormant seed
amplitude relative to the parent, \(F_j\) is the required active
three-dimensional concentration, and \(M_j\) is the separation/collar
window.  These quantities are not free after (0.1) is chosen.

## 1. C132: material-phase determinant and the \(A=32\) corrector ledger

Let \(X(t,a)\) be the flow map of a smooth incompressible velocity on the
one-cell interval and \(F(t,a)=D_aX(t,a)\).  Define three material phases by

\[
                 \Phi_r(t,X(t,a))=a_r,\qquad r=1,2,3.       \tag{1.1}
\]

Then

\[
 (\nabla_x\Phi_1,\nabla_x\Phi_2,\nabla_x\Phi_3)=F^{-T},
 \qquad
 \det(\nabla_x\Phi_1,\nabla_x\Phi_2,\nabla_x\Phi_3)=1.       \tag{1.2}
\]

Thus smooth incompressible transport cannot make the three phases
algebraically dependent.  It can make their norms and the associated
endpoint chart badly conditioned; (1.2) alone gives no \(j\)-uniform
analytic bound.

The previously selected Gevrey-2 corrector resolution can be frozen at

\[
                         K_j=n^{32}.                         \tag{1.3}
\]

Indeed, if the order-\(m\) coefficients obey the **conditional** majorant

\[
                    \|V_m\|\le C^{m+1}(m!)^2               \tag{1.4}
\]

and

\[
              m_j=\left\lfloor\eta{j^2\over\log n}\right\rfloor,
                                                                    \tag{1.5}
\]

then Stirling's bound gives

\[
 C^{m_j+1}(m_j!)^2K_j^{-m_j}
       \le \exp(-c_\eta j^2)                                \tag{1.6}
\]

for all sufficiently large \(j\).  Any fixed polynomial factors from
\(F_j=n^{12}\), the \(n^2\) chart, or finitely many derivatives are absorbed
by (1.6).  The number 32 is therefore more than sufficient at the scalar
level; it does not prove (1.4) for the localized nonlinear hierarchy.

Equation (1.2) is **PROVED**.  Equation (1.6) is an **EXACT CONDITIONAL
IMPLICATION** from (1.4).

## 2. C133: the activation trichotomy

There are three arithmetically distinct orders of operation.  This section
is a scale gate: it assumes no WKB/Beltrami or propagated-response
cancellation and uses the displayed carrier-derivative and envelope-divisor
rules.  It does not derive those PDE estimates from the Navier--Stokes
evolution.  Conditional on that raw-feedback model, the comparisons below
are exact.

### 2.1 Route weakly at parent frequency, then activate

For a route lasting the logarithmic gain time \(O(\log q_j)\), the worst
uncancelled parent-carrier feedback relative to the parent is

\[
              \varepsilon_j^{\rm parent}
              \lesssim b_j\log q_j
              ={8\log n\over n^2}.                          \tag{2.1}
\]

This is summable.  An envelope term carrying the separation divisor is even
smaller:

\[
              \varepsilon_j^{\rm env}
              \lesssim {b_j\log q_j\over M_j}
              ={8\log n\over n^{11/2}}.                    \tag{2.2}
\]

Thus this branch has a summable raw scale and is the only candidate among
the three branches below.  Summability alone is not the stage-map tolerance:
if (2.1) returns to an active chart, it is much larger than the \(n^{-6}\)
pre-chart BAFL allowance.  The one-cell theorem must therefore place this
feedback in the intended dynamics, cancel it, or export it to a retained
wake with the stronger bound in (3.6).

### 2.2 Route the weak seed after moving it to child frequency

Without a WKB/Beltrami cancellation, one carrier derivative changes the
relative feedback to

\[
               b_jq_j\log q_j=8n^6\log n,                   \tag{2.3}
\]

which grows rather than tends to zero.

### 2.3 Amplify before routing

Even if an envelope divisor \(M_j^{-1}\) is retained, pre-amplification gives

\[
          {F_jb_j\log q_j\over M_j}
          =8n^{13/2}\log n,                                 \tag{2.4}
\]

again divergent.

Therefore the only order not already divergent in the uncancelled scalar
ledger is

\[
 \boxed{\text{weak parent-frequency route }\longrightarrow
        \text{ destination activation }\longrightarrow
        \text{ active focus}.}                              \tag{2.5}
\]

This proves a trichotomy of scale factors within the displayed model.  It
does **not** prove that the first branch is realized by an exact solution,
nor does it rule out (2.3) or (2.4) when an exact structural cancellation is
present.

## 3. C134: a finite Fourier circuit does not close the PDE

Use the equal-shell integer A2 roots already fixed by C114:

\[
 \mathcal H=\{\pm k_1,\pm k_2,\pm k_3\},\qquad
 k_1=(1,-1,0),\quad k_2=(0,1,-1),\quad k_3=(-1,0,1).        \tag{3.1}
\]

They satisfy \(k_1+k_2+k_3=0\) and \(|k_i|^2=2\).  The hexagon
contains its intended resonant triads, but it is not closed under all
quadratic sums.  For example, the difference shell contains

\[
                        k_1-k_2=(1,-2,1)\notin\mathcal H,
             \qquad |k_1-k_2|^2=6.                          \tag{3.2}
\]

For a Fourier solution the exact complement equation contains

\[
 \partial_t\widehat u(k)+\nu|k|^2\widehat u(k)
 =-iP_k\sum_{p+q=k}(\widehat u(p)\cdot q)\widehat u(q).     \tag{3.3}
\]

Consequently support combinatorics alone never makes the A2 ladder an exact
Navier--Stokes subsystem.  Some coefficients in (3.3) may vanish because of
polarization, Beltrami structure, or pairwise cancellation; those
cancellations must be proved for the transported, localized packet rather
than inferred from the finite circuit.

Let \(P_j(t)\) project onto the intended pump/ladder/seed channels and let
\(Q_j(t)=I-P_j(t)\) be the retained complement/wake state.  For a
divergence-free candidate trajectory \(v_j\), form its full unforced
residual

\[
 \mathcal S_j=\partial_tv_j+\mathbb P\nabla\!\cdot(v_j\otimes v_j)
                    -\nu\Delta v_j.                         \tag{3.4}
\]

After removing only terms already covered by a proved or explicitly
conditional paid ledger, denote the remainder by
\(\mathcal S_j^{\rm leak}\).  If the projections move, their time and
spatial commutators are included in this source.  In particular, the
difference-shell part of
\(Q_j\mathbb P\nabla\!\cdot(v_j\otimes v_j)\) cannot be omitted.

Let \(\mathcal U_j(t,s)\) be the full linearized evolution family about
\(v_j\), including the remaining focus gain but excluding the final
finite-dimensional chart inverse, and set

\[
 w_j(t)=\int_{t_j^-}^{t}\mathcal U_j(t,s)
                    \mathcal S_j^{\rm leak}(s)\,ds,\qquad
 \mathfrak L_j^{\rm act}=\sup_t\|P_j(t)w_j(t)\|_{X_j},\qquad
 \mathfrak L_j^{\rm wake}=\sup_t\|Q_j(t)w_j(t)\|_{X_j}.     \tag{3.5}
\]

This response definition retains cancellations inside Duhamel's integral.
An integral of propagator norms against \(\|\mathcal S_j^{\rm leak}\|\)
is a sufficient upper bound, but is not BAFL itself.  The missing
two-channel estimate is

\[
 \boxed{\mathfrak L_j^{\rm act}\le Cn^{-6},\qquad
        \mathfrak L_j^{\rm wake}\le Cn^{-4},\qquad
        L_j\mathfrak L_j^{\rm act}\le Cn^{-4}.}              \tag{3.6}
\]

We call (3.6), together with construction of the trajectory on which it is
measured, the **backward-weighted active-focus leakage estimate (BAFL)**.
The \(n^{-4}\) allowance belongs to a complement component which remains in
the retained wake.  It is not an allowed pre-chart active error.  Because
\(L_j\le n^2\), a component which returns to the active chart must gain two
extra powers, \(n^{-6}\), before chart inversion.  The backward weight is
essential: an unweighted small defect made before focusing can be multiplied
by as much as \(F_j=n^{12}\) before the endpoint.  An ordinary forward
energy estimate does not imply (3.6).

The non-closure (3.2)--(3.3) is **PROVED**.  BAFL is **OPEN** and is the named
PDE obstruction retained by this reconstruction.

## 4. C135: exact active-focus algebra

The parent cell volume is \(V_j\asymp\ell_j^3\), while the child cell volume
is

\[
                         V_{j+1}=q_j^{-3}V_j.                \tag{4.1}
\]

A dormant seed of amplitude \(b_ja_j\) has energy scale

\[
                         E_j^{\rm seed}=b_j^2E_j.            \tag{4.2}
\]

If an active, energy-transferring focus places that same leading energy in
the child volume, its amplitude must be

\[
 a_j^{\rm child}=b_ja_j\left({V_j\over V_{j+1}}\right)^{1/2}
                    =b_jq_j^{3/2}a_j
                    =n^{10}a_j=a_{j+1}.                     \tag{4.3}
\]

The child energy ratio is correspondingly

\[
 { (a_j^{\rm child})^2V_{j+1}\over a_j^2V_j}
        =b_j^2=n^{-4}={E_{j+1}\over E_j}.                   \tag{4.4}
\]

Thus the apparently special \(R^{3/2}\) target is forced by
three-dimensional volume scaling, with \(R=q_j\).  It is exactly compatible
with the factorial schedule.

The allowed exit-chart inverse loss also fits exactly:

\[
                         L_j\le n^2
              \quad\Longrightarrow\quad L_jb_j\le1.        \tag{4.5}
\]

Finally the stage target retains, rather than erases, a complement wake
increment with

\[
                    \|w_{j+1}^{\rm new}\|_{X_{j+1}}
                         \le Cn^{-4}.                        \tag{4.6}
\]

The series in (4.6) converges.  This is distinct from a mode which returns
to an active chart coordinate: that mode must be \(O(n^{-6})\) before the
\(O(n^2)\) chart inverse, as in (3.6).  Which norm \(X_j\) is stable under
scale change and nonlinear wake coupling is part of BAFL.

Equations (4.1)--(4.5) are **PROVED SCALE ALGEBRA**.  The active focus in
(4.3) and the wake estimate (4.6) are **CONDITIONAL TARGETS**.  Passive
incompressible material transport cannot prove (4.3), because it preserves
the volume of every material set.  The focus must be an active nonlinear
energy transfer from the decaying pump/reservoir.

## 5. C136: the one-cell closure theorem, stated without hiding the gap

Spatial gates can remove some local products exactly.  If two smooth fields
have disjoint supports at a fixed time, their cross tensors vanish
pointwise.  They do not become independent Navier--Stokes solutions:
the heat flow immediately develops tails, and the Leray kernel for
\(\mathbb P\nabla\cdot\) has size \(O(d^{-4})\) at separation \(d\).

Assume the one-source/one-target first-Picard estimate before chart inversion
is \(C\log q_j/M_j\).  The chart-weighted tail is then

\[
 L_j{\log q_j\over M_j}
 \le {8\log n\over n^{3/2}},                                \tag{5.1}
\]

which is summable.  This proves that the chosen separation is large enough
for that single interaction to lie inside the conditional exit chart with a
summable right-inverse control cost.  The chart must cancel this paid term;
the quantity in (5.1) is not admitted as an \(n^{-4}\) terminal active error.
It does not control the full Duhamel series or leakage created by the
correcting control.

The missing theorem can now be stated in one line.

> **One-cell BAFL closure theorem (OPEN).**  For every sufficiently large
> \(j\), construct an exact unforced Navier--Stokes evolution on the existing
> A2 one-cell geometry which starts from the \(j\)-th decaying pump, its
> dormant parent-frequency seed, and the admitted carried wake; realizes the
> C118/C119 depletion and C125 gain blocks; routes the weak seed before
> activation; transfers its energy into the \(q_j^{-3}\) child volume with
> the \(q_j^{3/2}\) amplitude concentration (4.3); has exit-chart inverse
> loss at most \(n^2\); retains a dark/exported complement wake bounded by
> \(Cn^{-4}\); and bounds every complement component returning to the active
> chart by \(Cn^{-6}\) before that chart inverse, as in (3.6).

This is one theorem about the full PDE trajectory, not a list of independent
controls that may be prescribed afterward.  If it holds with constants
uniform in \(j\), and its admitted-wake and entrance-chart classes are
closed under the stated endpoint errors, the factorial identities of C127
and the summable tails (4.6), (5.1), and (3.6) allow iteration.  Closure of
those state classes is part of the theorem, not a consequence of scalar
summability alone.  If it fails, the failure must occur
in BAFL: the active focus cannot be generated without an overweight
off-ladder mode, or its retained wake cannot meet the same backward-weighted
bound.

At present BAFL is **OPEN**.  The just-in-time ordering is selected by exact
arithmetic within the displayed raw-feedback model, not constructed; the
\(R^{3/2}\) focus is the exact same-leading-energy handoff scaling, not an
existence result; and only a one-source/one-target pressure tail has been
closed.  Mode leakage under the remaining focus is the single named
obstruction to the one-cell stage map on this skeleton.

## Claim boundary

* C132: material determinant exact; Gevrey residual conditional on its
  displayed coefficient bound.
* C133: exact scalar comparison conditional on the displayed raw-feedback
  rules; realization open.
* C134: finite-circuit non-closure exact; BAFL open.
* C135: \(q^{3/2}\), energy ratio, and chart compatibility exact; active
  transfer and wake norm open.
* C136: first separated tail arithmetic exact; full one-cell BAFL closure
  open.

The checker verifies these algebraic statements and deliberately exits with
an `OPEN` status line for the PDE theorem.
