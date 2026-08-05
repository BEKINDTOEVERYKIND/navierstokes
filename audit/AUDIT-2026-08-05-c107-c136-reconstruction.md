# Adversarial audit of the reconstructed C107--C136 one-cell checkpoint

Date: 2026-08-05  
Status: **independently re-derived within this repository; mixed
EXACT/CONDITIONAL/OPEN claim boundary; no stage map or blow-up theorem**

## 1. Provenance

The registry on the surviving branch reserved C107--C136 and named their
topics, but the corresponding in-session files were no longer reachable
from any local or remote Git ref. The six research notes and six checkers
audited here are therefore conservative reconstructions from:

- the surviving C50--C106 and C137--C139 artifacts;
- the auditor's transported-writer and wake-slaving notes;
- the topic ledger preserved in CLAIMS.md; and
- fresh exact Fourier, ODE, and scale derivations.

They are not claimed to be byte-for-byte recovery of the lost drafts. This
provenance distinction is now explicit in every registry-level statement.

## 2. Independent audit assignments

Three separate adversarial passes re-derived disjoint claim blocks, followed
by a fourth pass over the one-cell reduction and registry.

| Block | Audit result | Material corrections |
|---|---|---|
| C107--C117 | Clebsch--Piola, Weber, homochiral darkness, all cyclic helical interactions, exhaustive real four-mode convolution, and the second-Picard outputs were rechecked exactly. | Directional wake protection was separated from genuine profile slaving; the zero-mode hypothesis was added to the single-shell equivalence; the negative result was narrowed to the bare paired gate, not a larger ladder. |
| C118--C125 | Energy conservation, depletion heteroclinics, the weighted characteristic polynomial and rational enclosure, Beltrami pump, decaying-pump gain, and seed bootstrap were re-derived. | C124 was narrowed to a finite normal form. C125 was repaired from a linearized identity to the exact perturbation identity including \(P_H\mathcal N(v,v)\), and now concludes only retained-amplitude gain under its displayed bound. Support arithmetic is not used to assert a nonzero physical coefficient. |
| C126--C136 | Factorial exponents, logarithmic gain-time costs, localization implications, chart scaling, Gevrey remainder, activation ordering, focus scaling, and BAFL budgets were rechecked. | The logarithmic window was charged to time/heat/dissipation; the chart radius constant was exposed; C133 was made conditional on its raw-feedback model; the pressure tail was classified as a control cost; moving-projection/full-residual terms were retained; iteration was conditioned on closure of entrance and wake state classes. |
| One-cell reduction/registry | The general \(A_2\) polarization calculation, two-channel BAFL norm, fixed-point bootstrap, and C107--C136 status table were cross-checked against all six notes. | BAFL was split into pre-chart active \(O(n^{-6})\) and retained-wake \(O(n^{-4})\) channels; the contraction uses the paid residual through \(KZ_j\le1/8\); the derivative argument was limited to robust Lipschitz perturbative closures; all exact/conditional labels and equation references were synchronized. |

## 3. Exact positive layer

The audit accepts the following as exact under the hypotheses stated in the
individual notes:

1. The integer roots
   \[
   k_1=(1,-1,0),\quad k_2=(0,1,-1),\quad k_3=(-1,0,1)
   \]
   form an equal-shell \(A_2\) hexagon.
2. The dual-helicity C116 parent pair cancels its first reality-difference
   output while retaining a nonzero terminal child.
3. That child necessarily creates at least one of
   \[
   e_1=(2,-1,-1),\qquad e_2=(1,1,-2)
   \]
   in the next Picard interaction of the bare gate.
4. The abstract pump--leaf normal form has the C118/C119 depletion
   heteroclinics.
5. The weighted ladder polynomial and enclosure are
   \[
   \lambda^6-9\lambda^4+18\lambda^2-9,\qquad
   {633\over250}<\sigma_*<{2533\over1000}.
   \]
6. The homochiral \(A_2\) pump is an exact unforced Navier--Stokes
   heat-decaying background, while its finite retained mode has ideal gain
   \[
   G(t_*)={d\over\delta}(R-1-\log R).
   \]
7. The factorial schedule has the exact handoff algebra
   \[
   \ell_j=(j!)^{-8},\quad a_j=(j!)^{10},\quad
   q_j=n^8,\quad b_j=n^{-2},\quad F_j=n^{12},\quad
   b_jF_j=n^{10},\quad n=j+1,
   \]
   and all explicitly recorded scalar budgets, including logarithmic
   gain-time factors, are summable.

None of items 2--7 identifies an invariant localized Navier--Stokes ladder.

## 4. The exact named obstruction

For the fixed nondegenerate paired gate, stopping when the child reaches
\(b_j=n^{-2}\) gives a second-Picard output with nonzero leading coefficient
and natural size

\[
                     \Theta(b_j^2)=\Theta(n^{-4}).
\]

This matches the admitted retained-wake budget. It is too large to return
through a worst-case \(O(n^2)\) active chart. The single remaining package
is therefore **backward-weighted active-focus leakage (BAFL)**:

\[
 \boxed{
 \mathfrak L_j^{\rm act}\le Cn^{-6},\qquad
 \mathfrak L_j^{\rm wake}\le Cn^{-4},\qquad
 L_j\mathfrak L_j^{\rm act}\le Cn^{-4}.}
\]

The response is the full cancellation-sensitive Duhamel response about the
localized candidate, with all remaining focus gain, returned complement
modes, nonlinear wake feedback, and moving-projection effects included.
Construction of the trajectory on which this response is measured is part
of the BAFL package.

This unifies the C107 fork with C117: the named second-Picard mode must
remain on a dark/slaved or spatially exported wake branch, or its propagated
active return must gain an additional factor \(n^{-2}\).

## 5. What remains open

The audits do **not** establish:

- physical polarization coefficients realizing the abstract weighted
  hexagon and depletion block;
- a localized invariant/slaved ladder;
- the just-in-time routing trajectory;
- active \(q_j^{3/2}\) energy concentration;
- a uniform endpoint inverse with the required radius constant;
- BAFL for the full nonlinear Leray/heat evolution;
- closure of the admitted entrance and wake classes under iteration; or
- any Navier--Stokes stage map, singularity, or Millennium result.

Accordingly C136 is an OPEN theorem target, not a theorem.

## 6. Verification

The following dependency-free checks are the executable boundary of this
checkpoint:

- checks/transported_writer_reconcile_c107_c113.py
- checks/terminal_triad_hexagon_c114_c117.py
- checks/hexagon_depletion_eigenmode_c118_c120.py
- checks/unforced_decaying_pump_c121_c125.py
- checks/factorial_stage_schedule_c126_c131.py
- checks/active_focus_activation_c132_c136.py
- checks/one_cell_stage_map_obstruction.py

Each passes, each compiles, and repository-wide checker execution,
audit/independent_recheck.py, whitespace checking, and Markdown delimiter
scans pass at this checkpoint. A passing checker certifies only its stated
finite algebra or arithmetic; it does not certify BAFL.
