# Navier--Stokes frontier

**Branch:** `agent/aug2-integrated-transition`
**Registry frontier:** C191 (2026-08-27 checkpoint; C189 is the ingested
auditor cross-audit)
**Boot rule:** every research session starts from this file, `CLAIMS.md`, and
the referenced artifacts on the branch. Chat history is not a premise.
Every checkpoint commit must update this file in the same commit.

## Current architecture target

Define the exact unforced nonautonomous stage-renormalization skew-product

\[
 (X,\mu,n)\longmapsto
 \left(\mathcal R_{n,\mu}X,{q_n\over g_n}\mu,n+1\right),
 \qquad
 \mathcal R_{n,\mu}=\mathcal C_{{\rm exit},n}\circ
 S^{NS,\mu}_{T_n},
\]

on one complete structured state space containing the active packet, carried
wake, phase/modulation variables, the normalized-viscosity coordinate, and
the exit chart.  Construct a nonautonomous invariant trapping tube with one
expanding direction.  An autonomous inviscid-face fixed point is a further
conditional target only after proving a normalized limit
\(\mathcal R_{n,0}\to\mathcal R_{\infty,0}\).  C188 proves that a fixed
supercritical scaling \(g>q\) cannot have an autonomous fixed point at
positive normalized viscosity because \(\mu'=(q/g)\mu<\mu\); the actual
factorial map has no finite-\(n\) fixed point.
Leak, pressure, and wake channels are components of the single residual
\(\mathcal R_{n,\mu}(X_n)-X_{n+1}\); they are not independent modules to be
closed one at a time.

The unmatched leading object is the **unforced viscous physical-velocity
scale return (UVSR)**:

> a nonzero smooth boundaryless three-dimensional unforced Navier--Stokes
> orbit/profile whose renormalized endpoint is a smaller-copy structured
> state and whose physical velocity normalization realizes a same-energy
> focus \(F\) and net endpoint gain \(g=bF\), with the admitted wake
> included in the state.  In the equal-shape ledger
> \(F=q^{3/2}\) and \(g=q^\gamma\); C188 records the explicit correction
> for unequal normalized \(L^2\) shapes.

C188 replaces the hard-coded C127 net exponent \(\gamma=5/4\).  If the
bounded-profile C176 worst-case upper envelope
\(C(1+\Lambda)q^{-1}J^{7/2}/b^3\to0\) remains in the specification, its
exact polynomial window is \(7/6<\gamma<3/2\), with a sufficient
slowly-diverging logarithmic correction at the lower boundary.  This is
not a physical necessity theorem because C176 supplies no matching lower
bound.  Here and below \(C_{\rm col}\) absorbs the declared uniform bound
on \(1+\Lambda_j\).  These intervals describe normalized pure powers.  The
exact equal-shape all-sequence demand specifications are
\[
 g=q\rho,\quad 1<\rho=bq^{1/2}<q^{1/2},\quad
 {\log\rho\over\log q}\longrightarrow0
 \qquad\hbox{(direct class)},
\]
and, with \(K_*=C_{\rm col}C_J^{7/2}\),
\[
 g=K_*^{1/3}q^{7/6}(\log q)^{7/6}\omega,\quad
 \omega\longrightarrow\infty,\quad
 {\log\omega\over\log q}\longrightarrow0,\quad b<1
 \qquad\hbox{(declared C176 envelope)}.
\]
Both sequence specifications must also satisfy the exact global energy
condition (2.0a) in C188.  The explicit boundary witnesses there do.
In particular \(q=n^4\), \(g=2q\), \(b=2n^{-2}\) doubles Reynolds at
every return while attaining direct polynomial order one; it is a scalar
witness and does not reuse C180's \(q=n^8\) shell theorem.
One convenient lower-demand scalar schedule which keeps C180's
proved \(q=n^8\) shell is
\[
 b=n^{-5/2},\qquad g=n^{19/2}=q^{19/16}.
\]
It requires C161's explicit respecification
\(H=n^{51/2}\), \(\varepsilon=n^{-28}\), and
\(J_{\rm split}=\lceil n^{5/2}\rceil\), as recorded and checked in C188.
If a direct full-residual UVSR certificate cancels or absorbs the separate
C176 collar term, the broader normalized pure-power window is
\(1<\gamma<3/2\), with the direct all-sequence boundary stated above.
Changing the equal-shape \(q^{3/2}\) power would require nonuniform profile
constants or an explicit energy-disposal/wake ledger.

The fixed-point theorem itself is not the novel component. Once a UVSR
approximate profile with a certified full residual exists, the intended
imports are ABC's full-operator spectral splitting and Duhamel contraction,
Elgindi's symmetry-adapted coercivity and compactness, and Chen--Hou's
coercive-bulk/validated-finite-rank trapping method. The detailed primary-
source map is recorded in
`research/2026-08-23-fixed-point-literature-map-c185.md`.

## C188 pre-registration and verdict

The following outcome rule was fixed before the C188 checker was run.

1. Optimize the net demand over a declared schedule class. If no schedule
   meets all imported scalar constraints, issue an architecture-level
   no-go and stop UVSR. If the class is nonempty, its exact exponent window
   becomes the UVSR specification.
2. Test that window against CKN/Lin epsilon regularity, the local energy
   inequality, and Tao's quantitative critical-norm rate. If those theorems
   make the corridor empty under their actual hypotheses, issue the
   architecture-level no-go. Otherwise retain the exact surviving window.

**Verdict:** nonempty at the landed scalar-ledger level.  The
bounded-profile legacy-collar worst-case-envelope window is
\((7/6,3/2)\); the direct-full-residual window is \((1,3/2)\).  These are
the normalized pure-power windows; the exact decorated boundary
specifications are displayed above.  No further
numerical power restriction follows from CKN/Tao without the open
singular-center, time-occupancy, pressure, and wake hypotheses.  The
proposed inference that every smooth intermediate stage center must remain
above one epsilon threshold is false: C188 gives explicit \(r^2,r^3,r^4\)
upper bounds there.  This does not control backward cylinders centered at
an eventual singular point.  No criterion was weakened after seeing the
result.

## C190 pre-registered PPRG verdict

The C186 rotating-gradient orbit has now completed outcome (b) of the
unchanged witness test.  Along its fixed origin, any two consecutive maps
between three distinct endpoint times can both have determinant one only
when the conserved vertical charge is \(m=0\); this includes the two
quarter episodes.  A genuine \(2\pi\) coefficient return also forces
\(m=0\).  For the two equal quarter episodes, in the canonical common
co-rotating orthonormal fiber frame,

\[
 \Phi_1=\Phi_2=I+{\pi\over2}\beta E_{12},\quad |\beta|\le1,\qquad
 (\Phi_i-I)^2=0,qquad
 \operatorname {tr}(\Phi_2\Phi_1)=2.
\]

Thus the ordered \(1/100\) boxes about \(I+E_{12}\) and \(I+E_{21}\)
cannot both hold, and the direct \(2+\delta\) branch fails for every
\(\delta>0\).  The genuine full-return block
\(S=I+2\pi\beta E_{12}\) obeys the explicit bound

\[
                         \|S^N\|_2\le1+{44\over7}N.
\]

This is a theorem-grade obstruction for the chosen orbit, not a failed
float search and not a no-go for the complete admissible Kelvin class.
Accordingly it feeds but does not fire the architecture trigger.  No
criterion was weakened and no successor gate was created.

## C191 C185-deficit verdict and class-scope correction

C185's allowed operator-norm floor has enough abstract scalar exponent for
the raw C182 power.  On C188's \(q=n^8\) schedule,

\[
 R_\Delta=\left\lceil{15\over8}\log q\right\rceil
 \quad\Longrightarrow\quad
 \|G_{R_\Delta T}\|_{2\to2}\ge q^{3/8},
\]

and the required inertial pump action is enclosed explicitly by

\[
 {45\over8}\log q<TR_\Delta
 <{57\over10}\log q+{76\over25}.
\]

If the same floor is also assigned to C188's formal preparation factor,
\[
R_*=\left\lceil{285\over16}\log q\right\rceil
\quad\Longrightarrow\quad
\|G_{R_*T}\|_{2\to2}\ge Hq^{3/8}.
\]
Conditionally on C176/C188's declared collar constants, this longer clock
has nonviscous normalized collar at most
\(C_{\rm col}856^{7/2}n^{-1/4}\) and viscous-collar factor at most
\(C_{\nu{\rm col}}\nu856^{3/2}(j!)^{-3/2}n^{437/28}\).  The explicit
C191 heat condition retains at least \(99/100\) of the pump and charges at
most a \(100/99\) factor in normalized time; on C188's physical scale,
\[
(j!)^{-35/2}TR_*\le t_*\le{100\over99}(j!)^{-35/2}TR_*.
\]

This does not close the physical endpoint.  The corpus supplies no lower
stage-action theorem guaranteeing those returns.  C182 is an
\(L^\infty\) upper bound on one packet class, while C185 is an unrestricted
supremal \(L^2\) operator lower bound with no common finite-\(q\) retained
band.  Treating \(q^{3/8}\) as a common scalar multiplier overruns child
energy by \(q^{3/4}=n^6\); dividing the entrance seed by \(q^{3/8}\) restores
the final energy and cancels the point multiplier.  Under uniform
rescaling, C147's coherent-\(q^3\)-packet writer diagnostic remains
divergent, \((57/2)n^8\log n\), but its transfer to C188/C161's displayed
\(q^2\)-source coefficient is not proved.  The exact reservoir
specification is therefore coherent fixed-final-energy
\(L^\infty/L^2\) concentration by \(q^{3/8}\) on one retained-band witness,
with window, viscosity, active-retention, depletion, and wake losses
charged explicitly.  Bare \(L^2\) growth receives zero focus credit.

C191 also corrects the proposed class dichotomy.  C152's background is
itself passive 2D3C:
\[
U=N\times\nabla f-\sqrt2fN=v_h+\Theta n,\qquad
n=N/\sqrt3,\quad\Theta=-\sqrt6f,\quad v_h\cdot\nabla\Theta=0.
\]
C159/C183/C184/C185 therefore already supplies an \(m\ne0\) exact return
and growth inside the broad C179/C183 class.  Universal secular no-return
is false on that class.  Calling this member the active pump rather than an
auxiliary reservoir is not a mathematical class exclusion, and no narrower
class is invented here.  The two requested forms are accepted program
verdict forms, not an exhaustive theorem about the landed class.  Form (i)
is unavailable; no third form is introduced.

## Current positive target

1. **Final accepted PPRG witness form.** Produce one non-fixed-point base
   trajectory with genuinely incommensurate frequencies and either a
   certified \(m\ne0\) covector return or certified Lyapunov growth in the
   physically transported frame.  The same witness must retain the stage
   band and prove coherent fixed-final-energy concentration by at least
   \(q^{3/8}\), with explicit window, viscosity, finite-frequency,
   active-retention, depletion, and wake constants.  A failed search is
   reported as a failed search, not converted into a class no-go or a new
   gate.
2. **Same-witness viscous conversion.** C185's inviscid operator-norm floor
   may be used only in its C189-approved form.  The existing
   \(\mathfrak d,\mathfrak b,\mathfrak n\ge99/100\) and
   \(\varepsilon_{\rm ff}\le1/100\) bookkeeping remains a necessary
   three-loss audit, but C191 shows it is not sufficient unless the same
   retained-band vector also proves the normalized concentration endpoint.
   Its operator-norm-only floor would be
   \[
    {6\over5}\left({99\over100}\right)^3-{1\over100}
    ={2885897\over2500000}>{23\over20},
   \]
   which is not a concentration estimate.
   No \(L^2\) exponent may be subtracted from the \(q^{3/8}\) focus demand
   before that common-witness theorem.
3. **UVSR profile search.** Only after items 1--2, specify the complete
   Banach state, modulation slice, and viscosity coordinate for
   \({\cal R}_{n,\mu}\), then search for one approximate nonautonomous
   trapping orbit by minimizing the norm of the full stage residual (not
   selected leak projections).  An inviscid-face fixed-point search is
   permitted only after certifying a normalized autonomous limit.  Use C188's
   \(7/6<\gamma<3/2\) bounded-profile legacy-envelope specification unless
   the full residual itself certifies removal of that collar comparison.

## Proved frontier

- C121: the homochiral \(A_2\) pump is an exact unforced heat-decaying
  Navier--Stokes background.
- C152/C159: one exact zero-drift periodic orbit has a rigorously certified
  Kelvin cone with one-period multiplier strictly larger than
  \(e^{1/5}\); before this checkpoint it was registered only as a
  finite-dimensional principal-cocycle statement.
- C161/C176/C179/C180: the required support cardinality, correlated tube,
  exact passive 2D3C background, and factorial-shell tight-star arithmetic
  exist, but no physical velocity endpoint follows from them.
- C181: a fixed-cone static vertical shear has an explicit uniform bound.
- C182: under its stated entrance regularity hypotheses, fixed-time endpoint
  growth is at most \(bq^{9/8}\), short of \(bq^{3/2}\).
- C183: the exact passive 2D3C Lagrangian gauge gives the return-resonance and
  linear/quadratic covector-drift formulas; a common-frame
  \(O(q^{-1/2})\) perturbation gives no extra \(\sqrt q\).
- C184: degree-\(q\) stationary polynomial palettes pay an exponential
  off-line cost and cannot meet a polynomial reservoir budget.
- C185/C189: C159 plus the abstract-level Shvydkoy spectral inclusion gives
  the robust infinite-dimensional operator-norm estimate
  \(\|G_{nT}\|_{L^2\to L^2}\ge e^{n/5}>(6/5)^n\).  The
  essential-spectral-radius form is citation-held until the paper-body
  Theorem 4.1 check lands.  C189 independently confirms C186 and C187 in
  full and C185 with precisely this split.
- C186: conditionally, candidate unipotent polarization blocks do not
  algebraically force a common flag. The exact pair
  \(I+E_{12},I+E_{21}\) has per-episode exponent \(>12/25\), robustly
  \(>9/20\) in its entrywise \(1/100\) boxes. For separately established
  square-zero generators, the exact decision scalar is
  \(\tau=\operatorname{tr}(N_1N_2)\). A bounded exact passive 2D3C flow can
  rotate its scalar gradient through transverse lines, but this does not
  determine the Kelvin polarization flags.
- C187: ordinary finite-stage linearized Navier--Stokes Duhamel continuity
  holds in \(H^3\) with the explicit constant
  \(44\sqrt{2T/\nu}\exp(7744V^2T/\nu)\). It is not scale-uniform.
- C188: equal-normalized-shape same-energy transfer into a \(q^{-3}\) child
  has ledger focus \(q^{3/2}\); unequal shapes contribute the explicit
  square root of their \(L^2\)-constant ratio.  The net gain has sharp
  polynomial infimum \(7/6\) within the declared bounded-profile C176
  worst-case envelope, and infimum \(1\) for a direct full-residual UVSR.
  The direct floor has the exact scalar witness
  \(q=n^4,g=2q,b=2n^{-2}\), with total energy at most \(16/3\) in its
  displayed normalization.
  The exact fixed-\(q=n^8\) scalar schedule \(b=n^{-5/2}\) lowers \(5/4\)
  to \(19/16\), with C161 seed/split normalization respecified explicitly.
  At a smooth stage center, the standard scaled energy/cubic/enstrophy
  quantities are bounded by explicit constants times \(r^2,r^3,r^4\); no
  numerical CKN/Tao power follows from the scalar ledger alone.  Also
  \(\mu'=(q/g)\mu\), and the factorial schedule is an augmented
  nonautonomous map, not an autonomous fixed point.
- C190: on the explicit C186 rotating-gradient orbit, the determinant-one
  gate for any two consecutive episodes forces \(m=0\).  In the canonical
  common co-rotating frame both maps equal
  \(I+(\pi/2)\beta E_{12}\), \(|\beta|\le1\); their product trace is
  exactly two, their
  square-zero generators share a flag, and the full-return powers satisfy
  \(\|S^N\|_2\le1+(44/7)N\).  This is the pre-registered orbit-specific
  obstruction, not the complete-class PPRG verdict.
- C191: the C185 floor supplies the raw scalar \(q^{3/8}\) after
  \(R_\Delta=\left\lceil(15/8)\log q\right\rceil\) returns and is compatible with the
  declared scalar collar/heat powers after an explicit logarithmic-window
  allocation.  No actual lower stage coverage is proved.  The landed
  unrestricted \(L^2\) operator norm cannot be multiplied by C182's
  \(L^\infty\) upper bound; common scalar reuse overruns energy by
  \(q^{3/4}\), fixed-energy normalization cancels the gain, and C185
  supplies no C125 relative-return cancellation.  Also, the
  C159/C185 returning growing covector lies in the broad passive-2D3C
  class, falsifying universal \(m\ne0\) secular no-return on that class.

## Open

- UVSR: no exact or rigorously trapped unforced viscous scale-return profile
  with physical velocity focus is known.
- PPRG realization: the broad landed passive-2D3C class contains both
  C190's secular orbit and C159/C185's exact returning growing orbit, so a
  universal secular-lock theorem is false there.  The sole unresolved
  accepted witness form is one genuinely incommensurate non-fixed orbit
  with an \(m\ne0\) return or transported-frame Lyapunov growth.  It must
  carry the same-witness retained-band, viscosity, fixed-final-energy
  \(q^{3/8}\) concentration, C125, depletion, and wake estimates.
- A complete state space and normalization for \(\mathcal R\), including the
  retained wake and normalized-viscosity coordinate, have not been fixed.
- No normalized convergence \(\mathcal R_{n,0}\to\mathcal R_{\infty,0}\)
  has been proved, so there is not yet an autonomous inviscid-face fixed-
  point problem for the actual factorial schedule.
- No coercive bulk estimate or spectral gap has yet been proved for the
  proposed renormalized Navier--Stokes operator.
- **UVSR standing terminal obligation:** terminal singular-center tracking
  has not been proved.  C188's CKN/Lin/Tao clearance is confined to the
  scalar ledger and smooth intermediate stage centers.  Any terminal UVSR
  object must still control backward cylinders at the proposed singular
  center, including time occupancy, pressure/wake contributions, the local
  energy inequality, and the relevant critical-norm comparison.  No
  residual-stable forced epsilon criterion or torus transfer of Tao's
  quantitative theorem has been proved.
- C125's relative return, the old BAFL split, and RIGM remain unproved, but
  are to be absorbed into one full residual/trapping inequality rather than
  promoted into further architecture gates.

## Pre-registered kill criteria

These verdicts are fixed before the next computation.

1. **PPRG architecture-trigger scope.** The proposed universal
   \(m\ne0\) secular-lock branch is unavailable on the broad C179/C183
   class because C159/C185 is a landed counterexample.  Architectural role
   labels do not remove it from that mathematical class, and no narrower
   reservoir class is introduced to rescue the statement.  Failure to find
   the remaining incommensurate witness is not a theorem and does not fire
   the architecture trigger.  No third verdict form or successor gate is
   permitted.
2. **Witness boundary.** An abstract matrix pair or unrestricted
   operator-norm exponent does not validate PPRG.  Promotion requires one
   exact non-fixed incommensurate passive 2D3C orbit, its physically
   transported return/growth, and on the same witness a finite-frequency
   PDE estimate proving fixed-final-energy concentration by \(q^{3/8}\)
   after every explicit loss.
3. **Fixed-time regular branch.** Any candidate satisfying C182's entrance
   hypotheses on only \(O(1)\) normalized time is rejected as a source of
   \(bq^{3/2}\) focus.
4. **Stationary palette branch.** A candidate whose outer selected
   coefficient forces super-polynomial collateral norm, as in C184, is
   rejected unless the collateral modes are explicitly part of the useful
   endpoint state.
5. **No sixth gate.** If the sole remaining accepted witness search fails,
   record a failed search.  Do not rename the residue, narrow the class ad
   hoc, or create another gate; the architecture trigger remains unfired.
6. **Anti-formalism checkpoint rule.** A checkpoint must deliver an
   explicit-constant estimate, a strict corridor narrowing, or a completed
   dichotomy branch. State-space/operator definitions alone are not a
   checkpoint. Numerical residual minima count only with a full interval
   certificate.

## Audit obligations

- **Resolved by C187:** C136's (5.1) is folded into the single full
  residual/trapping hypothesis. The explicit \(H^3\) substitute is a
  finite-stage continuity theorem only; no uniform structured-state
  constant is claimed.
- **Withdrawn by C187:** the session-only infinite-ladder spectral
  enclosures \(0.66855\ldots\) and \(2.63707\ldots\) have no landed
  operator, normalization, interval certificate, tail resolvent bound, or
  checker. They carry zero evidentiary weight. C120's certified finite
  \(6\times6\) enclosure is the only relevant landed eigenvalue interval.
- Every load-bearing inequality added from this checkpoint onward has an
  explicit numerical constant. Bare \(O(\cdot)\) and \(o(\cdot)\) may
  describe non-load-bearing orientation only.
- No obstruction may be the sole advance of two consecutive checkpoints.
- **Resolved by C188:** the demand-side scalar schedule optimization is
  nonempty and the CKN/Tao “every smooth stage center” premise is
  withdrawn.  The bounded-profile C176 worst-case-envelope infimum is
  \(7/6\), not the current \(5/4\); the equal-shape ledger multiplier is
  \(q^{3/2}\), with an explicit physical shape-factor correction.
- **Standing UVSR obligation after C188:** the CKN/Lin/Tao conclusion is
  scalar-ledger-scoped.  It does not clear backward-cylinder regularity at
  a terminal singular center.  Terminal center tracking, time occupancy,
  pressure/wake terms, the local energy inequality, and critical-norm
  comparison remain attached to UVSR and may not be dropped at profile
  certification.
- **C185 citation restriction from C189:** until the Shvydkoy Theorem 4.1
  paper-body check lands, every downstream use cites only
  \(\|G_{nT}\|_{L^2\to L^2}\ge e^{n/5}>(6/5)^n\).  The
  \(r_{\rm ess}\) form remains source-held and is not a premise.
- **Completed by C190 for the chosen orbit:** the PPRG two-box/trace test
  landed as pre-registered outcome (b), with exact square-zero and trace
  calculations.  C191's class-scope correction does not alter that
  orbit-specific theorem.  A partial float search is not a checkpoint.
- **Resolved by C191:** the C185 floor has sufficient abstract scalar
  exponent in an explicitly long logarithmic action allocation, but the
  direct C182 reconciliation is negative on the landed statements.  There
  is no lower stage-coverage theorem, no common retained-band concentration
  witness or viscous multiplier, fixed-energy scalar reuse cancels, and
  C125 remains open.  Do not credit
  C185's \(L^2\) exponent toward the raw \(q^{3/8}\) concentration demand
  before a same-witness theorem.
- **Class correction from C191:** universal \(m\ne0\) secular no-return is
  false on the broad C179/C183 passive-2D3C class because C159/C185 is a
  returning growing member.  The incommensurate transported-frame witness
  is the only unresolved accepted form.  Failed search is not a class
  theorem, and no third formulation is permitted.
