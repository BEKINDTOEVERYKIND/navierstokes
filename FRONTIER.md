# Navier--Stokes frontier

**Branch:** `agent/aug2-integrated-transition`
**Registry frontier:** C190 (2026-08-25 checkpoint; C189 is the ingested
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

## Current positive target

1. **PPRG complete-class decision.** C190 excludes the chosen C186 orbit
   exactly but leaves the broader admissible passive 2D3C class open.  Any
   positive realization must still give two consecutive maps in one
   common C183 frame, either in C186's respective entrywise \(1/100\)
   boxes or with
   \(|\operatorname{tr}(\Phi_2\Phi_1)|\ge2+\delta\) for an explicit
   \(\delta>0\), and it must include a witness-specific finite-frequency
   error budget preserving a stated positive exponent.  The scalar
   \(\tau=\operatorname{tr}(N_1N_2)\) remains available only after
   \((\Phi_i-I)^2=0\) is proved for the actual maps.  A negative verdict
   still requires the complete-class theorem in the kill criteria below.
2. **Viscous operator bridge.** Starting from C185's robust inviscid
   operator-norm lower bound \(\|G_T\|>6/5\), certify three separate
   multiplicative operator/norm-
   amplitude factors:
   viscous damping \(\mathfrak d\), retained band fraction
   \(\mathfrak b\), and physical normalization \(\mathfrak n\).
   The pre-registered pass threshold is
   \[
     \mathfrak d,\mathfrak b,\mathfrak n\ge{99\over100},
     \qquad \varepsilon_{\rm ff}\le{1\over100}.
   \]
   It gives the strict explicit normalized gain
   \[
    G_{\rm visc}>{6\over5}\left({99\over100}\right)^3-{1\over100}
    ={2885897\over2500000}>{23\over20}.
   \]
   Each factor and the additive finite-frequency error must be proved for
   the same witness and stage; they may not be merged into an unnamed
   constant.
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

## Open

- UVSR: no exact or rigorously trapped unforced viscous scale-return profile
  with physical velocity focus is known.
- PPRG realization: abstract noncommuting polarization growth is not enough.
  C190 proves that the chosen C186 rotating-gradient orbit cannot realize
  the two-block witness.  The complete admissible passive 2D3C class is
  still undecided; any positive pair must arise from one exact unforced
  evolution in a common frame and carry finite-frequency control.
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

1. **PPRG architecture trigger.** C186 keeps only abstract matrix algebra
   alive, and C190 excludes only its explicit rotating-gradient candidate.
   PPRG is negative only after a rigorous theorem for the complete
   admissible Kelvin class: for example, exact Euler/passive integrability
   forces a common invariant flag together with an explicit growth upper
   bound too small for C182's \(q^{3/8}\) deficit, or an exhaustive cocycle
   enclosure gives the same upper bound. Excluding the canonical two-box
   witness, or failing a non-exhaustive search, is not such a theorem. Once
   the full negative dichotomy lands, do not create a successor gate. The
   next deliverable must be a standalone architecture-level no-go assembled
   from C181--C186, followed by a re-derivation of the base stack from the
   surviving constraints.
2. **Witness boundary.** An abstract pair of noncommuting unipotents keeps
   PPRG algebraically alive but does not validate the architecture. PPRG is
   promoted only after those blocks are realized by one exact passive 2D3C
   orbit and yield a finite-frequency PDE estimate with explicit constants.
3. **Fixed-time regular branch.** Any candidate satisfying C182's entrance
   hypotheses on only \(O(1)\) normalized time is rejected as a source of
   \(bq^{3/2}\) focus.
4. **Stationary palette branch.** A candidate whose outer selected
   coefficient forces super-polynomial collateral norm, as in C184, is
   rejected unless the collateral modes are explicitly part of the useful
   endpoint state.
5. **No sixth gate.** A negative PPRG verdict ends the modular
   scalar-ledger architecture. Renaming its residue is not progress.
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
  calculations.  The complete-class success and architecture-negative
  criteria remain unchanged.  A later candidate may not relax them, and a
  partial float search is not a checkpoint.
