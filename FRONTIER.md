# Navier--Stokes frontier

**Branch:** `agent/aug2-integrated-transition`
**Registry frontier:** C187 (2026-08-23 checkpoint)
**Boot rule:** every research session starts from this file, `CLAIMS.md`, and
the referenced artifacts on the branch. Chat history is not a premise.
Every checkpoint commit must update this file in the same commit.

## Current architecture target

Define the exact unforced stage-renormalization operator

\[
   \mathcal R=\mathcal C_{\rm exit}\circ S^{NS}_{T}
\]

on one complete structured state space containing the active packet, carried
wake, phase/modulation variables, and the exit chart. Construct either a
fixed point or an invariant trapping set with one expanding direction.
Leak, pressure, and wake channels are components of the single residual
\(\mathcal R(X)-X\); they are not independent modules to be closed one at a
time.

The unmatched leading object is the **unforced viscous physical-velocity
scale return (UVSR)**:

> a nonzero smooth boundaryless three-dimensional unforced Navier--Stokes
> orbit/profile whose renormalized endpoint is a smaller-copy structured
> state and whose physical velocity normalization realizes C135's
> \(bq^{3/2}\) concentration, with the admitted wake included in the state.

The fixed-point theorem itself is not the novel component. Once a UVSR
approximate profile with a certified full residual exists, the intended
imports are ABC's full-operator spectral splitting and Duhamel contraction,
Elgindi's symmetry-adapted coercivity and compactness, and Chen--Hou's
coercive-bulk/validated-finite-rank trapping method. The detailed primary-
source map is recorded in
`research/2026-08-23-fixed-point-literature-map-c185.md`.

## Current positive target

1. **PPRG admissibility test.** On one explicit smooth bounded unforced
   passive 2D3C orbit, derive two consecutive return-fiber polarization
   maps \(\Phi_1,\Phi_2\in SL(2,\mathbb R)\) in a common C183 return
   frame. Success is pre-registered as either entrywise enclosure of
   \(\Phi_1,\Phi_2\) in C186's respective \(1/100\) boxes about the
   canonical \(U,V\), or a direct enclosure
   \(|\operatorname{tr}(\Phi_2\Phi_1)|\ge2+\delta\) with explicit
   \(\delta>0\). It must also include a witness-specific finite-frequency
   error budget which preserves a stated positive exponent. The scalar
   \(\tau=\operatorname{tr}(N_1N_2)\) is available only if
   \((\Phi_i-I)^2=0\) is separately proved; C183 does not prove that.
2. **UVSR profile search.** Specify the complete Banach state and modulation
   slice for \({\cal R}\), then search for one approximate fixed point or
   trapping orbit by minimizing the norm of the full residual (not selected
   leak projections).
3. **Viscous operator bridge.** Starting from C185's inviscid essential
   growth, prove one explicit finite-stage positive lower bound for the
   actual viscous linearized evolution at the UVSR stage scales.

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
- C185: C159 plus Shvydkoy's theorem gives
  \(r_{\rm ess}(G_T)\ge e^{1/5}>6/5\) for the actual
  infinite-dimensional velocity-form linearized Euler group and
  \(\|G_{nT}\|\ge e^{n/5}\).
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

## Open

- UVSR: no exact or rigorously trapped unforced viscous scale-return profile
  with physical velocity focus is known.
- PPRG realization: abstract noncommuting polarization growth is not enough;
  the two episode blocks must arise from one smooth bounded passive scalar
  and the exact unforced 2D3C evolution, with finite-frequency control.
- A complete state space and normalization for \(\mathcal R\), including the
  retained wake, have not been fixed.
- No coercive bulk estimate or spectral gap has yet been proved for the
  proposed renormalized Navier--Stokes operator.
- C125's relative return, the old BAFL split, and RIGM remain unproved, but
  are to be absorbed into one full residual/trapping inequality rather than
  promoted into further architecture gates.

## Pre-registered kill criteria

These verdicts are fixed before the next computation.

1. **PPRG architecture trigger.** C186 keeps only abstract matrix algebra
   alive. PPRG is negative only after a rigorous theorem for the complete
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
