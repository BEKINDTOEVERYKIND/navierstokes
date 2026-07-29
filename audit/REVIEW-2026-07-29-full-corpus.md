# Full-corpus review, 2026-07-29 (auditing model)

Scope: all sixteen `research/*.md` notes present at commit `c961046`, reviewed
by four independent readers (one line-by-line pass on the two central notes;
three parallel deep-read passes on the remainder, each instructed to
recompute load-bearing algebra adversarially). This file contains findings
and suggestions only; no file under `research/` or `checks/` has been
modified. Statuses are mirrored in `CLAIMS.md`.

## Overall assessment

The corpus is of high quality. Every note carries an accurate claim-boundary
block; no note claims a singularity or a Millennium result; the distinction
"refutes the bridge/proof, not the theorem" is applied consistently in the
external audits. Across roughly 8,300 lines, **no hard sign, exponent, or
algebra error was found** by any of the four readers on any checkable
numbered equation. Items independently recomputed and confirmed include:
the Kelvin amplifier system and its explicit flow map with exact
circulation provenance (kelvin-reynolds §1, §5); the strict-drain cone
R:S=-kappa|S|^2 with PSD two-ray factorization (§2); the single-shear
darkness and pressure-quadrupole/amplification incompatibility (§3, §6);
the two-endpoint cokernel correction (§9); the cylinder pressure operator
multiplier (iκ-2γ)(iκ+1-2γ)-l(l+1) with no integer resonance in 1<γ<3/2
(log-periodic §3); the five zero-flux exponents (§5); the quadratic
no-linear-bifurcation observation (§7.5); the det=-384 lattice rank
computation; the sideband obstruction (2.3) (separate check file); the
CMZ-retrofit polynomial bounds; the helicity covering factor det(O)·r; the
Waleffe-triad Lyapunov classification; the Laurent null-pole steps 1-4; the
Gevrey/Stirling flatness ledgers in all three carrier schedules; and the
1<γ<3/2 window consistency across every note that uses it.

The sixteen notes have converged on **one open object stated in three
equivalent coordinate systems**: (i) projected endpoint submersion Γ =
Π U(S,0)|_Z onto the child constraints with tame right inverse
(kelvin-reynolds §10); (ii) invariant-graph resolvent problem
(λ_pq I − L)Z_pq = G_pq (forward-corrector §5); (iii) periodic cylinder
wake on S²×T_L plus a flux-carrying annular transition (log-periodic §10).
That convergence, plus the §9 cokernel self-correction (which fixed the
earlier two-endpoint right-inverse formulation), is what a healthy program
looks like.

## Flags (ordered by importance; none is a found error)

1. **Affine-germ pressure gauge** (kelvin-reynolds §2.1). For exactly
   affine fields the trace-free part of the pressure Hessian is a free
   gauge, so any trace-free drift S' is attainable; "S'+qκS=0" is a
   convention at germ level, and the drain statement acquires content only
   in the localized problem where pressure is fixed by decay. The note
   half-acknowledges this ("only an exact polynomial germ"); suggest the
   localized statement be treated as the claim of record.
2. **One underived pivot in the Chen–Hou bridge audit**: s ≍ C_ω(τ)
   (chen-hou-boundary §2). All of β = 2.920561 and ν s^{1−2β} hangs on it;
   currently supported only by a BKM self-consistency loop. Needs a
   derivation from Chen–Hou (2.6)–(2.12)/(2.23). Similarly, the Chen
   Euler-II audit needs c > 1/2 across the *entire* constructed α-range,
   not only near α = 1/3.
3. **"Standard energy identity" overstatement** (chen-euler §3): the
   finite-Dirichlet-energy Liouville theorem for stationary 3D NS is open
   in general; the note's conclusion survives because it restricts to
   fast-decay profiles, and should say so where it says "standard".
   Also the B = 2/5 pin is drift-sign-convention dependent (conclusion
   unaffected either way).
4. **Two textual linchpins in the Shahmurov refutations** need verbatim
   anchoring to the target preprints: (a) that the cutoff scale λ is
   isotropic (else C = λ²b may misread the paper's blowup quantity);
   (b) that eq. (5)'s 3/r coefficient is applied to the incompressibility
   of the *physical* meridional velocity with no intervening reweighting.
   Both refutations are airtight conditional on those readings.
5. **Three flatness schedules coexist** (M_j ≍ j against ε_j and Λ_j^{-1};
   M_j ≍ κj against exponential ρ_j; M_j ≍ j²/log j against polynomial
   K_j^{-1}). Each is internally correct; a one-table reconciliation of
   which truncation pairs with which small parameter would prevent silent
   mixing. Related: a merely-Gevrey expansion *in j* of the O(1/j)
   carrier modulation gives only e^{-cj}, not e^{-cj²}
   (polynomial-carrier §4) — that caveat deserves prominence.
6. **Five-ray minimality sub-claim** (localized-stress §3) is
   under-justified (a bijective coefficient map can carry a positive
   vector to the isotropic class in principle); it is hedged and
   non-load-bearing — suggest downgrading to an explicit conjecture.
7. **Near-circularity to keep visible**: the invariant-graph machinery
   presupposes the inviscid fixed cell F_{0,0}(P)=P and tameness of
   λ_pq I − DF_{0,0}(P) — i.e., the hyperbolic structure of the object
   whose existence is the open problem. Both notes flag this; it should
   stay flagged in any future condensation.
8. **Chat-claimed artifacts still absent** (CLAIMS C20–C26: affine no-go,
   O(r^{7/2}) corrector, r^{-√2} tail, strain-action bound, Chae–Wolf
   closure, Zenodo refutation, separated-wake estimate). Standing request:
   push the v5 checkpoint files so these can enter audit at all.

## One constructive disagreement: the cylinder-wake search is computable now

log-periodic §10 states "No GPU calculation is warranted until either an
explicit homogeneous base for (8.1) or a finite-dimensional symmetry
reduction of (2.8) is found." I disagree with the strength of this. The
nonperturbative periodic-shape alternative the note itself names is a
well-posed fixed-point search on a **compact** domain S²×T_L:

- unknowns: (v, w) in spherical harmonics × log-Fourier modes, quadratic
  nonlinearity, amplitude normalized (free by scaling), zero-flux
  constraints holding automatically per §5;
- pressure eliminated through the operator A_γ **whose invertibility on
  the whole window the note itself proves** (§3) — so the reduced system
  is clean and polynomial;
- method: spectral Newton–Krylov from structured random seeds (odd
  symmetry imposed per kelvin-reynolds §4), with arclength continuation
  in (γ, L) — the standard technology for nonperturbative equilibria and
  periodic orbits in shear flows, which does not require a bifurcation
  base;
- decision value, either way: a nonconstant solution makes the
  cylinder-wake theorem target concrete (and gives the annular-transition
  problem an explicit boundary state); systematic Newton failure across
  seeds, resolutions, and (γ, L) sweeps is evidence *for* the Liouville
  alternative the note says "would decisively kill this outer-wake
  route."

This is not a broad DNS and does not measure a scalar gain; it is a
structured existence probe for the exact object the program now needs,
with a pre-registrable either-way outcome — the same standard the program
applies elsewhere. Proposed division of labor as before: I build the
cylinder solver independently (audit namespace), the research model
audits it, Richard runs the continuation sweeps.

## Priority queue going forward

1. Cylinder-system solver spec + implementation (above), pre-registered.
2. Flag 2 (s ≍ C_ω derivation) and flag 4 (Shahmurov verbatim anchors) —
   external-facing correctness of the strongest refutations.
3. Line-level audit of the §1 compactness no-go (return-cell-no-go-map)
   — structure already reviewed, no gaps found at sketch level.
4. C20–C26 upon artifact arrival.
