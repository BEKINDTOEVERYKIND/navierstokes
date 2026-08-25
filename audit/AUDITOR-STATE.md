# Auditor state (handoff file)

Purpose: lets ANY fresh auditor session resume with zero loss. Boot rule:
read this file, `FRONTIER.md`, `CLAIMS.md`, and the artifacts they cite.
Chat history is not a premise. Updated: 2026-08-25, checkpoint 0c5bf3b +
branch `audit/c185-c187-cross-audit`.

## Standing constraints (from R, verbatim where it matters)
- Never modify `research/` or `checks/` — auditor namespace is `audit/`,
  `CLAIMS.md` (shared registry), `docs/`, `engines/`, `gate/`, `data/`.
- "the goal here is for us to find novel results, not to verify the public
  results"; R values only clear breakthroughs, not incremental progress.
- Registry protocol: no claim is established until the party that did not
  produce it re-derives it from artifacts in this repository
  (VERIFIED-CROSS). Self-review does not upgrade status.
- Deliver working-notes updates as deltas, not whole files.

## Where the program stands
- Authoritative statement: `FRONTIER.md` at 0c5bf3b. Architecture is the
  stage-renormalization fixed point/trapping; the single unmatched object
  is UVSR; positive targets are the PPRG SL(2) witness test (pre-registered
  boxes/trace criterion), the UVSR profile search, and the viscous bridge
  from C185.
- Latest cross-audit: `audit/AUDIT-2026-08-25-c185-c187.md` (C188).
  C186, C187 VERIFIED-CROSS in full. C185 VERIFIED-CROSS with a split:
  norm-growth half robust from the abstract-level Shvydkoy inclusion;
  r_ess half awaits paper-body Thm 4.1 (math-ph/0412019). C159 premise
  independently corroborated (float re-implementation, wide margins).
  Prior auditor flags ((5.1) uniformity; session ladder numbers) are
  RESOLVED and closed by C187.

## Open auditor queue (priority order)
1. Shvydkoy Theorem 4.1 paper-body check (closes the C185 split).
2. PPRG witness artifacts (Φ₁, Φ₂ enclosures) when they land — verify the
   interval enclosures independently; the pre-registered criteria are in
   FRONTIER.md §1 and C186 §2 (1/100 boxes or |tr| ≥ 2+δ).
3. UVSR corridor: check the required concentration against CKN/ε-regularity
   and quantitative critical-norm bounds — NOTE: the entire corpus has zero
   engagement with this constraint theory; it is the referee's first
   objection and has been recommended to the research model as the
   demand-side computation to run before the profile hunt.
4. Earlier standing items: external review of the countable-stagnation
   Liouville theorem (branch 2026-08-01); auditor claims C102–C106 and
   C137–C139 (transported writer, wake ledger/floor, wake slaving) remain
   VERIFIED-SELF and open for their cross-audit.

## Numbering discipline
Research model's registry is canonical once landed (currently through
C187; auditor audit = C188). On any collision, the landed number wins and
the auditor renumbers.

## Operational notes
- This session's git proxy blocks pushes (repo not in the session's source
  set); reads work. Deliverable flow: commit locally → `git bundle` →
  hand to R → applied/pushed from an authorized clone (the research
  model's session pushes fine). A fresh session with the repo attached as
  a source at creation restores direct push.
- GPU (Colab A100) is available via R for validated-numerics work; prefer
  interval/outward-rounded certificates over floats for anything
  load-bearing.
