# Claims registry

Three-party protocol (Richard; the model writing `research/`; Claude Fable 5
writing `audit/`, `docs/`, `engines/`, `gate/`, `data/`): **no claim is
treated as established until the other model has independently re-derived or
reproduced it from artifacts in this repository.** Self-review does not
upgrade status. This rule exists because cross-audit has already reversed
four self-audited conclusions in this collaboration (fidelity-sliver metric,
scale-step loophole, norm convention, float32 C1/C0 inflation — see
`docs/` Addenda G.1–G.3).

Statuses: `VERIFIED-CROSS` (independently re-derived/reproduced by the other
party), `VERIFIED-SELF` (one party only), `ARITH-CHECKED` (stated arithmetic
verified; derivation not yet audited), `UNAUDITED`, `NO-ARTIFACT` (claimed in
chat, files not yet in repo), `REFUTED`, `OPEN` (research target, not a
claim).

None of the entries below is a solution of the Millennium problem, and no
entry asserts a Navier–Stokes singularity.

## Established numerical results (adversarial DNS program, sessions 1–14)

| ID | Claim | Artifacts | Status |
|----|-------|-----------|--------|
| C1 | Dyadic shape-return plateau F* ≈ 0.11 (corrected sliver-proof metric, 48³/64³ stable) | `docs/` Add. G; `engines/return_map_opt.py` | VERIFIED-CROSS |
| C2 | Second-octave endpoint collapse: Q3(shell 2) = 0.4488 at N=384, per-octave 0.617; integrator-independent | `data/stage5/shell2_*` | VERIFIED-CROSS |
| C3 | One-octave scale-fair endpoint crossing at r0=10^6: Q3_oct = 1.00513/1.00769/1.00897/1.00940 at N=128/192/256/320 (their engine, f64, dt-converged) and 1.005491 at 96³ (independent c128 engine) | `data/stage5/`; Add. G.3 | VERIFIED-CROSS |
| C4 | Fixed frontier seed saturates under r0 alone: Q3_oct → 0.96494 (Richardson-flat by r0=10^4) | `data/stage5/high_seed_r0_scan/` | VERIFIED-CROSS |
| C5 | The crossing does not iterate; failure is shape/dephasing (F*), not endpoint mass | Add. G.3 §caveats | VERIFIED-CROSS |

## Palasek-embedding gate (Addenda H, H.1)

| ID | Claim | Artifacts | Status |
|----|-------|-----------|--------|
| C6 | Lagged-intermittency certificate: window 1<b<(1+√13)/4, 2b<β<1+3/(2b) nonempty; α=1+3/(2b) forced; envelope volume ≡ Palasek's fractal-support hypothesis; first-order cancellation fails and second-order clears **on the entire window** | `gate/bundle/parameter_certificate.py`; `audit/window_generalization.py` | VERIFIED-CROSS |
| C6a | Leakage-deficit exponent (b−1)(β−1) | `gate/bundle/GATE_AUDIT.md` §4 | UNAUDITED (asserted, not derived; alternative accounting (3/2)(1−1/b) gives same dichotomy) |
| C7 | Downshifter sweep: ρ²-scaled leak/parent = 1.42/2.88/2.81 at ρ=4/5/8, strength 0.045–0.071, winners revalidated in c128 to 5×10⁻⁷; **not converged** (bests at cutoff); verdict ambiguous-leaning-positive under the pre-registered rule | `gate/results/`; `audit/independent_recheck.py` | VERIFIED-CROSS |
| C8 | Winner mechanism: single-sign helicity 0.79→0.94 rising with ρ (quasi-Beltrami), radii K,K±1; residual leakage 78% in the 2K band at ρ=8 | Add. H.1 | VERIFIED-CROSS |
| C9 | ρ=16 discriminator (pre-registered H.1) | Add. H.1 | OPEN — de-prioritized: `research/2026-07-29-return-cell-no-go-map.md` §1–§2 (if they survive audit) remove fixed finite/helical cells as terminal templates; the equal-radius transfer cell is not a recurrent amplifier. Still informative for transfer-cell scaling only. |

## Repo research notes of 2026-07-29 (the other model)

| ID | Claim | Artifacts | Status |
|----|-------|-----------|--------|
| C10 | Palasek formal full-Laplacian window (2.1): 1<b<a/2, 2b<β<a≤5/2; dormancy flatness θ=β(b−1)/b² | `research/...forced-spectral-recurrence.md` §2; `checks/forced_recurrence_ledger.py` | VERIFIED-CROSS (consistent with H audit and the paper's stated hypotheses) |
| C11 | CMZ locking q/R=d, threshold d<(22−8√7)/9 | §3 | UNAUDITED (citation real: arXiv 2407.06776; exponents to be checked against the paper) |
| C12 | CMZ-retrofit principal-heat no-go: (3.2)+(3.3) ⇒ polynomial K_k, contradicting required exp(cN_k^θ) gain; WKB order cannot repair | §3, (3.1)–(3.5) | UNAUDITED — priority audit target (short algebra, conditional on CMZ support/cutoff scalings) |
| C13 | Bounded-ratio ladder rows (5.4) and window (5.5); dyadic design point r=32, μ=2048, q=2 → E×1/4, D×1/8, Re×2, R×1/16, L³↑ | §5; ledger script | ARITH-CHECKED (design point exact; rows consistent with slender-vortex scaling E~S²a⁴R, D~rate×τ, ‖u‖₃~Sa(a²R)^{1/3}); derivations pending |
| C14 | Asymptotically-finite active cell no-go via Kishimoto–Yoneda compactness (conditional on lattice compatibility + strong normalized bounds) | `research/...return-cell-no-go-map.md` §1 | UNAUDITED — priority audit target |
| C15 | Pairwise reality-sideband obstruction (2.3): |p|≠|q|, P_{p−q}T_−=0 ⇒ T_+=0 | §2 | UNAUDITED (planned: direct sympy check) |
| C16 | Same-eigenvalue Beltrami darkness P[(B·∇)V+(V·∇)B]=0; cross-eigenvalue = (λ−μ)P(B×V) | §4, (4.1)–(4.2) | VERIFIED-CROSS (re-derived: ∇(B·V) identity; in-session) |
| C17 | Additive-closure rank trichotomy; det[p,q1,q2]=−384 for the tested chain ⇒ volume-filling; single-diagonal branch is rank-2 (2D3C) | §3 | VERIFIED-CROSS for the determinant and rank facts (in-session); trichotomy is standard |
| C18 | Growing-ratio schedule exponents (5.1)–(5.2); chain (5.3) forces r_j→∞ | §5 | UNAUDITED |
| C19 | Renormalized unstable-manifold return theorem (five properties) | both notes §6 | OPEN (the program's stated missing theorem) |

## Claimed in chat, not yet in repo (request: push checkpoint v5 files)

| ID | Claim | Status |
|----|-------|--------|
| C20 | Exact affine no-go: c∉{1/5,2/5} ⇒ W=0 (first-moment ∫W⊗W=αI + relative energy α=0) | NO-ARTIFACT |
| C21 | Nonlinear-core compact corrector bound ‖W‖₂=O(r^{7/2}); kills β<2 concentrating corrections | NO-ARTIFACT |
| C22 | Noncompact escape tail U~r^{−√2}, energy growth R^{3−2√2}, order-R⁵ transition ledger | NO-ARTIFACT |
| C23 | Circulation strain-action lower bound 2log(a_n/a_{n+1})−O(1); kills bounded-turnover resets | NO-ARTIFACT |
| C24 | Chae–Wolf closes the finite-energy bounded-gradient β=2/5 periodic Euler-profile route | NO-ARTIFACT (citation real: arXiv 1706.02020) |
| C25 | Refutation of the Zenodo July-2026 helical quasi-trapping flux theorem (six-mode counterexample) | NO-ARTIFACT |
| C26 | Separated-wake pressure estimate (support-gap, cross-stress, periodic-image caveats) | NO-ARTIFACT |

## Audit queue (priority order)

1. C12 (CMZ-retrofit no-go) — short, load-bearing for the strategic pivot.
2. C14 (finite-cell no-go) — load-bearing for abandoning the gate extrapolation.
3. C15 (sideband obstruction) — mechanical sympy check.
4. C13 ledger derivations; C18 schedule.
5. C20–C26 upon artifact arrival.
