# Toward the localized transported writer: gauge rigidity, wake ledger, floor, and the surviving window

Date: 2026-08-03. Author: auditor (Claude). Status: **attempt**, not a proof.
Verification scripts: `audit/wake_ledger.py`, `audit/affine_transport.py`,
`audit/window_search.py` ( all identities below marked [S] are exact sympy checks, all
numbers marked [N] are outputs of the scripts).

## 0. Claim boundary

This note attacks the decisive missing PDE result identified in
`research/2026-08-01-localized-chirp-scale-ledger.md` §6 (the five-item
localized transported geometric-optics lemma). It does **not** prove that
lemma and does not construct a Navier–Stokes singularity. What it proves,
at the same frozen-ledger level of rigor as the existing notes:

1. **Gauge rigidity** (Lemma 1): the common carrier phase is not a lever —
   it drops out of the leading low beat exactly. The collar wake cannot be
   phase-gauged away.
2. **Exact wake ledger** (Lemma 2): for the flattened chirp × characteristic
   collar, the entire non-child low output at leading order consists of
   exactly three terms with explicit coefficients and disjoint routing
   directions.
3. **Wake floor and a second derivation of the nested-torus condition (12)**
   (Lemma 3): any coherent transverse collar has an irreducible DC wake, and
   wake-smallness *forces* `b(2α−4) > 1`. Condition (12) now has two
   independent origins (curvature coherence and wake routing).
4. **Exact affine transport** (Lemma 4): on an affine parent the writer's
   function classes (quadratic chirps, affine jets) are *exactly*
   transport-invariant, and the leading beat identity is form-invariant
   under frame advection: `L₀ = 2i (q·W) W`, `W = t×(w×t)`. Item 2 of the
   five-item lemma is thereby reduced to curvature error, which is
   controlled precisely by (12).
5. **The window theorem** (§5): the full set of constraints — including two
   new ones ((L2), (L3)) — has a nonempty interior with balanced margins
   ≈ 0.046 at `(α,β,b) = (2.454, 2.407, 1.173)` [N]. A natural stronger
   reading of the wake tolerance ((L2′)) is *strictly infeasible* on the
   whole window; which reading is correct is now the sharpest open fork.

Throughout: parent frequency `N`, strain `g = N^β`, carrier `Λ = N^{β/2}`,
child `Q = N^b`, `ε = b−1`, `Δ = α−β`, `μ = β−2b`, window
`1 < b`, `2b < β < α < 5/2`. Child block: minor `r_Q = N^{−b}`, major
`R_Q = N^{b(4−2α)}`. Gain action `G = (log Λ)²`. Define

```
m12 := b(2α−4) − 1        (margin of the nested-torus condition (12)),
```

so `R_Q = N^{−(1+m12)}` and `R_Q/r_N = N^{−m12}`.

## 1. Lemma 1 (gauge rigidity of the low beat)

Take the exact-curl two-carrier block of
`research/2026-08-01-generalized-phase-localization.md` §1, but with a
general *common* phase `Θ` in addition to the relative phase `Ψ`:

```
φ_± = Λ t·x + Θ(x) ± Ψ(x)/2,      U_± = (1/iΛ) ∇×( c e^{iφ_±} ),   c = w×t.
```

**Lemma 1.** For arbitrary smooth `Θ, Ψ` (no size assumptions), the
demodulated low beat `e^{−iΨ}[B(U₊,Ū₋)+B(Ū₋,U₊)]` has `1/Λ`-expansion whose
order-0 term is exactly `2i d w`, `d = w·∇Ψ` — **independent of Θ**. `Θ`
first appears at order `1/Λ`, and there only through bilinear terms
`∇Ψ·∇Θ` and Hessians of `Θ`. [S: wake_ledger.py (A)]

Consequences. (a) There is no O(1) common-phase lever: any attempt to cancel
the collar wake by re-splitting the two carrier phases fails at leading
order. The wake floor of Lemma 3 is therefore structural, not an artifact of
the symmetric split `±Ψ/2`. (b) Conversely, `Θ` is *harmlessly free* at
leading order: a common phase of size `|∇Θ| ≲ Q`, `|∇²Θ| ≲ Q²` perturbs
only the `O(Q²/Λ)` layer — available for curvature compensation without
touching the child.

## 2. Lemma 2 (exact wake decomposition)

Take the flattened chirp × characteristic collar,

```
Ψ = ψ(s) χ(τ, ζ),    s = r·x,  τ = t·x,  ζ = (h−Hr)·x,   w = r + Hh,
```

with `ψ' = 2as` on the core `|s| ≤ R_c`, `ψ'` ramping to 0 over the s-flank
of length `L_s`, and `χ` a transverse collar of lengths `(L_t, L_ζ)`.

**Lemma 2.** Exactly (all [S: wake_ledger.py (B)]):

```
d  = D_w Ψ = ψ'(s) χ            (the characteristic collar keeps d clean),
q_t        = ψ χ_τ              (the forced high-leak coefficient),
```

and the leading low output, after removing the exact gradient
`2∇e^{iΨ}` (pressure), is `2i[(w·q)w − q]e^{iΨ}` with

```
(w·q)w − q  =  H ψ' χ · h        (child jet)
             − ψ χ_τ · t         (t-wake)
             − ψ χ_ζ · (h−Hr)    (ζ-wake).
```

With `χ ≡ 1` this is the pure child `Hψ' h`; the `O(Λ^{−2})` remainder is
`(v·∇)v/(2Λ²)` as in their (2.1), of relative size `1/(GΛ²R_c²)` on the
core — subsumed by (L3) below. The three non-child terms have *disjoint
supports and explicit directions*: the s-flank modifies the jet on
`R_c < |s| < R_c+L_s` with strain-reversed coefficient `≤ (R_c/L_s)`·(core
strain); the t- and ζ-wakes live on the transverse collars with
coefficients `max|ψ|/L_t`, `max|ψ|/L_ζ`.

## 3. Lemma 3 (wake floor; (12) is forced twice)

**Lemma 3.** (a) *DC floor.* For any collar profile `χ` descending from 1
to 0 over length `L`, `∫ ∂_τ χ dτ = −1`; no oscillatory or staircase
profile changes the zero-mode [N: window_search.py D]. Since `max|ψ| ≥ aR_c²`
whenever the core writes the full jet, the wake-to-child coefficient ratio
obeys the pointwise floor

```
wake/child  ≥  max|ψ| / (L · 2aH R_c)  ≥  R_c / (2H L).
```

(b) *Coherence cap.* A coherent collar must sit inside the region where the
parent strain is a single frozen matrix to relative accuracy `1/G` over the
gain window; this caps `L ≤ r_N/G²`. With `R_c = R_Q`:

```
wake/child  ≥  (G²/2H) · R_Q/r_N  =  (G²/2H) · N^{−m12}.
```

(c) *Second derivation of (12).* Wake-smallness therefore *requires*
`m12 > 0`, i.e. `b(2α−4) > 1` — the nested-torus condition — independently
of its original motivation (fitting the child inside the parent minor
radius / affine-coherence of the curvature error). Two unrelated failure
modes pin the same inequality; on this route (12) is not optional.

All three routing directions (s-flank, t-collar, ζ-collar) enjoy the same
uniform suppression `δ := G² N^{−m12}` when their lengths are taken maximal.

*Loophole recorded (open).* The cap (b) assumes the collar writes
coherently. A deliberately decoherent collar (carrier phases drifting O(1)
beyond `ℓ_coh`) might time-average part of the wake; nothing here rules
that in or out. This is flagged, not used.

## 4. Lemma 4 (exact affine transport)

On an affine parent `U_p = Mx` (`tr M = 0`, `M` arbitrary, including
nonnormal):

**Lemma 4.** (a) *Phases transport exactly and algebraically.*
`φ = k(t)·x + ½ xᵀP(t)x` with `k' = −Mᵀk`, `P' = −MᵀP − PM` satisfies
`∂_tφ + (Mx)·∇φ = 0` [S]. The quadratic-chirp class is exactly invariant —
no WKB error from transport of the phase.
(b) *Affine jets are exactly invariant.* For `u = A(t)x`, the linearized
pressure closes in the quadratic class (`−∇(½xᵀSx) = −Sx` for symmetric S),
so only the antisymmetric part of `A' + AM + MA` is constrained and the
affine class propagates exactly [S].
(c) *Polarizations follow Kelvin dynamics.* `a' = −Ma + 2(k·Ma)k/|k|² −
ν|k|²a` preserves `a·k = 0` to machine precision including Jordan-type `M`
[N: ≤ 2.5e−15 over 3 strain times].
(d) *The beat identity is form-invariant under frame advection.* For a
general instantaneous frame (`t, w` arbitrary, not orthonormal — as produced
by advection `k ↦ e^{−Mᵀt}k`), the order-0 low beat of the curl block is
**exactly**

```
L₀ = 2i (q·W) W,        W := t×(w×t) = |t|² w − (t·w) t     [S].
```

For orthonormal frames this is the known `2i d w`; under advection the
writer keeps writing along the projected polarization `W(t)` with coupling
`q·W`. Combined with (a)–(c): on an affine parent the *entire* frozen
writer transports exactly at leading order with computable time-dependent
coefficients; the only errors are the `O(Q/Λ)` polarization corrections
(Kelvin vs frozen `w`) and the *curvature* of the true Gavrilov parent,
`∇²U_p ≠ 0`, whose accumulated effect over the gain window is
`G·R_Q/r_N = G N^{−m12}` — again controlled exactly by (12).

Design consequence: because the affine flow maps affine jets to affine jets
and quadratic chirps to quadratic chirps, the writer should be *designed in
Lagrangian coordinates*: choose `Ψ₀` so that its advected image at the
effective writing time (the last `O(1/g)` of the window, which dominates
the `∫e^{2(1−ϑ)gt}` accumulation) is the desired child phase. This is an
O(1) redesign inside an invariant class, not a perturbation.

## 5. The window theorem and the fork

New exponent conditions derived above, joined to the existing ledger:

```
(1)    1 < b,  2b < β < α < 5/2
(12)   m12 = b(2α−4) − 1 > 0                       [forced twice: §3, §4]
(17a)  5 − 2α > εΔ          (17b)  μ > 2εΔ         [their bare-chirp ledger]
(L2)   m12 > 2b(b−1)Δ                              [wake ≤ ρ_Q² tolerance
                                                    of the NEXT parent]
(L3)   β/2 > 1 + m12                               [Λ R_Q > 1: carrier
                                                    oscillates inside core]
```

**Theorem (window, numerical).** The joint window (1)∧(12)∧(17)∧(L2)∧(L3)
is nonempty; the max-min-margin point is

```
(α, β, b) = (2.4541, 2.4072, 1.1726),   min margin 0.0459,
```

with all nine margins in `[0.046, 0.173]` [N: window_search.py B]. The
review candidates (2.45, 2.35, 1.15) and (2.49, 2.45, 1.21) both lie inside
(min margins 0.0005 and 0.0100).

**Theorem (fork, one branch closed).** The stronger "paranoid" tolerance —
wake strain measured against the *current* parent at its `ρ_N²`
invariant-breaking accuracy —

```
(L2′)  m12 > (b−1)(2α−β−1)
```

is **strictly infeasible** on the whole window: since `α ≤ 5/2` gives
`m12 ≤ b−1`, and `α > 2`, `β < α` give `(b−1)(2α−β−1) > (b−1)(α−1) > b−1`,
the margin is negative everywhere (numerical sup `−7·10⁻⁵ → 0⁻` only at the
degenerate corner `b→1`) [N: window_search.py A].

So the route lives or dies on a single question:

> **Fork.** Is the collar wake, relative to the child it accompanies, a
> *structured co-growing perturbation* (it is generated by the same beat,
> hence proportional to the child at every instant, entering the capture
> ODE as a coefficient shift of relative size δ = G²N^{−m12}) — in which
> case (L2) suffices, `Σ_k δ_k < ∞` by double-exponential growth of `N_k`,
> and the window above survives? Or does it enter the *next* stage as an
> unstructured parent-equation residual subject to the r²-amplified `ρ²`
> tolerance in the invariant-breaking direction — in which case the
> relevant condition is (L2) (satisfiable, margins above), or its paranoid
> strengthening (L2′) (infeasible, route dead)?

The distinction is precisely the one their affine-core note (R3)/(R11)
draws between coefficient perturbations (benign) and additive forcing
(malignant). The wake is *not* arbitrary additive error: it is proportional
to the child throughout the window (both are written by the same carriers),
so it never threatens the `e^{−G}` seed with O(1) noise — the seed-survival
trap of (R11) does not apply to it. My assessment: the co-growing/(L2)
reading is correct, but this is exactly the point a transported version of
their capture-robustness estimates must settle. It is a finite, well-posed
ODE-perturbation question, not a PDE question.

**Finite-N honesty.** At the optimum, the floor factor `G²N^{−m12} < 1`
requires `N > N* ≈ 10^247` [N: window_search.py C]. For the analysis this
is harmless (choose `N₀ ≥ N*`; the cascade `N_k = N_{k−1}^b` then keeps all
stages valid, and the initial datum remains smooth). But no simulation will
ever see this regime — consistent with their closing line that no GPU
search resolves these requirements, and a warning against numerical
"validation" of the wake bounds at accessible N.

## 6. What this advances, and what remains

Against the five-item lemma of their §6:

1. *Phase localization with support ratio (4)* — *advanced*: localization
   achieved by phase (evading their §2 exact envelope obstruction and §3
   tube no-go, which needed small geometric intersections), at the exact
   cost of the three wake terms of Lemma 2, uniformly suppressed by
   `δ = G²N^{−m12}` (Lemma 3), with (12) now forced.
2. *Invariance under the curved Gavrilov parent* — *advanced*: exactly
   invariant on the affine part (Lemma 4); remaining error is curvature
   only, `G N^{−m12}`, plus `O(Q/Λ)` polarization corrections. What is NOT
   done: the transported BAS/Floquet gain estimate on the rotating flank
   (the polarization must stay near the growing Kelvin branch through the
   full gain window on the actual Gavrilov frame dynamics).
3. *Robustness of viscous capture under those errors* — *reduced* to the
   Fork above: a finite ODE-perturbation question (co-growing coefficient
   shift vs invariant-breaking forcing) against their existing capture
   robustness series.
4. *Conversion of the written jet into the full next Gavrilov block* —
   *not advanced* here, beyond noting the written object
   `Sym∇(γ s h) = (γ/2)(h⊗r + r⊗h)` has exactly the `(μ,0,−μ)` one-colour
   strain form their multicolour note §5 requires, and that the
   strain-reversed s-flank is pushed to distance `L_s` with the same δ
   suppression.
5. *Summability through the cascade* — only the trivial part: `Σ_k
   G_k²N_k^{−m12} < ∞` for any `m12 > 0`.

Also still open from their notes: the separated-envelope Sobolev estimate
(their §6 last paragraph), and the decoherent-collar loophole recorded in
Lemma 3.

## 7. Claim registrations

Registered in `CLAIMS.md` as **C102–C106** (the range C50–C101 is reserved:
the research model's own ledger has grown into that range in-session, but
those artifacts have not yet landed in this repository and remain
NO-ARTIFACT under the registry protocol):

- C102: Lemma 1 (gauge rigidity), status VERIFIED-SELF (sympy, this note).
- C103: Lemma 2 (exact wake decomposition), VERIFIED-SELF.
- C104: Lemma 3 (wake floor; (12) forced by routing), VERIFIED-SELF
  (algebra) + the coherence cap is a modeling step, flagged.
- C105: Lemma 4 (exact affine transport; form-invariant beat
  `L₀ = 2i(q·W)W`), VERIFIED-SELF.
- C106: window theorem + fork ((L2′) infeasible; (L2) window nonempty,
  optimum (2.4541, 2.4072, 1.1726), min margin 0.046; N* ≈ 10^247),
  VERIFIED-SELF (numerical + the analytic (L2′) chain).

Cross-model verification requested on: the coherence cap `L ≤ r_N/G²`
(Lemma 3b — the one modeling assumption in the floor), the Fork resolution
via transported capture robustness, and independent re-derivation of
`L₀ = 2i(q·W)W`.
