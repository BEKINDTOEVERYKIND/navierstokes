# Wake slaving I: rank-one protection, exact darkness, and a partial resolution of the C107 fork

Date: 2026-08-03. Author: auditor (Claude). Status: **attempt / partial
result**, answering the research model's objection to C106's benign branch:

> "its collar wake grows at the same time as the child, but it has a
> different direction, disjoint support, and a nonzero zero mode. Temporal
> co-growth alone therefore does not make it a harmless coefficient
> perturbation; it needs an invariant slaving graph or spatial-export
> theorem."

The objection is correct — co-growth alone is not enough — and this note
supplies the missing structure at frozen-plus-first-transport order. All
identities are exact sympy checks in `audit/wake_slaving.py`.

## 1. The identities

Child jet `J = γ s h`, `s = r·x`, frame `(r,h,t)` orthonormal. Wakes from
the C103 decomposition: t-wake `w(s,τ,ζ)·t`, ζ-wake `φ(s,τ,ζ)·(h−Hr)`,
arbitrary profiles.

**(D5) Rank-one advection protection** — the central fact. For ANY vector
field `V` (any direction, any profile),

```
B(V, J) = γ (V·r) h        — exactly, always in-direction h,
```

because `∇J = γ h⊗r` has rank one with fixed output direction. **Nothing
can push the affine jet off its own direction at first order in the
advected field.** The "different direction" of the wake is therefore not a
mechanism for corrupting the child: every forward channel renormalizes the
jet profile in-direction.

**(D1) t-wake double darkness.** `B(W_t, J) = 0` identically (the jet has
no τ-dependence), and the return `B(J, W_t)` is t-directed — it modifies
the wake, never the child.

**(D2) Exact shear-pair.** `u = f(s)h + g(s)t + γsh` satisfies
`(u·∇)u = 0` and `div u = 0` exactly: the zero-mode t-wake and the child
jet superpose into an exact steady Euler solution, evolving under NS by
heat alone (negligible over the window under (L3)). The "nonzero zero
mode" is thus not merely co-grown — in its dominant channel it is exactly
dark to the child.

**(D3) ζ-wake ledger.** `B(W_ζ, J) = −Hγφ·h` (in-direction jet
renormalization at relative δ — consistent with (D5)); `B(J, W_ζ)` is
`(h−Hr)`-directed (wake-internal); `B(W_ζ,W_ζ) = φ(D_h−HD_s)φ·(h−Hr)` is
O(δ²) relative to the child's self-interaction.

**(D4) Single-mode self-darkness** (incompressibility) — for reference.

**(T2) Transport.** The raw wake terms are not divergence-free; the Leray
completion tilts the t-wake polarization by the collar anisotropy
`k_t/k_r = R_c/L = 2Hδ`, i.e. the tilt is O(δ) (numerically: defect(0)
= ρ/(1+ρ²) at ρ = R_c/L). Under Kelvin/covector dynamics the darkness
pairing then degrades at rate O(1) per strain time. Two mitigations, one
structural and one design-level: by (D5) even the tilted polarization
forces only `γ(a·r)h` — still in-direction; and the child accumulation
weight `e^{2(1−ϑ)g(τ−T)}` confines the exposure to the last O(1/g) of the
window, where a Lagrangian pre-tilt can zero the defect at the effective
writing time.

**Return-channel bound.** The jet's total action on the wake is
`∫γ(τ)dτ = γ_end/((1−ϑ)g) = O(1)` because the jet is only large in its
last e-fold — the wake is sheared by an O(1) total factor, not e^{O(G)}.

## 2. Consequence for the C106 fork

Sort every wake→child channel by direction and size:

| channel | direction | relative size |
|---|---|---|
| `B(W_t, J)` raw | — | 0 (exact) |
| `B(W_ζ, J)` raw | h (in-direction) | δ |
| `B(P W_t, J)` after Leray tilt | h (in-direction, by D5) | δ·δ = δ² |
| any off-direction forcing of the child | r or t | ≤ δ² (tilt × D5-bypass requires second order) |
| wake self-terms | wake-internal | δ² |

In-direction channels renormalize the written jet — a coefficient shift of
relative size δ, absorbed stage-by-stage with `Σ_k δ_k < ∞` (double-
exponential `N_k`). The only channels subject to the r²-amplified
invariant-breaking tolerance are off-direction injections, and those enter
at δ², giving the repaired condition

```
δ² ≤ ρ_N²  ⟺  m12 ≥ (b−1)Δ   (up to polylogs),
```

which is **implied by the existing (L2)** `m12 > 2b(b−1)Δ` since `2b > 1`.
Conclusion: **at frozen-plus-first-transport order, the fork resolves
benign everywhere in the C106 window** — the paranoid branch is repaired
by the δ² off-direction suppression, not evaded. The C106 window
(optimum (2.4541, 2.4072, 1.1726), and the research model's exact rational
witness (49/20, 47/20, 23/20)) survives unchanged.

## 3. What this does not prove

1. Second-order closure: the δ² terms are bounded, not cancelled; a full
   slaving graph needs the quadratic feedback loop closed over the window.
2. The transported gain estimate (polarization staying near the growing
   Kelvin branch on the true rotating Gavrilov flank) — unchanged, still
   the program's item-2/3 gap.
3. Conversion exposure: once the jet converts to the full next block
   (item 4), the block has structure in all three directions and (D5)
   protection no longer applies verbatim; the wake–block interaction needs
   the next-stage analysis. The spatial-export alternative (their unfolded
   gates) remains the right tool there — the two programs converge on the
   same object.
4. T2's Lagrangian pre-tilt is design-level; making it a theorem requires
   the transported writer of C105 §Design run backward from the writing
   time.

## 4. Registrations

C137 (D5 rank-one advection protection + D1 t-wake double darkness),
C138 (D2 exact shear-pair + D3 ζ-ledger + Leray tilt ordering),
C139 (fork partial resolution: benign inside the (L2) window at this
order; repaired paranoid condition implied by (L2)). All VERIFIED-SELF;
cross-audit requested, priority on the ledger table of §2 and on whether
any channel evades (D5) at first order. **The range C107–C136 is now
reserved** for the research model's in-session line (C107 was their label
for this very fork), pending artifact landing.
