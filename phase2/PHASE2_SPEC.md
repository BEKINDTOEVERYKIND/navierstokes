# Phase 2 — The affine-core continuation: full specification

**Object.** A stationary Euler similarity profile `U` on R³ solving, after Leray projection,

    L_γ U + P(U·∇U) = 0,   div U = 0,   L_γ = (1−γ) + γ y·∇        (P2.1)

for some γ in the Type-II window, connecting a prescribed **inner** structure (an exceptional
fixed-line affine jet) to a prescribed **outer** decay law (the matched Euler tail). A verified
solution is the missing core of a Clay-negative certificate (v3 §5.2 items 1–2); a verified
*nonexistence* across the window is the sharpest available positive-direction rigidity evidence.
Both outcomes are wins. All formulas below were verified symbolically in-session.

## 1. The window (exact statements)

- CIV Theorem 2.1 (arXiv:2602.17570): finite kinetic energy + the stated gradient-rate hypothesis
  force **γ ≥ 2/5** (non-strict as published).
- A nonzero matched far-field tail with locally finite energy forces **γ > 2/5 strictly**
  (the shell energy integral is logarithmic at γ = 2/5 — v4 §8.4, independently confirmed).
- Endpoint L³ regularity forces the critical norm to diverge: **γ < 1/2**.
- Practical scan window: **γ ∈ [0.41, 0.49]**, with γ = 0.45 as the first working point.

## 2. Inner data: the canonical affine family (verified)

Every vortical affine exact solution has, after a rigid rotation, the two-parameter form

    M(λ,α) = [[1, 0, 0],
              [0, λ, −α],
              [0, α, −1−λ]],       trace M = 0,  M + M² symmetric automatically.   (P2.2)

The similarity field (M + γI)y has a **fixed line** iff

    α² = (λ+γ)(λ+1−γ),                                                            (P2.3)

and then spec(M+γI) = {0, 2γ−1, 1+γ} automatically — the spectrum any exceptional
fixed-set profile must exhibit (v3 §7.4, eq. 11). Free inner parameters: **γ, λ** (α fixed by
P2.3 up to sign; overall rotation pinned by convention). U₀(y) = M y solves (P2.1) exactly
(with P₀ = −y·(M+M²)y/2): this is the **anchor solution** — infinite energy, wrong far field.

## 3. Outer data: the matched-tail expansion (verified)

Matching to a finite nonzero outer flow forces the leading exterior behavior

    U ~ U₁ = R^{d₁} Φ(ŷ),   d₁ = 1 − 1/γ ∈ (−3/2, −1),   div U₁ = 0,             (P2.4)

with Φ a free divergence-compatible angular profile (this is the outer shooting data — a
*function*, not a number). The exterior series U ~ Σ U_n, deg U_n = 1 − n/γ, is **nonresonant**:
L_γ U_n = (1−n)U_n and deg[(U_i·∇)U_j] = d_{i+j}, so

    U_n = (1/(n−1)) P Σ_{i+j=n} (U_i·∇)U_j,   n ≥ 2.                              (P2.5)

Robin boundary condition at the truncation radius R_out (replaces any artificial Dirichlet/energy BC):

    R ∂_R U = (1 − 1/γ) U + O(R^{1−2/γ}).                                          (P2.6)

## 4. Formulation and homotopy

Unknowns: the field U on R³ (or B_{R_out} with P2.6), plus (γ, λ) and the tail profile Φ.
Structure: a nonlinear eigenvalue/matching problem — the fixed-line jet contributes a strongly
unstable local direction (eigenvalue 1+γ ≈ 1.45), the tail a weakly stable one (2γ−1 ≈ −0.1),
so expect a finite-codimension matching selecting (γ, λ, Φ) jointly. There is no known solution
to continue from; use the **anchor homotopy**:

    Solve on B_R with boundary data  (1−θ)·(M y)|_{∂B_R} + θ·(tail form P2.4/P2.6),
    continue θ: 0 → 1 at fixed moderate R, then continue R upward re-imposing P2.6.

θ = 0 is the exact affine anchor (machine-zero residual — the code sanity gate), so Newton
starts converged and continuation discovers where (if anywhere) the branch fails.

## 5. Three solver routes

**A. Newton–Krylov collocation (reference).** Compactified radius ρ = R/(R₀+R) with mapped
Chebyshev (N_R = 64–128); angular expansion in vector spherical harmonics, ℓ ≤ L = 24–48,
divergence-free (toroidal/poloidal) basis in the frame with the fixed line along e₃; pressure
eliminated by projection. Jacobian-free GMRES, preconditioned by the linearization about the
affine jet; pin the rotation-about-fixed-line symmetry. Continuation in (θ, R, γ, λ) with
pseudo-arclength. Cost: workstation/CPU-scale; best precision.

**B. Stabilized τ-relaxation (exploratory).** Integrate the time-dependent similarity flow
(eq. 2 of the note) with ν = 0 from blended data, suppressing the known 1+γ unstable direction
by projective stabilization; cheap dynamics to locate basins before Newton.

**C. PINN least-squares (Colab-native — recommended first).** Ansatz
U = χ(R)·M y + (1−χ)·R^{1−1/γ}·N_Φ(ŷ) + N_c(y) with a smooth partition χ, learned angular
tail N_Φ, learned correction N_c (SIREN/Fourier-feature MLP, float64); loss = |P2.1 residual|²
at collocation points (radially stratified, ρ-uniform) + div penalty + jet-match penalty near
the fixed line + tail-form penalty at large R. Optimize Adam → L-BFGS → Gauss–Newton
(the DeepMind high-precision recipe, arXiv:2509.14185). Scan (γ, λ) on a coarse grid.

## 6. Validation anchors (all pre-registered)

1. θ = 0 anchor: discrete residual of the exact affine jet ≤ 10× machine epsilon (starter script checks this).
2. Fixed-line spectrum: computed DV eigenvalues at the fixed set = {0, 2γ−1, 1+γ} within tolerance.
3. Measured far-field exponent = 1 − 1/γ within tolerance at R_out/2.
4. Refinement: residual and (γ, λ) estimates stable under N_R, L (or network width/points) doubling.
5. Negative control: with outgoing-geometry constraints imposed, **no** solution should appear at
   γ ≥ 1/2-side scans (CIV exclusion) — if one does, the discretization is broken.

## 7. Outputs and decision tree

- **Converged candidate** (residual → 0 under refinement, anchors pass): freeze (γ, λ, Φ), compute
  the linearized similarity spectrum (v3 §5.2 item 3), count unstable directions, and open the
  certification conversation (HWY-style rigorous numerics). This would be a genuine Type-II candidate.
- **Clean nonexistence signals** (residual floors scaling like a power of resolution, Newton basin
  collapse at a consistent θ*(γ,λ) < 1, obstruction eigenvalue crossing): map θ* over the window and
  formulate the quantitative rigidity conjecture — the positive-direction Liouville target aligned
  with Seregin's ancient-limit class.
- **Inconclusive**: report conditioning numbers and the specific failing anchor. No silent failure.

## 8. Budget and reproducibility

Route C: hours per (γ, λ) point on an A100 at float64; a 5×5 window scan in ~2–4 GPU-days.
Route A: 1–2 weeks of focused development. All seeds, scripts, and grids ship with results —
both directions of this collaboration have now been bitten once by unattached artifacts.

Starter code: `phase2_starter.py` — verifies P2.2/P2.3 and the exact anchor symbolically for
general (γ, λ), extracts the fixed-line spectrum, checks a candidate tail's eigenrelation and
divergence, and assembles the unprojected n = 2 exterior source term (Leray projection of
homogeneous fields is the first real task for the Colab session).
