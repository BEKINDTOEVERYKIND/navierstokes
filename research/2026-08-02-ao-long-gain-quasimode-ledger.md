# Albritton--Ożański edge deficit and a long-gain quasimode ledger

**Date:** 2026-08-02
**Status:** source-anchored asymptotic and conditional semigroup lemma;
self-derived, not cross-audited
**Scope:** linear carrier survival through a polynomial action window.  This
note does not construct the finite-curvature Gavrilov mode, the nonlinear
transition, or a Navier--Stokes singularity.

## 1. The missing AO constant is not arbitrary

Here `AO` means the ring modes of Albritton and Ożański.  In their notation,
Theorem 1.1 and equations (1.28), (1.31) give, for every fixed Weber index
`m`,

\[
 \omega_{m,n}=n\Lambda _0+i\sqrt{b_0}
 +(1-i)c_m n^{-1/2}+O(n^{-1+\delta}),
 \qquad c_m>0.                                      \tag{1.1}
\]

The time dependence is `exp(-i omega t)`, so the growth rate is

\[
 \lambda_{m,n}=\operatorname{Im}\omega_{m,n}
 =\sqrt{b_0}-c_m n^{-1/2}+O(n^{-1+\delta}).          \tag{1.2}
\]

The leading local centrifugal/WKB rate at the concentrating ring is
`sqrt(b_0)`.  Consequently the ratio to that **local leading rate** obeys

\[
 {\lambda_{m,n}\over\sqrt{b_0}}
 =1-{c_m\over\sqrt{b_0}}n^{-1/2}+O(n^{-1+\delta})
 \longrightarrow1.                                  \tag{1.3}
\]

Thus the numerical ratio itself is not the obstruction in C50: relative to
the same leading edge, every fixed threshold below one, including `7/9`, is
eventually cleared.

For a resonant helical bicharacteristic with radial covector `ell`, the
column BAS reduction gives

\[
 \sigma^2(\ell;r)
 =b(r){n^2/r^2+\alpha^2\over
          \ell^2+n^2/r^2+\alpha^2}
 \le b(r).                                           \tag{1.4}
\]

Equality holds at `ell=0`.  This identifies `sqrt(b(r_0))` with the local
resonant BAS rate at the AO concentration radius.  For a generic Assumption-A
profile it does **not** prove that `r_0` maximizes `b(r)` over the whole
column.  The explicit Batchelor profile in
`2026-08-02-ao-batchelor-global-bas-certificate.md` now supplies that
global `b` certificate in the selected fixed-`beta` sector.  The full
reduction in `2026-08-02-ao-batchelor-full-bas-cocycle.md` proves that
non-resonant exponents vanish, but also exhibits a slightly faster exact
resonance at a different helical ratio.  Thus this first explicit profile
does not identify the unrestricted edge needed on the `j^2` window.  The
repaired construction in
`2026-08-02-ao-batchelor-full-edge-matched-profile.md` simultaneously
maximizes in radius and helical ratio and closes this principal-symbol
profile-selection issue.

There is a crucial qualification.  Albritton--Ożański explicitly note after
their equation (3.32) that the semigroup/resolvent constant can depend on
the Fourier mode `(n,alpha)`.  They therefore cannot identify the supremum
of the mode spectral bounds with their global bound.  Equation (1.3) must
not be promoted to

\[
 \lambda_{m,n}/\Lambda_{\rm full}\to1
\]

without a uniform sector estimate.  The unfinished thin-Gavrilov-ring
calculation has exactly the same obligation.

Primary source: D. Albritton and W. Ożański, *Linear and nonlinear
instability of vortex columns*, arXiv:2310.20674v3, equations (1.28),
(1.31), and the warning following (3.32).

## 2. Finite-action quasimode estimate

The nonlinear-instability exit time in C50 is only logarithmic in frequency.
The cascade needs a longer action, of order `G_j ~ j^2`, to amplify a seed of
size `exp(-G_j)`.  The relevant estimate is still elementary, but it must
retain the edge deficit.

Let `S_M(t)=exp(tL_M)` act on a Banach space `Z`.  Assume

\[
 \|S_M(t)\|_{Z\to Z}
 \le C(1+t)^d e^{\Lambda_Mt},                       \tag{2.1}
\]

with `C,d` independent of `M`, and let

\[
 \|\phi_M\|_Z=1,\qquad
 \|(L_M-\lambda_M)\phi_M\|_Z\le\varepsilon_M,
 \qquad \Delta_M:=\Lambda_M-\lambda_M\ge0.          \tag{2.2}
\]

Variation of constants gives

\[
 \begin{aligned}
 \|e^{-\lambda_Mt}S_M(t)\phi_M-\phi_M\|_Z
 &\le C\varepsilon_M
   \int_0^t(1+t-s)^d e^{\Delta_M(t-s)}\,ds\\
 &\le C\varepsilon_M t(1+t)^d e^{\Delta_Mt}.
                                                               \tag{2.3}
 \end{aligned}
\]

This is a relative error.  Multiplying the incoming mode by `exp(-G)` also
multiplies the defect by `exp(-G)`; unlike an arbitrary additive forcing, it
cannot cancel the exponentially small seed unless the relative error becomes
order one.

Assume, for some `q,r>0`,

\[
 \varepsilon_M\lesssim M^{-q},\qquad
 \Delta_M\lesssim M^{-r},                             \tag{2.4}
\]

and choose

\[
 G_j\asymp j^g,\qquad M_j=j^A,\qquad
 T_j\asymp G_j/\lambda_{M_j},qquad
 \inf_j\lambda_{M_j}>0.                               \tag{2.5}
\]

Then (2.3) becomes

\[
 E_j\lesssim
 j^{-Aq+g(d+1)}
 \exp\!\left(Cj^{g-Ar}\right).                        \tag{2.6}
\]

Hence

\[
 \boxed{Ar\ge g,\qquad Aq>g(d+1)}                    \tag{2.7}
\]

is sufficient for `E_j -> 0`.  Equality in the first inequality leaves the
edge exponential bounded by `exp(C)`; decay still comes from the second
inequality.  The errors are summable if the second strict
inequality is strengthened to

\[
                 Aq>g(d+1)+1.                         \tag{2.8}
\]

The polynomial prefactor in (2.1) is not cosmetic.  A generic estimate
`C_delta exp((Lambda+delta)t)` with an uncontrolled `C_delta` does not imply
(2.1) uniformly on the growing window.

## 3. The AO/curvature design point

For the AO asymptotic,

\[
 r={1\over2}.                                         \tag{3.1}
\]

The recovered thin-torus calculation reported a normalized curvature
residual `O(M^-1)+exp(-cM)`, corresponding conditionally to

\[
 q=1.                                                 \tag{3.2}
\]

For the cascade action `g=2` and a uniformly diagonalizable Floquet sector
(`d=0`), (2.7) reads

\[
 {A\over2}\ge2,\qquad A>2.
\]

Thus

\[
                         \boxed{A\ge4}                \tag{3.3}
\]

controls both effects and makes the relative errors summable: at the
endpoint the estimate is `E_j=O(j^-2)`.  At this
design point the edge deficit, not the `O(M^-1)` curvature residual, is the
limiting exponent.  The separate Gevrey-2 polynomial-carrier ledger C41
imposes the strict condition `A>2 sigma=4`, so the integrated cascade still
selects `A>4` even though this long-gain estimate includes the endpoint.

If (2.1) carries a factor `(1+t)^d`, the sufficient condition becomes

\[
 A\ge4,\qquad A>2(d+1)                               \tag{3.4}
\]

for convergence, and

\[
 A\ge4,\qquad A>2(d+1)+1                             \tag{3.5}
\]

for absolute summability.

## 4. Connection to viscous capture

The robust affine-capture theorem on the remote continuation branch begins
at a fixed incoming section and distinguishes multiplicative carrier error
from arbitrary additive error.  Estimate (2.3) supplies exactly the former
kind during the gain interval, provided the finite-curvature block admits
(2.1)--(2.4).  Once the carrier reaches the fixed section, the capture
theorem only needs a small `L^1` perturbation on a compact transition.  If
the residual remains proportional to the decaying carrier during cleanup,
the terminal spectral gap preserves an exponentially small exit packet.

This removes one reason to demand an exact isolated AO eigenvalue: a
first-order, edge-matched quasimode can survive the full `j^2` gain window at
`A=4`; the independent Gevrey-2 ledger asks for `A>4` in the complete
cascade.

It does **not** remove the following obligations.

1. The principal-symbol base-profile gate is now closed by
   `2026-08-02-locked-pitch-gavrilov-edge-profile.md`: an explicit compact
   pressure modulation in the actual `W_G=V_G/sqrt(2)` family has a unique
   full BAS edge.  Its q-free transfer of the AO eigenmode proof is
   self-audited and requires independent verification; a direct local Weber
   quasimode is the fallback.
2. For that compatible profile, construct the divergence-free
   finite-curvature ring packet and prove the normalized `q=1` residual in
   the norm used by (2.1).  Curvature is now the first geometric mismatch,
   but the claimed residual still has to be derived.
3. Prove a Fourier/microlocal bound (2.1) for the compact locked-pitch
   profile, including a mode-uniform polynomial prefactor and the
   `Delta_M=O(M^-1/2)` comparison to the **full relevant** growth edge.  The
   generalized C66 BAS argument gives `d=0` only at principal-symbol level;
   AO's published PDE theorem does not provide the required uniformity.
4. Pass from the complex AO harmonic to a real physical perturbation with a
   uniform lower norm throughout the gain window.  A single straight-column
   harmonic should inherit this from helical translation symmetry, but a
   curved Floquet mode with sidebands requires a non-cancellation proof.
5. Couple at least the two real phases needed to write the child and retain
   their sideband, collar, pressure, and viscous wakes.
6. Close the nonlinear endpoint map and the exact/terminally-flat
   Navier--Stokes residual.  Linear carrier survival alone does not do so.

The finite-frequency transverse-Weber/global-resolvent theorem is therefore
replaced by a sharper target, not declared solved: a uniform polynomial-
prefactor sector propagator plus a first-order finite-curvature quasimode.
