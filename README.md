# Navier–Stokes breakdown research checkpoint

This repository records a theorem-oriented attempt to resolve the
three-dimensional Navier–Stokes Millennium problem. It is **not** a claimed
solution.

The current target is Clay alternative (D): construct smooth periodic initial
data and a smooth, rapidly time-decaying force for which no global smooth
solution exists. The official formulation permits such a force.

The July 2026 checkpoint now makes a sharper strategic distinction:

- a direct Palasek–Córdoba–Martínez-Zoroa–Zheng material-stretching retrofit
  fails for the ordinary Laplacian because the required dormant-seed gain
  creates a super-polynomial active frequency and therefore a fatal principal
  heat cost;
- bounded-ratio **global spectral recurrence** passes the scalar exponent
  ledger but fails as a clean return: finite active templates degenerate,
  scaled torus copies conflict with helicity/sublattice invariance, and
  volume-filling copies eventually lose to viscosity;
- a localized, non-precompact wake cascade remains formally compatible in
  the sharp window
  \(\ell_j=r^{-j}\), \(a_j=\ell_j^{-\gamma}\),
  \(1<\gamma<3/2\).

The remaining prize-level theorem is a localized Kelvin--Reynolds return
cell.  An exact affine Kelvin wave proves that true Euler strain can realize
the required amplitude-frequency law \(A\sim k^\gamma\), with only
\(O(\nu/(a_j\ell_j))\) full-Laplacian damping per bounded-ratio stage.  What
is not proved is the hard half: localize the affine reservoir, drain the old
parent through a multiwave Reynolds stress, create the smaller parent and
fresh seed, close pressure/wakes/helicity/circulation, and obtain a tame
all-order viscous return map.

See
[research/2026-07-29-forced-spectral-recurrence.md](research/2026-07-29-forced-spectral-recurrence.md)
and the follow-up
[research/2026-07-29-return-cell-no-go-map.md](research/2026-07-29-return-cell-no-go-map.md),
which rules out asymptotically finite active templates and narrows the
surviving target to a recurrent leading infinite/localized tail;
[research/2026-07-29-global-scaled-return-obstruction.md](research/2026-07-29-global-scaled-return-obstruction.md),
which closes clean volume-filling scaled copies; and
[research/2026-07-29-flat-force-selfsimilar-audit.md](research/2026-07-29-flat-force-selfsimilar-audit.md),
which gives the exact surviving exponent window, the Kelvin amplifier, and
the conditional all-order flat-force theorem.  A
separate global-regularity audit,
[research/2026-07-29-quadratic-lyapunov-no-go.md](research/2026-07-29-quadratic-lyapunov-no-go.md),
proves that no positive Fourier-diagonal quadratic Lyapunov functional
stronger than energy can work for the full equation.

To reproduce the scalar recurrence ledger, install no extra dependencies
and run:

```bash
python checks/forced_recurrence_ledger.py
```

Primary references are linked in the research note. All “pass” statements in
the verifier concern algebraic consistency only; none asserts existence of the
required Navier–Stokes solution.

## Audit workspace (added by the auditing model)

`CLAIMS.md` is the registry of every substantive claim in this program with
its verification status. Protocol: nothing is "established" until the party
that did not produce it re-derives or reproduces it from artifacts here.

- `docs/` — the working notes (source + rendered) and per-session deltas;
  Addenda A–H.1 record the full audited history including the corrected
  metric-gaming episodes.
- `engines/` — the independent solvers and optimizers (complex128 replay
  engine, v2.3 standard-norm optimizer, r0-scan).
- `gate/` — the downshifter-gate bundle, its A100 results, and the
  independent float64 validator.
- `audit/` — symbolic and numerical audit scripts (certificate window
  generalization, metric recomputation, null baselines).
- `data/stage5/` — the endpoint-crossing evidence base (Addendum G.3).
- `phase2/` — the affine-core continuation spec (open since Addendum E).
