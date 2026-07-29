# Navier–Stokes breakdown research checkpoint

This repository records a theorem-oriented attempt to resolve the
three-dimensional Navier–Stokes Millennium problem. It is **not** a claimed
solution.

The current target is Clay alternative (D): construct smooth periodic initial
data and a smooth, rapidly time-decaying force for which no global smooth
solution exists. The official formulation permits such a force.

The July 2026 checkpoint makes one strategic distinction:

- a direct Palasek–Córdoba–Martínez-Zoroa–Zheng material-stretching retrofit
  fails for the ordinary Laplacian because the required dormant-seed gain
  creates a super-polynomial active frequency and therefore a fatal principal
  heat cost;
- a bounded-ratio **spectral recurrence** is still exponent-compatible,
  because a genuine unstable mode can grow in amplitude without suffering the
  same material frequency inflation.

The remaining prize-level theorem is a compact renormalized return cell:
amplify a smooth dormant child, transfer energy into it, and return—after
translation, rotation, scaling, and phase—to the same parent/unstable-mode
class with a strictly super-parabolic strain gain.

See
[research/2026-07-29-forced-spectral-recurrence.md](research/2026-07-29-forced-spectral-recurrence.md)
and the follow-up
[research/2026-07-29-return-cell-no-go-map.md](research/2026-07-29-return-cell-no-go-map.md),
which rules out asymptotically finite active templates and narrows the
surviving target to a recurrent leading infinite/localized tail,
and run:

```bash
python checks/forced_recurrence_ledger.py
```

Primary references are linked in the research note. All “pass” statements in
the verifier concern algebraic consistency only; none asserts existence of the
required Navier–Stokes solution.
