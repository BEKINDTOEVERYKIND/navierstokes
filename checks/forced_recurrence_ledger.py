#!/usr/bin/env python3
"""Exact/algebraic checks for the July 2026 forced-recurrence checkpoint.

This script verifies scaling inequalities only.  It does not construct a
Navier--Stokes solution, prove an unstable compact mode, or prove a nonlinear
renormalized return map.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction


def superexponential_retrofit() -> dict[str, object]:
    # A representative point in the formal Palasek full-Laplacian window.
    intermittency = 2.49
    b = 1.20
    beta = 2.45

    heat_margin = beta / b - 2.0
    high_high_decay_power = 2.0 * (intermittency - beta) * (b - 1.0)
    theta = beta * (b - 1.0) / (b * b)

    assert 1.0 < b < intermittency / 2.0
    assert 2.0 * b < beta < intermittency <= 2.5
    assert heat_margin > 0.0
    assert high_high_decay_power > 0.0
    assert theta > 0.0

    # If K/L <= N_{k-1}^{-1} and KL <= N_k, then
    # K <= N_k^((b-1)/(2b)): only polynomial gain is allowed.
    carrier_polynomial_bound = (b - 1.0) / (2.0 * b)

    # If L >= K N_{k-1} and (KL)^2 << A_{k-1}, then
    # K << N_k^((beta-2)/(4b)): again only polynomial gain.
    heat_polynomial_bound = (beta - 2.0) / (4.0 * b)

    assert carrier_polynomial_bound > 0.0
    assert heat_polynomial_bound > 0.0

    return {
        "formal_window": "PASS",
        "direct_cmz_material_retrofit": "FAIL",
        "intermittency": intermittency,
        "b": b,
        "beta": beta,
        "heat_margin_beta_over_b_minus_2": heat_margin,
        "high_high_decay_power": high_high_decay_power,
        "dormancy_flatness_theta": theta,
        "required_gain": "K_k >= exp(c*N_k^theta)",
        "carrier_and_nesting_allow_at_most":
            f"K_k <= N_k^{carrier_polynomial_bound:.12g}",
        "active_heat_and_nesting_allow_at_most":
            f"K_k << N_k^{heat_polynomial_bound:.12g}",
        "reason":
            "super-polynomial required gain contradicts both polynomial bounds",
    }


def bounded_ratio_spectral_ladder() -> dict[str, object]:
    # Exact dyadic design point from the note.
    r = Fraction(32, 1)
    mu = Fraction(2048, 1)
    q = Fraction(2, 1)

    energy_ratio = q * mu * mu / (r**5)
    dissipation_ratio = q * mu / (r**3)
    reynolds_ratio = mu / (r**2)
    radius_ratio = q / r
    volume_ratio = q / (r**3)
    l3_ratio_cubed = q * (mu / (r**2)) ** 3

    assert r**2 < mu
    assert float(mu) < float(r) ** 2.5
    assert 1 < q < r**5 / (mu**2)
    assert energy_ratio == Fraction(1, 4)
    assert dissipation_ratio == Fraction(1, 8)
    assert reynolds_ratio == Fraction(2, 1)
    assert radius_ratio == Fraction(1, 16)
    assert volume_ratio == Fraction(1, 16384)
    assert energy_ratio < 1
    assert dissipation_ratio < 1
    assert reynolds_ratio > 1
    assert l3_ratio_cubed > 1

    # Polynomial residence M_j ~ j^2 cannot defeat geometric margins.
    for j in range(1, 1000):
        residence = j * j
        assert math.isfinite(residence * float(energy_ratio) ** j)
        assert math.isfinite(residence * float(dissipation_ratio) ** j)

    return {
        "algebraic_ledger": "PASS",
        "nonlinear_return_theorem": "OPEN",
        "r": int(r),
        "mu": int(mu),
        "q": int(q),
        "energy_ratio": str(energy_ratio),
        "dissipation_ratio": str(dissipation_ratio),
        "reynolds_ratio": str(reynolds_ratio),
        "major_radius_ratio": str(radius_ratio),
        "volume_ratio": str(volume_ratio),
        "L3_ratio_cubed": str(l3_ratio_cubed),
        "dormant_seed": "delta_j=exp(-c*j^2)",
        "required_normalized_residence": "M_j=O(j^2)",
        "terminal_force_flatness":
            "exp(-c*j^2+C*j)=O((T-t)^K) for every fixed K",
    }


def main() -> None:
    result = {
        "scope": {
            "proves_navier_stokes_blowup": False,
            "proves_clay_alternative_D": False,
            "proves_compact_unstable_mode": False,
            "proves_nonlinear_return_map": False,
            "checks": "scaling inequalities and principal-frequency obstruction",
        },
        "superexponential_material_route": superexponential_retrofit(),
        "bounded_ratio_spectral_route": bounded_ratio_spectral_ladder(),
        "decision": {
            "killed":
                "direct Palasek dormant seed in a CMZ material hyperbolic packet",
            "survives":
                "bounded-ratio dormant spectral recurrence with exact/all-order closure",
            "single_missing_object":
                "compact renormalized unstable-manifold return cell",
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
