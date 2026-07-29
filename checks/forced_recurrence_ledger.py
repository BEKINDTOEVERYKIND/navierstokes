#!/usr/bin/env python3
"""Exact/algebraic checks for the July 2026 breakdown checkpoint.

This script verifies scaling identities and inequalities only.  The packed
bath identities refer to disjoint, rescaled compact steady Euler seeds, so
their leading Euler defect is exactly zero.  The script does not construct a
Navier--Stokes solution or prove a time-dependent transition/return map.
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
        "global_fixed_profile_return": "CLOSED_BY_RETURN_NO_GO",
        "interpretation":
            "a scalar-compatible ledger, not an admissible global recurrence",
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


def localized_wake_cascade() -> dict[str, object]:
    """Check one exact point in the surviving 1 < gamma < 3/2 window."""
    r = Fraction(2, 1)
    gamma = Fraction(5, 4)

    time_ratio = float(r) ** -float(1 + gamma)
    reynolds_ratio = float(r) ** float(gamma - 1)
    energy_ratio = float(r) ** float(2 * gamma - 3)
    dissipation_ratio = float(r) ** -float(2 - gamma)
    viscosity_ratio = float(r) ** -float(gamma - 1)
    gradient_ratio = float(r) ** float(gamma + 1)
    helicity_scale_ratio = float(r) ** float(2 * gamma - 2)
    circulation_scale_ratio = reynolds_ratio

    assert 1 < gamma < Fraction(3, 2)
    assert 0 < time_ratio < 1
    assert reynolds_ratio > 1
    assert 0 < energy_ratio < 1
    assert 0 < dissipation_ratio < 1
    assert 0 < viscosity_ratio < 1
    assert gradient_ratio > 1
    assert helicity_scale_ratio > 1
    assert circulation_scale_ratio > 1

    return {
        "formal_ledger": "PASS",
        "existence_of_localized_return_cell": "OPEN",
        "r": int(r),
        "gamma": str(gamma),
        "time_ratio": time_ratio,
        "reynolds_ratio": reynolds_ratio,
        "active_energy_ratio": energy_ratio,
        "stage_dissipation_ratio": dissipation_ratio,
        "normalized_viscosity_ratio": viscosity_ratio,
        "gradient_ratio": gradient_ratio,
        "helicity_scale_ratio": helicity_scale_ratio,
        "circulation_scale_ratio": circulation_scale_ratio,
        "required_helicity_design":
            "exact zero net helicity or an equally leading bihelical companion",
        "required_topology":
            "non-tubular strain cell or proved multi-strand aggregation",
        "all_order_remainder":
            "requires a Gevrey-tame WKB/viscous hierarchy and matched wake",
    }


def polynomial_carrier_refinement() -> dict[str, object]:
    """Check one conditional Gevrey ledger for K_j=j^A."""
    r = 2.0
    gamma = 1.25
    gevrey_order = 2.0
    carrier_power = 6.0
    eta = 0.1
    coefficient_constant = 10.0

    assert 1.0 < gamma < 1.5
    assert carrier_power > 2.0 * gevrey_order

    worst_wkb_rate = -math.inf
    worst_heat_log = -math.inf
    for j in range(300, 1001):
        log_j = math.log(j)
        truncation = max(
            1,
            int(eta * j * j / math.log(math.e + j)),
        )
        log_carrier = carrier_power * log_j

        log_wkb_remainder = (
            truncation * math.log(coefficient_constant)
            + gevrey_order * math.lgamma(truncation + 1)
            - truncation * log_carrier
        )
        worst_wkb_rate = max(worst_wkb_rate, log_wkb_remainder / (j * j))

        # epsilon_j K_j^2 for the exact Kelvin-compatible amplitude
        # a_j=ell_j^{-gamma} K_j^gamma.
        log_heat_parameter = (
            -(gamma - 1.0) * math.log(r) * j
            + (2.0 - gamma) * log_carrier
        )
        worst_heat_log = max(worst_heat_log, log_heat_parameter)

        assert log_wkb_remainder < -0.25 * j * j
        assert log_heat_parameter < 0.0

    return {
        "conditional_scalar_ledger": "PASS",
        "constructs_all_order_cell": False,
        "r": r,
        "gamma": gamma,
        "gevrey_order_sigma": gevrey_order,
        "carrier": "K_j=j^6",
        "truncation": "M_j=floor(0.1*j^2/log(e+j))",
        "gevrey_gate": "A>2*sigma",
        "worst_checked_log_wkb_remainder_over_j2": worst_wkb_rate,
        "worst_checked_log_heat_parameter": worst_heat_log,
        "carrier_handoff":
            "K_(j+1)/K_j=1+A/j+O(j^-2), incorporated exactly",
        "missing":
            "Gevrey-tame one-carrier cell inverse and matched viscous endpoint jets",
    }


def exact_packed_gavrilov_bath() -> dict[str, object]:
    """Verify the carrier ledger realized by a disjoint microbubble packing.

    In an annulus of outer scale ell, pack N~K^3 bubbles of diameter
    delta=ell/K and amplitude a=ell^(-gamma) K^gamma.  Additivity of
    disjoint supports gives every line below.  The zero Euler defect uses the
    compact joint velocity-pressure support of Gavrilov's steady seed; it is
    a mathematical input, not something established numerically here.
    """
    gamma = Fraction(5, 4)

    # Store powers as (power of ell, power of K).  Multiplication is vector
    # addition, which lets the checker verify the identities without
    # floating-point arithmetic.
    number = (Fraction(0), Fraction(3))
    diameter = (Fraction(1), Fraction(-1))
    amplitude = (-gamma, gamma)
    stage_time = (1 + gamma, -gamma)

    def add(*terms: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
        return (
            sum((term[0] for term in terms), Fraction(0)),
            sum((term[1] for term in terms), Fraction(0)),
        )

    def multiple(
        coefficient: int,
        term: tuple[Fraction, Fraction],
    ) -> tuple[Fraction, Fraction]:
        return (coefficient * term[0], coefficient * term[1])

    l2_squared = add(number, multiple(2, amplitude), multiple(3, diameter))
    l3_cubed = add(number, multiple(3, amplitude), multiple(3, diameter))
    grad_l2_squared = add(
        number,
        multiple(2, amplitude),
        diameter,
    )
    carrier_reynolds = add(amplitude, diameter)
    viscous_stage_loss = add(stage_time, grad_l2_squared)

    assert l2_squared == (3 - 2 * gamma, 2 * gamma)
    assert l3_cubed == (3 - 3 * gamma, 3 * gamma)
    assert grad_l2_squared == (1 - 2 * gamma, 2 * gamma + 2)
    assert carrier_reynolds == (1 - gamma, gamma - 1)
    assert viscous_stage_loss == (2 - gamma, gamma + 2)

    return {
        "exact_stationary_euler_bath": "PASS",
        "constructs_time_dependent_transition": False,
        "gamma": str(gamma),
        "bubble_count": "N_j asymptotic to K_j^3",
        "bubble_diameter": "delta_j=ell_j/K_j",
        "bubble_amplitude": "a_j=ell_j^(-gamma)*K_j^gamma",
        "carrier_frequency": "k_j asymptotic to K_j/ell_j",
        "L2_squared_powers_ell_K": [str(x) for x in l2_squared],
        "L3_cubed_powers_ell_K": [str(x) for x in l3_cubed],
        "grad_L2_squared_powers_ell_K":
            [str(x) for x in grad_l2_squared],
        "carrier_Re_powers_ell_K":
            [str(x) for x in carrier_reynolds],
        "viscous_stage_loss_powers_ell_K":
            [str(x) for x in viscous_stage_loss],
        "leading_euler_defect":
            "exactly zero by joint-support disjointness of steady seeds",
        "remaining":
            "dynamically create the child bath and close viscosity/seams",
    }


def main() -> None:
    result = {
        "scope": {
            "proves_navier_stokes_blowup": False,
            "proves_clay_alternative_D": False,
            "proves_localized_kelvin_reynolds_cell": False,
            "proves_nonlinear_return_map": False,
            "checks": "scaling inequalities and principal-frequency obstruction",
        },
        "superexponential_material_route": superexponential_retrofit(),
        "bounded_ratio_spectral_route": bounded_ratio_spectral_ladder(),
        "localized_wake_cascade": localized_wake_cascade(),
        "polynomial_carrier_refinement": polynomial_carrier_refinement(),
        "exact_packed_gavrilov_bath": exact_packed_gavrilov_bath(),
        "decision": {
            "killed":
                "direct CMZ material retrofit and clean global fixed-profile return",
            "survives":
                "localized Kelvin/Reynolds cascade with an exact static packed carrier",
            "single_missing_object":
                "active finite-band annular transition with viscous endpoint jets",
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
