#!/usr/bin/env python3
"""Exact arithmetic/scaling checks for nonlinear periodic wake carry.

This script checks only identities in
research/2026-08-03-nonlinear-periodic-wake-carry.md.  It does not verify
the missing weighted nonlinear endpoint estimate.
"""

from __future__ import annotations

from fractions import Fraction
from math import exp, isclose


def term_ledger() -> None:
    """Powers are (wake amplitude b, frequency N, viscosity mu)."""

    dt_z = (1, 0, 0)
    transport = (1, 1, 0)
    stretch = (1, 0, 0)
    self_advection = (2, 1, 0)
    heat = (1, 2, 1)

    assert dt_z == stretch
    assert transport == (1, 1, 0)
    assert self_advection == (2, 1, 0)
    assert heat == (1, 2, 1)

    # Substituting N=1, K, and M*K gives mu, mu*K^2, mu*(M*K)^2.
    assert heat[2] == 1 and heat[1] == 2
    print("one-turnover linear/nonlinear/heat power ledger: PASS")


def heat_no_go() -> None:
    stage_length = Fraction(3, 2)
    values = [exp(-float(stage_length) * 10.0 ** (-n)) for n in range(1, 8)]
    assert all(0.0 < value < 1.0 for value in values)
    assert all(values[n] < values[n + 1] for n in range(len(values) - 1))
    assert isclose(values[-1], 1.0, rel_tol=0.0, abs_tol=2e-7)
    print("lowest-mode heat multiplier tends to one as Re tends to infinity: PASS")


def exponentially_grown_source() -> None:
    growth = 0.7
    for gain_length in (1.0, 4.0, 25.0, 100.0):
        integral = (1.0 - exp(-2.0 * growth * gain_length)) / (2.0 * growth)
        assert 0.0 < integral <= 1.0 / (2.0 * growth)
    print("terminal exponential-ramp quadratic source is uniformly integrable: PASS")


def dilation_identity() -> None:
    # D_j f(y)=q^{-gamma} f(y/r).  In the |y|^(alpha+m) derivative
    # seminorm, r^(-m) from differentiation cancels r^m from the weight.
    r = Fraction(5, 1)
    gamma = Fraction(5, 4)
    alpha = Fraction(1, 1)
    k_ratio = Fraction(11, 10)  # K_{j+1}/K_j

    # Use rational powers only through the fourth power (gamma=5/4).
    factor_fourth = (r ** (4 * alpha)) / ((r * k_ratio) ** (4 * gamma))
    expected_fourth = (r ** (4 * (alpha - gamma))) / (k_ratio ** (4 * gamma))
    assert factor_fourth == expected_fourth
    assert factor_fourth < 1

    # The factor is independent of derivative order m.
    for derivative_order in range(8):
        weight_power = alpha + derivative_order
        combined_r_power = weight_power - derivative_order
        assert combined_r_power == alpha

    print("weighted terminal-dilation identity and subcritical contraction: PASS")


def critical_telescoping() -> None:
    # At alpha=gamma, each step is (K_j/K_{j+1})^gamma.
    # Raising to the fourth power avoids floating arithmetic for gamma=5/4.
    gamma_times_four = 5
    carrier_power = 8
    i, j = 3, 29

    product_fourth = Fraction(1, 1)
    for n in range(i, j):
        k_n = n**carrier_power
        k_next = (n + 1) ** carrier_power
        product_fourth *= Fraction(k_n, k_next) ** gamma_times_four

    expected_fourth = Fraction(i, j) ** (carrier_power * gamma_times_four)
    assert product_fourth == expected_fourth
    print("critical-weight carrier factors telescope and are not uniform: PASS")


def l2_dilation() -> None:
    gamma = Fraction(5, 4)
    exponent = Fraction(3, 2) - gamma
    assert exponent == Fraction(1, 4) > 0
    print("renormalized L2 shell exponent 3/2-gamma is positive: PASS")


def nonlinear_contraction_gate() -> None:
    r = 16.0
    gamma = 1.25
    alpha = 0.75
    endpoint_lipschitz = 2.0
    carried_lipschitz = endpoint_lipschitz * r ** (alpha - gamma)
    assert isclose(carried_lipschitz, 0.5)
    assert carried_lipschitz < 1.0

    # No fixed r can beat an endpoint bound exp(C*j^2).
    fixed_scale_gain = r ** (alpha - gamma)
    long_window = [fixed_scale_gain * exp(0.1 * j * j) for j in range(1, 20)]
    assert long_window[-1] > 1.0
    assert all(long_window[n] < long_window[n + 1] for n in range(len(long_window) - 1))
    print("scale-shift contraction criterion and long-window obstruction: PASS")


def main() -> None:
    term_ledger()
    heat_no_go()
    exponentially_grown_source()
    dilation_identity()
    critical_telescoping()
    l2_dilation()
    nonlinear_contraction_gate()
    print("all nonlinear periodic wake-carry checks passed")


if __name__ == "__main__":
    main()
