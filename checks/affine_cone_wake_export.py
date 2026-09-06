#!/usr/bin/env python3
"""Exact algebra/scaling checks for affine-cone wake export.

This script checks identities in
research/2026-08-03-affine-cone-wake-export.md.  It does not verify the
Gaussian kernel estimate, nonlinear cone stability, or a localized
finite-energy export construction.
"""

from __future__ import annotations

from fractions import Fraction
from math import exp, isclose, log


def trace_and_characteristics() -> None:
    alpha = Fraction(1, 1)
    beta = Fraction(5, 4)
    gamma = alpha + beta
    rates = (-alpha, -beta, gamma)

    assert sum(rates, Fraction(0, 1)) == 0

    # Logarithmic determinant of exp(tS) is t*trace(S).
    time = Fraction(7, 3)
    log_det = time * sum(rates, Fraction(0, 1))
    assert log_det == 0

    # Cone slopes x_i/x_3 decay at rates alpha+gamma and beta+gamma.
    assert rates[0] - rates[2] == -(alpha + gamma)
    assert rates[1] - rates[2] == -(beta + gamma)
    print("trace-free characteristic and cone-narrowing ledger: PASS")


def l1_biot_savart_exponents() -> None:
    gamma = Fraction(9, 4)

    # L1 vorticity costs +gamma.  The m-th velocity derivative of the
    # Biot-Savart kernel at exported distance costs -(2+m)*gamma.
    for derivative_order in range(9):
        total = gamma - (2 + derivative_order) * gamma
        expected = -(1 + derivative_order) * gamma
        assert total == expected
        assert total < 0
    print("Biot-Savart distance beats affine L1 stretching: PASS")


def viscous_covariance_identities() -> None:
    alpha = 1.0
    beta = 1.25
    gamma = alpha + beta

    for time in (1.0e-4, 0.03, 0.7, 5.0):
        a1 = (exp(2.0 * alpha * time) - 1.0) / (2.0 * alpha)
        a2 = (exp(2.0 * beta * time) - 1.0) / (2.0 * beta)
        a3 = (1.0 - exp(-2.0 * gamma * time)) / (2.0 * gamma)

        b1 = (1.0 - exp(-2.0 * alpha * time)) / (2.0 * alpha)
        b2 = (1.0 - exp(-2.0 * beta * time)) / (2.0 * beta)
        b3 = (exp(2.0 * gamma * time) - 1.0) / (2.0 * gamma)

        assert isclose(exp(-2.0 * alpha * time) * a1, b1, rel_tol=1e-11)
        assert isclose(exp(-2.0 * beta * time) * a2, b2, rel_tol=1e-11)
        assert isclose(exp(2.0 * gamma * time) * a3, b3, rel_tol=1e-11)

        # In the outgoing coordinate, physical variance divided by the
        # square of the affine expansion is uniformly bounded.
        normalized_outgoing_variance = exp(-2.0 * gamma * time) * 2.0 * b3
        expected = (1.0 - exp(-2.0 * gamma * time)) / gamma
        assert isclose(normalized_outgoing_variance, expected, rel_tol=1e-11)
        assert normalized_outgoing_variance < 1.0 / gamma

        assert 2.0 * b1 < 1.0 / alpha
        assert 2.0 * b2 < 1.0 / beta

    print("anisotropic material/physical Gaussian covariance identities: PASS")


def window_uniformity() -> None:
    gamma = 2.25
    radius = 8.0

    for derivative_order in range(6):
        values = [
            radius ** (-2 - derivative_order)
            * exp(-(1 + derivative_order) * gamma * time)
            for time in (0.0, 1.0, 10.0, 100.0)
        ]
        assert all(values[i] >= values[i + 1] for i in range(len(values) - 1))
        assert values[0] == radius ** (-2 - derivative_order)

    print("incoming-wake endpoint bound is independent of gain-window length: PASS")


def noncompact_shell_summability() -> None:
    # A velocity tail D^{-delta} has shell vorticity L1 size D^{2-delta}.
    # The core Biot-Savart weight D^{-2-m} therefore leaves D^{-delta-m}.
    delta = Fraction(5, 4)
    shell_ratio = Fraction(4, 1)

    for derivative_order in range(6):
        vorticity_l1_power = 2 - delta
        kernel_power = -2 - derivative_order
        contribution_power = vorticity_l1_power + kernel_power
        assert contribution_power == -delta - derivative_order < 0

        ratio = float(shell_ratio) ** float(contribution_power)
        assert 0.0 < ratio < 1.0

    print("noncompact dyadic wake is summable in the core-weighted norm: PASS")


def nonlinear_bootstrap_form() -> None:
    # If the nonlinear self-strain integral is A, the vorticity bound loses
    # exp(A), but it does not acquire the gain-window length separately.
    self_strain_integral = Fraction(3, 5)
    assert self_strain_integral >= 0

    # Material outgoing displacement remains positive under a weighted
    # perturbation smaller than half the initial separation.
    initial_separation = Fraction(10, 1)
    weighted_displacement = Fraction(4, 1)
    assert weighted_displacement < initial_separation / 2
    print("conditional nonlinear cone-bootstrap inequalities: PASS")


def finite_energy_displacement_gate() -> None:
    # Endpoint inequality:
    # V*D <= initial first moment + sqrt(V)*integrated L2 norm.
    volume = Fraction(1, 64)
    sqrt_volume = Fraction(1, 8)
    distance = Fraction(64, 1)
    initial_first_moment = Fraction(0, 1)
    integrated_l2 = Fraction(8, 1)
    lhs = volume * distance
    rhs = initial_first_moment + sqrt_volume * integrated_l2
    assert sqrt_volume * sqrt_volume == volume
    assert lhs == rhs

    # Solving sqrt(V)*R*exp(gamma*T) <= E*T gives the squared-volume
    # exponent -2*gamma*T and polynomial factor T^2.
    gamma = Fraction(9, 4)
    exponential_power_before_squaring = -gamma
    exponential_power_after_squaring = 2 * exponential_power_before_squaring
    time_power_before_squaring = 1
    time_power_after_squaring = 2 * time_power_before_squaring
    assert exponential_power_after_squaring == -2 * gamma
    assert time_power_after_squaring == 2

    # For G_j=c*j^2, a geometric-in-j volume exp(-C*j) is eventually much
    # larger than the necessary poly(j)*exp(-2*gamma*c*j^2) upper bound.
    c_gain = 0.5
    geometric_rate = 6.0
    comparisons = []
    for stage in (5, 10, 20, 40):
        log_geometric_volume = -geometric_rate * stage
        log_required_upper_bound = (
            4.0 * log(stage) - 2.0 * float(gamma) * c_gain * stage * stage
        )
        comparisons.append(log_geometric_volume > log_required_upper_bound)
    assert comparisons[-2:] == [True, True]
    print("finite-energy material-displacement and j^2-window gate: PASS")


def main() -> None:
    trace_and_characteristics()
    l1_biot_savart_exponents()
    viscous_covariance_identities()
    window_uniformity()
    noncompact_shell_summability()
    nonlinear_bootstrap_form()
    finite_energy_displacement_gate()
    print("all affine-cone wake-export checks passed")


if __name__ == "__main__":
    main()
