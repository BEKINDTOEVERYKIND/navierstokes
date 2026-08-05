#!/usr/bin/env python3
"""Arithmetic/combinatorial checks for recovered claims C132--C136.

No output of this script asserts existence of the open BAFL Navier--Stokes
stage map.
"""

from fractions import Fraction
import itertools
import math


Vector = tuple[int, int, int]


def add(p: Vector, q: Vector) -> Vector:
    return (p[0] + q[0], p[1] + q[1], p[2] + q[2])


def neg(p: Vector) -> Vector:
    return (-p[0], -p[1], -p[2])


def focus_and_activation_checks() -> None:
    for n in range(2, 100):
        q = n**8
        b = Fraction(1, n**2)
        focus = n**12
        gain = n**10
        chart = n**2

        assert focus * focus == q**3
        assert b * focus == gain
        assert chart * b == 1

        # Energy conservation through ideal focus and child-volume scaling.
        child_energy_ratio = (b * focus) ** 2 * Fraction(1, q**3)
        assert child_energy_ratio == b**2 == Fraction(1, n**4)

        logq = math.log(q)
        separation = n**3.5
        parent_route = float(b) * logq
        envelope_route = float(b) * logq / separation
        child_route = float(b) * q * logq
        preamplified_route = float(b) * focus * logq / separation
        pressure_after_chart = chart * logq / separation

        assert math.isclose(parent_route, 8.0 * math.log(n) / n**2,
                            rel_tol=2e-15)
        assert math.isclose(envelope_route, 8.0 * math.log(n) / n**5.5,
                            rel_tol=2e-15)
        assert math.isclose(child_route, 8.0 * n**6 * math.log(n),
                            rel_tol=2e-15)
        assert math.isclose(preamplified_route,
                            8.0 * n**6.5 * math.log(n), rel_tol=2e-15)
        assert math.isclose(pressure_after_chart,
                            8.0 * math.log(n) / n**1.5,
                            rel_tol=2e-15)

    # Passing route tends to zero; the two wrong orders grow.
    n0, n1 = 100, 1000
    assert 8.0 * math.log(n1) / n1**2 < 8.0 * math.log(n0) / n0**2
    assert 8.0 * n1**6 * math.log(n1) > 8.0 * n0**6 * math.log(n0)
    assert 8.0 * n1**6.5 * math.log(n1) > 8.0 * n0**6.5 * math.log(n0)


def a2_nonclosure_check() -> None:
    h1 = (1, -1, 0)
    h2 = (0, 1, -1)
    h3 = (-1, 0, 1)
    hexagon = {h1, h2, h3, neg(h1), neg(h2), neg(h3)}
    assert add(add(h1, h2), h3) == (0, 0, 0)
    assert sum(x * x for x in h1) == 2
    assert sum(x * x for x in h2) == 2
    assert sum(x * x for x in h3) == 2
    assert (0, 0, 0) not in hexagon
    assert add(h1, neg(h2)) == (1, -2, 1)
    assert add(h1, neg(h2)) not in hexagon

    outside_sums = {
        add(p, q)
        for p, q in itertools.product(hexagon, repeat=2)
        if add(p, q) != (0, 0, 0) and add(p, q) not in hexagon
    }
    assert outside_sums
    assert (1, -2, 1) in outside_sums
    assert sum(x * x for x in (1, -2, 1)) == 6
    # This is a support statement only; polarization coefficients may cancel.


def material_and_gevrey_checks() -> None:
    # Diagonal determinant test for F^{-T}; the identity itself is general.
    for a, b in ((2.0, 3.0), (0.5, 7.0), (11.0, 0.25)):
        c = 1.0 / (a * b)  # det F = 1
        det_f_inverse_transpose = (1.0 / a) * (1.0 / b) * (1.0 / c)
        assert math.isclose(det_f_inverse_transpose, 1.0, rel_tol=1e-14)

    # Logarithmic Gevrey-2 remainder with K=n^32 and
    # m=floor(eta*j^2/log(n)).  A fixed C and n^20 prefactor are included.
    eta = 0.1
    coefficient_constant = 10.0
    for j in (50, 100, 200, 500):
        n = j + 1
        m = max(1, int(eta * j * j / math.log(n)))
        log_remainder = (
            (m + 1) * math.log(coefficient_constant)
            + 2.0 * math.lgamma(m + 1)
            - 32.0 * m * math.log(n)
            + 20.0 * math.log(n)
        )
        assert log_remainder < -j * j


def backward_weight_check() -> None:
    # Scalar Duhamel identity: z'=gamma*z+f, z(0)=0.  It demonstrates why
    # terminal leakage is weighted by the remaining gain.
    gamma = 3.0
    terminal = 2.0
    forcing = 0.25
    exact = forcing * (math.exp(gamma * terminal) - 1.0) / gamma
    weighted_integral = (
        forcing * math.exp(gamma * terminal)
        * (1.0 - math.exp(-gamma * terminal)) / gamma
    )
    unweighted_integral = forcing * terminal
    assert math.isclose(exact, weighted_integral, rel_tol=1e-14)
    assert exact > 30.0 * unweighted_integral


def summable_target_checks() -> None:
    # Exact exponents: retained wake n^-4 and chart-weighted separated tail
    # 8 log(n)n^-3/2 are summable.  The integral tail formula is positive and
    # tends to zero.
    for n in (10, 100, 10_000):
        pressure_integral_tail = 8.0 * (2.0 * math.log(n) + 4.0) / math.sqrt(n)
        assert pressure_integral_tail > 0.0
    assert 8.0 * (2.0 * math.log(10_000) + 4.0) / 100.0 < 1.8

    wake_partial = sum(1.0 / n**4 for n in range(2, 100_000))
    assert wake_partial < 0.083

    # A retained dark/exported complement may be n^-4.  A component that
    # returns through an n^2 chart must be n^-6 before chart inversion.
    for n in range(2, 100):
        chart_loss = n**2
        active_prechart = Fraction(1, n**6)
        retained_wake = Fraction(1, n**4)
        assert chart_loss * active_prechart == retained_wake


def main() -> None:
    focus_and_activation_checks()
    a2_nonclosure_check()
    material_and_gevrey_checks()
    backward_weight_check()
    summable_target_checks()
    print("PASS: C132--C136 material, focus, activation, and leakage algebra")
    print("OPEN: backward-weighted active-focus leakage (BAFL) stage theorem")


if __name__ == "__main__":
    main()
