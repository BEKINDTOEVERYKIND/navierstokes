#!/usr/bin/env python3
"""Exact/rational ledgers for C121--C125.

Checks the A2 Beltrami algebra, the decaying-pump gain and its rational
optimality inequalities, the nonlinear finite-normal-form energy/comparison
ledger, and the same-order centre/second-sideband support outputs.

This checker deliberately does not label the six-leaf block invariant under
the full PDE.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import exp, log


Vector = tuple[int, int, int]

N: Vector = (1, 1, 1)
ROOTS: tuple[Vector, ...] = (
    (1, -1, 0),
    (0, 1, -1),
    (-1, 0, 1),
)


def add(left: Vector, right: Vector) -> Vector:
    return tuple(x + y for x, y in zip(left, right))  # type: ignore[return-value]


def sub(left: Vector, right: Vector) -> Vector:
    return tuple(x - y for x, y in zip(left, right))  # type: ignore[return-value]


def scale(number: int, vector: Vector) -> Vector:
    return tuple(number * entry for entry in vector)  # type: ignore[return-value]


def dot(left: Vector, right: Vector) -> int:
    return sum(x * y for x, y in zip(left, right))


def cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def gain(t: float, sigma: float, pump: float, delta: float, damping: float) -> float:
    return sigma * pump * (1.0 - exp(-delta * t)) / delta - damping * t


def check_beltrami_algebra() -> None:
    # h=t+i*sqrt(2)*n.  Store the real t and the coefficient n of
    # i*sqrt(2).  The identity i r x h=sqrt(2)h is equivalent to
    # r x t=2n and r x n=-t.
    for root in ROOTS:
        tangent = cross(N, root)
        assert dot(root, tangent) == 0
        assert dot(root, N) == 0
        assert cross(root, tangent) == scale(2, N)
        assert cross(root, N) == scale(-1, tangent)
        # For h=a+i*sqrt(2)*b, the positive-helicity condition is
        # k x a=2b and -k x b=a.  Check both the positive mode
        # (k,a,b)=(r,t,n) and its real-field conjugate (-r,t,-n).
        assert scale(-1, cross(root, N)) == tangent
        negative_root = scale(-1, root)
        negative_imaginary = scale(-1, N)
        assert cross(negative_root, tangent) == scale(
            2, negative_imaginary
        )
        assert scale(
            -1, cross(negative_root, negative_imaginary)
        ) == tangent
        # |root|^2=2 gives Delta=-2 on the unit root shell.
        assert dot(root, root) == 2

    for viscosity in (F(1, 10), F(7, 13)):
        for k_value in (1, 5, 29):
            delta = 2 * viscosity * k_value * k_value
            # P'= -delta P for the exact heat-decaying pump.
            pump = F(11, 7)
            assert -delta * pump == -2 * viscosity * k_value * k_value * pump


def check_linear_gain_and_optimum() -> None:
    cases = (
        # sigma, P0, delta, d
        (2.0, 8.0, 0.5, 3.0),
        (2.532, 20.0, 0.2, 4.0),
        (1.25, 50.0, 0.75, 2.5),
    )
    for sigma, pump, delta, damping in cases:
        ratio = sigma * pump / damping
        assert ratio > 1.0
        t_star = log(ratio) / delta
        exponent = gain(t_star, sigma, pump, delta, damping)
        closed = damping * (ratio - 1.0 - log(ratio)) / delta
        assert abs(exponent - closed) < 1e-12
        assert exponent > gain(0.9 * t_star, sigma, pump, delta, damping)
        assert exponent > gain(1.1 * t_star, sigma, pump, delta, damping)
        # The derivative vanishes exactly up to floating evaluation.
        assert abs(sigma * pump * exp(-delta * t_star) - damping) < 1e-12

    # Rational positivity certificate R-1-log R>0 is represented by the
    # integral of (s-1)/s from 1 to R.  The elementary derivative is
    # positive for every rational sample R>1.
    for ratio in (F(1001, 1000), F(3, 2), F(2), F(10)):
        assert (ratio - 1) / ratio > 0
        assert float(ratio - 1) - log(float(ratio)) > 0


def check_shell_and_reynolds_ledger() -> None:
    for viscosity in (F(1, 100), F(3, 1000)):
        for k_value in (2, 9, 40):
            for m_value in (1, 7, 31):
                centre = scale(m_value, N)
                pump_decay = 2 * viscosity * k_value * k_value
                leaf_squared = 3 * m_value * m_value + 2 * k_value * k_value
                leaf_damping = viscosity * leaf_squared
                assert leaf_damping / pump_decay == F(
                    leaf_squared, 2 * k_value * k_value
                )
                for root in ROOTS:
                    for sign in (-1, 1):
                        leaf = add(centre, scale(sign * k_value, root))
                        assert dot(leaf, leaf) == leaf_squared


def check_nonlinear_normal_form_energy_and_seed_bound() -> None:
    # Exact cancellation in p'=-delta*p-2*sigma*sum(a_i^2),
    # a_i'=(sigma*p-d)*a_i, E=p^2+2*sum(a_i^2).
    samples = (
        (F(2), F(3), F(1, 5), (F(1, 7),)),
        (F(7, 3), F(5, 2), F(2, 9), (F(1, 8), F(2, 11), F(3, 13))),
    )
    for sigma, damping, delta, leaves in samples:
        pump = F(11, 5)
        leaf_square_sum = sum((leaf * leaf for leaf in leaves), F(0))
        pump_dot = -delta * pump - 2 * sigma * leaf_square_sum
        leaf_dots = tuple((sigma * pump - damping) * leaf for leaf in leaves)
        energy_dot = 2 * pump * pump_dot + 4 * sum(
            (leaf * leaf_dot for leaf, leaf_dot in zip(leaves, leaf_dots)), F(0)
        )
        expected = -2 * delta * pump * pump - 4 * damping * leaf_square_sum
        assert energy_dot == expected

    # A concrete strict nonlinear bootstrap gate from (5.6).
    sigma = 2.0
    pump = 8.0
    delta = 0.5
    damping = 3.0
    ratio = sigma * pump / damping
    t_star = log(ratio) / delta
    g_star = damping * (ratio - 1.0 - log(ratio)) / delta
    eta = 0.01
    seed = 1e-10
    maximum_leaf = seed * exp(g_star)
    assert 2 * sigma * t_star * maximum_leaf**2 < damping / (2 * sigma)
    assert sigma**2 * maximum_leaf**2 * t_star**2 < eta
    assert exp(g_star - eta) > 1.0


def check_same_order_support_outputs() -> None:
    for k_value in (1, 8, 37):
        for m_value in (1, 5, 23):
            centre = scale(m_value, N)
            first_shell = {
                add(centre, scale(sign * k_value, root))
                for root in ROOTS
                for sign in (-1, 1)
            }
            for root in ROOTS:
                leaf = add(centre, scale(k_value, root))
                centre_output = sub(leaf, scale(k_value, root))
                outer_output = add(leaf, scale(k_value, root))
                assert centre_output == centre
                assert outer_output == add(centre, scale(2 * k_value, root))
                assert centre_output not in first_shell
                assert outer_output not in first_shell

                # A desired retained edge is also exactly one pump addition.
                other_roots = [candidate for candidate in ROOTS if candidate != root]
                retained = add(leaf, scale(k_value, other_roots[0]))
                assert retained in first_shell


def main() -> None:
    check_beltrami_algebra()
    check_linear_gain_and_optimum()
    check_shell_and_reynolds_ledger()
    check_nonlinear_normal_form_energy_and_seed_bound()
    check_same_order_support_outputs()
    print("C121-C125 unforced decaying-pump checker: PASS")
    print("  exact background: heat-decaying A2 Beltrami pump")
    print("  exact ideal gain: (d/delta)*(R-1-log R)")
    print("  unresolved PDE gate: same-order off-ladder/nonlinear Duhamel norm")


if __name__ == "__main__":
    main()
