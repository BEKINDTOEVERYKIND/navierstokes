#!/usr/bin/env python3
"""Dependency-free exact ledgers for C178.

This checker verifies:

* the Fourier Duhamel multiplier and local Taylor preparation for a
  normalized polynomial ramp;
* the analytic-buffer terminal semigroup identity and its contraction
  bound, using independent deterministic quadrature;
* exact commutation of divergence with repeated Laplacians on a polynomial
  curl representative;
* monotone convergence and the Poisson-tail inequality used in C178;
* the C176 preparation, heat-number, C125, and degree-one exponents; and
* the tuned C142 factorial comparison.

It does not certify analytic unique continuation, a Gevrey tail, the full
A2 evolution family, C125, MCKC, LCE, BAFL, or an unforced stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import exp, factorial, log


Poly = dict[tuple[int, int, int], F]


def padd(a: Poly, b: Poly) -> Poly:
    result = dict(a)
    for monomial, coefficient in b.items():
        result[monomial] = result.get(monomial, F(0)) + coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def pscale(c: F, p: Poly) -> Poly:
    return {monomial: c * coefficient for monomial, coefficient in p.items()}


def derivative(p: Poly, axis: int) -> Poly:
    result: Poly = {}
    for monomial, coefficient in p.items():
        exponent = monomial[axis]
        if exponent == 0:
            continue
        reduced = list(monomial)
        reduced[axis] -= 1
        key = tuple(reduced)
        result[key] = result.get(key, F(0)) + coefficient * exponent
    return result


def laplacian(p: Poly) -> Poly:
    result: Poly = {}
    for axis in range(3):
        result = padd(result, derivative(derivative(p, axis), axis))
    return result


def divergence(vector: tuple[Poly, Poly, Poly]) -> Poly:
    result: Poly = {}
    for axis in range(3):
        result = padd(result, derivative(vector[axis], axis))
    return result


def vector_laplacian(vector: tuple[Poly, Poly, Poly]) -> tuple[Poly, Poly, Poly]:
    return tuple(laplacian(component) for component in vector)  # type: ignore[return-value]


def ramp_moment(order: int, horizon: F) -> F:
    """Moment for theta'(s)=2s/T^2 on [0,T]."""

    return F(2) * horizon**order / F(order + 2)


def truncated_backward_multiplier(
    x: float, viscosity: float, horizon: float, degree: int
) -> float:
    # Integral of (2s/T^2) sum_{k<=M}(nu*s*x)^k/k! ds.
    return sum(
        (viscosity * x) ** order
        * float(ramp_moment(order, F(str(horizon))))
        / factorial(order)
        for order in range(degree + 1)
    )


def exact_backward_multiplier(x: float, viscosity: float, horizon: float) -> float:
    # Stable quadrature-free formula for integral_0^T (2s/T^2)e^(nu*x*s) ds.
    rate = viscosity * x
    if rate == 0:
        return 1.0
    return 2.0 * (exp(rate * horizon) * (rate * horizon - 1.0) + 1.0) / (
        horizon**2 * rate**2
    )


def terminal_error_multiplier(
    x: float, viscosity: float, horizon: float, degree: int
) -> float:
    return exp(-viscosity * horizon * x) * (
        truncated_backward_multiplier(x, viscosity, horizon, degree)
        - exact_backward_multiplier(x, viscosity, horizon)
    )


def simpson_integral(function, left: float, right: float, panels: int = 4096) -> float:
    """Deterministic Simpson quadrature, with an even panel count."""

    assert panels > 0 and panels % 2 == 0
    step = (right - left) / panels
    total = function(left) + function(right)
    for index in range(1, panels):
        total += (4.0 if index % 2 else 2.0) * function(left + index * step)
    return total * step / 3.0


def check_analytic_buffer_identity() -> None:
    # theta'(s)=2s/T^2.  The prepared multiplier is evaluated two ways:
    # independently by quadrature and by the closed backward multiplier.
    horizon = 1.5
    viscosity = 0.17
    heat_time = viscosity * horizon

    for sigma in (heat_time, 1.4 * heat_time, 2.0 * heat_time):
        for x in (0.0, 0.1, 1.0, 7.0, 30.0):
            ramp = lambda s: 2.0 * s / horizon**2
            datum_quadrature = simpson_integral(
                lambda s: ramp(s) * exp(-(sigma - viscosity * s) * x),
                0.0,
                horizon,
            )
            datum_closed = exp(-sigma * x) * exact_backward_multiplier(
                x, viscosity, horizon
            )
            assert abs(datum_quadrature - datum_closed) < 2e-12

            buffered_profile = exp(-sigma * x)
            forced_terminal = simpson_integral(
                lambda s: ramp(s)
                * exp(-viscosity * (horizon - s) * x)
                * buffered_profile,
                0.0,
                horizon,
            )
            homogeneous_terminal = exp(-heat_time * x) * datum_quadrature
            assert abs(homogeneous_terminal - forced_terminal) < 2e-12

            # sigma>=nu*T makes every factor in the datum a forward heat
            # multiplier, so both L2 operator norms are at most one.
            assert 0.0 <= buffered_profile <= 1.0
            assert 0.0 <= datum_quadrature <= 1.0 + 2e-12


def check_duhamel_and_taylor_preparation() -> None:
    horizon = F(3, 2)
    assert ramp_moment(0, horizon) == 1
    assert ramp_moment(1, horizon) == 1
    assert ramp_moment(2, horizon) == F(9, 8)

    viscosity = 0.17
    t_float = float(horizon)
    for x in (0.0, 0.1, 1.0, 7.0, 30.0):
        exact = exact_backward_multiplier(x, viscosity, t_float)
        assert exact >= 1.0 - 1e-12
        errors = [
            abs(terminal_error_multiplier(x, viscosity, t_float, degree))
            for degree in (0, 1, 2, 4, 8, 16, 32)
        ]
        assert errors[-1] < 1e-11
        # Taylor partial sums are monotone for this nonnegative ramp.
        assert all(left >= right for left, right in zip(errors, errors[1:]))

        # Direct Duhamel identity: terminal homogeneous part minus the
        # forced response equals the displayed error multiplier.
        for degree in (0, 1, 3, 7):
            prepared = exp(-viscosity * t_float * x) * (
                truncated_backward_multiplier(x, viscosity, t_float, degree)
            )
            forced = exp(-viscosity * t_float * x) * exact
            assert abs((prepared - forced) - terminal_error_multiplier(
                x, viscosity, t_float, degree
            )) < 1e-13


def check_divergence_laplacian_commutation() -> None:
    # Treat psi as a compactly supported smooth test profile algebraically;
    # polynomial differentiation checks curl=(d_y psi,-d_x psi,0).
    psi: Poly = {
        (5, 4, 3): F(7, 3),
        (3, 6, 2): F(-5, 2),
        (2, 2, 7): F(11, 5),
    }
    zero: Poly = {}
    field = (derivative(psi, 1), pscale(F(-1), derivative(psi, 0)), zero)
    assert divergence(field) == zero

    current = field
    for _ in range(5):
        current = vector_laplacian(current)
        assert divergence(current) == zero

    # Exact linear combination of Laplacian powers remains divergence free.
    nu = F(2, 7)
    horizon = F(3, 2)
    prepared = ({}, {}, {})
    current = field
    components = [dict() for _ in range(3)]
    for order in range(4):
        coefficient = nu**order * ramp_moment(order, horizon) / factorial(order)
        for axis in range(3):
            components[axis] = padd(
                components[axis], pscale(coefficient, current[axis])
            )
        current = vector_laplacian(current)
    prepared = (components[0], components[1], components[2])
    assert divergence(prepared) == zero


def poisson_tail(y: float, cutoff: int) -> float:
    # Stable finite recurrence plus a rigorously negligible numerical tail
    # for the moderate test values below.
    term = exp(-y)
    total = 0.0
    for order in range(0, 400):
        if order > cutoff:
            total += term
        term *= y / (order + 1)
    assert term < 1e-50
    return total


def check_poisson_tail_bound() -> None:
    for y in (0.0, 0.01, 0.2, 1.0, 3.0, 10.0, 30.0):
        for cutoff in (0, 1, 2, 5, 12):
            lhs = poisson_tail(y, cutoff)
            rhs = y ** (cutoff + 1) / factorial(cutoff + 1)
            assert lhs <= rhs + 2e-14


def check_stage_exponents() -> None:
    # q=n^8, q^(-5/2)=n^-20, seed=n^-28.
    for n in (2, 3, 7, 20):
        q = n**8
        assert q**2 * n**4 == n**20
        source_without_lambda_j = F(1, n**20)
        seed = F(1, n**28)
        assert source_without_lambda_j / seed == n**8
        # Squared norms: prep energy n^-40 versus seed energy n^-56.
        assert source_without_lambda_j**2 / seed**2 == n**16

        # If Lambda*J^(7/2)<=n^2, total preparation is at most n^-18.
        prep_edge = F(1, n**18)
        assert prep_edge / seed == n**10

        # With alpha<=q^-1 and M=1, the low-frequency Taylor exponent is
        # q^-2 times a polylog factor: the pure n power is n^-16.
        low_power = F(1, q**2)
        assert low_power == F(1, n**16)
        assert prep_edge * low_power == F(1, n**34)
        # A chosen Gevrey tail n^-12 gives n^-30 at the edge.
        assert prep_edge * F(1, n**12) == F(1, n**30)


def check_heat_number_identity() -> None:
    # Choose the normalized tuned constants exactly:
    # log h=12 L, J=12 L, Lambda=mu*q^3/L.
    for n in (2, 3, 5, 11):
        q = F(n**8)
        logarithm_symbol = F(7, 5)  # cancels exactly
        mu = F(2, 13)
        log_h = 12 * logarithm_symbol
        j_slices = log_h
        lam_ratio = mu * q**3 / logarithm_symbol
        alpha = mu * q**2 * log_h / (lam_ratio * j_slices**2)
        assert alpha == F(1, 12) / q


def check_tuned_factorial_gate() -> None:
    # Squared version of mu*n^32*J^(7/2)/log(n), with J=12log(n).
    # Squaring removes the half exponent and preserves convergence to zero.
    values: list[float] = []
    for n in (20, 30, 40, 60, 80):
        stage_factorial = factorial(n - 1)
        mu = 1.0 / stage_factorial**2
        j_slices = 12.0 * log(n)
        ratio = mu * n**32 * j_slices**3.5 / log(n)
        values.append(ratio)
    assert all(left > right for left, right in zip(values, values[1:]))
    assert values[-1] < 1e-30


def main() -> None:
    check_duhamel_and_taylor_preparation()
    print("PASS C178: exact heat/Duhamel and local Taylor preparation")
    check_analytic_buffer_identity()
    print("PASS C178: exact analytic-buffer identity and contraction")
    check_divergence_laplacian_commutation()
    print("PASS C178: Laplacian preparation preserves divergence/support algebra")
    check_poisson_tail_bound()
    print("PASS C178: quantitative Poisson-tail bound")
    check_stage_exponents()
    check_heat_number_identity()
    print("PASS C178: C176 heat/preparation and n^-28 ledgers")
    check_tuned_factorial_gate()
    print("PASS C178: tuned C142 preparation is factorially below seed")
    print("OPEN: actual A2 evolution, pressure tails, C125, MCKC, LCE, BAFL")


if __name__ == "__main__":
    main()
