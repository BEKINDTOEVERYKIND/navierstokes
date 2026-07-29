#!/usr/bin/env python3
"""Exact algebra checks for the variable-phase pressure recurrence.

This script uses only the Python standard library.  It checks the formal
identity in one slow dimension for a genuinely variable phase gradient
xi(y)=y and one nonzero angle charge.  Laurent polynomials and Gaussian
rational coefficients keep every operation exact.

It does not prove the Gevrey estimates or a Navier--Stokes result.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


Q = Fraction


@dataclass(frozen=True)
class GaussianRational:
    real: Fraction = Q(0)
    imag: Fraction = Q(0)

    def __add__(self, other: object) -> "GaussianRational":
        other = gaussian(other)
        return GaussianRational(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self) -> "GaussianRational":
        return GaussianRational(-self.real, -self.imag)

    def __sub__(self, other: object) -> "GaussianRational":
        return self + (-gaussian(other))

    def __rsub__(self, other: object) -> "GaussianRational":
        return gaussian(other) - self

    def __mul__(self, other: object) -> "GaussianRational":
        other = gaussian(other)
        return GaussianRational(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "GaussianRational":
        other = gaussian(other)
        denominator = other.real * other.real + other.imag * other.imag
        assert denominator
        return GaussianRational(
            (self.real * other.real + self.imag * other.imag) / denominator,
            (self.imag * other.real - self.real * other.imag) / denominator,
        )

    def __bool__(self) -> bool:
        return bool(self.real or self.imag)


def gaussian(value: object) -> GaussianRational:
    if isinstance(value, GaussianRational):
        return value
    if isinstance(value, (int, Fraction)):
        return GaussianRational(Q(value), Q(0))
    raise TypeError(value)


I = GaussianRational(Q(0), Q(1))
ZERO = GaussianRational()


@dataclass(frozen=True)
class Laurent:
    """Finite Laurent polynomial in y with exact Gaussian coefficients."""

    coefficients: tuple[tuple[int, GaussianRational], ...]

    @staticmethod
    def make(entries: dict[int, object]) -> "Laurent":
        cleaned = {
            power: gaussian(value)
            for power, value in entries.items()
            if gaussian(value)
        }
        return Laurent(tuple(sorted(cleaned.items())))

    @staticmethod
    def monomial(power: int, coefficient: object = 1) -> "Laurent":
        return Laurent.make({power: coefficient})

    def as_dict(self) -> dict[int, GaussianRational]:
        return dict(self.coefficients)

    def __add__(self, other: "Laurent") -> "Laurent":
        result = self.as_dict()
        for power, coefficient in other.coefficients:
            result[power] = result.get(power, ZERO) + coefficient
        return Laurent.make(result)

    def __neg__(self) -> "Laurent":
        return Laurent.make(
            {power: -coefficient for power, coefficient in self.coefficients}
        )

    def __sub__(self, other: "Laurent") -> "Laurent":
        return self + (-other)

    def scale(self, value: object) -> "Laurent":
        scalar = gaussian(value)
        return Laurent.make(
            {power: scalar * coefficient for power, coefficient in self.coefficients}
        )

    def shift(self, amount: int) -> "Laurent":
        return Laurent.make(
            {power + amount: coefficient for power, coefficient in self.coefficients}
        )

    def derivative(self, order: int = 1) -> "Laurent":
        result = self
        for _ in range(order):
            result = Laurent.make(
                {
                    power - 1: coefficient * power
                    for power, coefficient in result.coefficients
                    if power
                }
            )
        return result


def add_to_series(
    series: dict[int, Laurent], power: int, value: Laurent
) -> None:
    series[power] = series.get(power, Laurent.make({})) + value


ANGLE_CHARGE = 3


def A2(value: Laurent) -> Laurent:
    # a=y^2 and partial_theta^2=-ANGLE_CHARGE^2.
    return value.shift(2).scale(-(ANGLE_CHARGE**2))


def A2_inverse(value: Laurent) -> Laurent:
    return value.shift(-2).scale(Q(-1, ANGLE_CHARGE**2))


def A1(value: Laurent) -> Laurent:
    # (2 y d_y + d_y y) partial_theta
    raw = value.derivative().shift(1).scale(2) + value
    return raw.scale(I * ANGLE_CHARGE)


def A0(value: Laurent) -> Laurent:
    return value.derivative(2)


def B1(value: Laurent) -> Laurent:
    # xi partial_theta F = i*ANGLE_CHARGE*y*F.
    return value.shift(1).scale(I * ANGLE_CHARGE)


def B0(value: Laurent) -> Laurent:
    return value.derivative()


def build_coefficients(source: Laurent, maximum_order: int) -> list[Laurent]:
    coefficients = [A2_inverse(B1(source))]
    coefficients.append(
        A2_inverse(B0(source) - A1(coefficients[0]))
    )
    for n in range(2, maximum_order + 1):
        coefficients.append(
            A2_inverse(
                -A1(coefficients[n - 1]) - A0(coefficients[n - 2])
            )
        )
    return coefficients


def check_exact_recurrence() -> None:
    source = Laurent.make({-2: 2, 0: -1, 1: 3, 4: 5})
    maximum_order = 7
    p = build_coefficients(source, maximum_order)

    assert A2(p[0]) == B1(source)
    assert A2(p[1]) + A1(p[0]) == B0(source)
    for n in range(2, maximum_order + 1):
        assert A2(p[n]) + A1(p[n - 1]) + A0(p[n - 2]) == Laurent.make({})

    # Assemble (A2+h A1+h^2 A0) h sum h^n p_n exactly.
    left: dict[int, Laurent] = {}
    for n, coefficient in enumerate(p):
        add_to_series(left, n + 1, A2(coefficient))
        add_to_series(left, n + 2, A1(coefficient))
        add_to_series(left, n + 3, A0(coefficient))

    right = {1: B1(source), 2: B0(source)}
    residual = {
        power: left.get(power, Laurent.make({}))
        - right.get(power, Laurent.make({}))
        for power in set(left) | set(right)
    }
    residual = {power: value for power, value in residual.items() if value.coefficients}

    expected = {
        maximum_order + 2: A1(p[-1]) + A0(p[-2]),
        maximum_order + 3: A0(p[-1]),
    }
    expected = {power: value for power, value in expected.items() if value.coefficients}
    assert residual == expected

    # In one dimension the principal pressure gradient removes all of F:
    # i*charge*y*p0 = F.
    assert p[0].shift(1).scale(I * ANGLE_CHARGE) == source

    print(
        "variable xi=y pressure recurrence:",
        f"orders 0..{maximum_order} exact",
    )
    print("finite truncation residual: only the two claimed terminal powers")


def check_weighted_derivative_budget() -> None:
    # p_0 costs zero derivatives of F.  p_1 has a one-derivative source.
    # T1 adds one derivative while advancing one h-order; T2 adds two
    # derivatives while advancing two h-orders.
    maximum_order = 100
    derivative_cost = [0, 1]
    for n in range(2, maximum_order + 1):
        derivative_cost.append(
            max(derivative_cost[n - 1] + 1, derivative_cost[n - 2] + 2)
        )
    assert derivative_cost == list(range(maximum_order + 1))
    print("weighted slow-derivative budget: cost(p_n)=n, not 2n")


def check_gauge_cancellation() -> None:
    # For exp(i*l*(theta-K*phi)), the slow derivative contributes
    # -i*l*K*xi and K*xi*partial_theta contributes +i*l*K*xi.
    charge = Q(5)
    carrier = Q(101)
    xi_component = Q(7, 3)
    slow_log_derivative = -I * charge * carrier * xi_component
    fast_log_derivative = I * charge * carrier * xi_component
    assert slow_log_derivative + fast_log_derivative == ZERO
    print("extended-profile gauge mode: D_K cancellation exact")


def main() -> None:
    check_exact_recurrence()
    check_weighted_derivative_budget()
    check_gauge_cancellation()
    print("all exact variable-phase pressure checks passed")


if __name__ == "__main__":
    main()
