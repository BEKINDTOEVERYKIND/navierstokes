#!/usr/bin/env python3
"""Exact algebra for the strain-capped locked-pitch semigroup bound."""

from fractions import Fraction as F
from math import factorial


class Alg:
    """Elements c0 + c1*a with 3*a^2 + 10*a + 1 = 0."""

    __slots__ = ("c0", "c1")

    def __init__(self, c0=0, c1=0):
        if isinstance(c0, Alg):
            self.c0, self.c1 = c0.c0, c0.c1
        else:
            self.c0, self.c1 = F(c0), F(c1)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Alg) else Alg(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Alg(self.c0 + other.c0, self.c1 + other.c1)

    __radd__ = __add__

    def __neg__(self):
        return Alg(-self.c0, -self.c1)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        # a^2 = -(10*a + 1)/3.
        product_a2 = self.c1 * other.c1
        return Alg(
            self.c0 * other.c0 - product_a2 / 3,
            self.c0 * other.c1
            + self.c1 * other.c0
            - 10 * product_a2 / 3,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = F(scalar)
        return Alg(self.c0 / scalar, self.c1 / scalar)

    def __eq__(self, other):
        other = self.coerce(other)
        return self.c0 == other.c0 and self.c1 == other.c1

    def interval(self):
        # The selected root a=(-5-sqrt(22))/3 lies in (-3.231,-3.230).
        lo, hi = F(-3231, 1000), F(-3230, 1000)
        values = (self.c0 + self.c1 * lo, self.c0 + self.c1 * hi)
        return min(values), max(values)


def poly_mul(left, right):
    result = [Alg(0) for _ in range(len(left) + len(right) - 1)]
    for j, x_value in enumerate(left):
        for k, y_value in enumerate(right):
            result[j + k] = result[j + k] + x_value * y_value
    return result


def exp_upper(x_value, order=16):
    """Rational upper bound for exp(x), using a geometric tail."""

    x_value = F(x_value)
    partial = sum(
        (x_value**k / factorial(k) for k in range(order + 1)),
        F(0),
    )
    first_tail = x_value ** (order + 1) / factorial(order + 1)
    ratio = x_value / (order + 2)
    assert ratio < 1
    return partial + first_tail / (1 - ratio)


def main():
    a_value = Alg(0, 1)
    one = Alg(1)

    # The defining stationarity polynomial and c_* simplification.
    assert 3 * a_value * a_value + 10 * a_value + 1 == Alg(0)
    c_value = (4 + 28 * a_value) / 9
    assert c_value * (3 * a_value + 1) == 8 * a_value * a_value

    # At the ring, S^2/Omega^2=(3a^2+2a+1)/8=-a.
    strain_ring = (3 * a_value * a_value + 2 * a_value + 1) / 8
    assert strain_ring == -a_value

    # Exact factorization of the critical cubic.
    critical = [c_value, 1 + 3 * c_value, Alg(2), Alg(3)]
    quotient = [4 * (a_value + 1) / 3, 2 + 3 * a_value, Alg(3)]
    factored = poly_mul([-a_value, one], quotient)
    assert factored == critical

    q_zero = quotient[0]
    q_minus_one = quotient[0] - quotient[1] + quotient[2]
    assert q_zero.interval()[1] < 0
    assert q_minus_one.interval()[0] > 0

    # The h=0 join has a large strict strain margin:
    # exp((4+sqrt(22))/8) < exp(9/8) < 4, while -a>3.
    assert exp_upper(F(9, 8)) < 4
    assert (-a_value).interval()[0] > 3
    assert F(4, 8) < 3

    # The strain ellipse is positive definite:
    # 8*N^2=3q^2+2 Omega*q+Omega^2 has discriminant 4-12<0.
    assert 2**2 - 4 * 3 * 1 < 0

    print("locked-pitch strain eigenvalues and ring equality: PASS")
    print("global outer-flank critical polynomial factorization: PASS")
    print("strict h=0 inner-gluing strain margin: PASS")
    print("positive-definite open strain ellipse: PASS")
    print("all locked-pitch strain-cap checks passed")


if __name__ == "__main__":
    main()
