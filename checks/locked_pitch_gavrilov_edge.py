#!/usr/bin/env python3
"""Dependency-free certificate for the locked-pitch Gavrilov BAS edge.

The analytic proof is in
``research/2026-08-02-locked-pitch-gavrilov-edge-profile.md``.  All
inequality decisions below use rational intervals; floating point is used
only when printing human-readable values.
"""

from fractions import Fraction as F
from math import comb

from ao_batchelor_full_edge_matched import exp_point
from ao_batchelor_global_bas import I, enclose_sqrt


def as_i(value):
    return value if isinstance(value, I) else I(value)


class Jet2:
    """Value, first derivative, and second derivative interval jet."""

    def __init__(self, value, first=0, second=0):
        self.value = as_i(value)
        self.first = as_i(first)
        self.second = as_i(second)

    def __add__(self, other):
        other = as_jet(other)
        return Jet2(
            self.value + other.value,
            self.first + other.first,
            self.second + other.second,
        )

    __radd__ = __add__

    def __neg__(self):
        return Jet2(-self.value, -self.first, -self.second)

    def __sub__(self, other):
        return self + (-as_jet(other))

    def __rsub__(self, other):
        return as_jet(other) - self

    def __mul__(self, other):
        other = as_jet(other)
        return Jet2(
            self.value * other.value,
            self.first * other.value + self.value * other.first,
            self.second * other.value
            + 2 * self.first * other.first
            + self.value * other.second,
        )

    __rmul__ = __mul__

    def reciprocal(self):
        value_squared = self.value * self.value
        return Jet2(
            1 / self.value,
            -self.first / value_squared,
            2 * self.first * self.first / (value_squared * self.value)
            - self.second / value_squared,
        )

    def __truediv__(self, other):
        return self * as_jet(other).reciprocal()

    def __rtruediv__(self, other):
        return as_jet(other) / self


def as_jet(value):
    return value if isinstance(value, Jet2) else Jet2(value)


def poly_add(left, right):
    size = max(len(left), len(right))
    return [
        (left[j] if j < len(left) else F(0))
        + (right[j] if j < len(right) else F(0))
        for j in range(size)
    ]


def poly_scale(poly, scalar):
    return [F(scalar) * coefficient for coefficient in poly]


def poly_mul(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for j, left_coefficient in enumerate(left):
        for k, right_coefficient in enumerate(right):
            result[j + k] += left_coefficient * right_coefficient
    return result


def poly_derivative(poly):
    return [F(j) * poly[j] for j in range(1, len(poly))]


def poly_shift(poly, amount):
    """Ascending coefficients of p(x + amount)."""
    amount = F(amount)
    return [
        sum(
            poly[k] * F(comb(k, j)) * amount ** (k - j)
            for k in range(j, len(poly))
        )
        for j in range(len(poly))
    ]


def poly_divmod(numerator, denominator):
    """Exact ascending-coefficient polynomial division."""
    numerator = list(map(F, numerator))
    denominator = list(map(F, denominator))
    while len(numerator) > 1 and numerator[-1] == 0:
        numerator.pop()
    quotient = [F(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        degree = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[degree] = coefficient
        for j, value in enumerate(denominator):
            numerator[degree + j] -= coefficient * value
        while len(numerator) > 1 and numerator[-1] == 0:
            numerator.pop()
    return quotient, numerator


def locked_f(h, sigma, helical_ratio):
    """The off-resonance b/Omega^2 function at fixed beta*r."""
    return (
        -2
        * helical_ratio
        * ((1 + h) * sigma + helical_ratio * (h + 2))
        / (1 + helical_ratio * helical_ratio)
    )


def main():
    sqrt_22 = enclose_sqrt(
        I(22), F(4690415759, 10**9), F(4690415760, 10**9)
    )
    sqrt_2 = enclose_sqrt(
        I(2), F(1414213562, 10**9), F(1414213563, 10**9)
    )
    h_star = (-5 - sqrt_22) / 3
    sigma = 1 / sqrt_2
    h_dot = 8 * h_star * h_star / (1 + 3 * h_star)
    helical_ratio = h_star / (sigma * (1 + h_star))

    assert h_star.subset(F(-323014, 100000), F(-323013, 100000))
    assert h_dot.subset(F(-96049, 10000), F(-96048, 10000))
    assert helical_ratio.subset(F(20483, 10000), F(20484, 10000))

    # The transverse stationarity polynomial and the exact edge value.
    stationarity_polynomial = [F(1), F(10), F(3)]
    assert 10**2 - 4 * 3 == 4 * 22  # roots (-5 +/- sqrt(22))/3
    stationarity = 3 * h_star * h_star + 10 * h_star + 1
    assert stationarity.lo < 0 < stationarity.hi
    numerator = 3 * h_star * h_star + 6 * h_star + 1
    denominator = 3 * h_star * h_star + 2 * h_star + 1
    edge = -2 * h_star * numerator / denominator
    assert (edge + h_star).lo < 0 < (edge + h_star).hi

    # Exact polynomial-remainder certificates for f(h_*)=-h_* and
    # f'(h_*)/f(h_*)=-(3h_*+1)/(4h_*).
    n_poly = [F(1), F(6), F(3)]
    d_poly = [F(1), F(2), F(3)]
    edge_identity = poly_add(poly_scale(n_poly, 2), poly_scale(d_poly, -1))
    _, edge_remainder = poly_divmod(edge_identity, stationarity_polynomial)
    assert edge_remainder == [F(0)]
    n_prime = poly_derivative(n_poly)
    d_prime = poly_derivative(d_poly)
    # Clear 4*h*N*D from log(f)' + (3h+1)/(4h).
    log_identity = poly_scale(poly_mul(n_poly, d_poly), 4)
    log_identity = poly_add(
        log_identity,
        [F(0)] + poly_scale(poly_mul(n_prime, d_poly), 4),
    )
    log_identity = poly_add(
        log_identity,
        [F(0)] + poly_scale(poly_mul(n_poly, d_prime), -4),
    )
    log_identity = poly_add(
        log_identity,
        poly_mul([F(1), F(3)], poly_mul(n_poly, d_poly)),
    )
    _, log_remainder = poly_divmod(log_identity, stationarity_polynomial)
    assert log_remainder == [F(0)]

    # Derive the numerator of (log F)'' on the steep lobe.  With x=-h,
    # A=3x^2-6x+1 and B=3x^2-2x+1, the common denominator is x^2 A^2 B^2.
    a_poly = [F(1), F(-6), F(3)]
    b_poly = [F(1), F(-2), F(3)]
    a_prime = poly_derivative(a_poly)
    b_prime = poly_derivative(b_poly)
    log_second_numerator = poly_scale(
        poly_mul(poly_mul(a_poly, a_poly), poly_mul(b_poly, b_poly)), -1
    )
    a_log_part = poly_mul(
        poly_add(poly_scale(a_poly, 6), poly_scale(poly_mul(a_prime, a_prime), -1)),
        poly_mul(b_poly, b_poly),
    )
    b_log_part = poly_mul(
        poly_add(poly_scale(b_poly, 6), poly_scale(poly_mul(b_prime, b_prime), -1)),
        poly_mul(a_poly, a_poly),
    )
    log_second_numerator = poly_add(
        log_second_numerator, [F(0), F(0)] + a_log_part
    )
    log_second_numerator = poly_add(
        log_second_numerator,
        poly_scale([F(0), F(0)] + b_log_part, -1),
    )
    positive_poly = [
        F(1), F(-16), F(132), F(-504), F(1110),
        F(-1152), F(612), F(-216), F(81),
    ]
    assert log_second_numerator == poly_scale(positive_poly, -1)

    # R(x)>0 for x>=1 follows immediately from the positive coefficients
    # of R(1+y).  The steep unstable lobe has x>1+sqrt(2/3)>1.
    shifted = poly_shift(positive_poly, 1)
    assert shifted == [
        F(48), F(224), F(672), F(1632), F(2640),
        F(2520), F(1368), F(432), F(81),
    ]
    assert all(coefficient > 0 for coefficient in shifted)

    # The mild lobe has x<1/5 and F<2x<2/5.  Its exponential weight is
    # below exp((4+sqrt(22))/8)<exp(9/8)<4, so it is below 8/5.  The
    # selected steep-lobe value is -h_*>3.
    _, exp_upper = exp_point(F(9, 8))
    assert exp_upper < 4
    assert -h_star.hi > 3
    assert F(8, 5) < 3

    # Exact second-order local geometry in t=log r.  The interval widths
    # account for sqrt(2) and sqrt(22); no float decides a sign.
    h_jet = Jet2(h_star, h_dot, 0)
    y_fixed = Jet2(helical_ratio, helical_ratio, helical_ratio)
    omega_squared = Jet2(
        1, 2 * h_star, 2 * h_dot + 4 * h_star * h_star
    )
    b_fixed = omega_squared * locked_f(h_jet, sigma, y_fixed)
    assert b_fixed.first.lo < 0 < b_fixed.first.hi
    assert b_fixed.second.subset(F(-193), F(-192))

    y_resonant = h_jet / ((1 + h_jet) * sigma)
    full_edge = omega_squared * locked_f(h_jet, sigma, y_resonant)
    assert full_edge.first.lo < 0 < full_edge.first.hi
    assert full_edge.second.subset(F(-207), F(-205))

    y_variable = Jet2(helical_ratio, 1, 0)
    transverse = locked_f(Jet2(h_star), sigma, y_variable)
    assert transverse.first.lo < 0 < transverse.first.hi
    assert transverse.second.hi < 0

    # Omega(0)=-1, Omega_t=Omega*h, Omega_tt=Omega*(h^2+h_dot).
    omega = Jet2(-1, -h_star, -(h_star * h_star + h_dot))
    radius = Jet2(1, 1, 1)
    phase = radius * sigma * helical_ratio * omega - omega
    assert phase.first.lo < 0 < phase.first.hi
    assert phase.second.subset(F(753, 100), F(754, 100))
    assert phase.value.subset(F(-449, 1000), F(-448, 1000))

    # Rational substitution checks for q=-Gamma'/W' and
    # Phi=2 V Gamma'/r^2 in the two Rayleigh coefficients.
    cases = (
        (F(7, 5), F(11, 6), F(-5, 4), F(9, 7), F(-13, 8)),
        (F(2, 3), F(5, 2), F(17, 9), F(-4, 5), F(7, 3)),
        (F(13, 7), F(8, 5), F(3, 11), F(19, 6), F(5, 12)),
    )
    for beta, radius_value, velocity, w_prime, gamma_prime in cases:
        q_value = -gamma_prime / w_prime
        phi_value = 2 * velocity * gamma_prime / (radius_value * radius_value)
        original_b = (
            beta
            * radius_value**2
            * (1 - beta * q_value)
            * phi_value
            / (q_value * (1 + beta**2 * radius_value**2))
        )
        regular_b = (
            -2
            * beta
            * velocity
            * (w_prime + beta * gamma_prime)
            / (1 + beta**2 * radius_value**2)
        )
        assert original_b == regular_b
        original_a_inner = (
            (beta * radius_value**2 + q_value)
            * w_prime
            / (radius_value * (1 + beta**2 * radius_value**2))
        )
        regular_a_inner = (
            beta * radius_value**2 * w_prime - gamma_prime
        ) / (radius_value * (1 + beta**2 * radius_value**2))
        assert original_a_inner == regular_a_inner

    print("locked-pitch BAS formula and joint stationarity: PASS")
    print("unique global steep-lobe edge and mild-lobe gap: PASS")
    print("local fixed-beta b maximum and local Lambda minimum geometry: PASS")
    print("q-free Rayleigh coefficient rational substitutions: PASS")
    print("3.23013 < squared full BAS edge < 3.23014")
    print("all locked-pitch Gavrilov edge checks passed")


if __name__ == "__main__":
    main()
