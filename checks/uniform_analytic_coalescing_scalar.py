#!/usr/bin/env python3
"""Exact ledgers for the uniform analytic coalescing scalar packet.

This checker uses only rational and Gaussian-rational arithmetic.  It checks
the centered Taylor coefficients, dilation powers, phase sign, and stage
exponents.  It does not replace the analytic WKB majorant proof.
"""

from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class QI:
    """Gaussian rational a+ib."""

    re: F
    im: F = F(0)

    @staticmethod
    def make(value):
        if isinstance(value, QI):
            return value
        return QI(F(value))

    def __add__(self, other):
        other = self.make(other)
        return QI(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-self.make(other))

    def __rsub__(self, other):
        return self.make(other) - self

    def __mul__(self, other):
        other = self.make(other)
        return QI(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    __rmul__ = __mul__

    def inverse(self):
        denominator = self.re * self.re + self.im * self.im
        assert denominator
        return QI(self.re / denominator, -self.im / denominator)

    def __truediv__(self, other):
        return self * self.make(other).inverse()

    def __rtruediv__(self, other):
        return self.make(other) / self

    def __pow__(self, power):
        assert power >= 0
        answer = QI(F(1))
        base = self
        while power:
            if power & 1:
                answer *= base
            base *= base
            power //= 2
        return answer


def pad(series, degree):
    return list(series) + [QI(F(0))] * (degree + 1 - len(series))


def series_add(left, right, degree):
    left = pad(left, degree)
    right = pad(right, degree)
    return [left[k] + right[k] for k in range(degree + 1)]


def series_scale(value, series, degree):
    value = QI.make(value)
    series = pad(series, degree)
    return [value * series[k] for k in range(degree + 1)]


def series_mul(left, right, degree):
    left = pad(left, degree)
    right = pad(right, degree)
    return [
        sum(
            (left[k] * right[n - k] for k in range(n + 1)),
            QI(F(0)),
        )
        for n in range(degree + 1)
    ]


def series_inverse(series, degree):
    series = pad(series, degree)
    answer = [series[0].inverse()]
    for n in range(1, degree + 1):
        tail = sum(
            (series[k] * answer[n - k] for k in range(1, n + 1)),
            QI(F(0)),
        )
        answer.append(-tail / series[0])
    return answer


def series_divide(left, right, degree):
    return series_mul(left, series_inverse(right, degree), degree)


def exponent_add(left, right):
    """Monomial exponents in (s, eta)."""
    return (left[0] + right[0], left[1] + right[1])


def exponent_scale(power, value):
    return (power * value[0], power * value[1])


def check_exact_taylor_coefficients():
    degree = 3
    s = F(1, 100)
    sigma = F(1)
    r_c = F(1) + sigma * s
    y = F(19, 10)

    # Analytic model with the same required ring jets:
    # A=r^2/(1+r^2), b=4+(r-1)^2,
    # Lambda=(3/2)(r-1)^2+(1/3)(r-1)^3.
    A_c = r_c * r_c / (1 + r_c * r_c)
    b_c = F(4) + (r_c - 1) ** 2
    eta2 = (b_c / y**2 - 1) / A_c
    assert eta2 > 0

    r = [QI(r_c), QI(s * eta2)]
    one = [QI(F(1))]
    r_squared = series_mul(r, r, degree)
    A = series_divide(r_squared, series_add(one, r_squared, degree), degree)

    q = series_add(r, [QI(-1)], degree)
    q_squared = series_mul(q, q, degree)
    q_cubed = series_mul(q_squared, q, degree)
    b = series_add([QI(F(4))], q_squared, degree)
    lam = series_add(
        series_scale(F(3, 2), q_squared, degree),
        series_scale(F(1, 3), q_cubed, degree),
        degree,
    )

    # gamma=s^-2(Lambda(r(X))-Lambda(r_c))-iy.
    gamma = list(lam)
    gamma[0] = QI(F(0), -y)
    for k in range(1, degree + 1):
        gamma[k] = gamma[k] / s**2

    gamma_squared = series_mul(gamma, gamma, degree)
    W = series_scale(
        1 / eta2,
        series_add(one, series_divide(b, gamma_squared, degree), degree),
        degree,
    )

    # Exact centering.
    assert W[0] == QI(-A_c)

    b_prime_c = 2 * (r_c - 1)
    lambda_prime_c = 3 * (r_c - 1) + (r_c - 1) ** 2
    d_s = -s * b_prime_c / y**2
    c_s = 2 * b_c * lambda_prime_c / (s * y**3)
    assert W[1] == QI(d_s, c_s)
    assert c_s > 0

    # K=-W/A.  For xi_0=-1, Im(phi''(0)) is positive and exact.
    K = series_scale(-1, series_divide(W, A, degree), degree)
    assert K[0] == QI(F(1))
    xi_0 = -1
    phi_second = QI(xi_0) * K[1] / 2
    assert phi_second.im == c_s / (2 * A_c)
    assert phi_second.im > 0

    # Principal bracket at (0,xi_0) is nonzero.
    bracket = 2 * A_c * xi_0 * c_s
    assert bracket != 0


def check_dilation_powers():
    # Exponents are powers of (s, eta).
    h = (F(1), F(-3))
    second_after_dilation = (F(2), F(-6))
    first_after_dilation = (F(3), F(-4))
    lower_potential = (F(2), F(-2))

    assert exponent_scale(2, h) == second_after_dilation
    assert exponent_add(exponent_scale(3, h), (F(0), F(5))) == (
        first_after_dilation
    )
    assert exponent_add(exponent_scale(2, h), (F(0), F(4))) == (
        lower_potential
    )


def action_power(carrier_power, gain_power):
    return (F(carrier_power) - 3 * F(gain_power)) / 2


def check_stage_ledgers():
    g = F(2)
    assert action_power(8, g) == 1
    assert action_power(10, g) == 2
    assert (action_power(F(13, 2), g) > 0)  # strict A>6
    assert (action_power(6, g) == 0)
    # Beating arbitrary fixed-order endpoint losses exp(C_N*j) requires
    # action power > 1: A>3g+2.  A=8 is the critical, insufficient edge.
    assert action_power(8, g) == 1
    assert action_power(F(17, 2), g) > 1

    # Long-gain curvature and straight-viscosity powers.  If
    # epsilon=j^-E and mu=j^-M, they decay precisely under
    # E>A+g and M>2A+g.
    A = F(8)
    E = F(11)
    M = F(19)
    curvature_decay = E - A - g
    viscosity_decay = M - 2 * A - g
    assert curvature_decay == 1
    assert viscosity_decay == 1

    # The cascade specialization epsilon=ell^beta and
    # mu=nu*ell^(gamma-1)*p^(-gamma) has positive geometric powers
    # under 0<beta<(gamma-1)/2.  Polynomial powers cannot defeat them.
    similarity_gamma = F(5, 4)
    beta = F(1, 10)
    assert 0 < beta < (similarity_gamma - 1) / 2
    straight_viscous_ell_power = similarity_gamma - 1
    curved_viscous_ell_power = beta + similarity_gamma - 1
    straight_viscous_j_power = A * (2 - similarity_gamma) + g
    assert straight_viscous_ell_power > 0
    assert curved_viscous_ell_power > straight_viscous_ell_power
    assert straight_viscous_j_power == 8


def main():
    check_exact_taylor_coefficients()
    check_dilation_powers()
    check_stage_ledgers()
    print("exact centered scalar Taylor coefficients: PASS")
    print("nonzero bracket and decaying Gaussian phase sign: PASS")
    print("fixed-X dilation powers: PASS")
    print("long-gain, endpoint-flatness, and A=8/A=10 rates: PASS")
    print("curvature, viscosity, and geometric stage ledgers: PASS")
    print("all uniform coalescing scalar checks passed")


if __name__ == "__main__":
    main()
