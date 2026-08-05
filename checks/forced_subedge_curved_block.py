#!/usr/bin/env python3
"""Exact algebra checks for the recovered forced/curved spectral block.

The analytic arguments are recorded in
``research/2026-08-03-recovered-subedge-curved-block.md``.  This checker
uses Gaussian rationals throughout; no floating-point comparison decides a
claim.
"""

from fractions import Fraction as F
from itertools import permutations


class C:
    """A Gaussian rational number."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        if isinstance(re, C):
            self.re, self.im = re.re, re.im
        else:
            self.re, self.im = F(re), F(im)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, C) else C(value)

    def __add__(self, other):
        other = self.coerce(other)
        return C(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return C(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    __rmul__ = __mul__

    def conjugate(self):
        return C(self.re, -self.im)

    def __truediv__(self, other):
        other = self.coerce(other)
        norm = other.re * other.re + other.im * other.im
        if norm == 0:
            raise ZeroDivisionError
        return self * other.conjugate() * C(1 / norm)

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, exponent):
        if exponent == 2:
            return self * self
        if exponent == -1:
            return C(1) / self
        raise ValueError("only powers 2 and -1 are needed")

    def __eq__(self, other):
        other = self.coerce(other)
        return self.re == other.re and self.im == other.im

    def __repr__(self):
        return f"C({self.re}, {self.im})"


I = C(0, 1)
ZERO = C(0)
ONE = C(1)


def determinant(matrix):
    size = len(matrix)
    answer = ZERO
    for perm in permutations(range(size)):
        inversions = sum(
            perm[j] > perm[k]
            for j in range(size)
            for k in range(j + 1, size)
        )
        term = ONE if inversions % 2 == 0 else -ONE
        for row, column in enumerate(perm):
            term *= matrix[row][column]
        answer += term
    return answer


def matmul(left, right):
    return [
        [
            sum(
                (C.coerce(left[i][k]) * C.coerce(right[k][j])
                 for k in range(len(right))),
                ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def check_schur_factorization():
    cases = (
        # gamma, r, beta, eta, Omega, Gamma'/r, W'
        (C(F(3, 2), F(-5, 4)), F(7, 5), F(4, 3), F(2, 5),
         F(5, 6), F(-7, 4), F(9, 8)),
        (C(F(-2, 3), F(-11, 7)), F(9, 4), F(3, 5), F(-8, 7),
         F(-4, 9), F(13, 6), F(-5, 3)),
        (C(F(5, 8), F(7, 6)), F(11, 6), F(7, 9), F(6, 5),
         F(8, 7), F(2, 3), F(5, 11)),
    )

    for gamma, radius, beta, eta, omega, a_theta, w_prime in cases:
        matrix = [
            [I * gamma, -2 * omega, 0, I * eta],
            [a_theta, I * gamma, 0, -I / radius],
            [w_prime, 0, I * gamma, I * beta],
            [eta, -1 / radius, beta, 0],
        ]
        complement = [
            [I * gamma, 0, -I / radius],
            [0, I * gamma, I * beta],
            [-1 / radius, beta, 0],
        ]

        d_value = 1 + beta * beta * radius * radius
        complement_expected = gamma * d_value / (radius * radius)
        assert determinant(complement) == complement_expected

        gamma_prime = radius * a_theta
        velocity = radius * omega
        omega_prime = (a_theta - 2 * omega) / radius
        lambda_prime = beta * w_prime - omega_prime
        b_value = (
            -2
            * beta
            * velocity
            * (w_prime + beta * gamma_prime)
            / d_value
        )
        radial_schur = (
            I
            * gamma
            * (
                1
                + radius * radius * eta * eta / d_value
                + b_value / (gamma * gamma)
            )
            - radius * radius * eta * lambda_prime / d_value
        )
        assert determinant(matrix) == complement_expected * radial_schur


def check_forced_pressure_correction():
    cases = (
        (F(7, 5), F(4, 3), 11, F(5, 6), F(-3, 7)),
        (F(9, 4), F(2, 5), 17, F(-8, 9), F(11, 6)),
        (F(5, 3), F(7, 8), 23, F(13, 10), F(2, 9)),
    )
    gamma = C(F(3, 4), F(-5, 7))
    for radius, beta, n_value, f_theta, f_z in cases:
        d_value = 1 + beta * beta * radius * radius
        q_force = I * radius * (f_theta - beta * radius * f_z) / (
            n_value * d_value
        )
        # This is exactly the force-only part of the Fourier divergence
        # equation after the theta and z momentum equations are solved.
        divergence_force = (
            -I * n_value * n_value * d_value * q_force
            / (radius * radius * gamma)
            - n_value * f_theta / (radius * gamma)
            + beta * n_value * f_z / gamma
        )
        assert divergence_force == ZERO


def check_hodge_symbol():
    cases = (
        (F(7, 5), F(6, 5), (F(2), F(-3), F(5, 2))),
        (F(9, 4), F(11, 10), (F(-4, 3), F(7, 5), F(8, 7))),
        (F(5, 3), F(13, 9), (F(9, 8), F(2, 11), F(-5, 6))),
    )
    for radius, h_value, xi in cases:
        inverse_metric = (F(1), 1 / (radius * radius), 1 / (h_value * h_value))
        sharp = tuple(inverse_metric[j] * xi[j] for j in range(3))
        q_value = sum((xi[j] * sharp[j] for j in range(3)), F(0))
        projector = [
            [
                C(F(int(i == j)) - sharp[i] * xi[j] / q_value)
                for j in range(3)
            ]
            for i in range(3)
        ]
        assert matmul(projector, projector) == projector
        for column in range(3):
            annihilated = sum(
                (C(xi[row]) * projector[row][column] for row in range(3)),
                ZERO,
            )
            assert annihilated == ZERO

        # Exact factorization of the only curved metric coefficient.
        epsilon = h_value - 1
        assert 1 / (h_value * h_value) - 1 == (
            -epsilon * (2 + epsilon) / (h_value * h_value)
        )


def check_characteristic_and_bracket():
    # A=1, b=5, y=1 makes eta_0=2 exactly.
    a_value = F(1)
    b_value = F(5)
    y_value = F(1)
    eta = F(2)
    kappa = F(3)
    sigma = F(2)
    assert a_value * eta * eta + 1 - b_value / (y_value * y_value) == 0

    derivative_eta_real = 2 * a_value * eta
    derivative_y_imag = 2 * b_value * kappa * sigma / (y_value**3)
    bracket = derivative_eta_real * derivative_y_imag
    expected = (
        4 * a_value * eta * b_value * kappa * sigma / (y_value**3)
    )
    assert bracket == expected == 240


def main():
    check_schur_factorization()
    check_forced_pressure_correction()
    check_hodge_symbol()
    check_characteristic_and_bracket()
    print("forced pressure correction and divergence cancellation: PASS")
    print("4x4 velocity-pressure determinant and scalar Schur factor: PASS")
    print("thin-torus principal Hodge projector identities: PASS")
    print("fixed-subedge characteristic and Poisson bracket: PASS")
    print("all forced/curved block checks passed")


if __name__ == "__main__":
    main()
