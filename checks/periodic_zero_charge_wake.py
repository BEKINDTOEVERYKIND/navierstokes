#!/usr/bin/env python3
"""Dependency-free exact checks for periodic zero-charge routing."""

from fractions import Fraction as F
from math import factorial


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), F(0))


def real_symbol(k_value, f_value):
    """Return A with Fourier stress equal to -i*A."""

    modulus = dot(k_value, k_value)
    pairing = dot(k_value, f_value)
    return [
        [
            (
                k_value[i] * f_value[j]
                + f_value[i] * k_value[j]
                - (pairing if i == j else 0)
            )
            / modulus
            for j in range(3)
        ]
        for i in range(3)
    ]


def main():
    samples = (
        ((1, 0, 0), (0, 2, -3)),
        ((1, -2, 3), (4, 5, -1)),
        ((7, 5, -4), (-3, 2, 9)),
        ((11, -13, 17), (19, 23, -29)),
    )

    for raw_k, raw_f in samples:
        k_value = tuple(F(value) for value in raw_k)
        f_value = tuple(F(value) for value in raw_f)
        stress = real_symbol(k_value, f_value)
        assert all(stress[i][j] == stress[j][i] for i in range(3) for j in range(3))

        # Since Rhat=-i*A, multiplication by i*k in div R gives A*k.
        divergence = tuple(
            sum((stress[i][j] * k_value[j] for j in range(3)), F(0))
            for i in range(3)
        )
        assert divergence == f_value

        pairing = dot(k_value, f_value)
        modulus = dot(k_value, k_value)
        assert sum((stress[i][i] for i in range(3)), F(0)) == -pairing / modulus

        # The multiplier is exactly homogeneous of order minus one.
        for dilation in (2, 5, 19):
            dilated = real_symbol(
                tuple(F(dilation) * value for value in k_value), f_value
            )
            assert all(
                dilated[i][j] == stress[i][j] / dilation
                for i in range(3)
                for j in range(3)
            )

    # A solenoidal datum gives an exactly trace-free stress.
    solenoidal_k = (F(2), F(-3), F(5))
    solenoidal_f = (F(3), F(2), F(0))
    assert dot(solenoidal_k, solenoidal_f) == 0
    solenoidal_stress = real_symbol(solenoidal_k, solenoidal_f)
    assert sum((solenoidal_stress[i][i] for i in range(3)), F(0)) == 0

    # Taylor coefficients of the one-mode Duhamel solution solve
    # z' + a*z = F exactly: c_n=F*(-a)^(n-1)/n!, n>=1.
    damping = F(7, 3)
    forcing = F(11, 5)
    order = 16
    coefficients = [F(0)] + [
        forcing * (-damping) ** (n - 1) / factorial(n)
        for n in range(1, order + 1)
    ]
    assert coefficients[0] == 0
    assert coefficients[1] == forcing
    for power in range(1, order):
        assert (power + 1) * coefficients[power + 1] + damping * coefficients[power] == 0

    print("periodic symmetric anti-divergence symbol: PASS")
    print("divergence, symmetry, and trace-free solenoidal case: PASS")
    print("exact order-minus-one scaling: PASS")
    print("single-mode heat-wake coefficient recursion: PASS")
    print("all periodic zero-charge wake checks passed")


if __name__ == "__main__":
    main()
