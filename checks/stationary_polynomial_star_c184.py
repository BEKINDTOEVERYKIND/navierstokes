#!/usr/bin/env python3
"""Exact checks for C184's stationary-polynomial collateral obstruction.

All coefficient and Gram calculations use Fraction arithmetic.  The
script checks the algebraic ledgers; it does not certify a finite-frequency
localized Navier--Stokes propagator.
"""

from __future__ import annotations

from fractions import Fraction as F


Index = tuple[int, int]
Laurent = dict[Index, F]

ALPHA = F(2, 5)
BETA = F(1, 2)


def add_poly(left: Laurent, right: Laurent) -> Laurent:
    out = dict(left)
    for index, coefficient in right.items():
        out[index] = out.get(index, F(0)) + coefficient
        if out[index] == 0:
            del out[index]
    return out


def scale_poly(scalar: F, value: Laurent) -> Laurent:
    return {
        index: scalar * coefficient
        for index, coefficient in value.items()
        if scalar * coefficient
    }


def convolve(left: Laurent, right: Laurent) -> Laurent:
    out: Laurent = {}
    for (j, k), left_coefficient in left.items():
        for (r, s), right_coefficient in right.items():
            index = (j + r, k + s)
            out[index] = out.get(index, F(0)) + left_coefficient * right_coefficient
    return {index: coefficient for index, coefficient in out.items() if coefficient}


def convolution_powers(maximum: int) -> list[Laurent]:
    measure: Laurent = {
        (1, 0): BETA,
        (-1, 0): BETA,
        (0, 1): BETA,
        (0, -1): BETA,
        (1, 1): ALPHA,
        (-1, -1): ALPHA,
    }
    powers: list[Laurent] = [{(0, 0): F(1)}]
    for _ in range(maximum):
        powers.append(convolve(powers[-1], measure))
    return powers


def hex_radius(index: Index) -> int:
    j, k = index
    return max(abs(j), abs(k), abs(j - k))


def exact_coefficient_space() -> None:
    powers = convolution_powers(10)
    for degree, coefficients in enumerate(powers):
        assert all(hex_radius(index) <= degree for index in coefficients)
        for (j, k), coefficient in coefficients.items():
            assert coefficient == coefficients.get((k, j), F(0))
            assert coefficient == coefficients.get((-j, -k), F(0))
            assert coefficient == coefficients.get((-k, -j), F(0))

        for radius in range(1, 11):
            if degree < radius:
                assert coefficients.get((radius, radius), F(0)) == 0
        if degree > 0:
            assert coefficients[(degree, degree)] == ALPHA**degree
            assert coefficients[(degree, 0)] == BETA**degree
            assert coefficients[(0, degree)] == BETA**degree

    # The nonzero triangular diagonal proves linear independence of
    # C_0,...,C_q without a floating rank computation.
    for degree in range(0, 11):
        pivot = (degree, degree) if degree else (0, 0)
        assert powers[degree][pivot] != 0
        assert all(powers[lower].get(pivot, F(0)) == 0 for lower in range(degree))


def diagonal_interpolation() -> None:
    degree = 10
    first_radius = 6
    powers = convolution_powers(degree)
    derivative_at_zero = F(-7, 5)
    desired = {
        radius: F(3 * radius - 7, 11 * radius + 5)
        for radius in range(first_radius, degree + 1)
    }
    polynomial = {0: F(0), 1: derivative_at_zero}
    for radius in range(degree, first_radius - 1, -1):
        already = sum(
            polynomial.get(power, F(0))
            * powers[power].get((radius, radius), F(0))
            for power in range(radius + 1, degree + 1)
        )
        polynomial[radius] = (desired[radius] - already) / (ALPHA**radius)
    for power in range(2, first_radius):
        polynomial[power] = F(0)

    fourier: Laurent = {}
    for power, coefficient in polynomial.items():
        fourier = add_poly(fourier, scale_poly(coefficient, powers[power]))
    assert polynomial[0] == 0
    assert polynomial[1] == derivative_at_zero
    for radius, target in desired.items():
        assert fourier[(radius, radius)] == target
        assert fourier[(-radius, -radius)] == target

    # H(t)-c*t has zero constant and linear terms: it is exactly t^2 S(t).
    remainder = dict(polynomial)
    remainder[1] -= derivative_at_zero
    assert remainder[0] == 0 and remainder[1] == 0


def outer_face_collateral() -> None:
    powers = convolution_powers(12)
    for degree in range(2, 13):
        # The fixed c*f jet has support radius one and cannot cancel this
        # outer face when degree >= 2.
        assert all(
            powers[1].get((degree, second), F(0)) == 0
            for second in range(degree + 1)
        )
        outer_target = F(3, 7 * degree)
        leading = outer_target / (ALPHA**degree)
        face = []
        for second in range(degree + 1):
            coefficient = leading * powers[degree][(degree, second)]
            expected = (
                leading
                * F(combination(degree, second))
                * BETA ** (degree - second)
                * ALPHA**second
            )
            assert coefficient == expected
            face.append(coefficient)
        assert face[-1] == outer_target
        collateral = face[:-1]
        collateral_l1 = sum((abs(value) for value in collateral), F(0))
        exact_l1 = abs(outer_target) * ((F(9, 4) ** degree) - 1)
        assert collateral_l1 == exact_l1
        collateral_l2_squared = sum((value * value for value in collateral), F(0))
        assert F(degree) * collateral_l2_squared >= collateral_l1**2

        # Every physical face frequency degree*r1+second*r2 has
        # |g|^2=2(degree^2+second^2-degree*second)>=3 degree^2/2.
        for second, coefficient in enumerate(collateral):
            wavevector_squared = 2 * (
                degree * degree + second * second - degree * second
            )
            assert F(wavevector_squared) >= F(3, 2) * degree * degree
            assert coefficient != 0


def combination(n: int, k: int) -> int:
    k = min(k, n - k)
    output = 1
    for value in range(1, k + 1):
        output = output * (n - k + value) // value
    return output


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(scalar, vector):
    return tuple(scalar * value for value in vector)


def project(vector, covector):
    return add(vector, scale(-dot(vector, covector) / dot(covector, covector), covector))


def edge(source, gate, normal, amplitude):
    charge = dot(normal, source)
    raw = add(scale(charge, amplitude), scale(dot(amplitude, gate), normal))
    return project(raw, add(source, gate))


def selected_tangent_gram() -> None:
    # Exact separated model cone: p=(0,q,2q), n=e3,
    # g=(+/-r,0,0), q/2<=r<=q.  The reality pair cancels the Gram
    # cross-term.  With tau=1/(q sqrt(Q)), the generalized Gram relative
    # to the input basis a1=(1,0,0), a2=(0,-2,1) has eigenvalues
    # in [7/3,7/2] and exactly 4, respectively.  Work with tau^2 to
    # avoid irrational arithmetic.
    normal = (F(0), F(0), F(1))
    for q in (8, 12, 20, 32):
        radii = tuple(range(q // 2, q + 1))
        signed_count = 2 * len(radii)
        tau_squared = F(1, q * q * signed_count)
        source = (F(0), F(q), F(2 * q))
        basis = ((F(1), F(0), F(0)), (F(0), F(-2), F(1)))
        assert all(dot(source, amplitude) == 0 for amplitude in basis)
        gram = [[F(0), F(0)], [F(0), F(0)]]
        for radius in radii:
            for sign in (-1, 1):
                gate = (F(sign * radius), F(0), F(0))
                outputs = tuple(edge(source, gate, normal, amplitude) for amplitude in basis)
                for row in range(2):
                    for column in range(2):
                        gram[row][column] += tau_squared * dot(outputs[row], outputs[column])
        assert gram[0][1] == gram[1][0] == 0
        first_generalized = gram[0][0] / dot(basis[0], basis[0])
        second_generalized = gram[1][1] / dot(basis[1], basis[1])
        assert F(7, 3) <= first_generalized <= F(7, 2)
        assert second_generalized == 4


def factorial_does_not_beat_collateral() -> None:
    # The dimensionless initial-rate ledger is
    # (j!)^-2 (j+1)^-2 [(9/4)^((j+1)^8)-1].
    # Exact rational comparisons already become enormous and then increase.
    # Use the exact elementary bounds
    #   (9/4)^q-1 > 2^q,
    #   (j!)^2 (j+1)^2 <= (j+1)^(2j+2) <= 2^(2j^2+2j).
    # It is enough to check the resulting integer exponent tends upward;
    # no enormous rational power need be materialized.
    previous_exponent = None
    for stage in range(2, 101):
        q = (stage + 1) ** 8
        assert stage + 1 <= 2**stage
        exponent = q - 2 * stage * stage - 2 * stage
        assert exponent > 0
        if previous_exponent is not None:
            assert exponent > previous_exponent
        previous_exponent = exponent
    assert previous_exponent is not None and previous_exponent > 10**15


def main() -> None:
    exact_coefficient_space()
    diagonal_interpolation()
    outer_face_collateral()
    selected_tangent_gram()
    factorial_does_not_beat_collateral()
    print("C184 stationary-polynomial collateral checks passed")


if __name__ == "__main__":
    main()
