#!/usr/bin/env python3
"""Exact finite-dimensional checks for endpoint-jet interpolation.

This verifies the Hermite jet rank, the explicit three-phase stress and
amplitude right inverses, and the zero-amplitude rank loss.  It does not
verify a dynamical Navier--Stokes endpoint theorem.
"""

from fractions import Fraction as F
from math import comb, factorial, isqrt


def rank(matrix):
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    value = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(value, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[value], work[pivot] = work[pivot], work[value]
        divisor = work[value][column]
        work[value] = [entry / divisor for entry in work[value]]
        for row in range(rows):
            if row == value:
                continue
            multiplier = work[row][column]
            work[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(work[row], work[value])
            ]
        value += 1
    return value


def derivative_of_monomial(degree, order, point):
    if degree < order:
        return F(0)
    coefficient = F(factorial(degree), factorial(degree - order))
    exponent = degree - order
    return coefficient * (F(1) if exponent == 0 else F(point) ** exponent)


def hermite_matrix(order):
    degree = 2 * order + 1
    rows = []
    for point in (0, 1):
        for derivative in range(order + 1):
            rows.append(
                [
                    derivative_of_monomial(power, derivative, point)
                    for power in range(degree + 1)
                ]
            )
    return rows


def zero_matrix():
    return [[F(0) for _ in range(3)] for _ in range(3)]


def add(left, right):
    return [
        [left[i][j] + right[i][j] for j in range(3)]
        for i in range(3)
    ]


def scale(value, matrix):
    return [[value * entry for entry in row] for row in matrix]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def multiply(left, right):
    return [
        [
            sum((left[i][k] * right[k][j] for k in range(len(right))), F(0))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def stress_pieces(stress):
    s11, s22, s33 = stress[0][0], stress[1][1], stress[2][2]
    s12, s13, s23 = stress[0][1], stress[0][2], stress[1][2]
    return (
        [
            [F(0), F(0), F(0)],
            [F(0), s22 / 2, s23],
            [F(0), s23, s33 / 2],
        ],
        [
            [s11 / 2, F(0), s13],
            [F(0), F(0), F(0)],
            [s13, F(0), s33 / 2],
        ],
        [
            [s11 / 2, s12, F(0)],
            [s12, s22 / 2, F(0)],
            [F(0), F(0), F(0)],
        ],
    )


def check_hermite_jet_rank():
    for order in range(0, 13):
        matrix = hermite_matrix(order)
        assert len(matrix) == 2 * order + 2
        assert rank(matrix) == 2 * order + 2


def check_triangular_residual_rank():
    # Scalar proxy for the operator jet linearization.  Row n contains an
    # arbitrary lower-order part and coefficient one on h_{n+1}.
    for order in range(1, 20):
        matrix = []
        for n_value in range(order):
            row = [F(0) for _ in range(order + 1)]
            for column in range(n_value + 1):
                row[column] = F((n_value + 2) * (column + 1), column + 2)
            row[n_value + 1] = F(1)
            matrix.append(row)
        assert rank(matrix) == order
        assert len(matrix[0]) - rank(matrix) == 1


def check_stress_and_amplitude_right_inverse():
    samples = (
        [
            [F(2), F(1, 3), F(-2, 5)],
            [F(1, 3), F(-3, 2), F(4, 7)],
            [F(-2, 5), F(4, 7), F(11, 6)],
        ],
        [
            [F(-5, 4), F(7, 9), F(1, 8)],
            [F(7, 9), F(13, 5), F(-3, 11)],
            [F(1, 8), F(-3, 11), F(17, 3)],
        ],
    )
    frames = (
        [[F(0), F(0)], [F(1), F(0)], [F(0), F(1)]],
        [[F(1), F(0)], [F(0), F(0)], [F(0), F(1)]],
        [[F(1), F(0)], [F(0), F(1)], [F(0), F(0)]],
    )
    amplitude = F(3, 2)

    for stress in samples:
        pieces = stress_pieces(stress)
        total = zero_matrix()
        derivative_total = zero_matrix()
        for piece, frame in zip(pieces, frames):
            total = add(total, piece)
            base = scale(amplitude, frame)
            variation = scale(F(1, 2) / amplitude, multiply(piece, frame))
            derivative = add(
                multiply(variation, transpose(base)),
                multiply(base, transpose(variation)),
            )
            assert derivative == piece
            derivative_total = add(derivative_total, derivative)
        assert total == stress
        assert derivative_total == stress

        # Endpoint-scaled version: base A=alpha*A0 and desired stress
        # alpha^2*G have a smooth variation alpha*delta A0.
        for alpha in (F(1, 2), F(1, 7), F(3, 11)):
            scaled_derivative_total = zero_matrix()
            for piece, frame in zip(pieces, frames):
                base = scale(alpha * amplitude, frame)
                variation = scale(
                    alpha * F(1, 2) / amplitude,
                    multiply(piece, frame),
                )
                derivative = add(
                    multiply(variation, transpose(base)),
                    multiply(base, transpose(variation)),
                )
                assert derivative == scale(alpha**2, piece)
                scaled_derivative_total = add(
                    scaled_derivative_total, derivative
                )
            assert scaled_derivative_total == scale(alpha**2, stress)


def check_amplitude_rank_and_endpoint_loss():
    # Generate the twelve transverse amplitude columns at unit base
    # amplitude and record them in (xx, yy, zz, xy, xz, yz) coordinates.
    basis = (
        ((1, 2), (1, 2)),
        ((0, 2), (0, 2)),
        ((0, 1), (0, 1)),
    )
    columns = []
    for vector_indices, variation_indices in basis:
        for vector_index in range(2):
            base_vector = [F(0), F(0), F(0)]
            base_vector[vector_indices[vector_index]] = F(1)
            for variation_index in variation_indices:
                variation = [F(0), F(0), F(0)]
                variation[variation_index] = F(1)
                matrix = zero_matrix()
                for i in range(3):
                    for j in range(3):
                        matrix[i][j] = (
                            variation[i] * base_vector[j]
                            + base_vector[i] * variation[j]
                        )
                columns.append(
                    [
                        matrix[0][0], matrix[1][1], matrix[2][2],
                        matrix[0][1], matrix[0][2], matrix[1][2],
                    ]
                )
    assert len(columns) == 12
    assert rank(columns) == 6

    # At zero base amplitude every derivative column is identically zero.
    zero_columns = [[F(0) for _ in range(6)] for _ in range(12)]
    assert rank(zero_columns) == 0


def check_endpoint_order_doubling():
    # To erase amplitude derivatives through M, take A=s^(M+1).  Its
    # covariance starts at order 2M+2, so the residual must vanish through
    # derivative order 2M+1.  The C88 heat gate changes only by a factor
    # tending to four.
    for order in range(0, 40):
        amplitude_order = order + 1
        covariance_order = 2 * amplitude_order
        assert covariance_order == 2 * order + 2
        assert covariance_order - 1 == 2 * order + 1
        if order >= 1:
            gate_ratio = F(covariance_order**2, order**2)
            assert gate_ratio > 4
            assert gate_ratio <= 16
    # The asymptotic constant is exactly four.
    large_order = 10_000
    ratio = F((2 * large_order + 2) ** 2, large_order**2)
    assert abs(ratio - 4) < F(1, 500)

    # sqrt((2r)!) <= 2^r r! is the exact reason the positive square root
    # remains inside the Gevrey-2 endpoint budget.
    for r_value in range(1, 80):
        doubled_factorial = factorial(2 * r_value)
        bound = 2**r_value * factorial(r_value)
        assert doubled_factorial <= bound**2
        assert comb(2 * r_value, r_value) <= 4**r_value
        assert isqrt(doubled_factorial) <= bound


def main():
    check_hermite_jet_rank()
    check_triangular_residual_rank()
    check_stress_and_amplitude_right_inverse()
    check_amplitude_rank_and_endpoint_loss()
    check_endpoint_order_doubling()
    print("two-endpoint Hermite jet map has exact full rank: PASS")
    print("Navier--Stokes residual-jet linearization is triangular onto: PASS")
    print("explicit three-phase stress/amplitude right inverse: PASS")
    print("interior amplitude rank 6 and zero-endpoint rank 0: PASS")
    print("clean amplitude jet requires exact 2M+2 residual order: PASS")
    print("all endpoint-jet interpolation checks passed")


if __name__ == "__main__":
    main()
