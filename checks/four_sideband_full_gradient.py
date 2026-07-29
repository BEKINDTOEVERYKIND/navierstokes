#!/usr/bin/env python3
"""Exact full-gradient chart for four common-carrier sidebands.

The preferred one-carrier construction uses k=K e3, a=e1 and partners
r=q-Ke3.  The exact child map is onto q-perp for every listed q and every
positive integer K.  This checker verifies that four such child planes
span all of sl(3), and gives an exact periodic synthesis of the affine
Kelvin strain diag(-5/4, 9/4, -1).

Only the Python standard library is used.
"""

from __future__ import annotations

from fractions import Fraction as F


Vector = tuple[F, F, F]

QS: tuple[Vector, ...] = (
    (F(-45), F(-36), F(20)),
    (F(-4), F(-5), F(9)),
    (F(1), F(1), F(1)),
    (F(1), F(2), F(3)),
)

# Coordinates of a trace-free matrix with M33 omitted.
MATRIX_COORDINATES = (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
)


def transverse_basis(q: Vector) -> tuple[Vector, Vector]:
    alpha, beta, delta = q
    return (
        (-beta, alpha, F(0)),
        (-delta, F(0), alpha),
    )


def gradient_column(c: Vector, q: Vector) -> list[F]:
    return [c[i] * q[j] for i, j in MATRIX_COORDINATES]


def determinant(matrix: list[list[F]]) -> F:
    work = [row[:] for row in matrix]
    out = F(1)
    for column in range(len(work)):
        pivot = next(
            row
            for row in range(column, len(work))
            if work[row][column] != 0
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            out = -out
        pivot_value = work[column][column]
        out *= pivot_value
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] / pivot_value
            for entry in range(column + 1, len(work)):
                work[row][entry] -= multiplier * work[column][entry]
    return out


def solve(matrix: list[list[F]], target: list[F]) -> list[F]:
    size = len(matrix)
    work = [matrix[row][:] + [target[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            row
            for row in range(column, size)
            if work[row][column] != 0
        )
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = work[row][column]
            work[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return [work[row][-1] for row in range(size)]


def child_block_determinant(q: Vector, k_value: int) -> F:
    """Determinant of the exact partner-to-child map in natural bases."""

    alpha, _, delta = q
    q_squared = sum((entry * entry for entry in q), F(0))
    return alpha * alpha * (q_squared - F(2 * k_value) * delta) / q_squared


def main() -> None:
    relation = (-15, 47, -669, 182)
    assert all(
        sum(
            F(relation[index]) * QS[index][coordinate]
            for index in range(4)
        )
        == 0
        for coordinate in range(3)
    )

    columns = [
        gradient_column(c, q)
        for q in QS
        for c in transverse_basis(q)
    ]
    matrix = [[columns[column][row] for column in range(8)] for row in range(8)]
    chart_determinant = determinant(matrix)
    assert chart_determinant == F(15_451_090_200)

    # These are all exceptional K=|q|^2/(2 q3) values.  Their
    # nontrivial denominators prove, without a finite search cutoff, that
    # every child block is nonsingular for every integer carrier K.
    exceptional_carriers = tuple(
        sum((entry * entry for entry in q), F(0)) / (2 * q[2])
        for q in QS
    )
    assert exceptional_carriers == (
        F(3721, 40),
        F(61, 9),
        F(3, 2),
        F(7, 3),
    )
    assert all(value.denominator != 1 for value in exceptional_carriers)

    gamma = F(5, 4)
    target = [
        -gamma,
        F(0),
        F(0),
        F(0),
        gamma + 1,
        F(0),
        F(0),
        F(0),
    ]
    coefficients = solve(matrix, target)
    expected = [
        F(-169, 73367),
        F(-17027, 9537710),
        F(9753, 62440),
        F(193829, 2435160),
        F(9693, 6580),
        F(144689, 42770),
        F(145119, 209620),
        F(-163724, 157215),
    ]
    assert coefficients == expected

    full = [[F(0) for _ in range(3)] for _ in range(3)]
    index = 0
    for q in QS:
        for c in transverse_basis(q):
            for i in range(3):
                for j in range(3):
                    full[i][j] += coefficients[index] * c[i] * q[j]
            index += 1
    assert full == [
        [-gamma, F(0), F(0)],
        [F(0), gamma + 1, F(0)],
        [F(0), F(0), F(-1)],
    ]

    print("four-sideband common-carrier gradient chart: exact rank 8")
    print("chart determinant:", chart_determinant)
    print("finite-K child blocks: rank 2 for every positive integer K")
    print("exact affine target: diag(-5/4, 9/4, -1)")
    print("slow quotient: -15 q1+47 q2-669 q3+182 q4=0")


if __name__ == "__main__":
    main()
