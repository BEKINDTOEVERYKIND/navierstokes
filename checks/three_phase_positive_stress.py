#!/usr/bin/env python3
"""Exact checks for corrected two- and three-phase positive-stress charts.

The script checks fixed-direction ranks, the variable-kernel submersion,
and rational positive decompositions.  It does not construct material
phases or evolve Euler/Navier--Stokes.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def rank(matrix: list[list[Fraction]]) -> int:
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


def check_two_phase_rank_loss() -> None:
    # Symmetric coordinates: xx, yy, zz, xy, xz, yz.
    s_e1 = [
        [Q(0), Q(1), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(1), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(0), Q(1)],
    ]
    s_e2 = [
        [Q(1), Q(0), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(1), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(1), Q(0)],
    ]
    combined = s_e1 + s_e2
    assert rank(combined) == 5

    # Trace-free coordinates: xx-zz, yy-zz, xy, xz, yz.
    trace_free = [
        [row[0] - row[2], row[1] - row[2], row[3], row[4], row[5]]
        for row in combined
    ]
    assert rank(trace_free) == 4
    print("two orthogonal phases: symmetric rank 5, trace-free rank 4")


def check_variable_two_phase_submersion() -> None:
    # Baseline A=diag(0,q,a), B=diag(q,0,q-a), with q=2,a=1.
    # Coordinates are xx, yy, zz, xy, xz, yz.  Five block variations plus
    # a rotation of ker(A) from e1 toward e2 span all of Sym_3.
    columns = [
        [Q(1), Q(0), Q(0), Q(0), Q(0), Q(0)],   # B_xx
        [Q(0), Q(1), Q(0), Q(0), Q(0), Q(0)],   # A_yy
        [Q(0), Q(0), Q(1), Q(0), Q(0), Q(0)],   # A_zz
        [Q(0), Q(0), Q(0), Q(-2), Q(0), Q(0)],  # delta k=e2
        [Q(0), Q(0), Q(0), Q(0), Q(1), Q(0)],   # B_xz
        [Q(0), Q(0), Q(0), Q(0), Q(0), Q(1)],   # A_yz
    ]
    assert rank(columns) == 6
    print("two variable phases: full symmetric rank 6")


def decompose(matrix: list[list[Fraction]]):
    q11, q22, q33 = matrix[0][0], matrix[1][1], matrix[2][2]
    q12, q13, q23 = matrix[0][1], matrix[0][2], matrix[1][2]
    r1 = [
        [Q(0), Q(0), Q(0)],
        [Q(0), q22 / 2, q23],
        [Q(0), q23, q33 / 2],
    ]
    r2 = [
        [q11 / 2, Q(0), q13],
        [Q(0), Q(0), Q(0)],
        [q13, Q(0), q33 / 2],
    ]
    r3 = [
        [q11 / 2, q12, Q(0)],
        [q12, q22 / 2, Q(0)],
        [Q(0), Q(0), Q(0)],
    ]
    return r1, r2, r3


def positive_two_by_two(block: list[list[Fraction]]) -> bool:
    return block[0][0] > 0 and (
        block[0][0] * block[1][1] - block[0][1] ** 2 > 0
    )


def check_three_phase_decomposition() -> None:
    matrices = [
        [
            [Q(1), Q(1, 10), Q(-1, 8)],
            [Q(1, 10), Q(9, 10), Q(1, 12)],
            [Q(-1, 8), Q(1, 12), Q(11, 10)],
        ],
        [
            [Q(4), Q(-1, 2), Q(1, 3)],
            [Q(-1, 2), Q(9, 2), Q(-1, 4)],
            [Q(1, 3), Q(-1, 4), Q(15, 4)],
        ],
    ]
    for matrix in matrices:
        pieces = decompose(matrix)
        total = [
            [sum(piece[i][j] for piece in pieces) for j in range(3)]
            for i in range(3)
        ]
        assert total == matrix
        assert all(pieces[i][i][j] == 0 for i in range(3) for j in range(3))
        assert positive_two_by_two(
            [[pieces[0][1][1], pieces[0][1][2]],
             [pieces[0][2][1], pieces[0][2][2]]]
        )
        assert positive_two_by_two(
            [[pieces[1][0][0], pieces[1][0][2]],
             [pieces[1][2][0], pieces[1][2][2]]]
        )
        assert positive_two_by_two(
            [[pieces[2][0][0], pieces[2][0][1]],
             [pieces[2][1][0], pieces[2][1][1]]]
        )
    print("explicit three-phase positive decompositions: exact")


def check_full_three_phase_rank() -> None:
    # Add S(e3) to the two spaces above and recover all six symmetric
    # coordinates, hence all five trace-free coordinates after projection.
    columns = [
        [Q(0), Q(1), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(1), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(0), Q(1)],
        [Q(1), Q(0), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(1), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(1), Q(0)],
        [Q(1), Q(0), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(1), Q(0), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(1), Q(0), Q(0)],
    ]
    assert rank(columns) == 6
    trace_free = [
        [row[0] - row[2], row[1] - row[2], row[3], row[4], row[5]]
        for row in columns
    ]
    assert rank(trace_free) == 5
    print("three orthogonal phases: symmetric rank 6, trace-free rank 5")


def main() -> None:
    check_two_phase_rank_loss()
    check_variable_two_phase_submersion()
    check_three_phase_decomposition()
    check_full_three_phase_rank()
    print("all corrected positive-stress checks passed")


if __name__ == "__main__":
    main()
