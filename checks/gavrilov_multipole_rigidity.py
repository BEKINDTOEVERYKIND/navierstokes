#!/usr/bin/env python3
"""Exact algebra checks for the Gavrilov multipole/dynamic audit.

The checks cover finite-dimensional consequences of compact steady Euler:

* vanishing harmonic quadratic moments force isotropic covariance;
* rotating/scaling an isotropic covariance gives only the isotropic ray,
  with zero trace-free orbit rank; and
* a vector that is both longitudinal (a pressure gradient) and transverse
  (a divergence-free time tangent) at nonzero Fourier frequency is zero.

They do not construct Gavrilov's flow or a Navier--Stokes transition.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations, product


Q = Fraction


def matrix_rank(rows: list[list[Fraction]]) -> int:
    work = [row[:] for row in rows]
    if not work:
        return 0
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            multiplier = work[row][column]
            if multiplier:
                work[row] = [
                    entry - multiplier * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[rank])
                ]
        rank += 1
        if rank == column_count:
            break
    return rank


def matmul(left, right):
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(3))
            for j in range(3)
        ]
        for i in range(3)
    ]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def symmetric_coordinates(matrix):
    return [
        matrix[0][0],
        matrix[1][1],
        matrix[2][2],
        matrix[0][1],
        matrix[0][2],
        matrix[1][2],
    ]


def signed_permutation_rotations():
    rotations = []
    for permutation in permutations(range(3)):
        inversion_count = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        permutation_sign = -1 if inversion_count % 2 else 1
        for signs in product((-1, 1), repeat=3):
            if permutation_sign * signs[0] * signs[1] * signs[2] != 1:
                continue
            matrix = [[Q(0) for _ in range(3)] for _ in range(3)]
            for row, column in enumerate(permutation):
                matrix[row][column] = Q(signs[row])
            rotations.append(matrix)
    assert len(rotations) == 24
    return rotations


def check_quadratic_virial_chart() -> None:
    # Symmetric coordinates are C11,C22,C33,C12,C13,C23.
    # These are contractions with Hessians of
    # x^2-y^2, x^2-z^2, xy, xz, yz (irrelevant scalar factors omitted).
    harmonic_quadratic_rows = [
        [Q(1), Q(-1), Q(0), Q(0), Q(0), Q(0)],
        [Q(1), Q(0), Q(-1), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(1), Q(0), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(1), Q(0)],
        [Q(0), Q(0), Q(0), Q(0), Q(0), Q(1)],
    ]
    assert matrix_rank(harmonic_quadratic_rows) == 5
    isotropic = [Q(1), Q(1), Q(1), Q(0), Q(0), Q(0)]
    assert all(
        sum(entry * value for entry, value in zip(row, isotropic)) == 0
        for row in harmonic_quadratic_rows
    )
    print("quadratic pressure-dark covariance space: isotropic ray exactly")


def check_rotation_scaling_orbit() -> None:
    identity = [
        [Q(1), Q(0), Q(0)],
        [Q(0), Q(1), Q(0)],
        [Q(0), Q(0), Q(1)],
    ]
    orbit = []
    amplitude = Q(3, 2)
    spatial_scale = Q(5, 3)
    scalar = amplitude * amplitude * spatial_scale**3
    for rotation in signed_permutation_rotations():
        covariance = matmul(matmul(rotation, identity), transpose(rotation))
        covariance = [
            [scalar * entry for entry in row]
            for row in covariance
        ]
        assert covariance == [
            [scalar if i == j else Q(0) for j in range(3)]
            for i in range(3)
        ]
        orbit.append(symmetric_coordinates(covariance))
    assert matrix_rank(orbit) == 1

    trace_free_orbit = [
        [
            row[0] - (row[0] + row[1] + row[2]) / 3,
            row[1] - (row[0] + row[1] + row[2]) / 3,
            row[3],
            row[4],
            row[5],
        ]
        for row in orbit
    ]
    assert matrix_rank(trace_free_orbit) == 0
    print("rotated/scaled covariance orbit: rank 1, trace-free rank 0")


def check_gradient_transverse_intersection() -> None:
    # A compact/decaying time tangent to divergence-free fields is
    # transverse to k.  A pressure gradient is longitudinal.  Their
    # intersection is zero for every nonzero k.
    for wavevector in ((1, 0, 0), (1, 2, 3), (-2, 5, 7)):
        norm_squared = sum(Q(value * value) for value in wavevector)
        assert norm_squared
        # If r=lambda*k and k dot r=0, then lambda*|k|^2=0.
        coefficient_of_lambda = norm_squared
        assert coefficient_of_lambda != 0
    print("divergence-free/gradient time-tangent intersection: zero")


def main() -> None:
    check_quadratic_virial_chart()
    check_rotation_scaling_orbit()
    check_gradient_transverse_intersection()
    print("all exact Gavrilov multipole-rigidity checks passed")


if __name__ == "__main__":
    main()
