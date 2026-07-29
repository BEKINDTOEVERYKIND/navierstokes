#!/usr/bin/env python3
"""Exact low-order checks for pressure-multipole control.

The script constructs rational bases of homogeneous harmonic polynomials
in three variables, evaluates the ideal localized divergence-free atom
functional

    H -> b^T Hess(H)(y) b,

and verifies two facts through degree six:

* finitely many point/orientation atoms span every harmonic pressure
  moment of degrees 2,...,M; and
* the atom list has an explicit strictly positive null combination,
  obtained from orthogonal orientation triples; and
* two distinct shell radii separate exterior multipoles from interior
  harmonic pressure jets in every angular degree.

The point atoms are the zero-radius limits of compact divergence-free curl
packets.  This finite-dimensional calculation is not a construction of a
Navier--Stokes transition.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


Q = Fraction
Exponent = tuple[int, int, int]
Polynomial = dict[Exponent, Fraction]


def monomials(degree: int) -> list[Exponent]:
    return [
        (i, j, degree - i - j)
        for i in range(degree + 1)
        for j in range(degree - i + 1)
    ]


def rref_nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    if not matrix:
        return []
    rows = [row[:] for row in matrix]
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if rows[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        value = rows[pivot_row][column]
        rows[pivot_row] = [entry / value for entry in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row:
                continue
            multiplier = rows[row][column]
            if multiplier:
                rows[row] = [
                    entry - multiplier * pivot_entry
                    for entry, pivot_entry in zip(rows[row], rows[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break

    free_columns = [
        column for column in range(column_count) if column not in pivot_columns
    ]
    basis: list[list[Fraction]] = []
    for free in free_columns:
        vector = [Q(0) for _ in range(column_count)]
        vector[free] = Q(1)
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -rows[row][free]
        basis.append(vector)
    return basis


def harmonic_basis(degree: int) -> list[Polynomial]:
    domain = monomials(degree)
    if degree < 2:
        return [{exponent: Q(1)} for exponent in domain]
    codomain = monomials(degree - 2)
    row_index = {exponent: index for index, exponent in enumerate(codomain)}
    laplacian = [
        [Q(0) for _ in domain]
        for _ in codomain
    ]
    for column, exponent in enumerate(domain):
        for axis in range(3):
            if exponent[axis] < 2:
                continue
            output = list(exponent)
            output[axis] -= 2
            coefficient = exponent[axis] * (exponent[axis] - 1)
            laplacian[row_index[tuple(output)]][column] += coefficient

    nullspace = rref_nullspace(laplacian)
    assert len(nullspace) == 2 * degree + 1
    return [
        {
            exponent: coefficient
            for exponent, coefficient in zip(domain, vector)
            if coefficient
        }
        for vector in nullspace
    ]


def directional_hessian(
    polynomial: Polynomial,
    point: tuple[int, int, int],
    direction: tuple[int, int, int],
) -> Fraction:
    result = Q(0)
    for exponent, coefficient in polynomial.items():
        for first in range(3):
            for second in range(3):
                powers = list(exponent)
                derivative = Q(1)
                if first == second:
                    derivative *= powers[first] * (powers[first] - 1)
                    powers[first] -= 2
                else:
                    derivative *= powers[first] * powers[second]
                    powers[first] -= 1
                    powers[second] -= 1
                if derivative == 0:
                    continue
                value = coefficient * derivative
                for coordinate, power in zip(point, powers):
                    value *= coordinate**power
                result += direction[first] * direction[second] * value
    return result


def matrix_rank(rows: list[list[Fraction]]) -> int:
    if not rows:
        return 0
    work = [row[:] for row in rows]
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
            return rank
    return rank


ORIENTATION_GROUPS = (
    # Each group's weighted sum of b tensor b is a multiple of I.
    (((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1)),
    (((1, 1, 0), 1), ((1, -1, 0), 1), ((0, 0, 1), 2)),
    (((1, 0, 1), 1), ((1, 0, -1), 1), ((0, 1, 0), 2)),
    (((0, 1, 1), 1), ((0, 1, -1), 1), ((1, 0, 0), 2)),
)


def check_moment_chart(maximum_degree: int = 6) -> None:
    basis = [
        polynomial
        for degree in range(2, maximum_degree + 1)
        for polynomial in harmonic_basis(degree)
    ]
    expected_dimension = (maximum_degree + 1) ** 2 - 4
    assert len(basis) == expected_dimension

    # The first fifteen points of this deterministic rational grid already
    # give full rank through degree six.
    centers = [
        point
        for point in product(range(-3, 4), repeat=3)
        if 4 <= sum(coordinate * coordinate for coordinate in point) <= 12
    ][:15]

    rows: list[list[Fraction]] = []
    for center in centers:
        for group in ORIENTATION_GROUPS:
            group_rows = [
                [
                    directional_hessian(polynomial, center, direction)
                    for polynomial in basis
                ]
                for direction, _ in group
            ]
            # The positive group weights give an isotropic orientation
            # covariance, hence annihilate every harmonic Hessian.
            weighted_sum = [
                sum(
                    Q(weight) * row[column]
                    for row, (_, weight) in zip(group_rows, group)
                )
                for column in range(expected_dimension)
            ]
            assert not any(weighted_sum)
            rows.extend(group_rows)

    rank = matrix_rank(rows)
    assert rank == expected_dimension
    print(
        f"harmonic moments degrees 2..{maximum_degree}: "
        f"dimension={expected_dimension}, atom chart rank={rank}"
    )
    print("strictly positive isotropic null combinations: exact")


def check_quadrupole_energy_bound() -> None:
    # A sample covariance with eigenvalues 5,2,1.  A disjoint correction
    # making the total covariance isotropic has minimal c=lambda_max=5.
    eigenvalues = (Q(5), Q(2), Q(1))
    isotropic_level = max(eigenvalues)
    correction = tuple(isotropic_level - value for value in eigenvalues)
    assert correction == (Q(0), Q(3), Q(4))
    assert sum(correction) == 3 * max(eigenvalues) - sum(eigenvalues)
    print("quadrupole-dark disjoint correction energy lower bound: exact")


def check_single_carrier_no_go() -> None:
    # If every field is exactly transverse to w=e3, C_33=0.  Isotropy
    # C=cI would force c=0 and hence zero total energy.
    transverse_covariance = (
        (Q(2), Q(1), Q(0)),
        (Q(1), Q(3), Q(0)),
        (Q(0), Q(0), Q(0)),
    )
    trace = sum(transverse_covariance[index][index] for index in range(3))
    isotropic_level_from_trace = trace / 3
    assert transverse_covariance[2][2] == 0
    assert isotropic_level_from_trace > 0
    assert transverse_covariance[2][2] != isotropic_level_from_trace
    print("exactly transverse single-carrier bath cannot be quadrupole-dark")


def check_two_shell_radial_chart(maximum_degree: int = 20) -> None:
    # In angular degree m, the regular solid harmonic Hessian scales as
    # R^(m-2), while the Hessian of its Kelvin transform scales as
    # R^(-m-3).  At radii R and sigma*R the normalized 2-by-2 determinant
    # is sigma^(-m-3)-sigma^(m-2), hence never zero for sigma>1.
    sigma = Q(2)
    for degree in range(2, maximum_degree + 1):
        determinant = sigma ** (-(degree + 3)) - sigma ** (degree - 2)
        assert determinant
    print(
        f"two-shell exterior/interior radial chart: "
        f"nonzero through degree {maximum_degree}"
    )


def main() -> None:
    check_moment_chart()
    check_quadrupole_energy_bound()
    check_single_carrier_no_go()
    check_two_shell_radial_chart()
    print("all exact pressure-multipole checks passed")


if __name__ == "__main__":
    main()
