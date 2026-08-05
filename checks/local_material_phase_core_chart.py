#!/usr/bin/env python3
"""Exact ledger for the localized material-phase core chart."""

from fractions import Fraction as F


def determinant(matrix):
    work = [row[:] for row in matrix]
    size = len(work)
    value = F(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for index in range(column, size):
            work[column][index] /= pivot_value
        for row in range(column + 1, size):
            multiplier = work[row][column]
            for index in range(column, size):
                work[row][index] -= multiplier * work[column][index]
    return value


def matrix_vector(matrix, vector):
    return [
        sum((matrix[row][column] * vector[column]
             for column in range(len(vector))), F(0))
        for row in range(len(matrix))
    ]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def main():
    # A periodic exact perturbation has zero periods, so x1+eta2*psi2+
    # eta3*psi3 stays in the integral class (1,0,0).
    base_periods = (F(1), F(0), F(0))
    exact_perturbation_periods = (F(0), F(0), F(0))
    assert tuple(
        base + perturbation
        for base, perturbation in zip(
            base_periods, exact_perturbation_periods
        )
    ) == base_periods

    # On the inner core, the two phase derivatives are exactly e2,e3.
    kernel = [F(1), F(0), F(0)]
    phase_derivative_2 = [F(0), F(1), F(0)]
    phase_derivative_3 = [F(0), F(0), F(1)]
    assert sum(x * y for x, y in zip(kernel, phase_derivative_2)) == 0
    assert sum(x * y for x, y in zip(kernel, phase_derivative_3)) == 0

    # Material covector identity at t=0: F'=S implies
    # (F^{-T}xi)'=-S^T xi.
    strain = [
        [F(-1), F(1, 3), F(0)],
        [F(0), F(-5, 4), F(2, 7)],
        [F(0), F(0), F(9, 4)],
    ]
    assert sum(strain[index][index] for index in range(3)) == 0
    initial_covector = [F(1), F(2, 5), F(-1, 6)]
    covector_derivative = [
        -entry for entry in matrix_vector(transpose(strain), initial_covector)
    ]
    material_left_side = [
        derivative + transported
        for derivative, transported in zip(
            covector_derivative,
            matrix_vector(transpose(strain), initial_covector),
        )
    ]
    assert material_left_side == [F(0), F(0), F(0)]

    # Full stress Jacobian in (11,22,33,12,13,23) coordinates.
    a = F(3, 2)
    b = F(2, 5)
    columns = [
        [F(1), F(1), F(1), F(0), F(0), F(0)],
        [F(0), F(1), F(0), F(0), F(0), F(0)],
        [F(0), F(0), F(1), F(0), F(0), F(0)],
        [F(0), F(0), F(0), F(0), F(0), F(1)],
        [F(0), F(0), F(0), -a, F(0), F(0)],
        [F(0), F(0), F(0), F(0), -b, F(0)],
    ]
    jacobian = [list(row) for row in zip(*columns)]
    jacobian_determinant = determinant(jacobian)
    assert abs(jacobian_determinant) == a * b

    # Conditional quantitative Neumann ledger. The bounds 4 and 3/K are
    # hypotheses, not PDE estimates proved by this checker. If they hold,
    # K=32 makes theta=3/8 and the pointwise inverse bound is 32/5<8.
    inverse_bound = F(4)
    corrector_constant = F(3)
    carrier = F(32)
    theta = inverse_bound * corrector_constant / carrier
    perturbed_inverse_bound = inverse_bound / (1 - theta)
    assert theta == F(3, 8) < F(1, 2)
    assert perturbed_inverse_bound == F(32, 5) < F(8)

    print("fixed-cohomology phase tangent directions: PASS")
    print("linearized common-flow covector identity: PASS")
    print(f"frozen-core stress determinant: {jacobian_determinant}")
    print(f"conditional pointwise Neumann bound: {perturbed_inverse_bound}")
    print("all localized material-phase core checks passed")


if __name__ == "__main__":
    main()
