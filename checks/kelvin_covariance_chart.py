#!/usr/bin/env python3
"""Exact ledger for finite-time Kelvin transport of the rank-two chart."""

from fractions import Fraction as F


class Dual:
    """Fraction-valued first jet in six control variables."""

    def __init__(self, value, derivative=None):
        self.value = F(value)
        self.derivative = (
            [F(0)] * 6 if derivative is None else list(derivative)
        )

    @classmethod
    def variable(cls, value, index):
        derivative = [F(0)] * 6
        derivative[index] = F(1)
        return cls(value, derivative)

    @staticmethod
    def lift(other):
        return other if isinstance(other, Dual) else Dual(other)

    def __add__(self, other):
        other = self.lift(other)
        return Dual(
            self.value + other.value,
            [left + right for left, right in zip(
                self.derivative, other.derivative
            )],
        )

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-entry for entry in self.derivative])

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        return Dual(
            self.value * other.value,
            [
                left * other.value + self.value * right
                for left, right in zip(self.derivative, other.derivative)
            ],
        )

    __rmul__ = __mul__

    def reciprocal(self):
        return Dual(
            1 / self.value,
            [-entry / (self.value * self.value)
             for entry in self.derivative],
        )

    def __truediv__(self, other):
        return self * self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other) / self


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


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), Dual(0))


def cross(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def matrix_vector(matrix, vector):
    return [
        sum((entry * component for entry, component in zip(row, vector)),
            Dual(0))
        for row in matrix
    ]


def kelvin_amplitude(deformation, inverse_transpose, covector, amplitude):
    final_covector = matrix_vector(inverse_transpose, covector)
    initial_vorticity = cross(covector, amplitude)
    final_vorticity = matrix_vector(deformation, initial_vorticity)
    norm_squared = dot(final_covector, final_covector)
    return [
        component / norm_squared
        for component in cross(final_vorticity, final_covector)
    ]


def main():
    # A genuinely non-diagonal volume-preserving deformation.
    deformation = [
        [F(1), F(1), F(0)],
        [F(0), F(1), F(1)],
        [F(0), F(0), F(1)],
    ]
    inverse_transpose = [
        [F(1), F(0), F(0)],
        [F(-1), F(1), F(0)],
        [F(1), F(-1), F(1)],
    ]

    rho = Dual.variable(F(0), 0)
    block_22 = Dual.variable(F(3, 2), 1)
    block_33 = Dual.variable(F(2, 5), 2)
    block_23 = Dual.variable(F(0), 3)
    eta_2 = Dual.variable(F(0), 4)
    eta_3 = Dual.variable(F(0), 5)

    covector = [Dual(1), eta_2, eta_3]
    transverse_2 = [-eta_2, Dual(1), Dual(0)]
    transverse_3 = [-eta_3, Dual(0), Dual(1)]
    assert dot(covector, transverse_2).value == 0
    assert dot(covector, transverse_3).value == 0

    amplitude_2 = kelvin_amplitude(
        deformation, inverse_transpose, covector, transverse_2
    )
    amplitude_3 = kelvin_amplitude(
        deformation, inverse_transpose, covector, transverse_3
    )
    final_covector = matrix_vector(inverse_transpose, covector)
    assert dot(final_covector, amplitude_2).value == 0
    assert dot(final_covector, amplitude_3).value == 0

    # The velocity law is neither the vector pushforward F*a nor the
    # covector law F^{-T}*a.  The vorticity k x a, rather than a itself,
    # is what F pushes forward.
    assert [entry.value for entry in amplitude_2] == [
        F(2, 3), F(1, 3), F(-1, 3)
    ]
    naive_vector_push = matrix_vector(deformation, transverse_2)
    naive_covector_push = matrix_vector(inverse_transpose, transverse_2)
    assert [entry.value for entry in amplitude_2] != [
        entry.value for entry in naive_vector_push
    ]
    assert [entry.value for entry in amplitude_2] != [
        entry.value for entry in naive_covector_push
    ]
    pushed_vorticity = matrix_vector(
        deformation, cross(covector, transverse_2)
    )
    recovered_vorticity = cross(final_covector, amplitude_2)
    assert [entry.value for entry in recovered_vorticity] == [
        entry.value for entry in pushed_vorticity
    ]
    unit_covariance = [
        [
            amplitude_2[row].value * amplitude_2[column].value
            + amplitude_3[row].value * amplitude_3[column].value
            for column in range(3)
        ]
        for row in range(3)
    ]
    assert unit_covariance == [
        [F(5, 9), F(1, 9), F(-4, 9)],
        [F(1, 9), F(2, 9), F(1, 9)],
        [F(-4, 9), F(1, 9), F(5, 9)],
    ]
    assert [
        sum(unit_covariance[row][column] * base_component
            for column, base_component in enumerate((F(1), F(-1), F(1))))
        for row in range(3)
    ] == [F(0), F(0), F(0)]

    covariance = [[Dual(0) for _ in range(3)] for _ in range(3)]
    for row in range(3):
        for column in range(3):
            covariance[row][column] = (
                block_22 * amplitude_2[row] * amplitude_2[column]
                + block_33 * amplitude_3[row] * amplitude_3[column]
                + block_23 * (
                    amplitude_2[row] * amplitude_3[column]
                    + amplitude_3[row] * amplitude_2[column]
                )
            )

    stress = [[-entry for entry in row] for row in covariance]
    for index in range(3):
        stress[index][index] += rho
    outputs = [
        stress[0][0], stress[1][1], stress[2][2],
        stress[0][1], stress[0][2], stress[1][2],
    ]
    jacobian = [entry.derivative for entry in outputs]
    jacobian_determinant = determinant(jacobian)

    base_final_covector = [F(1), F(-1), F(1)]
    norm_squared = sum(entry * entry for entry in base_final_covector)
    expected_magnitude = F(3, 2) * F(2, 5) / norm_squared**4
    assert abs(jacobian_determinant) == expected_magnitude
    assert expected_magnitude == F(1, 135)

    # For S=diag(-1,-5/4,9/4), n=e1 gives log s=t and hence
    # log |det J| = log(ab)-8t.  This records the exact conditioning
    # exponent without introducing floating-point exponentials.
    strain_rates = (F(-1), F(-5, 4), F(9, 4))
    log_covector_rate = -strain_rates[0]
    log_jacobian_rate = -8 * log_covector_rate
    assert log_covector_rate == 1
    assert log_jacobian_rate == -8

    # Unequal positive damping of the two labelled harmonics still gives
    # an invertible covariance-block map.  At B=diag(a,b), its derivative
    # determinant is delta1*delta2*c, with c strictly between the dampings.
    sqrt_a = F(2)
    sqrt_b = F(3)
    delta_1 = F(1, 2)
    delta_2 = F(1, 3)
    off_diagonal_factor = (
        delta_1 * sqrt_a + delta_2 * sqrt_b
    ) / (sqrt_a + sqrt_b)
    viscous_block_determinant = (
        delta_1 * delta_2 * off_diagonal_factor
    )
    assert off_diagonal_factor == F(2, 5)
    assert viscous_block_determinant == F(1, 15) > 0

    print("non-diagonal Kelvin incompressibility: PASS")
    print(f"finite-time stress Jacobian determinant: {jacobian_determinant}")
    print(f"predicted determinant magnitude: {expected_magnitude}")
    print(f"diagonal affine log-Jacobian rate: {log_jacobian_rate}")
    print(f"unequal-damping block determinant: {viscous_block_determinant}")
    print("all Kelvin covariance-chart checks passed")


if __name__ == "__main__":
    main()
