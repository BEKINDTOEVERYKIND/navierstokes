#!/usr/bin/env python3
"""Exact arithmetic for C118--C120.

Checks:

* the A2 root and equal-leaf-shell identities;
* the six-cycle incidence graph and its first-step off-ladder support
  outputs;
* the unweighted hexagon eigenspace decomposition;
* the weighted hexagon characteristic polynomial and a rational enclosure
  of its simple top growth eigenvalue; and
* the pump/leaf energy and one-/three-leaf heteroclinic identities.

Only the finite-dimensional statements are checked.  No invariant Fourier
subspace of the full Euler or Navier--Stokes equation is asserted.
"""

from __future__ import annotations

from fractions import Fraction as F


Vector = tuple[int, int, int]
Matrix3 = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]

ROOTS: tuple[Vector, ...] = (
    (1, -1, 0),
    (0, 1, -1),
    (-1, 0, 1),
)
B0: Matrix3 = (
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
)
BSTAR: Matrix3 = (
    (0, 1, 1),
    (1, 0, 2),
    (1, 1, 0),
)


def add(left: Vector, right: Vector) -> Vector:
    return tuple(x + y for x, y in zip(left, right))  # type: ignore[return-value]


def sub(left: Vector, right: Vector) -> Vector:
    return tuple(x - y for x, y in zip(left, right))  # type: ignore[return-value]


def scale(number: int, vector: Vector) -> Vector:
    return tuple(number * entry for entry in vector)  # type: ignore[return-value]


def dot(left: Vector, right: Vector) -> int:
    return sum(x * y for x, y in zip(left, right))


def mat_vec(matrix: Matrix3, vector: tuple[F, F, F]) -> tuple[F, F, F]:
    return tuple(
        sum((F(entry) * value for entry, value in zip(row, vector)), F(0))
        for row in matrix
    )  # type: ignore[return-value]


def transpose(matrix: Matrix3) -> Matrix3:
    return tuple(zip(*matrix))  # type: ignore[return-value]


def mat_mul(left: Matrix3, right: Matrix3) -> Matrix3:
    right_t = transpose(right)
    return tuple(
        tuple(sum(x * y for x, y in zip(row, column)) for column in right_t)
        for row in left
    )  # type: ignore[return-value]


def determinant(matrix: Matrix3) -> int:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def p(value: F) -> F:
    return value**3 - 9 * value**2 + 18 * value - 9


def polynomial_product(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply ascending-order integer coefficient tuples."""

    output = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            output[left_degree + right_degree] += left_value * right_value
    return tuple(output)


def check_a2_geometry() -> None:
    assert add(add(ROOTS[0], ROOTS[1]), ROOTS[2]) == (0, 0, 0)
    for index, root in enumerate(ROOTS):
        assert dot(root, root) == 2
        for other in ROOTS[index + 1 :]:
            assert dot(root, other) == -1

    root_hexagon = {root for r in ROOTS for root in (r, scale(-1, r))}
    for k_value in (1, 7, 31):
        for m_value in (1, 4, 19):
            centre = (m_value, m_value, m_value)
            positive = [add(centre, scale(k_value, root)) for root in ROOTS]
            negative = [sub(centre, scale(k_value, root)) for root in ROOTS]
            expected_shell = 3 * m_value * m_value + 2 * k_value * k_value
            assert all(dot(mode, mode) == expected_shell for mode in positive + negative)

            incidence = []
            for left in positive:
                row = []
                for right in negative:
                    difference = sub(left, right)
                    row.append(int(any(difference == scale(k_value, r) for r in root_hexagon)))
                incidence.append(tuple(row))
            assert tuple(incidence) == B0

            # Pump-support addition also permits the centre and an outer
            # second sideband after one opposite/same-root step.
            for root, leaf in zip(ROOTS, positive):
                assert sub(leaf, scale(k_value, root)) == centre
                assert add(leaf, scale(k_value, root)) == add(
                    centre, scale(2 * k_value, root)
                )
                assert centre not in positive + negative
                assert add(centre, scale(2 * k_value, root)) not in positive + negative


def check_unweighted_eigenspaces() -> None:
    one = (F(1), F(1), F(1))
    sum_zero = ((F(1), F(-1), F(0)), (F(1), F(0), F(-1)))
    assert mat_vec(B0, one) == tuple(2 * x for x in one)
    for vector in sum_zero:
        assert mat_vec(B0, vector) == tuple(-x for x in vector)
    assert determinant(
        (
            tuple(int(x) for x in one),
            tuple(int(x) for x in sum_zero[0]),
            tuple(int(x) for x in sum_zero[1]),
        )
    ) != 0

    # The six eigenvectors of [[0,B0],[B0,0]] are obtained by pairing
    # (v,v) and (v,-v).  Their eigenvalues give the exact polynomial.
    eigenvalues = (F(2), F(-2), F(1), F(1), F(-1), F(-1))
    assert sum(eigenvalues, F(0)) == 0
    assert sum(value**2 for value in eigenvalues) == 12
    # Ascending coefficients:
    # (x^2-4)(x^2-1)^2 = x^6-6x^4+9x^2-4.
    polynomial = polynomial_product(
        (-4, 0, 1), polynomial_product((-1, 0, 1), (-1, 0, 1))
    )
    assert polynomial == (-4, 0, 9, 0, -6, 0, 1)


def check_weighted_polynomial_and_enclosure() -> None:
    astar = mat_mul(BSTAR, transpose(BSTAR))
    assert astar == ((2, 2, 1), (2, 5, 1), (1, 1, 2))
    assert determinant(BSTAR) == 3
    assert determinant(astar) == 9

    # For a 3x3 matrix, the characteristic coefficients are trace,
    # sum of principal 2x2 minors, and determinant.
    trace = sum(astar[index][index] for index in range(3))
    principal_two = (
        astar[0][0] * astar[1][1] - astar[0][1] * astar[1][0]
        + astar[0][0] * astar[2][2] - astar[0][2] * astar[2][0]
        + astar[1][1] * astar[2][2] - astar[1][2] * astar[2][1]
    )
    assert (trace, principal_two, determinant(astar)) == (9, 18, 9)
    # p(lambda^2) in ascending coefficients.
    p_coefficients = (-9, 18, -9, 1)
    ladder_coefficients = [0] * 7
    for degree, coefficient in enumerate(p_coefficients):
        ladder_coefficients[2 * degree] = coefficient
    assert tuple(ladder_coefficients) == (-9, 0, 18, 0, -9, 0, 1)

    lower = F(633, 250)
    upper = F(2533, 1000)
    assert p(lower * lower) < 0 < p(upper * upper)
    assert p(F(5)) < 0
    # p'(5)>0 and p''(mu)>0 for mu>=5, so p is strictly increasing there.
    assert 3 * F(5) ** 2 - 18 * F(5) + 18 > 0
    assert 6 * F(5) - 18 > 0

    # Polynomial eigenvector identity.  The last residual is -p(mu).
    for mu in (F(6), lower * lower, upper * upper):
        xi = (mu - 3, mu, mu * mu - 7 * mu + 6)
        image = mat_vec(astar, xi)
        residual = tuple(image[i] - mu * xi[i] for i in range(3))
        assert residual[:2] == (0, 0)
        assert residual[2] == -p(mu)


def check_depletion_identities() -> None:
    # General reduced identities at rational test values.  S2 denotes
    # sech(z)^2 and T=tanh(z), so S2=1-T^2 exactly.
    for sigma in (F(1), F(2), F(7, 3)):
        for leaf_count in (1, 3, 8):
            for amplitude in (F(2), F(11, 5)):
                for tanh_value in (F(-3, 5), F(0), F(4, 5)):
                    sech_squared = 1 - tanh_value * tanh_value
                    pump = -amplitude * tanh_value
                    leaf_square = amplitude * amplitude * sech_squared / (2 * leaf_count)

                    pump_derivative = -sigma * amplitude * amplitude * sech_squared
                    leaf_log_derivative = -sigma * amplitude * tanh_value
                    assert pump_derivative == -2 * sigma * leaf_count * leaf_square
                    assert leaf_log_derivative == sigma * pump
                    assert pump * pump + 2 * leaf_count * leaf_square == amplitude * amplitude

    # The B0 coordinate normalization x=y=a*1 has sigma=2 and three
    # equal coordinate leaves: p'=-12a^2, a'=2pa, E=p^2+6a^2.
    p_value, a_value = F(7, 5), F(2, 9)
    p_dot = -12 * a_value * a_value
    a_dot = 2 * p_value * a_value
    energy_dot = 2 * p_value * p_dot + 12 * a_value * a_dot
    assert energy_dot == 0


def main() -> None:
    check_a2_geometry()
    check_unweighted_eigenspaces()
    check_weighted_polynomial_and_enclosure()
    check_depletion_identities()
    print("C118-C120 A2 hexagon/depletion checker: PASS")
    print("  weighted ladder polynomial: lambda^6-9 lambda^4+18 lambda^2-9")
    print("  top growth enclosure: 633/250 < sigma_* < 2533/1000")
    print("  boundary: support admits centre/second-sideband outputs in one step")


if __name__ == "__main__":
    main()
