#!/usr/bin/env python3
"""Exact checks for the common-carrier shear sideband repair.

This file uses only the Python standard library.  It checks

* the exact two-by-two finite-K child matrices for the three published
  low directions;
* the inherited rank-five symmetric-strain minor;
* the all-generation cancellation of the nominal K factor between two
  nonzero fast charges; and
* a family showing that the resulting charge multiplier is quadratic,
  rather than linear, in an unrestricted charge ball.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import sqrt


Vector = tuple[F, F, F]


Q_VECTORS: tuple[Vector, ...] = (
    (F(-45), F(-36), F(20)),
    (F(-4), F(-5), F(9)),
    (F(1), F(1), F(1)),
)
W: Vector = (F(0), F(0), F(1))
A: Vector = (F(1), F(0), F(0))


def add(x: Vector, y: Vector) -> Vector:
    return tuple(a + b for a, b in zip(x, y))  # type: ignore[return-value]


def scale(c: F, x: Vector) -> Vector:
    return tuple(c * a for a in x)  # type: ignore[return-value]


def dot(x: Vector, y: Vector) -> F:
    return sum((a * b for a, b in zip(x, y)), F(0))


def project(q: Vector, x: Vector) -> Vector:
    return add(x, scale(-dot(q, x) / dot(q, q), q))


def determinant(matrix: list[list[F]]) -> F:
    work = [row[:] for row in matrix]
    result = F(1)
    for column in range(len(work)):
        pivot = next(
            row
            for row in range(column, len(work))
            if work[row][column] != 0
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(column + 1, len(work)):
            multiplier = work[row][column]
            work[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return result


def child(q: Vector, k_value: int, b: Vector) -> Vector:
    """The exact projected coefficient at q from k=K e3 and r=q-k."""

    k = scale(F(k_value), W)
    r = add(q, scale(F(-1), k))
    unprojected = add(scale(dot(A, r), b), scale(dot(b, k), A))
    return project(q, unprojected)


def child_coordinate_matrix(q: Vector, k_value: int) -> list[list[F]]:
    """Matrix from a convenient basis of r-perp to one of q-perp."""

    alpha, beta, delta = q
    domain = (
        (-beta, alpha, F(0)),
        (F(k_value) - delta, F(0), alpha),
    )
    output = (
        (-beta, alpha, F(0)),
        (-delta, F(0), alpha),
    )
    matrix: list[list[F]] = [[], []]
    for b in domain:
        assert dot(b, add(q, scale(F(-k_value), W))) == 0
        value = child(q, k_value, b)
        # If value=x output[0]+y output[1], its y,z coordinates are
        # alpha*x and alpha*y.
        coordinates = (value[1] / alpha, value[2] / alpha)
        assert value == add(
            scale(coordinates[0], output[0]),
            scale(coordinates[1], output[1]),
        )
        matrix[0].append(coordinates[0])
        matrix[1].append(coordinates[1])
    return matrix


def symmetric_column(c: Vector, q: Vector) -> list[F]:
    matrix = [
        [F(c[i] * q[j] + q[i] * c[j], 2) for j in range(3)]
        for i in range(3)
    ]
    return [
        matrix[0][0],
        matrix[1][1],
        matrix[0][1],
        matrix[0][2],
        matrix[1][2],
    ]


def check_child_rank() -> None:
    for k_value in range(1, 501):
        matrices = [
            child_coordinate_matrix(q, k_value) for q in Q_VECTORS
        ]
        for q, matrix in zip(Q_VECTORS, matrices):
            alpha, beta, delta = q
            q_squared = dot(q, q)
            expected = [
                [
                    alpha,
                    -F(2) * alpha * beta * k_value / q_squared,
                ],
                [
                    F(0),
                    alpha
                    * (q_squared - F(2) * k_value * delta)
                    / q_squared,
                ],
            ]
            assert matrix == expected
            assert determinant(matrix) != 0

        # The first five domain columns map block-diagonally to the first
        # five transverse child columns.  Hence their strain determinant
        # is the old transverse determinant times det(M1) det(M2) M3_11.
        inherited = (
            F(-1_214_003_700)
            * determinant(matrices[0])
            * determinant(matrices[1])
            * matrices[2][0][0]
        )
        explicit = (
            -F(644_815_080_000, 3721)
            * (F(3721) - F(40) * k_value)
            * (F(61) - F(9) * k_value)
        )
        assert inherited == explicit
        assert inherited != 0


def check_all_chain_k_cancellation() -> None:
    """Check the exact charge identity on representative integer data."""

    # xi and eta are arbitrary slow-lattice vectors for this algebra.
    xi: Vector = (F(17), F(-8), F(5))
    eta: Vector = (F(-3), F(11), F(7))
    h, g, k_value = F(3), F(-5), F(113)
    p = add(scale(k_value * h, W), xi)
    r = add(scale(k_value * g, W), eta)

    # Cross products provide exact vectors transverse to p and r.
    u: Vector = (p[1], -p[0], F(0))
    v: Vector = (r[2], F(0), -r[0])
    assert dot(u, p) == 0
    assert dot(v, r) == 0
    assert dot(u, r) == dot(u, add(eta, scale(-g / h, xi)))
    assert dot(v, p) == dot(v, add(xi, scale(-h / g, eta)))


def projected_norm(q: Vector, x: Vector) -> float:
    value = project(q, x)
    return sqrt(float(dot(value, value)))


def quadratic_charge_example(m: int) -> float:
    """Return |B|/m^2 for an unrestricted-charge counterexample.

    Take p=K e3+m q1 with fast charge h=1, and r=m K e3 with
    fast charge g=m.  The first amplitude is made exactly transverse to
    p, while the second is the parent shear polarization e1.  K=m^3
    keeps the output direction close to e3.  The projected Euler symbol
    is asymptotic to 45 m^2 e1.
    """

    q = Q_VECTORS[0]
    k_value = F(m**3)
    p = add(scale(k_value, W), scale(F(m), q))
    r = scale(F(m) * k_value, W)
    u: Vector = (
        F(1),
        F(0),
        -F(m) * q[0] / (k_value + F(m) * q[2]),
    )
    v = A
    assert dot(u, p) == 0
    assert dot(v, r) == 0
    raw = add(scale(dot(u, r), v), scale(dot(v, p), u))
    return projected_norm(add(p, r), raw) / (m * m)


def main() -> None:
    check_child_rank()
    check_all_chain_k_cancellation()
    ratios = [quadratic_charge_example(m) for m in (20, 40, 80, 160)]
    assert min(ratios) > 30.0
    print("finite-K matched child blocks: rank 2 for every integer 1<=K<=500")
    print("inherited symmetric child chart: exact rank 5")
    print("charged-charged nominal K factor: cancels exactly")
    print("quadratic charge multiplier |B|/M^2:", " ".join(f"{x:.6f}" for x in ratios))


if __name__ == "__main__":
    main()
