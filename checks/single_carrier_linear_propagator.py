#!/usr/bin/env python3
"""Checks for the fixed-star single-carrier propagator note.

Only the Python standard library is used.  The script checks:

* the exact finite-K two-column child determinants in two complementary
  bases;
* the exact K-free charged interaction identities;
* the distinction between an O(M) fixed-star derivative and the sharp
  O(M^2) unrestricted bilinear multiplier; and
* the length-two phase-translation chain for a one-direction shear.

These are finite-dimensional algebra checks, not a PDE existence proof.
"""

from __future__ import annotations

from fractions import Fraction as F


Vector = tuple[F, F, F]


def add(x: Vector, y: Vector) -> Vector:
    return tuple(a + b for a, b in zip(x, y))  # type: ignore[return-value]


def scale(c: F, x: Vector) -> Vector:
    return tuple(c * a for a in x)  # type: ignore[return-value]


def dot(x: Vector, y: Vector) -> F:
    return sum((a * b for a, b in zip(x, y)), F(0))


def cross(x: Vector, y: Vector) -> Vector:
    return (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )


def child_triple_determinant(q: Vector, k: int) -> F:
    """Compute det[V1,V2,q] for the normalized domain basis in the note."""

    alpha, beta, delta = q
    denominator = F(k) - delta
    b1 = (F(1), F(0), alpha / denominator)
    b2 = (F(0), F(1), beta / denominator)
    e1: Vector = (F(1), F(0), F(0))

    def raw(b: Vector) -> Vector:
        return add(scale(alpha, b), scale(dot(b, q), e1))

    return dot(cross(raw(b1), raw(b2)), q)


def project(q: Vector, x: Vector) -> Vector:
    return add(x, scale(-dot(q, x) / dot(q, q), q))


def child_coordinate_determinant(q: Vector, k: int) -> F:
    """Compute the child determinant in bases valid for every K."""

    alpha, beta, delta = q
    q_squared = dot(q, q)
    r = (alpha, beta, delta - F(k))
    d1: Vector = (-beta, alpha, F(0))
    d2: Vector = (F(k) - delta, F(0), alpha)
    c1: Vector = (-beta, alpha, F(0))
    c2: Vector = (-delta, F(0), alpha)
    e1: Vector = (F(1), F(0), F(0))

    def child(b: Vector) -> Vector:
        assert dot(b, r) == 0
        raw = add(scale(alpha, b), scale(dot(b, q), e1))
        return project(q, raw)

    columns: list[tuple[F, F]] = []
    for b in (d1, d2):
        value = child(b)
        # value=x*c1+y*c2, whose y,z entries are alpha*x,alpha*y.
        coordinates = (value[1] / alpha, value[2] / alpha)
        assert value == add(
            scale(coordinates[0], c1),
            scale(coordinates[1], c2),
        )
        columns.append(coordinates)

    determinant = (
        columns[0][0] * columns[1][1]
        - columns[1][0] * columns[0][1]
    )
    expected = alpha * alpha * (q_squared - F(2 * k) * delta) / q_squared
    assert determinant == expected
    return determinant


def check_child_determinant() -> None:
    q_vectors: tuple[Vector, ...] = (
        (F(-45), F(-36), F(20)),
        (F(-4), F(-5), F(9)),
        (F(1), F(1), F(1)),
    )
    for q in q_vectors:
        alpha, _, delta = q
        q_squared = dot(q, q)
        for k in (2, 7, 9, 20, 31, 127):
            coordinate_value = child_coordinate_determinant(q, k)
            assert coordinate_value != 0
            if F(k) != delta:
                expected = (
                    alpha * alpha * (F(2 * k) * delta - q_squared)
                    / (F(k) - delta)
                )
                assert child_triple_determinant(q, k) == expected


def check_k_free_identity() -> None:
    w: Vector = (F(0), F(0), F(1))
    xi: Vector = (F(17), F(-8), F(5))
    eta: Vector = (F(-3), F(11), F(7))
    h, g, k = F(3), F(-5), F(113)
    p = add(scale(k * h, w), xi)
    r = add(scale(k * g, w), eta)

    # Cross products furnish exact transverse amplitudes.
    a = cross(p, (F(1), F(2), F(3)))
    b = cross(r, (F(2), F(-1), F(4)))
    assert dot(a, p) == 0
    assert dot(b, r) == 0
    assert dot(a, r) == dot(a, add(eta, scale(-g / h, xi)))
    assert dot(b, p) == dot(b, add(xi, scale(-h / g, eta)))


def multiplier_bound(
    h: int,
    xi_size: int,
    g: int,
    eta_size: int,
) -> F:
    """Scalar upper multiplier in the exact K-free estimate."""

    return (
        F(eta_size)
        + abs(F(g, h)) * xi_size
        + F(xi_size)
        + abs(F(h, g)) * eta_size
    )


def check_linear_versus_quadratic_degree() -> None:
    # A fixed background label (g,eta)=(1,1) gives O(M).
    linear_ratios: list[F] = []
    # Allowing both inputs to grow gives O(M^2).
    quadratic_ratios: list[F] = []
    for m in (20, 40, 80, 160):
        fixed = multiplier_bound(1, m, 1, 1)
        unrestricted = multiplier_bound(1, m, m, 0)
        linear_ratios.append(fixed / m)
        quadratic_ratios.append(unrestricted / (m * m))

    assert max(linear_ratios) < F(4)
    assert min(linear_ratios) > F(2)
    assert min(quadratic_ratios) >= F(1)
    assert max(quadratic_ratios) < F(2)


def check_nilpotent_phase_chain() -> None:
    # Work coefficientwise with W_p=a, p=K e3.  A constant e3 mode
    # produces K a at p.  A shear-polarized p mode has zero interaction
    # with both +p and -p background modes, so applying L again is zero.
    k = F(97)
    p: Vector = (F(0), F(0), k)
    minus_p = scale(F(-1), p)
    a: Vector = (F(1), F(0), F(0))
    e3: Vector = (F(0), F(0), F(1))

    first = scale(dot(e3, p), a)
    assert first == scale(k, a)

    # Possible second targets are 0 and 2p.  Every raw coefficient
    # contains a dot product of a with p or -p.
    assert dot(a, p) == 0
    assert dot(a, minus_p) == 0


def main() -> None:
    check_child_determinant()
    check_k_free_identity()
    check_linear_versus_quadratic_degree()
    check_nilpotent_phase_chain()
    print("single-carrier child determinant: exact")
    print("charged interaction K-cancellation: exact")
    print("fixed-star derivative O(M), unrestricted bilinear map O(M^2)")
    print("unmodulated shear phase chain: nilpotent of length two")


if __name__ == "__main__":
    main()
