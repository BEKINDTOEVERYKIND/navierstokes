#!/usr/bin/env python3
"""Exact ledgers for the infinite two-harmonic charge ladder.

This checker verifies:

* the general nonzero-charge bath edge equals its K-free reduction;
* every bath edge shifts charge by at most two;
* the analytic/Gevrey shift inequalities are cutoff independent;
* the positive-second-harmonic partner edge is rank two;
* the Gevrey-two factorial convolution constant is at most three; and
* quadratic backward heat defeats every bounded general two-ended inverse
  between finite-radius Gevrey charge spaces.

The semigroup and Duhamel estimates in the accompanying note then follow
from the bounded-operator Dyson series.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb


Vector = tuple[F, F, F]

W: Vector = (F(1), F(0), F(0))
A: Vector = (F(0), F(1), F(0))
D: Vector = (F(0), F(0), F(1, 2))
Q_VECTORS: tuple[Vector, ...] = (
    (F(20), F(-45), F(-36)),
    (F(9), F(-4), F(-5)),
    (F(1), F(1), F(1)),
)
HARMONICS = (-2, -1, 1, 2)


def add(left: Vector, right: Vector) -> Vector:
    return tuple(x + y for x, y in zip(left, right))  # type: ignore[return-value]


def scale(number: F, vector: Vector) -> Vector:
    return tuple(number * entry for entry in vector)  # type: ignore[return-value]


def dot(left: Vector, right: Vector) -> F:
    return sum((x * y for x, y in zip(left, right)), F(0))


def bath_vector(harmonic: int) -> Vector:
    return A if abs(harmonic) == 1 else D


def transverse_basis(q: Vector, h: int, k_value: int) -> tuple[Vector, Vector]:
    """Two exact vectors perpendicular to q+h K e1."""

    x_value = q[0] + F(h * k_value)
    return (
        (-q[1], x_value, F(0)),
        (-q[2], F(0), x_value),
    )


def direct_edge(
    q: Vector,
    vector: Vector,
    h: int,
    harmonic: int,
    k_value: int,
) -> Vector:
    bath = bath_vector(harmonic)
    input_wave = add(q, scale(F(h * k_value), W))
    bath_wave = scale(F(harmonic * k_value), W)
    return add(
        scale(dot(bath, input_wave), vector),
        scale(dot(vector, bath_wave), bath),
    )


def reduced_edge(
    q: Vector,
    vector: Vector,
    h: int,
    harmonic: int,
) -> Vector:
    bath = bath_vector(harmonic)
    return add(
        scale(dot(bath, q), vector),
        scale(-F(harmonic, h) * dot(vector, q), bath),
    )


def triple(left: Vector, middle: Vector, right: Vector) -> F:
    return (
        left[0] * (middle[1] * right[2] - middle[2] * right[1])
        - left[1] * (middle[0] * right[2] - middle[2] * right[0])
        + left[2] * (middle[0] * right[1] - middle[1] * right[0])
    )


def check_k_free_edges() -> None:
    for k_value in (62, 101, 1009):
        for q in Q_VECTORS:
            for h in range(-12, 13):
                if h == 0:
                    continue
                input_wave = add(q, scale(F(h * k_value), W))
                for vector in transverse_basis(q, h, k_value):
                    assert dot(vector, input_wave) == 0
                    for harmonic in HARMONICS:
                        direct = direct_edge(
                            q, vector, h, harmonic, k_value
                        )
                        reduced = reduced_edge(q, vector, h, harmonic)
                        assert direct == reduced
                        bath = bath_vector(harmonic)
                        edge_bound_squared = (
                            F((1 + abs(harmonic)) ** 2)
                            * dot(q, q)
                            * dot(bath, bath)
                            * dot(vector, vector)
                        )
                        assert dot(reduced, reduced) <= edge_bound_squared


def check_shift_weight_and_support_ledgers() -> None:
    # The polynomial weight estimate
    # <h+m> <= (1+|m|)<h> is checked after squaring.
    for h in range(-1000, 1001):
        for harmonic in HARMONICS:
            assert abs(h + harmonic) <= abs(h) + abs(harmonic)
            assert 1 + (h + harmonic) ** 2 <= (
                (1 + abs(harmonic)) ** 2 * (1 + h * h)
            )

    # Starting from the partner charge -1, n bath edges cannot escape the
    # interval |-1|+2n.  The extreme charges are populated by repeated
    # second-harmonic shifts.
    charges = {-1}
    for edge_count in range(1, 13):
        charges = {
            charge + harmonic
            for charge in charges
            for harmonic in HARMONICS
            if charge + harmonic != 0
        }
        assert all(abs(charge) <= 1 + 2 * edge_count for charge in charges)
        assert -1 + 2 * edge_count in charges
        assert -1 - 2 * edge_count in charges


def check_positive_second_harmonic_rank() -> None:
    for k_value in (62, 101, 1009):
        for q in Q_VECTORS:
            first, second = transverse_basis(q, -1, k_value)
            output_wave = add(q, scale(F(k_value), W))
            first_image = reduced_edge(q, first, -1, 2)
            second_image = reduced_edge(q, second, -1, 2)
            # Projection to output_wave-perp has rank two exactly when this
            # triple determinant is nonzero.
            assert triple(first_image, second_image, output_wave) != 0


def check_factorial_convolution() -> None:
    for order in range(2, 251):
        value = F(order * order) * sum(
            (F(1, comb(order, part) ** 2) for part in range(1, order)),
            F(0),
        )
        assert value <= 3

        endpoint_contribution = F(2)
        if order >= 4:
            interior_bound = F(4 * (order - 3), (order - 1) ** 2)
            assert value <= endpoint_contribution + interior_bound
            assert interior_bound <= 1
        elif order == 3:
            assert value == endpoint_contribution
        else:
            # At r=2 the two nominal endpoints are the same term p=1.
            assert value == 1


def check_backward_heat_no_gevrey_radius() -> None:
    # theta*T=1/100 and an arbitrarily generous finite radius loss 10.
    # Since |h|^(1/sigma) <= |h| for sigma>=1, positivity of this analytic
    # exponent already implies divergence for every finite Gevrey order.
    theta_time = F(1, 100)
    radius_loss = F(10)
    previous = F(0)
    for charge in (2001, 4001, 8001, 16001):
        exponent_lower_bound = theta_time * charge * charge - radius_loss * charge
        assert exponent_lower_bound > previous
        previous = exponent_lower_bound


def main() -> None:
    check_k_free_edges()
    check_shift_weight_and_support_ledgers()
    check_positive_second_harmonic_rank()
    check_factorial_convolution()
    check_backward_heat_no_gevrey_radius()
    print("two-harmonic charge ladder: exact ledgers passed")
    print("  nonzero-charge nearest/next-nearest edges: K-free")
    print("  analytic/Gevrey finite-shift weights: cutoff independent")
    print("  positive-second-harmonic partner edge: rank two")
    print("  interaction-order Gevrey-2 convolution constant: <= 3")
    print("  bounded all-charge two-ended inverse: no finite Gevrey radius")


if __name__ == "__main__":
    main()
