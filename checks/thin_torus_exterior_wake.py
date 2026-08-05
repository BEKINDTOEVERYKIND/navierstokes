#!/usr/bin/env python3
"""Exact exponent and tensor checks for the thin-torus exterior wake.

This checker verifies the monomial scaling, the orthogonal-axis
quadrupole cancellation, and the packed multipole gains in
research/2026-08-02-thin-torus-exterior-wake-uniformity.md.  It does not
prove the Newton-kernel, Bogovskii-extension, or higher-jet estimates.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction

# Exponent order: (nu, v, ell, K, eps, T, d).
Power = tuple[Fraction, ...]


def power(*entries: int | Fraction) -> Power:
    return tuple(Q(entry) for entry in entries)


def add(*powers: Power) -> Power:
    return tuple(sum(column, Q(0)) for column in zip(*powers))


def scale(multiplier: int | Fraction, powers: Power) -> Power:
    return tuple(Q(multiplier) * entry for entry in powers)


def check_thin_source_and_pressure_scaling() -> None:
    viscosity = power(1, 0, 0, 0, 0, 0, 0)
    velocity = power(0, 1, 0, 0, 0, 0, 0)
    laplacian_velocity = power(0, 1, -2, 0, -2, 0, 0)
    tube_volume = power(0, 0, 3, 0, 2, 0, 0)

    source_l1 = add(viscosity, velocity, laplacian_velocity, tube_volume)
    assert source_l1 == power(1, 2, 1, 0, 0, 0, 0)

    # One more derivative of S and the square root of tube volume give
    # the L2 scale of the order-one pressure multiplier.
    source_pointwise = add(viscosity, velocity, laplacian_velocity)
    one_derivative = power(0, 0, -1, 0, -1, 0, 0)
    sqrt_volume = scale(Q(1, 2), tube_volume)
    pressure_gradient_l2 = add(source_pointwise, one_derivative, sqrt_volume)
    assert pressure_gradient_l2 == power(1, 2, Q(-3, 2), 0, -2, 0, 0)

    far_kernel = power(0, 0, 0, 0, 0, 0, -4)
    far_velocity = add(source_l1, far_kernel)
    assert far_velocity == power(1, 2, 1, 0, 0, 0, -4)

    turnover_squared = power(0, -2, 2, 0, 0, 2, 0)
    normalized_far = add(far_velocity, turnover_squared)
    assert normalized_far == power(1, 0, 3, 0, 0, 2, -4)

    exterior_l2 = power(1, 2, Q(-3, 2), 0, 0, 0, 0)
    normalized_exterior_l2 = add(exterior_l2, turnover_squared)
    assert normalized_exterior_l2 == power(1, 0, Q(1, 2), 0, 0, 2, 0)
    print("thin source and exterior wake are aspect-uniform: exact")


def check_orthogonal_axis_cancellation() -> None:
    # An axisymmetric covariance is gp*I + D*n*n.  Summing over an
    # orthonormal frame produces the scalar matrix (3*gp+D)*I.
    gp = Q(7, 5)
    anisotropy = Q(11, 6)
    identity = [[Q(int(i == j)) for j in range(3)] for i in range(3)]
    axes = [
        [Q(1), Q(0), Q(0)],
        [Q(0), Q(1), Q(0)],
        [Q(0), Q(0), Q(1)],
    ]

    total = [[Q(0) for _ in range(3)] for _ in range(3)]
    for axis in axes:
        for i in range(3):
            for j in range(3):
                total[i][j] += gp * identity[i][j] + anisotropy * axis[i] * axis[j]

    scalar = 3 * gp + anisotropy
    assert total == [[scalar * identity[i][j] for j in range(3)] for i in range(3)]

    trace = sum(total[i][i] for i in range(3))
    trace_free = [
        [total[i][j] - Q(int(i == j)) * trace / 3 for j in range(3)]
        for i in range(3)
    ]
    assert trace_free == [[Q(0) for _ in range(3)] for _ in range(3)]

    # Central inversion preserves a second-order tensor moment and reverses
    # a first spatial moment.
    first_moment = [Q(2), Q(-3), Q(5)]
    inverted = [-entry for entry in first_moment]
    assert [left + right for left, right in zip(first_moment, inverted)] == [Q(0)] * 3
    print("three-axis quadrupole and inversion-pair dipole cancellations: exact")


def packed_stage_power(first_surviving_moment: int) -> Power:
    q = first_surviving_moment
    # tau^2 * || |y|^q S ||_1 for one microtorus:
    # nu*T^2*ell^2*delta^(q+1), delta=ell/K.
    one = power(1, 0, q + 3, -(q + 1), 0, 2, 0)
    population = power(0, 0, 0, 3, 0, 0, 0)
    kernel_at_macro_distance = power(0, 0, -(q + 4), 0, 0, 0, 0)
    return add(one, population, kernel_at_macro_distance)


def check_packed_wake_and_moment_gains() -> None:
    assert packed_stage_power(0) == power(1, 0, -1, 2, 0, 2, 0)
    assert packed_stage_power(1) == power(1, 0, -1, 1, 0, 2, 0)
    assert packed_stage_power(2) == power(1, 0, -1, 0, 0, 2, 0)

    # A single microtorus evaluated one micro-spacing away.
    # nu*T^2*ell^2*delta * delta^-4 = nu*T^2*K^3/ell.
    one_moment = power(1, 0, 3, -1, 0, 2, 0)
    micro_distance_kernel = power(0, 0, -4, 4, 0, 0, 0)
    neighbor = add(one_moment, micro_distance_kernel)
    assert neighbor == power(1, 0, -1, 3, 0, 2, 0)

    # Exterior L2 triangle bound:
    # one tail is nu*T^2*ell^2*delta^-3/2, then multiply by K^3.
    one_l2 = power(1, 0, Q(1, 2), Q(3, 2), 0, 2, 0)
    all_l2 = add(one_l2, power(0, 0, 0, 3, 0, 0, 0))
    assert all_l2 == power(1, 0, Q(1, 2), Q(9, 2), 0, 2, 0)
    print("packed K^2 outer, K^3 neighbor, and clustered K^0 powers: exact")


def check_reynolds_identities() -> None:
    # v=ell^-gamma*K^gamma.
    for gamma in (Q(11, 10), Q(5, 4), Q(7, 5)):
        inverse_velocity = power(0, -1, 0, 0, 0, 0, 0)
        outer = power(1, 0, -1, 2, 0, 2, 0)
        neighbor = power(1, 0, -1, 3, 0, 2, 0)

        # Substitute v^-1=ell^gamma*K^-gamma.
        substitution = power(0, 1, gamma, -gamma, 0, 0, 0)
        assert add(outer, inverse_velocity, substitution) == power(
            1, 0, gamma - 1, 2 - gamma, 0, 2, 0
        )
        assert add(neighbor, inverse_velocity, substitution) == power(
            1, 0, gamma - 1, 3 - gamma, 0, 2, 0
        )

    # Theta=nu*T*K^2/(v*ell*eps^2).
    theta = power(1, -1, -1, 2, -2, 1, 0)
    outer_from_theta = add(theta, power(0, 0, 0, 0, 2, 1, 0))
    neighbor_from_theta = add(theta, power(0, 0, 0, 1, 2, 1, 0))
    assert outer_from_theta == power(1, -1, -1, 2, 0, 2, 0)
    assert neighbor_from_theta == power(1, -1, -1, 3, 0, 2, 0)
    print("thin viscous parameter and wake-ratio identities: exact")


def main() -> None:
    check_thin_source_and_pressure_scaling()
    check_orthogonal_axis_cancellation()
    check_packed_wake_and_moment_gains()
    check_reynolds_identities()
    print("all thin-torus exterior-wake checks passed")


if __name__ == "__main__":
    main()
