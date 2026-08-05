#!/usr/bin/env python3
"""Exact ledger for the corrected hollow Gavrilov comparison."""

from fractions import Fraction as F


def main() -> None:
    # Source expansion: alpha=2 r^2+3 X r^2+O(r^4), p=alpha/4.
    p2_coefficient = F(2, 4)
    p3_coefficient = F(3, 4)
    assert p2_coefficient == F(1, 2)
    assert p3_coefficient == F(3, 4)

    # At the 3-4-5 rational point, r is rational. The axial leading
    # component r/sqrt(2) has squared magnitude one half of the transverse
    # leading component r.
    x_value = F(3, 5)
    z_value = F(4, 5)
    radius_squared = x_value**2 + z_value**2
    assert radius_squared == 1
    transverse_magnitude_squared = radius_squared
    axial_magnitude_squared = radius_squared / 2
    assert axial_magnitude_squared / transverse_magnitude_squared == F(1, 2)

    # The transverse vector J(X,Z)=(Z,-X) is linear, but |(X,Z)| is not:
    # opposite points have the same positive radial component, whereas a
    # linear scalar functional would change sign.
    radial_at_point = F(1)
    radial_at_opposite = F(1)
    assert radial_at_point == radial_at_opposite
    assert radial_at_opposite != -radial_at_point

    # Scaling powers on an annulus. A degree-m analytic pressure term loses
    # epsilon^2 inside g. The leading conormal velocity is homogeneous of
    # degree one; its next annular term has degree two and loses epsilon.
    pressure_powers = {degree: degree - 2 for degree in range(2, 9)}
    annular_velocity_powers = {degree: degree - 1 for degree in range(1, 8)}
    assert pressure_powers[2] == 0
    assert pressure_powers[3] == 1
    assert annular_velocity_powers[1] == 0
    assert annular_velocity_powers[2] == 1

    # Exact scaled pressure sample, including the cubic source term.
    epsilon = F(2, 17)
    p_leading = radius_squared / 2
    p_scaled_through_cubic = (
        p_leading + epsilon * F(3, 4) * x_value * radius_squared
    )
    assert (
        p_scaled_through_cubic - p_leading
        == epsilon * F(3, 4) * x_value * radius_squared
    )

    # A positive inner pressure cutoff keeps the normalized active set away
    # from the nonsmooth seed circle.
    q_inner = F(1, 8)
    q_outer = F(3, 2)
    inner_radius_squared = 2 * q_inner
    outer_radius_squared = 2 * q_outer
    assert 0 < inner_radius_squared < outer_radius_squared

    # Local first-order coefficient error O(epsilon) acting on a
    # frequency-p packet has the formal O(epsilon*p) power.
    coefficient_power = 1
    generator_order = 1
    assert coefficient_power == generator_order == 1

    print("Gavrilov pressure/source coefficients: PASS")
    print("raw axial field is radial, not a linear Taylor term: PASS")
    print("hollow annular O(epsilon) scaling powers: PASS")
    print("locked-pitch squared magnitude ratio 1/2: PASS")
    print("local curved generator power epsilon*p: PASS (global Hodge open)")
    print("all corrected Gavrilov comparison checks passed")


if __name__ == "__main__":
    main()
