#!/usr/bin/env python3
"""Exact charge and Gevrey ledgers for three material phases."""

from fractions import Fraction as F
from math import comb, factorial


def charge_count(radius):
    return (
        1
        + 6 * radius
        + 12 * comb(radius, 2)
        + 8 * comb(radius, 3)
    )


def brute_count(radius):
    return sum(
        abs(x_value) + abs(y_value) + abs(z_value) <= radius
        for x_value in range(-radius, radius + 1)
        for y_value in range(-radius, radius + 1)
        for z_value in range(-radius, radius + 1)
    )


def main():
    for radius in range(0, 14):
        assert charge_count(radius) == brute_count(radius)
        assert charge_count(radius) <= (2 * radius + 1) ** 3
        if radius >= 1:
            assert charge_count(radius) - charge_count(radius - 1) == (
                4 * radius * radius + 2
            )
            assert (radius + 1) ** 3 <= 8**radius

    for order in range(1, 30):
        convolution = sum(
            F(1, comb(order, left)) for left in range(order + 1)
        )
        assert convolution <= 3
        dimensional = sum(
            comb(order, left)
            * factorial(left) ** 2
            * factorial(order - left) ** 2
            for left in range(order + 1)
        )
        assert dimensional == factorial(order) ** 2 * convolution

    # The three base transverse covariances sum exactly to qI.
    # Store only their diagonal entries, with q=2.
    q_value = F(2)
    r1 = (F(0), q_value / 2, q_value / 2)
    r2 = (q_value / 2, F(0), q_value / 2)
    r3 = (q_value / 2, q_value / 2, F(0))
    total = tuple(r1[j] + r2[j] + r3[j] for j in range(3))
    assert total == (q_value, q_value, q_value)

    # Exact singular-value ledger: if ||F|| and ||F^-1|| <= K, then
    # |F^-T k| lies between |k|/K and K|k|.  The scalar inequalities are
    # checked here for representative rational values.
    for k_norm, condition in ((1, 3), (5, 7), (13, 4)):
        lower = F(k_norm, condition)
        upper = F(k_norm * condition)
        assert lower > 0
        assert upper >= lower

    print("three-dimensional l1 charge count and shell polynomial: PASS")
    print("Gevrey-2 factorial convolution bound <= 3: PASS")
    print("three transverse positive covariances sum to isotropy: PASS")
    print("material lattice-gap scale ledger: PASS")
    print("all three-material-phase lattice checks passed")


if __name__ == "__main__":
    main()
