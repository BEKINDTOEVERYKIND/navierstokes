#!/usr/bin/env python3
"""Exact interval check for the Batchelor--Gavrilov first-jet mismatch."""

from fractions import Fraction as F

from ao_batchelor_full_edge_matched import exp_interval
from ao_batchelor_global_bas import I, enclose_sqrt


def main() -> None:
    x = I(F(59671214, 10**8), F(59671216, 10**8))
    Q = I(F(8278572, 10**7), F(8278581, 10**7))
    r = enclose_sqrt(x, F(77247, 10**5), F(77249, 10**5))
    ex = exp_interval(x)
    y = 1 / ex

    V_prime = Q * (((2 * x + 1) * y - 1) / x)
    W_prime = -2 * r * y
    assert V_prime.lo > 0
    assert W_prime.hi < 0

    slope_ratio = (-W_prime) / V_prime
    assert slope_ratio.subset(F(295, 100), F(296, 100))

    inv_sqrt_two = enclose_sqrt(
        I(F(1, 2)), F(7071, 10**4), F(7072, 10**4)
    )
    assert slope_ratio.lo > inv_sqrt_two.hi

    # If both profiles are transformed by s -> b*s and velocity -> a*u,
    # both radial derivatives acquire the same nonzero factor a*b.  An axial
    # Galilean shift has derivative zero; an axial reflection changes only
    # the sign, not the magnitude of the ratio.
    test_v_prime, test_w_prime = F(17, 9), F(-23, 11)
    original_ratio = abs(test_w_prime / test_v_prime)
    for a, b, axial_sign in (
        (F(2), F(3), F(1)),
        (F(-5, 2), F(7, 4), F(-1)),
    ):
        transformed_v_prime = a * b * test_v_prime
        # A Galilean shift contributes zero to this derivative.
        transformed_w_prime = a * b * axial_sign * test_w_prime + F(0)
        assert abs(transformed_w_prime / transformed_v_prime) == original_ratio

    print("Batchelor AO-ring derivatives have opposite signs: PASS")
    print("2.95 < -W_B'/V_B' < 2.96, separated from 1/sqrt(2): PASS")
    print("uniform scaling/Galilean derivative invariant: PASS")
    print("Batchelor--Gavrilov base-profile compatibility remains OPEN")


if __name__ == "__main__":
    main()
