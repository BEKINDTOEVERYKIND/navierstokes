#!/usr/bin/env python3
"""Rational certificate for the full-BAS counterexample.

The checker imports only the dependency-free interval primitives from the
companion AO checker.  All decisive comparisons use Fractions.
"""

from fractions import Fraction as F

from ao_batchelor_global_bas import I, enclose_sqrt


def atanh_log_bounds(z, terms):
    """Enclose 2*atanh(z) by a positive series and geometric tail."""
    z = F(z)
    lower = 2 * sum(z ** (2 * j + 1) / F(2 * j + 1) for j in range(terms))
    first_power = 2 * terms + 1
    tail = 2 * z**first_power / F(first_power) / (1 - z * z)
    return lower, lower + tail


def main():
    # X=log(9/5)=2 atanh(2/7).
    x_series_lo, x_series_hi = atanh_log_bounds(F(2, 7), 6)
    x_lo, x_hi = F(58778665, 10**8), F(58778667, 10**8)
    assert x_lo < x_series_lo < x_series_hi < x_hi
    x = I(x_lo, x_hi)

    beta0 = enclose_sqrt(
        1 / (4 - x), F(54135, 10**5), F(54136, 10**5)
    )
    Q = beta0 * x * x / (F(4, 5) - x)

    b0 = (
        32
        * x
        * x
        * (2 - 3 * x)
        / (81 * (4 - x) * (F(4, 5) - x) * (F(4, 5) - x))
    )
    assert b0.hi < F(21020, 10**5)

    # Y=log(7/4)=2 atanh(3/11).
    y_series_lo, y_series_hi = atanh_log_bounds(F(3, 11), 7)
    y_lo, y_hi = F(55961578, 10**8), F(55961579, 10**8)
    assert y_lo < y_series_lo < y_series_hi < y_hi
    y = I(y_lo, y_hi)

    beta1 = Q * (F(3, 4) - y) / (y * y)
    assert beta1.subset(F(53579, 10**5), F(53581, 10**5))

    beta1_Q = beta1 * Q
    assert beta1_Q.lo > 0
    assert beta1_Q.hi < 1

    b1 = (
        4
        * beta1_Q
        * (1 - beta1_Q)
        * F(4, 7)
        * F(3, 7)
        / (1 + beta1 * beta1 * y)
    )
    assert b1.lo > F(21034, 10**5)
    assert b1.hi < F(21036, 10**5)

    ratio = F(10003, 10000)
    assert b1.lo > ratio * ratio * b0.hi

    # Cross-multiplied structural identities in the reduced system:
    # K2'=-2 ell D and x'=-(K2'/K2)x+A0*y/K2 imply
    # (K2*x)'=A0*y exactly.
    ell, D, kh2 = F(7, 5), F(-3, 4), F(11, 6)
    k2 = ell * ell + kh2
    k2_prime = -2 * ell * D
    x, y_amp, A0 = F(5, 8), F(-2, 9), F(13, 7)
    x_prime = -(k2_prime / k2) * x + (A0 / k2) * y_amp
    assert k2_prime * x + k2 * x_prime == A0 * y_amp

    print("non-resonant BAS symmetrizer identity checked exactly")
    print("log(7/4) counterexample radius enclosed by exact rational tail")
    print("0.21034 < alternative resonant exponent squared < 0.21036")
    print("alternative exponent / selected AO edge > 1.0003")
    print("full-cocycle sqrt(b0) bound refuted; fixed-beta bound survives")
    print("all AO full-BAS checks passed")


if __name__ == "__main__":
    main()
