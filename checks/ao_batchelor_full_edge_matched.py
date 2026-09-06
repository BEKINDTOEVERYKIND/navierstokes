#!/usr/bin/env python3
"""Rational certificate for the full-edge-matched Batchelor profile."""

import math
from fractions import Fraction as F

from ao_batchelor_global_bas import I, enclose_sqrt


def exp_point(value, terms=35):
    """Rigorous positive Taylor enclosure for exp(value), 0<=value<2."""
    value = F(value)
    assert 0 <= value < 2
    term = F(1)
    lower = term
    for j in range(1, terms + 1):
        term *= value / F(j)
        lower += term
    first_omitted = term * value / F(terms + 1)
    ratio_bound = value / F(terms + 2)
    upper = lower + first_omitted / (1 - ratio_bound)
    return lower, upper


def exp_interval(interval):
    lo, _ = exp_point(interval.lo)
    _, hi = exp_point(interval.hi)
    return I(lo, hi)


def J_at(value):
    value = F(value)
    ex = I(*exp_point(value))
    e2x = I(*exp_point(2 * value))
    return (2 * value + 1) * e2x - (7 * value + 2) * ex + (5 * value + 1)


def main():
    x_lo = F(59671214, 10**8)
    x_hi = F(59671216, 10**8)
    left = J_at(x_lo)
    right = J_at(x_hi)
    assert left.hi < 0 < right.lo

    # The analytic proof shows J/x^2 is strictly increasing.  These exact
    # coefficient checks cover its base and induction step.
    a3 = 2**3 * (3 + 1) - (7 * 3 + 2)
    assert a3 == 9
    for n in range(3, 100):
        an = 2**n * (n + 1) - (7 * n + 2)
        an1 = 2 ** (n + 1) * (n + 2) - (7 * (n + 1) + 2)
        assert an > 0
        assert an1 - 2 * an == 2 ** (n + 1) + 7 * n - 5 > 0

    x = I(x_lo, x_hi)
    ex = exp_interval(x)
    h = (ex - 1 - x) / (x * x)
    g = (2 - ex) / ((x + 1) * ex - (2 * x + 1))

    # J=0 is equivalent to g=h/(2+xh); the enclosure must contain zero.
    stationarity_gap = g - h / (2 + x * h)
    assert stationarity_gap.lo < 0 < stationarity_gap.hi

    beta = enclose_sqrt(
        g, F(5101679, 10**7), F(5101682, 10**7)
    )
    Q = beta / h
    assert Q.subset(F(8278572, 10**7), F(8278581, 10**7))

    beta_Q = beta * Q
    assert beta_Q.subset(F(4223461, 10**7), F(4223469, 10**7))
    assert beta_Q.hi < F(1, 2)
    assert h.lo > F(1, 2)  # beta/Q=h, hence beta>Q/2.
    assert Q.hi < 1

    y = 1 / ex
    bstar = (
        4
        * beta_Q
        * (1 - beta_Q)
        * y
        * (1 - y)
        / (1 + beta * beta * x)
    )
    assert bstar.subset(F(2090085, 10**7), F(2090092, 10**7))
    sqrt_bstar = enclose_sqrt(
        bstar, F(4571745, 10**7), F(4571753, 10**7)
    )

    print("unique positive joint-stationarity root bracketed")
    print("0.59671214 < x_* < 0.59671216")
    print("0.5101679 < beta_* < 0.5101682")
    print("0.8278572 < Q_* < 0.8278581")
    print("0.2090085 < full BAS edge squared < 0.2090092")
    print("0.4571745 < full BAS edge < 0.4571753")
    print("AO pitch and positivity gates certified")
    print("all full-edge-matched Batchelor checks passed")


if __name__ == "__main__":
    main()
