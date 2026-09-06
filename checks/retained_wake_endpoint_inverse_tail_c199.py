#!/usr/bin/env python3
"""Exact arithmetic support for C199's full-mode return inverse criterion.

Certifies the constants in the energy calculation, a noncommuting block
inverse with both low/high couplings present, and the all-mode heat example.
It does not certify a nonzero reference solution, actual stage chart, or
the missing finite-section singular-value enclosure.
"""

from fractions import Fraction as F
from math import factorial


def identity(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, factor):
    return [[factor * x for x in row] for row in a]


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def inverse(a):
    n = len(a)
    rows = [list(row) + unit for row, unit in zip(a, identity(n))]
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        value = rows[col][col]
        rows[col] = [x / value for x in rows[col]]
        for i in range(n):
            if i != col:
                value = rows[i][col]
                rows[i] = [x - value * y
                           for x, y in zip(rows[i], rows[col])]
    return [row[n:] for row in rows]


def frobenius_squared(a):
    return sum(x * x for row in a for x in row)


def young_constant():
    # Exact polynomial coefficient identity after multiplying by mu:
    # (mu^2/2)v^2 + 2 M0^2 w^2 - 2 mu M0 vw
    # = (mu v - 2 M0 w)^2/2.
    # Tuples are exponents of (mu, v, M0, w).
    terms = {(1, 1, 0, 0): F(1), (0, 0, 1, 1): F(-2)}
    square = {}
    for left, lc in terms.items():
        for right, rc in terms.items():
            powers = tuple(x + y for x, y in zip(left, right))
            square[powers] = square.get(powers, F(0)) + lc * rc / 2
    assert square == {
        (2, 2, 0, 0): F(1, 2),
        (1, 1, 1, 1): F(-2),
        (0, 0, 2, 2): F(2),
    }
    # Doubling y'/2 + mu f <= mu f/2 + 2 M0^2 y/mu
    # gives precisely y' + mu f <= 4 M0^2 y/mu.
    assert 2 * F(2) == 4


def full_coupling_inverse():
    # This nonsymmetric rational model has neither off-block coupling zero.
    l = [
        [F(3, 2), F(1, 4), F(1, 40)],
        [F(0), F(1, 4), F(1, 40)],
        [F(1, 2), F(0), F(1, 20)],
    ]
    unit = identity(3)
    p = [[F(1), F(0), F(0)],
         [F(0), F(1), F(0)],
         [F(0), F(0), F(0)]]
    q = add(unit, scale(p, -1))
    a = add(unit, scale(l, -1))
    low = [row[:2] for row in a[:2]]
    low_inv = inverse(low)
    embed = [row + [F(0)] for row in low_inv] + [[F(0)] * 3]
    b = add(unit, scale(mul(l, p), -1))
    woodbury = add(unit, mul(l, embed))
    assert mul(b, woodbury) == mul(woodbury, b) == unit
    assert woodbury == inverse(b)
    assert add(b, scale(mul(l, q), -1)) == a
    assert frobenius_squared(mul(mul(p, l), q)) > 0
    assert frobenius_squared(mul(mul(q, l), p)) > 0

    m, s, delta = F(2), F(3, 8), F(1, 16)
    assert frobenius_squared(l) <= m * m
    assert frobenius_squared(low_inv) == F(56, 9) <= 1 / (s * s)
    assert frobenius_squared(mul(l, q)) == F(3, 800) <= delta * delta
    lam = 1 + m / s
    assert lam == F(19, 3)
    assert delta * lam == F(19, 48) < 1
    bound = lam / (1 - delta * lam)
    assert bound == F(304, 29)
    assert frobenius_squared(inverse(a)) <= bound * bound


def all_mode_heat_inverse():
    # The Taylor tail from n=4 has consecutive ratios at most 1/5.
    partial = sum(F(1, factorial(n)) for n in range(4))
    assert partial == F(8, 3)
    tail_majorant = F(1, factorial(4)) / (1 - F(1, 5))
    assert partial + tail_majorant == F(87, 32) < F(11, 4)
    lower_e, upper_e = F(8, 3), F(11, 4)
    first_gap = 3 / upper_e - 1
    higher_gap = 1 - 3 / (lower_e * lower_e)
    assert first_gap == F(1, 11)
    assert higher_gap == F(37, 64) > first_gap
    assert 1 / first_gap == 11
    # Every nonzero integer frequency has |k|^2 = 1 or >=2;
    # the second bound decreases no further as |k|^2 increases.
    for squared_frequency in range(2, 101):
        assert 1 - 3 / lower_e**squared_frequency >= higher_gap


def conditional_cutoff_arithmetic():
    # A parameter sanity check, conditional on the separately needed s.
    mu, time, m0 = F(1, 16), F(1), F(1, 4)
    smoothing_squared = 1 / (mu * time) + 4 * m0 * m0 / (mu * mu)
    assert smoothing_squared == 80 < 9**2
    m, s, cutoff = F(1), F(1, 2), F(54)
    lam = 1 + m / s
    delta_upper = m * 9 / cutoff
    assert delta_upper * lam == F(1, 2)
    assert lam / (1 - delta_upper * lam) == 6


def main():
    young_constant()
    full_coupling_inverse()
    all_mode_heat_inverse()
    conditional_cutoff_arithmetic()
    print("C199 retained-wake return inverse constants: PASS")
    print("H^-1 dissipation coefficient: 4 M0^2 / mu")
    print("FULL LOW/HIGH COUPLING: exact inverse identity verified")
    print("ALL-MODE HEAT EXAMPLE: ||(I-3 exp(Delta))^-1|| <= 11")
    print("BOUNDARY: actual finite block and viscosity-uniform inverse open")


if __name__ == "__main__":
    main()
