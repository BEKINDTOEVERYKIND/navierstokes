#!/usr/bin/env python3
"""Exact arithmetic for C202; PDE identities are proved in its note."""
from fractions import Fraction as F
from math import comb, factorial


def main():
    # Coefficient and H1 energy estimates.
    assert 3 * F(15, 2) ** 2 < 13 ** 2
    assert 6 + F(13, 2) < 19
    assert 12 + F(13, 2) < 19
    assert 27 * (3 + 4 + 2) == 243
    assert 2 * 243 == 486
    # Squared energy integrating factor: e^(12t)e^(26t)=e^(38t).
    assert 2 * 19 - 12 == 26
    assert 12 + 26 == 2 * 19

    # Check actual factorial inequality on a range; its all-n proof is in text.
    # Fourth powers eliminate the fractional exponent and sqrt(n) bound.
    for n in range(4, 100):
        assert factorial(n - 1) ** 4 >= n ** (n - 1)
    assert 1 - F(3 * (2137 - 1), 8) == -800
    assert F(1250) * 8 == 10000
    assert F(152, 25) <= (1 - F(57, 400)) * 20000
    assert F(100, 99) < 2
    assert F(1, 10) ** 99 < F(1, 100)

    # Finite real viscous budget; no floating exponentials.
    assert F(19 * 152, 25) < 116
    assert 3 ** 116 < 10 ** 56
    assert 19 * F(57, 400) == F(1083, 400)
    exponent = 1 - 50 + F(1, 2) + F(1083, 400)
    assert exponent == -F(18317, 400)
    assert exponent < -45
    assert 486 * 10 ** 56 < 10 ** 59
    assert 59 - 45 * 10000 == -449941
    assert F(4, 10 ** 319) + F(1, 10 ** 449941) < F(1, 100)

    # Explicit integer carriers; use q=m^12 so both inherited eps are exact.
    for m in range(2, 20):
        q = m ** 12
        for inv_eps in (m, m ** 4):
            d = inv_eps
            M = q // (8 * d)
            assert F(1, 2 * q) + F(d * M, q) <= F(1, 4)
            assert d >= inv_eps
        d = m ** 4
        count = (2 * (q // (8 * d)) + 1) ** 3
        assert count * 4096 >= q ** 2
    assert F(1, 2 * 4) + F(1, 8) == F(1, 4)
    assert 8 ** 3 * 2 ** 3 == 4096

    # General smooth initial selector extension: C197's mixed majorant
    # dominates propagation of nonzero initial spatial derivatives.
    c = [1, 4, 20, 144, 1392, 16800]
    e = [F(1)]
    for r in range(1, 5):
        e.append(c[r] + 12 * sum(
            F(comb(r, j) * c[j], 3 * r - j + 1) * e[r-j]
            for j in range(1, r+1)))
    f = [(int(r == 0) + 2 * c[r]) * 1620
         + 48 * (c[r+1] + r*c[r]) * 216**3 * 286992
         for r in range(5)]
    g = []
    for r in range(5):
        g.append(12 * sum(
            F(comb(r,j)*c[j], 3*r-j+7)*g[r-j]
            for j in range(1,r+1)) + sum(
            F(comb(r,j),3*r-j+6)*f[j]*e[r-j]
            for j in range(r+1)))
    assert all(g[r] >= e[r] for r in range(5))
    assert all(c[r] >= factorial(r) for r in range(5))
    v = [1, 3, 14, 102, 984]
    W = sum(F(v[r] * 2**(r+1), factorial(r)) for r in range(5))
    assert W == 1654
    ledger = [5+5*W, 2+F(29,2)*W, 5*W]
    assert max(ledger) == 23985
    assert 262440 * max(ledger) < 10**10
    # C200 all-r jets summed through order four.
    assert 150 + 5*10**45 < 6*10**45
    assert 150 + 5*10**45 + 5*10**90 < 6*10**90
    assert 6*10**45 * 4**4 * 10**160 < 2*10**208
    assert 6*10**90 * 4**4 * 10**160 < 2*10**253
    assert 46*F(600,37) < 20000
    assert 320 + 5*46 - 616 == -66
    assert 319 + 5*14 - 45*10000 == -449611
    assert F(8,10**66) + F(4,10**449611) < F(1,100)
    # One fully explicit C-infinity plateau profile.
    assert 4 / F(8,3)**2 < 1
    assert 4**4 / F(8,3)**4 + 2*3**3 / F(8,3)**3 < 8
    assert 9+162 == 171
    assert 72+324+1296+5832 == 7524
    assert 4*171 == 684
    assert 16*7524 == 120384
    assert 3*120384+6*684**2 < 4*10**6
    assert 4*10**6*(F(8,10**66)+F(4,10**449611)) < F(1,10**58)
    print("C202 PASS: periodic patch, H1/viscous constants, actual-stage chain")


if __name__ == "__main__":
    main()
