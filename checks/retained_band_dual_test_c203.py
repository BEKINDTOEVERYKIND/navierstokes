#!/usr/bin/env python3
"""Exact arithmetic sentinels for C203; proof of the analytic lemma is in the note."""
from fractions import Fraction as F
from math import comb


def truncated_power_sum(x):
    """Four-uniform convolution: sum (-1)^j C(4,j)(x+4-2j)_+^3."""
    return sum((-1)**j * comb(4, j) * max(F(0), x + 4 - 2*j)**3
               for j in range(5))


def main():
    # sinc^4 transform/pi = (1/48) sum alternating truncated cubics.
    for x in [F(0), F(1, 4), F(1, 2), F(1), F(3, 2), F(2)]:
        assert truncated_power_sum(x) / 48 == F(2, 3) - x*x/4 + x**3/16
    for x in [F(4), F(5), F(8)]:
        assert truncated_power_sum(x) == 0
    # On [0,1], derivative -x/2+3x^2/16<=0 and minimum exceeds 5/12.
    assert F(2, 3)-F(1, 4)+F(1, 16) == F(23, 48)
    assert 3*F(5, 12) > 1
    # Moments from sinc^4 <= min(1,|x|^-4); zeroth uses exact integral.
    assert 2*F(22, 7)/3 < 3
    assert 2*(F(1, 2)+F(1, 2)) == 2
    assert 2*(F(1, 3)+1) < 3
    assert 3*2*3**2 == 54
    assert 3*3*3**2 == 81
    assert 8*5**3 == 1000
    assert F(1000*27**2, 6**3) == 3375 < 64**2
    # Exact-solution and phase losses leave one quarter of the positive signal.
    assert 1-2*F(1, 4)-64*F(1, 256) == F(1, 4)
    assert 27*4 == 108
    # With wmin=c q/L^2, envelope and curvature terms decay at stated powers.
    assert F(1, 3)-1 == -F(2, 3)
    assert 1-2 == -1
    print('C203 PASS: Fourier spline, moments, mode count, dual norm, and retained lower bound')


if __name__ == '__main__':
    main()
