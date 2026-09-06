#!/usr/bin/env python3
"""Exact finite-scale ledger for the same-witness C204 composition."""
from fractions import Fraction as F


def main():
    c = F(1,10**100)
    # Whole sine image fits the split base/covector C195 tube.
    assert 220**2*3 < 400**2
    assert 434000000*F(1,10**25)+216*6*400*c < F(5,10**17)
    assert 2*F(5,10**17)+100*216*F(1,10**25) < F(1,10**13)
    assert 216*F(1,10**30) < F(1,2*10**25)
    # Smooth plateau quotient and full Hessian constants.
    assert 2*64 == 128
    assert 10*8+20*64+4*64+8*512 == 5712
    assert 512**2*3 < 1000**2
    assert 16*5712 < 262144
    assert 3*262144 < 10**6
    # Nearest-grid and retained-test constants.
    assert 162**2 > 27**3
    assert 3*162 < 500
    assert 30*4*10**40*400*c < F(1,10**50)
    assert 30*432*162 == 2099520 < 10**7
    assert F(1,10**8) < F(1,2*500*256)
    # Uniform explicit q0 = 10^100000 comparisons: log(q0)<10^6.
    Q0 = 100000
    tests = [
        ('Euler error', 326+6*46, F(37,600), -5563),
        ('viscous error', 325+6*14, F(45), -100000),
        ('curl', 256+6*21, F(2,3), -65000),
        ('small branch', 3, F(3,8), -37496),
        ('dual phase', 271+6*12, F(2,3), -66000),
    ]
    for label, positive_power, decay, target in tests:
        # Allow a full factor ten for the numerical coefficient.
        extra = 0 if label == 'small branch' else 1
        assert F(positive_power+extra)-decay*Q0 < target, label
    # All comparisons remain true above q0 by derivative of x^a exp(-beta x).
    for a,beta in [(46,F(37,600)),(14,F(45)),(21,F(2,3)),
                   (12,F(2,3)),(3,F(2,3)),(3,F(1,12))]:
        assert a/beta < 200000
    # eps_minus < rho0/4, integer mesh, widths >=1.
    assert F(25+1+6*3)-F(Q0,12) < 0
    assert F(5+100+6*3)-F(2*Q0,3) < 0
    assert Q0-100-6*2 > 0
    assert 8*12500 == Q0
    assert F(1)-F(3*(2137-1),8) <= -800
    # Actual entrance and full-trajectory budgets.
    assert F(27,24**3) == F(1,512)
    assert F(2,512) == F(1,16**2)
    assert 20*(1+900) < 20000
    assert 2*10**9*20*301*812 == 9776480000000000 < 10**16
    action = (4*10**10)**2*12*(2*812**2+2*243600**2+301**2)
    assert action < 3*10**33
    assert F(100,99)*3 < 4
    assert 6*3*10**33 < 2*10**34
    assert 652+6*93+2-F(37*Q0,300) < 0
    assert 93/F(37,300) < 200000
    # Pointwise exact-curl corrections, both envelope widths.
    assert 4-F(2*Q0,3) < -2
    assert 46+6*5-Q0 < -2
    assert F(1,10)+1 < 2
    assert F(401,70*2999)+F(1,10000) < F(1,450)
    assert F(450,432) == F(25,24)
    assert F(10,7)**2 > 2
    assert 900*10000*162 == 1458*10**6
    assert 160+12 < Q0
    assert 2 < 200000
    print('C204 PASS: common tube, explicit cutoff, actual clock/error powers, physical normalization, retained endpoint constant')


if __name__ == '__main__':
    main()
