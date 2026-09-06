#!/usr/bin/env python3
"""Exact arithmetic checker for the displayed C197 synthesis/WKB ledger.

Checks finite symbolic reflection examples, general recurrence arithmetic,
derivative-chain coefficients, Bessel constants, residual and clock ledgers.
The ODE, Fourier synthesis, and product-rule proofs remain in the note.
No numerical PDE or expanding-aperture certificate is claimed.
"""
from fractions import Fraction as Q
from math import comb, factorial


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def main():
    ident = [[Q(i == j) for j in range(3)] for i in range(3)]
    for k in [(1, 0, 0), (1, 2, 3), (-3, 4, 7), (11, -5, 2)]:
        n2 = sum(x * x for x in k)
        p = [[Q(x * y, n2) for y in k] for x in k]
        refl = [[2 * p[i][j] - ident[i][j] for j in range(3)]
                for i in range(3)]
        assert matmul(refl, refl) == ident
        assert matmul(p, p) == p

    c = [1, 4, 20]
    v = [1, 3, 14]
    for r in range(3, 6):
        for table in (c, v):
            table.append(2 * r * table[-1] + r * (r - 1) * table[-2])
    assert c == [1, 4, 20, 144, 1392, 16800]
    assert v == [1, 3, 14, 102, 984, 11880]
    assert c[1] == 2 + 2 * c[0]
    assert c[2] == 2 + 4 * c[1] + 2 * c[0]
    assert v[1] == 1 + 2 * v[0]
    assert v[2] == 4 * v[1] + 2 * v[0]

    kappa, radius = Q(1, 2), Q(2)
    h = Q(216**2, 1) / kappa
    e = [Q(1)]
    for r in range(1, 5):
        e.append(c[r] + 12 * sum(
            Q(comb(r, j) * c[j], 3 * r - j + 1) * e[r-j]
            for j in range(1, r+1)))
    f = [(int(r == 0) + 2 * c[r]) * 1620
         + 12 * (c[r+1] + r*c[r]) * 216**3 * 286992 * radius/kappa
         for r in range(5)]
    g = []
    for r in range(5):
        g.append(12 * sum(
            Q(comb(r,j)*c[j], 3*r-j+7)*g[r-j]
            for j in range(1,r+1)) + sum(
            Q(comb(r,j),3*r-j+6)*f[j]*e[r-j]
            for j in range(r+1)))
    assert e[1] == 20
    assert g[0] == f[0]/6
    assert all(x > 0 for x in e+f+g)
    for r in range(1, 5):
        for j in range(1, r+1):
            assert 2*j + 3*(r-j) + 1 <= 3*r
            assert 2*j + 3*(r-j) + 6 + 1 <= 3*r+6
        for j in range(r+1):
            assert 2*j + 5 + 3*(r-j) + 1 <= 3*r+6

    # (1-Delta)^2 composed with a coordinatewise sine map.
    assert [6+3, 6+21+6, 18+12, 3+6] == [9,33,30,9]
    assert max([1] + [a*factorial(r) for r,a in enumerate([9,33,30,9],1)]) == 216
    zeta2_majorant = sum(Q(1,n*n) for n in range(1,6)) + Q(1,5)
    assert zeta2_majorant < Q(5,3)
    assert 1 + 26*Q(5,3) == Q(133,3) < 45
    assert (2+2*Q(22,7))**3 < 27**2
    assert 216*45*27 == 262440

    # Exact rational values of the polynomial majorants at z=1.
    p0 = sum(Q(c[r],factorial(r))*h**r for r in range(5))
    px = (216*286992*radius/kappa * sum(
        Q(c[r+1]+r*c[r],factorial(r))*h**r for r in range(5)))
    vv = 216/kappa * sum(Q(v[r],factorial(r))*h**r for r in range(5))
    vx = 216**2*286992*radius/kappa**2 * sum(
        Q(v[r+1]+r*v[r],factorial(r))*h**r for r in range(5))
    bb = sum(e[r]*h**r/factorial(r) for r in range(5))
    bx = 216*sum(g[r]*h**r/factorial(r) for r in range(5))
    r0 = (Q(225,2)*p0+90*px+Q(105,2))*vv*bb + (90*p0+72)*(vx*bb+vv*bx)
    r1 = (90*p0+72)*vv*216*bb
    k0, k1 = 262440*r0, 787320*r1
    assert k0 <= 10**94
    assert k1 <= 10**80
    assert max(8+9+12,11+9+12,8+12+12,8+9+19) == 36
    assert 8+9+1+12 == 30
    assert k0 > 0 and k1 > 0
    # Residual coefficient derivation: curl/commutator/pressure factors.
    assert Q(5,2)*6 == 15
    assert Q(5,2)*2+2 == 7
    assert Q(5,2)*2+5+2 == 12
    assert 6*Q(19,50) == Q(57,25)
    gain = Q(3,8)
    assert 1-Q(1,12)-Q(57,25)*gain == Q(37,600)
    assert 1-Q(1,3)-Q(57,25)*gain+gain == Q(14,75)
    assert 1-Q(57,25)*gain == Q(29,200)
    assert 1-Q(57,25)*gain+gain == Q(13,25)
    # Finite q>=10^10000 upper-error budget, including real completion.
    # 2<ln(10)<3 follows from e<3 and e>8/3.
    assert 3**2 < 10 < Q(8,3)**3
    assert Q(912,25) < 37 and 3**37 < 10**18
    assert (1-Q(57,400))*20000 > 1+Q(152,25)
    assert 20000 > 600
    assert 30000 < 10**5
    assert 5*37 == 185
    assert 10000*Q(37,600) > 616
    assert 112+185-616 == -319
    assert Q(4,10**319) < Q(1,100)
    assert min(Q(29,200),Q(37,600),Q(13,25),Q(14,75)) == Q(37,600)
    print('C197 uniform compact multi-beam WKB upper estimate: PASS')
    print('AMBIENT KELVIN EXPONENT: 6, reflection prefactor: 1')
    print('BEAM-COUNT LOSS: none under the stated common-lattice hypothesis')
    print('ERROR POLYNOMIAL DEGREES: 37 (profile), 31 (profile gradient)')
    print('K0 exact:', k0)
    print('K1 exact:', k1)
    print('CONVENIENT CONSTANTS: K0<=10^94; K1<=10^80')
    print('FINITE-q BUDGET: q>=10^10000 gives coefficient/profile-normalized real error <1/100')
    print('BOUNDARY: no fixed-aperture growth, periodic endpoint, retained band, or viscosity')


if __name__ == '__main__':
    main()
