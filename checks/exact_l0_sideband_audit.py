#!/usr/bin/env python3
"""Dependency-free algebra checks for exact_l0_sideband_audit.md."""

import cmath
import math
import random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def proj(k, v):
    kk = dot(k, k)
    return add(v, scale(-dot(k, v) / kk, k))


def norm(a):
    return math.sqrt(sum(abs(x) ** 2 for x in a))


def close(a, b, tol=2e-11):
    return norm(add(a, scale(-1, b))) <= tol * (1 + norm(a) + norm(b))


def interaction(k, u, ell, v):
    """Projected symmetric Euler vector, omitting the common factor -i."""
    out = add(scale(dot(u, ell), v), scale(dot(v, k), u))
    return proj(add(k, ell), out)


def main():
    er, et, eh = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    for _ in range(100):
        Lam = random.uniform(5.0, 20.0)
        Q = random.uniform(0.1, 0.8) * Lam
        delta = Q / (2 * Lam)
        H = -random.uniform(0.2, 6.0)
        q = scale(Q, er)
        kp = add(scale(Lam, et), scale(0.5, q))
        km = add(scale(Lam, et), scale(-0.5, q))
        pp = add(er, scale(-delta, et))
        pm = add(er, scale(delta, et))
        a = add(pp, scale(H, eh))
        b = add(pm, scale(H, eh))

        assert abs(dot(kp, a)) < 1e-10
        assert abs(dot(km, b)) < 1e-10

        # Central high sum is pure pressure.
        assert norm(interaction(kp, a, km, b)) < 1e-10

        # The low difference is the advertised nonzero e_h beat.
        low = interaction(kp, a, scale(-1, km), b)
        assert close(low, scale(2 * Q * H, eh))

        # Both real-child outward sidebands are nonzero and purely e_h.
        outp = interaction(q, eh, kp, a)
        outm = interaction(scale(-1, q), eh, km, b)
        assert close(outp, scale(Q, eh))
        assert close(outm, scale(-Q, eh))

        # Central nonlinear energy identity.
        x, y, Z = [random.uniform(-3, 3) for _ in range(3)]
        xdot = 0.0
        ydot = Q * Z * x
        Zdot = -2 * Q * x * y
        Edot = 4 * (1 + delta * delta) * x * xdot + 4 * y * ydot + 2 * Z * Zdot
        assert abs(Edot) < 1e-10

        # Exchange phases preserve the high-sum determinant.
        theta, phi = random.uniform(-3, 3), random.uniform(-3, 3)
        xp = cmath.exp(1j * theta) * x
        yp = cmath.exp(1j * theta) * y
        xm = cmath.exp(1j * (theta + phi)) * x
        ym = cmath.exp(1j * (theta + phi)) * y
        assert abs(xm * yp - xp * ym) < 1e-10

        # A child can quench one Bloch phase only by enhancing the opposite
        # phase.  The supremum of the inviscid radicand cannot decrease.
        g = random.uniform(0.2, 4.0)
        Zb = random.uniform(0.1, 3.0)
        rad0 = g * g + 2 * g * Q * Zb / H
        radpi = g * g - 2 * g * Q * Zb / H
        assert max(rad0, radpi) >= g * g

    # Equal-polarization, pressure-pure finite supports: brute-force the
    # support combinatorics for nonzero unit coefficients.  Singletons and
    # mirror pairs n+m=-1 are the only possibilities in this window; the
    # only pair that also writes the first child harmonic is {-1,0}.
    from itertools import combinations

    def high_sum_coeffs(support):
        out = {}
        for n in support:
            for m in support:
                s = n + m
                out[s] = out.get(s, 0) + (n - m) ** 2
        return out

    allowed = []
    sites = list(range(-4, 4))
    for size in range(1, 5):
        for support in combinations(sites, size):
            coeffs = high_sum_coeffs(support)
            if all(value == 0 for s, value in coeffs.items() if s != -1):
                allowed.append(support)
    assert all(len(s) <= 2 for s in allowed)
    assert all(len(s) == 1 or s[0] + s[1] == -1 for s in allowed)
    first_child = [s for s in allowed if len(s) == 2 and s[1] - s[0] == 1]
    assert first_child == [(-1, 0)]

    # Palasek's 2b<beta window gives Q/Lambda -> 0 for the deliberately
    # dissipative Lambda~N^(beta/2) carrier, while the first omitted
    # sign-flip order is not small after derivative leverage.
    for _ in range(100):
        bscale = random.uniform(1.001, 1.12)
        beta = random.uniform(2 * bscale + 1e-3, 2.49)
        alpha = random.uniform(beta + 1e-4, 2.5)
        eps_exp = bscale - beta / 2
        leak_exp = (bscale - 1) * (1 - 2 * (alpha - beta))
        assert eps_exp < 0
        assert leak_exp > 0

    # Conditional central turning formulas satisfy all three invariants.
    for _ in range(100):
        S = random.uniform(1.0, 1.2)
        H = -random.uniform(0.4, 8.0)
        C = S + H * H
        d = random.uniform(0.5, 3.0)
        r = random.uniform(0.01, 0.9)
        Q = r * d * math.sqrt(C)
        P0 = random.uniform(0.5, 4.0)
        ps = (r * r - 1) / (r * r + 1)
        P = ps * P0
        x = math.sqrt(2) * r * P0 / (math.sqrt(C) * (1 + r * r))
        y = 0.0
        Z = -2 * Q * H * P0 / (d * C * (1 + r * r))
        I1 = x * x + P * P / (2 * C)
        I10 = P0 * P0 / (2 * C)
        I2 = Z - Q * H * P / (d * C)
        I20 = -Q * H * P0 / (d * C)
        E = P * P + 2 * (S * x * x + y * y) + Z * Z
        assert abs(I1 - I10) < 1e-10
        assert abs(I2 - I20) < 1e-10
        assert abs(E - P0 * P0) < 1e-10

    print(
        "PASS: pressure cancellation, sideband/Bloch obstructions, "
        "scale ledger, energy, and turning formulas"
    )


if __name__ == "__main__":
    main()
