#!/usr/bin/env python3
"""Dependency-free component checks for generalized_phase_localization.md."""

from __future__ import annotations

import math
import random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def matvec(A, x):
    return tuple(sum(A[i][j] * x[j] for j in range(3)) for i in range(3))


def trace(A):
    return sum(A[i][i] for i in range(3))


def vnorm(a):
    return math.sqrt(sum(abs(x) ** 2 for x in a))


def close(name, a, b, tol=2e-10):
    if not isinstance(a, tuple):
        a, b = (a,), (b,)
    err = vnorm(sub(a, b))
    size = 1.0 + vnorm(a) + vnorm(b)
    if err > tol * size:
        raise AssertionError(f"{name}: error={err}, scale={size}\n{a}\n{b}")


def falling(n, k):
    out = 1
    for j in range(k):
        out *= n - j
    return out


# Monomials of the mixed cubic phase used in the symbolic derivation.
MONOMIALS = (
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (3, 0, 0), (2, 1, 0), (1, 1, 1),
    (1, 0, 2), (0, 3, 0), (0, 1, 2),
)


def derivative_value(coeffs, point, order):
    total = 0.0
    for coeff, powers in zip(coeffs, MONOMIALS):
        term = coeff
        for n, k, xx in zip(powers, order, point):
            if k > n:
                term = 0.0
                break
            term *= falling(n, k) * xx ** (n - k)
        total += term
    return total


def phase_jet(coeffs, point):
    q = tuple(
        derivative_value(coeffs, point, tuple(int(i == j) for i in range(3)))
        for j in range(3)
    )
    Hess = tuple(tuple(
        derivative_value(
            coeffs,
            point,
            tuple(int(i == j) + int(i == k) for i in range(3)),
        )
        for k in range(3)) for j in range(3))
    third = tuple(tuple(tuple(
        derivative_value(
            coeffs,
            point,
            tuple(int(i == j) + int(i == k) + int(i == ell)
                  for i in range(3)),
        )
        for ell in range(3)) for k in range(3)) for j in range(3))
    return q, Hess, third


def outer_phase_data(q, Hess, third, H, Lam):
    t = (0.0, 1.0, 0.0)
    w = (1.0, 0.0, H)
    c = cross(w, t)
    qt = dot(t, q)
    d = dot(w, q)
    gqt = matvec(Hess, t)
    gd = matvec(Hess, w)
    v = sub(scale(qt, w), scale(d, t))
    Jv = tuple(tuple(w[i] * gqt[j] - t[i] * gd[j]
                       for j in range(3)) for i in range(3))
    k = {
        sigma: add(scale(Lam, t), scale(0.5 * sigma, q))
        for sigma in (-1, 1)
    }
    a = {
        sigma: add(w, scale(sigma / (2.0 * Lam), v))
        for sigma in (-1, 1)
    }
    Ja = {
        sigma: tuple(tuple(sigma * Jv[i][j] / (2.0 * Lam)
                           for j in range(3)) for i in range(3))
        for sigma in (-1, 1)
    }
    gradn = tuple(
        sum(third[j][i][ell] * w[i] * w[ell]
            for i in range(3) for ell in range(3))
        for j in range(3)
    )
    return t, w, c, qt, d, gd, v, Jv, k, a, Ja, gradn


def B_phase(u, v, Jv, k_v):
    return add(matvec(Jv, u), scale(1j * dot(u, k_v), v))


def check_core_identities(q, Hess, third, H, Lam):
    (t, w, c, qt, d, gd, v, Jv, k, a, Ja, gradn) = outer_phase_data(
        q, Hess, third, H, Lam
    )

    for sigma in (-1, 1):
        a_cross = scale(1.0 / Lam, cross(k[sigma], c))
        close(f"curl polarization {sigma}", a[sigma], a_cross)
        close(f"phase transversality {sigma}", dot(a[sigma], k[sigma]), 0.0)
        close(f"a dot q {sigma}", dot(a[sigma], q), d)
        close(f"div a {sigma}", trace(Ja[sigma]), 0.0)

    Dvv = matvec(Jv, v)
    low = add(
        B_phase(a[1], a[-1], Ja[-1], scale(-1.0, k[-1])),
        B_phase(a[-1], a[1], Ja[1], k[1]),
    )
    low_expected = sub(scale(2j * d, w), scale(1.0 / (2 * Lam**2), Dvv))
    close("low identity", low, low_expected)

    high = add(
        B_phase(a[1], a[-1], Ja[-1], k[-1]),
        B_phase(a[-1], a[1], Ja[1], k[1]),
    )
    high_expected = sub(scale(1j * d / Lam, v),
                        scale(1.0 / (2 * Lam**2), Dvv))
    close("cross-high identity", high, high_expected)

    for sigma in (-1, 1):
        self_amp = B_phase(a[sigma], a[sigma], Ja[sigma], k[sigma])
        self_expected = add(
            scale(sigma / (2 * Lam), matvec(Jv, w)),
            scale(1.0 / (4 * Lam**2), Dvv),
        )
        close(f"self identity {sigma}", self_amp, self_expected)

    # Cross-high pressure decomposition.
    pi_high = -d**2 / (2 * Lam**2)
    grad_pi_high = scale(-d / Lam**2, gd)
    grad_phase_high = add(grad_pi_high, scale(2j * Lam * pi_high, t))
    R_high = sub(high, grad_phase_high)
    R_high_expected = add(
        scale(1j * d * qt / Lam, w),
        scale(1.0 / (2 * Lam**2), sub(scale(2 * d, gd), Dvv)),
    )
    close("cross-high pressure", R_high, R_high_expected)

    m = dot(w, matvec(Hess, t))
    n = dot(w, matvec(Hess, w))
    close("mixed derivatives", m, dot(t, gd))
    for sigma in (-1, 1):
        self_amp = B_phase(a[sigma], a[sigma], Ja[sigma], k[sigma])
        pi_self = sigma * 1j * n / (4 * Lam**2)
        grad_pi_self = scale(sigma * 1j / (4 * Lam**2), gradn)
        grad_phase_self = add(grad_pi_self, scale(2j * pi_self, k[sigma]))
        R_self = sub(self_amp, grad_phase_self)
        bracket = sub(add(Dvv, scale(n, q)), scale(sigma * 1j, gradn))
        R_self_expected = add(
            scale(sigma * m / (2 * Lam), w),
            scale(1.0 / (4 * Lam**2), bracket),
        )
        close(f"self pressure {sigma}", R_self, R_self_expected)


def check_exact_curl_divergence(q, Hess, third, H, Lam, point):
    del third
    t, w, c, _, _, _, _, _, k, a, Ja, _ = outer_phase_data(
        q, Hess, (((0.0,) * 3,) * 3,) * 3, H, Lam
    )
    x, y, z = point
    A = 1 + 2*x - 3*y + 5*z + 7*x*y + 11*y*z + 13*x*z + 17*x*x
    gA = (2 + 7*y + 13*z + 34*x, -3 + 7*x + 11*z, 5 + 11*y + 13*x)
    HA = ((34.0, 7.0, 13.0), (7.0, 0.0, 11.0), (13.0, 11.0, 0.0))

    for sigma in (-1, 1):
        bA = scale(-1j / Lam, cross(gA, c))
        U = add(scale(A, a[sigma]), bA)
        Jb = tuple(tuple(
            scale(-1j / Lam, cross(tuple(HA[i][j] for i in range(3)), c))[row]
            for j in range(3)) for row in range(3))
        JU = tuple(tuple(
            gA[j] * a[sigma][i] + A * Ja[sigma][i][j] + Jb[i][j]
            for j in range(3)) for i in range(3))
        div_phase = trace(JU) + 1j * dot(k[sigma], U)
        close(f"exact curl divergence {sigma}", div_phase, 0.0)


def check_quadratic_core(rng):
    t = (0.0, 1.0, 0.0)
    r = (1.0, 0.0, 0.0)
    h = (0.0, 0.0, 1.0)
    for _ in range(100):
        alpha = rng.uniform(-1.5, 1.5)
        x = rng.uniform(-1.0, 1.0)
        H = rng.uniform(-2.0, 2.0)
        Lam = rng.uniform(8.0, 30.0)
        q = (2 * alpha * x, 0.0, 0.0)
        Hess = ((2 * alpha, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        third = tuple(tuple(tuple(0.0 for _ in range(3))
                            for _ in range(3)) for _ in range(3))
        (_, w, _, _, d, gd, v, Jv, k, a, Ja, gradn) = outer_phase_data(
            q, Hess, third, H, Lam
        )
        Dvv = matvec(Jv, v)
        low = add(
            B_phase(a[1], a[-1], Ja[-1], scale(-1.0, k[-1])),
            B_phase(a[-1], a[1], Ja[1], k[1]),
        )
        close("quadratic low", low, scale(4j * alpha * x, w))
        close("quadratic child", sub(low, scale(2j, q)),
              scale(4j * alpha * H * x, h))

        high = add(
            B_phase(a[1], a[-1], Ja[-1], k[-1]),
            B_phase(a[-1], a[1], Ja[1], k[1]),
        )
        pi_high = -d**2 / (2 * Lam**2)
        grad_phase_high = add(scale(-d / Lam**2, gd),
                              scale(2j * Lam * pi_high, t))
        close("quadratic high remainder", sub(high, grad_phase_high),
              scale(4 * alpha**2 * x / Lam**2, r))

        n = dot(w, matvec(Hess, w))
        for sigma in (-1, 1):
            self_amp = B_phase(a[sigma], a[sigma], Ja[sigma], k[sigma])
            close("quadratic self before pressure", self_amp,
                  scale(-sigma * alpha / Lam, t))
            pi_self = sigma * 1j * n / (4 * Lam**2)
            grad_phase_self = add(scale(sigma * 1j / (4 * Lam**2), gradn),
                                  scale(2j * pi_self, k[sigma]))
            close("quadratic self remainder", sub(self_amp, grad_phase_self),
                  scale(alpha**2 * x / Lam**2, r))


def main():
    rng = random.Random(20260801)
    for _ in range(120):
        coeffs = [rng.uniform(-1.0, 1.0) for _ in MONOMIALS]
        point = tuple(rng.uniform(-0.8, 0.8) for _ in range(3))
        H = rng.uniform(-2.0, 2.0)
        Lam = rng.uniform(7.0, 40.0)
        q, Hess, third = phase_jet(coeffs, point)
        check_core_identities(q, Hess, third, H, Lam)
        check_exact_curl_divergence(q, Hess, third, H, Lam, point)
    check_quadratic_core(rng)
    print("PASS: generalized-phase curl, interactions, pressures, and quadratic core")


if __name__ == "__main__":
    main()
