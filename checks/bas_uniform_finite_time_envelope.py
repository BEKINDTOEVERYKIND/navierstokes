#!/usr/bin/env python3
"""Arithmetic checks for the full-column BAS finite-time envelope.

This is not a proof of the compactness argument in the companion note.  It
checks the exact normalization, the quadratic excess interval, its cubic
area majorant, and the geometric factor used in the oscillatory sector.
"""

from __future__ import annotations

import math
import random


def close(a: float, b: float, tol: float = 2e-11) -> None:
    assert abs(a - b) <= tol * (1.0 + abs(a) + abs(b)), (a, b)


def check_physical_normalization() -> None:
    rng = random.Random(2026080201)
    for _ in range(500):
        tau = rng.uniform(-8.0, 8.0)
        delta = rng.uniform(-3.0, 3.0)
        omega = rng.uniform(-2.0, 2.0)
        q = rng.uniform(-1.0, 1.0)
        p = math.sqrt(max(0.0, 1.0 - q * q))
        if rng.random() < 0.5:
            p = -p
        S = rng.uniform(-2.0, 2.0)
        T = (delta - p * S) / q if abs(q) > 0.05 else rng.uniform(-2.0, 2.0)
        delta = p * S + q * T
        A = 2.0 * omega * q
        c = p * T - q * (S + 2.0 * omega)
        g = 1.0 + tau * tau

        # Start with arbitrary x,y and compare the derivative obtained from
        # u=sqrt(g)x with the normalized matrix in equation (2.5).
        x = rng.uniform(-2.0, 2.0)
        y = rng.uniform(-2.0, 2.0)
        xdot = 2.0 * tau * delta / g * x + A / g * y
        gdot = -2.0 * tau * delta
        udot_chain = 0.5 * gdot / math.sqrt(g) * x + math.sqrt(g) * xdot
        u = math.sqrt(g) * x
        udot_matrix = delta * tau / g * u + A / math.sqrt(g) * y
        ydot_matrix = c / math.sqrt(g) * u
        close(udot_chain, udot_matrix)
        close(ydot_matrix, c * x)

    print("physical frequency normalization: PASS")


def check_excess_interval_and_area() -> None:
    rng = random.Random(2026080202)
    for _ in range(2000):
        lam = rng.uniform(0.1, 2.0)
        d = rng.uniform(1e-4, 3.0)
        mu = rng.uniform(1e-6, 5.0)
        tau = rng.uniform(-10.0, 10.0)
        g = 1.0 + tau * tau
        a = d * tau / g
        h2 = mu / g
        qplus = 0.5 * (a + math.sqrt(a * a + 4.0 * h2))
        F = lam * lam * g - lam * d * tau - mu
        close(
            F,
            lam * lam * (tau - d / (2.0 * lam)) ** 2
            - (mu + d * d / 4.0 - lam * lam),
        )
        assert (qplus > lam) == (F < 0.0)
        if F < 0.0:
            assert qplus - lam <= (-F) / (lam * g) * (1.0 + 2e-12)

        E = max(0.0, mu + d * d / 4.0 - lam * lam)
        if E == 0.0:
            continue
        center = d / (2.0 * lam)
        radius = math.sqrt(E) / lam
        n = 4000
        width = 2.0 * radius / n
        integral = 0.0
        for j in range(n):
            z = center - radius + (j + 0.5) * width
            gg = 1.0 + z * z
            aa = d * z / gg
            qq = 0.5 * (aa + math.sqrt(aa * aa + 4.0 * mu / gg))
            integral += max(0.0, qq - lam) * width / d
        majorant = 4.0 * E ** 1.5 / (3.0 * d * lam * lam)
        assert integral <= majorant * (1.0 + 2e-7) + 1e-10

    print("quadratic excess interval and cubic-area majorant: PASS")


def check_small_d_scaling() -> None:
    # If mu-lambda^2 <= Ld, the majorant is O(sqrt(d)); the d^2/4 term
    # only improves this at an exact edge point.
    lam = 0.47
    L = 1.3
    values = []
    for j in range(5, 25):
        d = 2.0 ** (-j)
        mu = lam * lam + L * d
        E = mu + d * d / 4.0 - lam * lam
        bound = 4.0 * E ** 1.5 / (3.0 * d * lam * lam)
        values.append(bound / math.sqrt(d))
    assert max(values) / min(values) < 1.02

    exact_edge = []
    for j in range(5, 20):
        d = 2.0 ** (-j)
        E = d * d / 4.0
        bound = 4.0 * E ** 1.5 / (3.0 * d * lam * lam)
        exact_edge.append(bound / (d * d))
    assert max(exact_edge) / min(exact_edge) < 1.0000001
    print("small-D excess is O(sqrt(D)); exact edge is O(D^2): PASS")


def check_geometric_factor() -> None:
    rng = random.Random(2026080203)
    for _ in range(5000):
        tau = rng.uniform(-100.0, 100.0)
        shift = rng.uniform(0.0, 30.0)
        ratio = math.sqrt(1.0 + (tau + shift) ** 2) / math.sqrt(
            1.0 + tau * tau
        )
        assert ratio <= 1.0 + shift + 2e-14
    print("oscillatory-sector geometric factor: PASS")


def check_batchelor_coefficients() -> None:
    rng = random.Random(2026080204)
    for _ in range(1000):
        r = 10.0 ** rng.uniform(-3.0, 1.2)
        Q = rng.uniform(0.2, 2.0)
        x = r * r
        y = math.exp(-x)
        omega = Q * (1.0 - y) / x
        S = 2.0 * Q * ((x + 1.0) * y - 1.0) / x
        close(S + 2.0 * omega, 2.0 * Q * y, 5e-10)

    # At the two noncompact ends the positive part of mu=Ac vanishes
    # uniformly: the axis limit is -4Q^2 q^2 and every factor decays at
    # infinity.  Sample the exact formulas as a regression check.
    Q = 0.88135
    for r in (1e-5, 10.0):
        x = r * r
        y = math.exp(-x)
        omega = Q * (-math.expm1(-x)) / x
        S = 2.0 * Q * ((x + 1.0) * y - 1.0) / x
        T = -2.0 * r * y
        max_positive = 0.0
        for j in range(2000):
            theta = 2.0 * math.pi * j / 2000.0
            p, q = math.cos(theta), math.sin(theta)
            A = 2.0 * omega * q
            c = p * T - q * (S + 2.0 * omega)
            max_positive = max(max_positive, A * c)
        assert max_positive < 1e-7
    print("Batchelor coefficient identities and endpoint decay: PASS")


if __name__ == "__main__":
    check_physical_normalization()
    check_excess_interval_and_area()
    check_small_d_scaling()
    check_geometric_factor()
    check_batchelor_coefficients()
