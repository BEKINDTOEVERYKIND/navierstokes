#!/usr/bin/env python3
"""Dependency-free algebra and exponent checks for localized_chirp_scale_ledger.md."""

import math
import random


def check_base_ledger():
    # Example used only for the bare residual window.
    b, beta, alpha = 1.05, 2.20, 2.35
    eps, Delta, mu = b - 1, alpha - beta, beta - 2 * b
    assert eps > 0 and Delta > 0 and mu > 0
    assert 2 * eps + mu + Delta < 0.5
    assert 5 - 2 * alpha > eps * Delta
    assert mu > 2 * eps * Delta
    assert abs(eps * Delta - 0.0075) < 1e-14
    assert abs(mu - 0.10) < 1e-14


def check_overlap_no_go():
    # If p>1/b and beta>2b, the best possible overlap margin is negative.
    for i in range(1, 1000):
        b = 1.0 + i / 1000.0
        upper = 3 * b - 2 * b * b - 1 / b
        assert upper < 1e-12

    # Random direct checks of the strict implication.
    rng = random.Random(815701)
    for _ in range(10000):
        b = rng.uniform(1.0001, 1.249)
        p = rng.uniform(1 / b + 1e-8, 0.9999)
        beta = rng.uniform(2 * b + 1e-8, 2.5)
        margin = b - p - beta * (b - 1)
        assert margin < 0


def check_characteristic_algebra():
    # In coordinates where a=(1,-delta), b=(1,+delta), eta=y+int delta ds
    # and xi=y-int delta ds, a.eta=0 and b.xi=0.
    rng = random.Random(277019)
    for _ in range(1000):
        delta = rng.uniform(-0.3, 0.3)
        grad_eta = (delta, 1.0)
        grad_xi = (-delta, 1.0)
        a = (1.0, -delta)
        b = (1.0, delta)
        assert abs(sum(x * y for x, y in zip(a, grad_eta))) < 1e-14
        assert abs(sum(x * y for x, y in zip(b, grad_xi))) < 1e-14

        # Exact-divergence envelopes F(eta),G(xi) leave the stated cross term.
        F, G = rng.uniform(0.2, 2), rng.uniform(0.2, 2)
        Fe, Gx = rng.uniform(-2, 2), rng.uniform(-2, 2)
        # a.grad G(xi)=-2 delta G_xi; b.grad F(eta)=2 delta F_eta.
        scalar_aG = -2 * delta * Gx
        scalar_bF = 2 * delta * Fe
        assert abs(F * scalar_aG + G * scalar_bF
                   - 2 * delta * (-F * Gx + G * Fe)) < 1e-14


def check_quadratic_chirp():
    rng = random.Random(409771)
    for _ in range(1000):
        R = rng.uniform(0.1, 2.0)
        chi = rng.uniform(1e-5, 0.1)
        a = chi / (R * R)
        s = rng.uniform(-R, R)
        exact = 2 * a * s * math.cos(a * s * s)
        affine = 2 * a * s
        if abs(affine) > 1e-12:
            rel = abs(exact / affine - 1)
            assert rel <= 0.51 * chi * chi


def main():
    check_base_ledger()
    check_overlap_no_go()
    check_characteristic_algebra()
    check_quadratic_chirp()
    print("PASS: localized-chirp ledger, tube no-go, and quadratic core")


if __name__ == "__main__":
    main()
