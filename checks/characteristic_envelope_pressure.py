#!/usr/bin/env python3
"""Exact arithmetic checks for the characteristic-envelope pressure ledger."""

import math
import random


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def mul(c, x):
    return tuple(c * a for a in x)


def norm(x):
    return math.sqrt(dot(x, x))


def proj(k, x):
    return add(x, mul(-dot(k, x) / dot(k, k), k))


def close(x, y, tol=2e-11):
    return norm(add(x, mul(-1.0, y))) < tol


def main():
    # Work in the orthonormal (e,v,t) coordinates.
    e = (1.0, 0.0, 0.0)
    v = (0.0, 1.0, 0.0)
    t = (0.0, 0.0, 1.0)

    for _ in range(200):
        H = random.uniform(-2.0, -0.2)
        W = math.sqrt(1.0 + H * H)
        Lam = random.uniform(10.0, 80.0)
        Q = random.uniform(0.1, 0.3 * Lam)
        delta = Q / (2.0 * Lam)
        a = add(mul(W, e), mul(-delta, t))
        b = add(mul(W, e), mul(delta, t))

        # Use a unit r in the (e,v)-plane with w.r=1.
        r = add(mul(1.0 / W, e), mul(abs(H) / W, v))
        assert abs(norm(r) - 1.0) < 1e-12
        k = add(mul(Lam, t), mul(Q / 2.0, r))
        ell = add(mul(Lam, t), mul(-Q / 2.0, r))
        K = add(k, ell)
        assert abs(dot(k, a)) < 1e-10
        assert abs(dot(ell, b)) < 1e-10

        phase = add(mul(dot(a, ell), b), mul(dot(b, k), a))
        claimed_phase = mul(-(Q * Q) / Lam, t)
        assert close(phase, claimed_phase)
        c = -Q * Q / (2.0 * Lam * Lam)
        assert close(phase, mul(c, K))
        assert abs(c + 2.0 * delta * delta) < 1e-13

        # Random first derivatives satisfying exact scalar-envelope
        # incompressibility. A_t and B_t determine A_e and B_e.
        A = random.uniform(0.2, 2.0)
        B = random.uniform(0.2, 2.0)
        At = random.uniform(-3.0, 3.0)
        Bt = random.uniform(-3.0, 3.0)
        Ae = delta * At / W
        Be = -delta * Bt / W
        Av = random.uniform(-2.0, 2.0)
        Bv = random.uniform(-2.0, 2.0)
        assert abs(W * Ae - delta * At) < 1e-12
        assert abs(W * Be + delta * Bt) < 1e-12

        agradB = W * Be - delta * Bt
        bgradA = W * Ae + delta * At
        F = add(mul(A * agradB, b), mul(B * bgradA, a))
        Pe = Ae * B + A * Be
        Pv = Av * B + A * Bv
        Pt = At * B + A * Bt
        F_claim = add(
            mul(2.0 * W * W * Pe, e),
            mul(-2.0 * delta * delta * Pt, t),
        )
        assert close(F, F_claim)

        gradP = (Pe, Pv, Pt)
        R = add(F, mul(-c, gradP))
        R_claim = add(
            mul(2.0 * (W * W + delta * delta) * Pe, e),
            mul(2.0 * delta * delta * Pv, v),
        )
        assert close(R, R_claim)
        assert abs(dot(R, K)) < 1e-10

        # The envelope-dark scalar is proportional to P_e, but says
        # nothing about the independent P_v pressure charge.
        dark = A * agradB + B * bgradA
        assert abs(dark - 2.0 * W * Pe) < 1e-10

        # Equal envelopes depending only on v obey both individual
        # divergence constraints and the old scalar dark equation, yet
        # retain the nonzero 2 delta^2 P_v v charge.
        fp = random.uniform(-2.0, 2.0)
        f = random.uniform(0.2, 2.0)
        equal_Pv = 2.0 * f * fp
        equal_R = mul(2.0 * delta * delta * equal_Pv, v)
        if abs(fp) > 1e-8:
            assert norm(equal_R) > 0.0

        # Characteristic coordinates are invariant along the appropriate
        # carrier velocities.
        alpha = delta / W
        grad_splus = add(t, mul(alpha, e))
        grad_sminus = add(t, mul(-alpha, e))
        assert abs(dot(a, grad_splus)) < 1e-12
        assert abs(dot(b, grad_sminus)) < 1e-12

        # Reparameterizing the longitudinal cutoff by its physical e-width
        # removes the apparent Lambda gain.
        L = random.uniform(0.3, 5.0)
        Re = L / alpha
        ratio_L = (W * W + delta * delta) / (
            2.0 * Lam * W * W * L
        )
        ratio_Re = (W * W + delta * delta) / (Q * W * Re)
        assert abs(ratio_L - ratio_Re) < 2e-12

        # Exact multiplier for a lateral v Fourier modulation.
        xi = random.uniform(-0.8 * Lam, 0.8 * Lam)
        if abs(xi) < 1e-8:
            xi = 0.1
        total = add(K, mul(xi, v))
        raw = mul(-(Q * Q) / Lam, t)
        exact = norm(proj(total, raw))
        formula = (
            (Q * Q / Lam)
            * abs(xi)
            / math.sqrt(4.0 * Lam * Lam + xi * xi)
        )
        assert abs(exact - formula) < 2e-11
        assert exact <= (
            (Q * Q / (2.0 * Lam * Lam)) * abs(xi) + 1e-12
        )

        # Equal longitudinal envelopes have a pressure-only raw high sum.
        # Their exact per-mode solenoidal correction is delta |xi|/|k+xi t|.
        xi_t = random.uniform(-0.4 * Lam, 0.4 * Lam)
        shifted_k = add(k, mul(xi_t, t))
        corrected_a = proj(shifted_k, a)
        correction = norm(add(corrected_a, mul(-1.0, a)))
        correction_formula = delta * abs(xi_t) / norm(shifted_k)
        assert abs(correction - correction_formula) < 2e-11
        assert correction <= 2.0 * delta * abs(xi_t) / Lam + 1e-12

        xi_e = random.uniform(-0.2 * Lam, 0.2 * Lam)
        xi_v = random.uniform(-0.2 * Lam, 0.2 * Lam)
        xi_t = random.uniform(-0.2 * Lam, 0.2 * Lam)
        xi_vec = (xi_e, xi_v, xi_t)
        assert norm(xi_vec) < 0.5 * Lam
        shifted_k = add(k, xi_vec)
        correction = norm(add(proj(shifted_k, a), mul(-1.0, a)))
        exact_formula = abs(W * xi_e - delta * xi_t) / norm(shifted_k)
        symbol_bound = 2.0 * (
            W * abs(xi_e) + delta * abs(xi_t)
        ) / Lam
        assert abs(correction - exact_formula) < 2e-11
        assert correction <= symbol_bound + 1e-12

        P_t = random.uniform(-3.0, 3.0)
        equal_F = mul(-delta * delta * P_t, t)
        equal_gradP = mul(P_t, t)
        equal_R = add(equal_F, mul(-c, equal_gradP))
        assert abs(norm(proj(K, equal_R))) < 1e-12

    print(
        "PASS: exact envelope-pressure criterion, residual, "
        "and core-buffer scales"
    )


if __name__ == "__main__":
    main()
