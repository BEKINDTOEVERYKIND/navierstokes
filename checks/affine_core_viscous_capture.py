#!/usr/bin/env python3
"""Dependency-free checks for affine_core_viscous_capture.md."""

from __future__ import annotations

import cmath
import math
import random


def close(a, b, tol=2e-11):
    scale = max(1.0, abs(a), abs(b))
    if abs(a - b) > tol * scale:
        raise AssertionError((a, b, abs(a - b)))


def rhs(z, alpha, rho):
    U, V, p, IU = z
    return (
        ((1.0 + rho) * p - 1.0 - alpha) * U + p * V,
        -p * U + (1.0 - (1.0 + rho) * p - alpha) * V,
        V * V - U * U,
        U * U,
    )


def perturbed_rhs(z, alpha, rho, e):
    """Four-dimensional system without imposing the affine invariant."""
    U, V, p, q = z
    eU, eV, ep, eq = e
    return (
        (q + rho * p - 1.0 - alpha) * U + q * V + eU,
        -q * U + (1.0 - alpha - rho * p - q) * V + eV,
        V * V - U * U + ep,
        V * V - U * U + eq,
    )


def check_algebra():
    rng = random.Random(731904)
    for _ in range(100):
        alpha = rng.uniform(0.05, 0.98)
        rho = rng.uniform(0.0, 0.2)
        U = rng.uniform(-1.0, 1.0)
        V = rng.uniform(-1.0, 1.0)
        p = rng.uniform(0.0, 1.0)
        dU, dV, dp, _ = rhs((U, V, p, 0.0), alpha, rho)
        Jdot = 2 * U * dU + 2 * V * dV + (2 * (1 + rho) * p - 2) * dp
        close(Jdot, -2 * alpha * (U * U + V * V))

        # K=UV-p^2/2 obeys K'=-2 alpha UV, independently of rho.
        Kdot = dU * V + U * dV - p * dp
        close(Kdot, -2 * alpha * U * V)


def check_robust_algebra():
    rng = random.Random(619427)
    alpha = 0.9
    for _ in range(200):
        rho = rng.uniform(0.0, 0.01)
        U = rng.uniform(-0.5, 0.5)
        V = rng.uniform(-0.5, 0.5)
        p = rng.uniform(0.13, 0.22)
        qchild = rng.uniform(0.13, 0.22)
        e = tuple(rng.uniform(-0.02, 0.02) for _ in range(4))
        dU, dV, dp, dq = perturbed_rhs(
            (U, V, p, qchild), alpha, rho, e
        )

        # The affine defect a=q-p is driven only by the error mismatch.
        close(dq - dp, e[3] - e[2])

        # Terminal-box carrier dissipation before additive forcing.
        dU0, dV0, _, _ = perturbed_rhs(
            (U, V, p, qchild), alpha, rho, (0.0, 0.0, 0.0, 0.0)
        )
        half_R2_dot = U * dU0 + V * dV0
        if half_R2_dot > -0.03 * (U * U + V * V) + 2e-14:
            raise AssertionError((rho, U, V, p, qchild, half_R2_dot))

        # If eq=ep and q=p, the perturbed exact identity has only the
        # displayed forcing terms on its right-hand side.
        eU, eV, ep, _ = e
        dU, dV, dp, dq = perturbed_rhs(
            (U, V, p, p), alpha, rho, (eU, eV, ep, ep)
        )
        identity_dot = (
            2 * U * dU
            + 2 * V * dV
            + (2 * (1 + rho) * p - 0.2) * dp
            + 3.6 * U * U
        )
        forcing = 2 * U * eU + 2 * V * eV + (
            2 * (1 + rho) * p - 0.2
        ) * ep
        close(identity_dot, forcing)

        # Physical-to-normalized invariant-breaking error conversion.
        r = rng.uniform(10.0, 100.0)
        dscale = rng.uniform(0.4, 3.0)
        P0 = rng.uniform(0.2, 2.0)
        EP = rng.uniform(-1e-3, 1e-3)
        EZ = rng.uniform(-1e-3, 1e-3)
        ep_phys = -r * r * EP / (dscale * P0 * P0)
        eq_phys = r * EZ / (dscale * P0 * P0)
        expected = r * (EZ + r * EP) / (dscale * P0 * P0)
        close(eq_phys - ep_phys, expected)

        # Conservative child-to-parent ratio on the terminal box.
        ratio_times_r = qchild / (1.0 - rho * p)
        if not (0.12 <= ratio_times_r <= 0.23):
            raise AssertionError((rho, p, qchild, ratio_times_r))


def rk4_step(z, h, alpha, rho):
    k1 = rhs(z, alpha, rho)
    z2 = tuple(z[i] + 0.5 * h * k1[i] for i in range(4))
    k2 = rhs(z2, alpha, rho)
    z3 = tuple(z[i] + 0.5 * h * k2[i] for i in range(4))
    k3 = rhs(z3, alpha, rho)
    z4 = tuple(z[i] + h * k3[i] for i in range(4))
    k4 = rhs(z4, alpha, rho)
    return tuple(
        z[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6
        for i in range(4)
    )


def check_capture_numerically():
    alpha = 0.9
    delta = 1.0 - alpha
    rho = 0.01
    v = 1.0e-6
    p = v * v / (2 * delta)
    U = v**3 / (4 * delta * (2 - alpha))
    z = (U, v, p, 0.0)

    h = 0.01
    # The seed needs about log(1/v)/delta dimensionless time to emerge.
    for _ in range(60000):
        z = rk4_step(z, h, alpha, rho)

    U, V, p, IU = z
    q = (2 * delta / ((1 + rho) * (2 - 3 * delta))) ** 2
    lower = (2 * delta - 4 * alpha * q / (1 - q)) / (1 + rho)
    upper = 2 * delta / (1 + rho)
    if not (lower < p < upper):
        raise AssertionError((lower, p, upper))
    if abs(U) + abs(V) > 1e-8:
        raise AssertionError((U, V))

    identity = (1 + rho) * p * p - 2 * delta * p + 4 * alpha * IU
    close(identity, 0.0, 3e-7)

    D = 1 - 2 * (1 + rho) * p + rho * (2 + rho) * p * p
    if math.sqrt(max(0.0, D)) >= alpha:
        raise AssertionError((p, D, alpha))

    # Integrated imbalance estimate behind the zero-mode chirp defect.
    q = (2 * delta / ((1 + rho) * (2 - 3 * delta))) ** 2
    imbalance_constant = 4 * math.sqrt(q) / (1 - q)
    if imbalance_constant >= 0.48:
        raise AssertionError(imbalance_constant)


def check_uniform_endpoint_interval():
    alpha = 0.9
    delta = 1.0 - alpha
    endpoints = []
    for r in (10.0, 20.0, 100.0, 1.0e6):
        rho = r**-2
        q = (2 * delta / ((1 + rho) * (2 - 3 * delta))) ** 2
        lower = (2 * delta - 4 * alpha * q / (1 - q)) / (1 + rho)
        upper = 2 * delta / (1 + rho)
        assert lower > 0.1489
        assert upper < 0.2
        endpoints.append(upper)
    # The old 0.1981 endpoint applies near r=10 but is not uniform as
    # r tends to infinity.
    assert endpoints[1] > 0.1981 and endpoints[-1] > 0.19999


def check_chirped_pair_algebra():
    rng = random.Random(904117)
    for sigma in (-1.0, 1.0):
        for _ in range(50):
            lam = rng.uniform(10.0, 100.0)
            psi1 = rng.uniform(-3.0, 3.0)
            psi2 = rng.uniform(-4.0, 4.0)
            d = psi1 / (2 * lam)
            dp = psi2 / (2 * lam)

            # Coordinates are (r,t,h), with a_sigma=(1,-sigma*d,H).
            H = rng.uniform(-2.0, 2.0)
            a = (1.0, -sigma * d, H)
            gradphi = (sigma * psi1 / 2, lam, 0.0)
            close(sum(a[i] * gradphi[i] for i in range(3)), 0.0)

            # The self ordered pair has only -sigma*d' in the t direction.
            self_t = -sigma * dp
            close(self_t, -sigma * psi2 / (2 * lam))

    # Cross high ordered-pair sum: derivative terms cancel and the two
    # phase terms leave -i psi'^2/Lambda in the t direction.
    for _ in range(100):
        lam = rng.uniform(10.0, 100.0)
        psi1 = rng.uniform(-3.0, 3.0)
        psi2 = rng.uniform(-4.0, 4.0)
        d = psi1 / (2 * lam)
        dp = psi2 / (2 * lam)
        # (+,-)
        t_pm = dp - 1j * psi1 * d
        # (-,+)
        t_mp = -dp - 1j * psi1 * d
        close(t_pm + t_mp, -1j * psi1 * psi1 / lam)

        # Leading Leray radial coefficient i F'/k.
        C = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
        Fprime = -1j * C * 2 * psi1 * psi2 / lam
        leading = 1j * Fprime / (2 * lam)
        close(leading, C * psi1 * psi2 / (lam * lam))

    # The self projection necessarily contains psi'''.
    lam = 37.0
    sigma = -1.0
    A = 1.2 - 0.4j
    psi1, psi2, psi3 = 0.7, -1.1, 2.3
    Gprime = (
        -sigma
        * A**2
        / (2 * lam)
        * (psi3 + 1j * sigma * psi1 * psi2)
    )
    leading = 1j * Gprime / (2 * lam)
    expected = A**2 * (psi1 * psi2 - 1j * sigma * psi3) / (4 * lam**2)
    close(leading, expected)


def main():
    check_algebra()
    check_robust_algebra()
    check_capture_numerically()
    check_uniform_endpoint_interval()
    check_chirped_pair_algebra()
    print("PASS: affine-core viscous capture and corrected chirped residual")


if __name__ == "__main__":
    main()
