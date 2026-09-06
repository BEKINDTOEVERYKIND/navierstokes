#!/usr/bin/env python3
"""Dependency-free checks for the Palasek--Gavrilov three-scale rectifier.

This is an algebra/scaling checker, not a Navier--Stokes proof.  It checks:
  * the exact P-A-Z heteroclinic and its dissipative energy identity;
  * endpoint transfer/depletion formulae;
  * the thin-torus intermittency and viscosity exponent ledger;
  * the reduced Leibovich--Stewartson sign on a cutoff's descending flank.
"""

import math
import random


def close(a: float, b: float, tol: float = 2e-11) -> None:
    if abs(a - b) > tol * (1.0 + abs(a) + abs(b)):
        raise AssertionError((a, b))


def check_rectifier(samples: int = 200) -> None:
    for _ in range(samples):
        d = 10.0 ** random.uniform(-1.0, 1.0)
        c = 10.0 ** random.uniform(-1.0, 2.0)
        p0 = 10.0 ** random.uniform(-1.0, 1.0)
        kappa = random.uniform(0.0, 0.95) * d * p0
        R = d * p0 - kappa
        t = random.uniform(-8.0, 8.0) / R

        r2 = d * d + c * c
        J = (c / d) * p0
        x = R * t
        th = math.tanh(x)
        sech = 1.0 / math.cosh(x)

        Q = -R * th
        A = R * sech / math.sqrt(r2)
        P = d * (Q + c * J + kappa) / r2
        Z = (d * d * J - c * Q - c * kappa) / r2

        Qdot = -R * R * sech * sech
        Adot = -R * A * th
        Pdot = d * Qdot / r2
        Zdot = -c * Qdot / r2

        close(Pdot, -d * A * A)
        close(Adot, (d * P - c * Z - kappa) * A)
        close(Zdot, c * A * A)
        close(2.0 * (P * Pdot + A * Adot + Z * Zdot), -2.0 * kappa * A * A)
        close(Q * Q + r2 * A * A, R * R)
        close(Z + (c / d) * P, J)

        p_plus = p0 - 2.0 * d * R / r2
        z_plus = 2.0 * c * R / r2
        close(p0 * p0 - p_plus * p_plus - z_plus * z_plus,
              4.0 * kappa * R / r2)

        rho = z_plus / p0
        e = kappa / (d * p0)
        theta = d / c
        close(rho, 2.0 * theta * (1.0 - e) / (1.0 + theta * theta))
        close(p_plus / p0, 1.0 - theta * rho)


def check_cutoff_discriminant() -> None:
    # For g(s)=exp[-1/(1-s^2)], h=s g'/g=-2s^2/(1-s^2)^2.
    # In the straight-tube limit the LS discriminant, after removing the
    # positive prefactor kappa^4 g^4, is
    # h[h(2+h)+(1+h)^2/2] = h(3h^2+6h+1)/2.
    for s in (0.80, 0.85, 0.90, 0.95):
        h = -2.0 * s * s / (1.0 - s * s) ** 2
        lhs = h * (h * (2.0 + h) + 0.5 * (1.0 + h) ** 2)
        rhs = 0.5 * h * (3.0 * h * h + 6.0 * h + 1.0)
        close(lhs, rhs)
        assert h < -2.0
        assert rhs < 0.0


def check_scale_ledger() -> None:
    alpha, beta, b = 2.4, 2.3, 1.1
    assert 2.0 * b < beta < alpha <= 2.5

    # r_N=N^-1, R_N=N^(4-2alpha):
    # R_N r_N^2=N^-2(alpha-1), delta=r_N/R_N=N^(2alpha-5).
    close((4.0 - 2.0 * alpha) - 2.0, -2.0 * (alpha - 1.0))
    assert 2.0 * alpha - 5.0 < 0.0

    # c_j~N^[b alpha-(b-1) beta] lies below N_child^alpha=N^(b alpha).
    close((b * alpha - (b - 1.0) * beta) - b * alpha,
          -(b - 1.0) * beta)

    # Carrier damping/growth ratio has exponent 2b-beta<0.
    assert 2.0 * b - beta < 0.0


if __name__ == "__main__":
    random.seed(20260801)
    check_rectifier()
    check_cutoff_discriminant()
    check_scale_ledger()
    print("PASS: rectifier orbit, energy ledger, cutoff sign, and scale exponents")
