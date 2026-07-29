#!/usr/bin/env python3
"""Exact rational check of the lagged-intermittency parameter window."""

from fractions import Fraction


b = Fraction(11, 10)
alpha = 1 + Fraction(3, 2) / b
beta = Fraction(9, 4)

margins = {
    "viscosity_beta_minus_2b": beta - 2 * b,
    "profile_alpha_minus_beta": alpha - beta,
    # A second-order correction in p/N, where N=p^b.
    "second_order_carrier_margin": (
        2 * (b - 1) - (b - 1) * (beta - 1)
    ),
    # A second-order correction in q/p, where p=q^b.
    "second_order_envelope_margin": (
        2 * (1 - 1 / b) - (b - 1) * (beta - 1)
    ),
}

print(f"b={b} ({float(b):.12g})")
print(f"alpha={alpha} ({float(alpha):.12g})")
print(f"beta={beta} ({float(beta):.12g})")
for name, value in margins.items():
    print(f"{name}={value} ({float(value):.12g})")
assert all(value > 0 for value in margins.values())
