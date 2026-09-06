#!/usr/bin/env python3
"""Exact exponent checks for the AO long-gain quasimode ledger.

This checks only the power inequalities in
research/2026-08-02-ao-long-gain-quasimode-ledger.md.  It does not prove a
uniform semigroup bound, construct a thin-ring quasimode, or identify the
full Euler growth edge.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def powers(A: Fraction, q: Fraction, r: Fraction, g: int, d: int):
    """Return the edge-action and polynomial relative-error exponents."""
    edge_action = Q(g) - A * r
    relative_error = Q(g * (d + 1)) - A * q
    return edge_action, relative_error


def admissible(A: Fraction, q: Fraction, r: Fraction, g: int, d: int):
    edge, error = powers(A, q, r, g, d)
    # Equality in the edge-action exponent leaves only a bounded factor
    # exp(C), so it is admissible.  Decay must still come from the
    # polynomial relative-error exponent.
    return edge <= 0 and error < 0


def summable(A: Fraction, q: Fraction, r: Fraction, g: int, d: int):
    edge, error = powers(A, q, r, g, d)
    return edge <= 0 and error < -1


def check_ao_design_point() -> None:
    q = Q(1)
    r = Q(1, 2)
    g = 2
    d = 0

    assert admissible(Q(4), q, r, g, d)
    assert summable(Q(4), q, r, g, d)
    assert admissible(Q(9, 2), q, r, g, d)
    assert summable(Q(9, 2), q, r, g, d)

    edge, error = powers(Q(5), q, r, g, d)
    assert edge == Q(-1, 2)
    assert error == Q(-3)
    print("AO q=1, r=1/2, G~j^2 long-gain gate: A>=4")


def check_polynomial_prefactors() -> None:
    q = Q(1)
    r = Q(1, 2)
    g = 2

    # d=1 needs A>4 for convergence but A>5 for summability.
    assert admissible(Q(9, 2), q, r, g, 1)
    assert not summable(Q(5), q, r, g, 1)
    assert summable(Q(11, 2), q, r, g, 1)

    # d=2 moves the polynomial-error condition above the edge condition.
    assert not admissible(Q(6), q, r, g, 2)
    assert admissible(Q(13, 2), q, r, g, 2)
    assert summable(Q(15, 2), q, r, g, 2)
    print("polynomial semigroup-prefactor exponents: exact")


def check_generic_equivalence() -> None:
    for g in range(1, 6):
        for d in range(5):
            for q_num in range(1, 6):
                for r_num in range(1, 6):
                    q = Q(q_num, 3)
                    r = Q(r_num, 4)
                    for A_num in range(1, 50):
                        A = Q(A_num, 4)
                        expected = A * r >= g and A * q > g * (d + 1)
                        assert admissible(A, q, r, g, d) == expected
    print("generic finite-action gate equivalence: exact")


if __name__ == "__main__":
    check_ao_design_point()
    check_polynomial_prefactors()
    check_generic_equivalence()
    print("all AO long-gain exponent checks passed")
