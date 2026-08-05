#!/usr/bin/env python3
"""Exact arithmetic checks for the embedded-quasimode threshold ledger.

This script verifies the rational inequalities in
research/2026-08-01-embedded-quasimode-instability-gate.md.  It does not
verify the Sobolev product estimate, construct an AO quasimode, or prove
nonlinear instability.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def sharper_threshold(beta: Fraction, q: int) -> Fraction:
    return (beta + q) / (beta + 2 * q)


def printed_fsv_threshold(beta: Fraction, q: int) -> Fraction:
    return (beta + 2 * q) / (beta + 3 * q)


def admissible_interval(
    beta: Fraction, q: int, chi: Fraction
) -> tuple[Fraction, Fraction]:
    assert chi > Q(1, 2)
    lower = beta * chi / (2 * chi - 1)
    upper = q * chi / (1 - chi)
    return lower, upper


def check_symbolic_equivalence() -> None:
    # Direct cross-multiplication of lower < upper in (5.2):
    # beta*chi/(2chi-1) < q*chi/(1-chi)
    # iff beta*(1-chi) < q*(2chi-1)
    # iff chi > (beta+q)/(beta+2q).
    beta = Q(5, 2)
    for q in range(1, 30):
        threshold = sharper_threshold(beta, q)
        for denominator in range(3, 80):
            for numerator in range(1, denominator):
                chi = Q(numerator, denominator)
                if chi <= Q(1, 2):
                    continue
                lower, upper = admissible_interval(beta, q, chi)
                assert (lower < upper) == (chi > threshold)
    print("frequency interval equivalence: exact")


def check_three_dimensional_table() -> None:
    beta = Q(5, 2)
    expected = {
        1: (Q(7, 9), Q(9, 11)),
        2: (Q(9, 13), Q(13, 17)),
        3: (Q(11, 17), Q(17, 23)),
        4: (Q(13, 21), Q(21, 29)),
    }
    for q, pair in expected.items():
        assert sharper_threshold(beta, q) == pair[0]
        assert printed_fsv_threshold(beta, q) == pair[1]
    print("three-dimensional threshold table: exact")


def check_four_fifths_example() -> None:
    beta = Q(5, 2)
    q = 1
    chi = Q(4, 5)
    s = Q(7, 2)
    eta = Q(15, 56)

    lower, upper = admissible_interval(beta, q, chi)
    assert lower < s < upper
    assert s * (1 / chi - 1) == Q(7, 8) < q
    assert eta < 1 - beta / s
    assert (1 + eta) * chi == Q(71, 70) > 1
    print("q=1, chi=4/5, s=7/2 witness: exact")


def check_lp_lever() -> None:
    # beta_p = 1 + 3/p in dimension three.  The p=infinity entry
    # is only the algebraic infimum over fixed finite p.
    expected = {
        2: (Q(5, 2), Q(7, 9)),
        3: (Q(2), Q(3, 4)),
        6: (Q(3, 2), Q(5, 7)),
    }
    for p, (beta, threshold) in expected.items():
        assert 1 + Q(3, p) == beta
        assert sharper_threshold(beta, 1) == threshold
    assert sharper_threshold(Q(1), 1) == Q(2, 3)
    print("three-dimensional Lp threshold lever: exact")


def check_required_order_inverse() -> None:
    beta = Q(5, 2)
    for chi in (Q(3, 5), Q(2, 3), Q(3, 4), Q(4, 5), Q(9, 10)):
        required = beta * (1 - chi) / (2 * chi - 1)
        for q in range(1, 20):
            assert (chi > sharper_threshold(beta, q)) == (q > required)
    print("inverse residual-order gate: exact")


def main() -> None:
    check_symbolic_equivalence()
    check_three_dimensional_table()
    check_four_fifths_example()
    check_lp_lever()
    check_required_order_inverse()
    print("all embedded-quasimode arithmetic checks passed")


if __name__ == "__main__":
    main()
