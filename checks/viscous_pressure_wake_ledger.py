#!/usr/bin/env python3
"""Exact scaling checks for the second viscous pressure-wake ledger.

The script checks monomial powers and finite geometric dominance.  It does
not verify Newton-kernel estimates or an all-order endpoint construction.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def add(*powers: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(entries) for entries in zip(*powers))


def scale(power: int, powers: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(power * entry for entry in powers)


def check_single_bubble_scaling() -> None:
    # Coordinates are powers of (nu, a, ell, K, d).
    viscosity = (1, 0, 0, 0, 0)
    amplitude_squared = (0, 2, 0, 0, 0)
    source_length = (0, 0, 1, 0, 0)
    kernel = (0, 0, 0, 0, -4)
    turnover_squared = (0, -2, 2, 0, 0)

    u2_far = add(viscosity, amplitude_squared, source_length, kernel)
    normalized = add(u2_far, turnover_squared)
    assert u2_far == (1, 2, 1, 0, -4)
    assert normalized == (1, 0, 3, 0, -4)

    at_adjacent_scale = add(normalized, (0, 0, -4, 0, 4))
    assert at_adjacent_scale == (1, 0, -1, 0, 0)
    print("single-bubble amplitude cancellation: exact")


def check_packed_scaling() -> None:
    # Per microbubble: nu * delta * ell^2 * d^-4, where
    # delta=ell/K.  Multiply by N=K^3.
    viscosity = (1, 0, 0, 0, 0)
    delta = (0, 0, 1, -1, 0)
    macro_time_length = (0, 0, 2, 0, 0)
    population = (0, 0, 0, 3, 0)
    kernel = (0, 0, 0, 0, -4)
    packed = add(viscosity, delta, macro_time_length, population, kernel)
    assert packed == (1, 0, 3, 2, -4)

    at_adjacent_scale = add(packed, (0, 0, -4, 0, 4))
    assert at_adjacent_scale == (1, 0, -1, 2, 0)
    print("packed K^3-bubble far-field loss K^2: exact")


def check_reynolds_identity() -> None:
    # a=ell^-gamma K^gamma and Re_carrier=a*ell/(nu*K).
    # K/Re_carrier = nu*ell^(gamma-1)*K^(2-gamma).
    for gamma in (Q(11, 10), Q(5, 4), Q(7, 5)):
        carrier_re = (-1, 1 - gamma, gamma - 1)  # (nu, ell, K)
        inverse = scale(-1, carrier_re)
        k_over_re = add(inverse, (0, 0, 1))
        assert k_over_re == (1, gamma - 1, 2 - gamma)
    print("K/Re_carrier identity: exact")


def check_geometric_dominance() -> None:
    # For ell_j=2^-j and K_j=(j+j0)^A, the last term controls the
    # preceding finite sum up to the pure geometric constant.  Since K is
    # increasing, term_j/term_(j+1) <= 2^-(1+m).
    for m in range(5):
        ratio = Q(1, 2 ** (1 + m))
        bound = 1 / (1 - ratio)
        for exponent in (1, 3, 7):
            for offset in (10, 100):
                terms = [
                    Q((j + offset) ** (2 * exponent))
                    * Q(2 ** ((1 + m) * j))
                    for j in range(12)
                ]
                assert sum(terms) <= bound * terms[-1]
    print("geometric-over-polynomial outer-wake dominance: exact")


def main() -> None:
    check_single_bubble_scaling()
    check_packed_scaling()
    check_reynolds_identity()
    check_geometric_dominance()
    print("all viscous pressure-wake scaling checks passed")


if __name__ == "__main__":
    main()
