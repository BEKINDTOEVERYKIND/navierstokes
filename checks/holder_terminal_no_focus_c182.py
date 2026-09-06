#!/usr/bin/env python3
"""Dependency-free checks for C182's finite-p conditional no-focus theorem.

The script checks exact Fourier pressure algebra and every C180/C161 power
ledger used by the note. The finite-p Mikhlin, transport-commutator,
Sobolev, and transport-diffusion estimates are analytic theorems, not
things a finite script proves.
"""

from __future__ import annotations

from fractions import Fraction as F


IV = tuple[int, int, int]


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def pressure_fourier_identity() -> None:
    """Check (1.2)--(1.3) on nontrivial incompressible Fourier pairs."""
    cases: tuple[tuple[IV, IV, IV, IV], ...] = (
        ((-2, -2, -2), (-2, 0, 2), (-2, -2, -1), (-2, 1, 2)),
        ((-2, -2, 0), (-2, 2, -2), (-2, -2, -2), (-2, 0, 2)),
        ((-2, -2, 1), (-2, 1, -2), (-2, -2, -2), (-2, 0, 2)),
        ((-2, -1, 0), (-1, 2, -2), (-2, -2, -2), (-2, 0, 2)),
    )
    for r, amplitude_u, s, amplitude_v in cases:
        assert dot(r, amplitude_u) == 0
        assert dot(s, amplitude_v) == 0
        k = add(r, s)
        assert dot(k, k) > 0

        direct_divergence = -(
            dot(k, amplitude_v) * dot(s, amplitude_u)
            + dot(k, amplitude_u) * dot(r, amplitude_v)
        )
        reduced_divergence = -2 * dot(s, amplitude_u) * dot(r, amplitude_v)
        assert direct_divergence == reduced_divergence
        assert reduced_divergence != 0

        # -Delta p equals the reduced divergence. The pressure divergence
        # Delta p then cancels the convection divergence exactly.
        pressure_coefficient = F(reduced_divergence, dot(k, k))
        pressure_divergence = -dot(k, k) * pressure_coefficient
        assert direct_divergence + pressure_divergence == 0

        # The pressure-gradient multiplier k_i k_j/|k|^2 has degree zero.
        for scale in (2, 7, 31):
            scaled_k = tuple(scale * x for x in k)
            for i in range(3):
                for j in range(3):
                    assert F(
                        scaled_k[i] * scaled_k[j], dot(scaled_k, scaled_k)
                    ) == F(k[i] * k[j], dot(k, k))


def c180_uniform_base_ledger() -> None:
    # Q belongs to [q,q+47], rho is O(1/(q sqrt(Q))), and hence
    # (Q rho)^2 is O(Q/q^2)=O(q^-1).
    previous = None
    for n in range(4, 80):
        q = n**8
        q_modes = q + 47
        squared_sum_without_fixed_constant = F(q_modes, q * q)
        assert squared_sum_without_fixed_constant <= F(2, q)
        if previous is not None:
            assert squared_sum_without_fixed_constant < previous
        previous = squared_sum_without_fixed_constant

    for heat in (F(0), F(1, 7), F(3, 4), F(1)):
        assert 0 <= heat <= 1
        assert heat * heat <= 1


def finite_p_power_ledger() -> None:
    # p=16, s=1/4 lies strictly above the 3/p point-evaluation threshold.
    p = 16
    s = F(1, 4)
    assert s > F(3, p)

    # L2/Linfinity interpolation: theta/2=1/16 gives theta=1/8.
    theta = F(1, 8)
    assert theta / 2 == F(1, 16)
    assert 1 - theta == F(7, 8)

    # In sixteenth powers:
    #   (b q^(7/8))^16 = b^16 q^14,
    #   (b q^(9/8))^16 = b^16 q^18,
    #   (b q^(3/2))^16 = b^16 q^24.
    for n in (4, 8, 16, 32):
        q = n**8
        b = F(1, n * n)
        j = n  # Eventually larger than every fixed multiple of log(n).

        interpolated_power = b**16 * q**14
        entrance_power = b**16 * q**18
        raw_target_power = b**16 * q**24
        taxed_target_power = raw_target_power / j**32

        assert entrance_power / interpolated_power == q**4
        assert entrance_power / raw_target_power == F(1, q**6)
        assert entrance_power / taxed_target_power == F(j**32, q**6)

        # With q=n^8 and the crude upper test J=n, the raw ratio is n^-3
        # and the J^2-taxed ratio is n^-1. Their sixteenth powers follow.
        assert F(1, q**6) == F(1, n**48)
        assert F(j**32, q**6) == F(1, n**16)

        raw_ratio = F(1, n**3)
        taxed_ratio = F(1, n)
        assert raw_ratio**16 == F(1, q**6)
        assert taxed_ratio**16 == F(j**32, q**6)

    # Normalized-Haar to unnormalized large-torus L16 conversion has a
    # K^(3/16) factor; its sixteenth power is K^3, not one.
    for parent_frequency in (2, 7, 31):
        assert parent_frequency**3 > 1

    # The rejected bounding-volume heuristic has L2 upper scale
    # b*J^2/sqrt(q), which is o(b) when J is polylogarithmic. Even J=n
    # gives the strict ratio n^-2 for q=n^8.
    for n in (4, 8, 16, 32):
        q = n**8
        j = n
        assert F(j * j, n**4) == F(1, n * n)
        assert F(j * j, n**4) < 1


def time_boundary() -> None:
    # Fixed terminal time gives a q-independent exponential constant.
    # A logarithmic interval instead yields a power of q.
    for q in (16, 81, 256, 625):
        fixed_time_factor = 7
        logarithmic_model_factor = q  # exp(log q)
        assert fixed_time_factor < logarithmic_model_factor


def main() -> None:
    pressure_fourier_identity()
    c180_uniform_base_ledger()
    finite_p_power_ledger()
    time_boundary()
    print("C182 finite-p conditional terminal no-focus checks passed")


if __name__ == "__main__":
    main()
