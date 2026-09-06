#!/usr/bin/env python3
"""Dependency-free exact arithmetic for C161.

This checks only the normalized Hilbert-space star, reality-symmetric
charge-translation cardinality, pure-normal gate darkness, and the
C146/C147 exponent ledger.  Fixed transfer-angle constants are stripped
from the exponent ledger.  It does not realize the star as a
Leray-projected Navier--Stokes block, identify normalized star time with a
physical gate amplitude, or turn a coefficient l1 scale into a physical
point-value lower bound.
"""

from fractions import Fraction as F
from math import isqrt


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def check_normalized_star():
    # q=n^8 is a perfect square.  The abstract calculation works for any
    # square q and therefore uses exact rational 1/sqrt(q) coefficients.
    for q_value in (4, 9, 16, 25, 64):
        square_root = isqrt(q_value)
        assert square_root * square_root == q_value
        coupling = F(1, square_root)

        # Bright vector B has q equal daughter coordinates 1/sqrt(q).
        bright_norm_squared = q_value * coupling * coupling
        assert bright_norm_squared == 1

        # Star generator: G e0=B, G B=-e0.  Its restriction to the
        # source/bright plane is the standard skew rotation generator.
        source_to_bright = q_value * coupling * coupling
        bright_to_source = -q_value * coupling * coupling
        assert source_to_bright == 1
        assert bright_to_source == -1

        # A representative daughter-dark vector e1-e2 is annihilated:
        # its daughter-coordinate sum is zero, so it has zero source
        # coupling.  More generally only the normalized sum couples.
        dark = (F(1), F(-1)) + (F(0),) * (q_value - 2)
        assert sum(dark) == 0
        assert coupling * sum(dark) == 0

        # Exact quarter-turn endpoint R=exp(pi G/2): e0 maps to B.  On a
        # daughter basis vector er, R er=er-B/sqrt(q)-e0/sqrt(q).
        # Check its norm and its orthogonality to R e0=B without invoking
        # transcendental arithmetic.
        daughter_part_norm_squared = (
            1 - F(2, q_value) + F(1, q_value)
        )
        source_part_norm_squared = F(1, q_value)
        assert daughter_part_norm_squared + source_part_norm_squared == 1
        daughter_sum = 1 - q_value * F(1, q_value)
        assert daughter_sum == 0

        # Collective normalization is essential: q couplings of size
        # theta/sqrt(q) have squared operator row norm theta^2.
        for theta in (F(1, 3), F(2, 7), F(5, 11)):
            assert q_value * (theta * coupling) ** 2 == theta**2


def check_charge_translation_cardinality():
    # A fixed-normal-charge q^2 source sheet translated by q pure-normal
    # shifts has q^3 distinct targets.  One shift serves every source.
    for q_value in (2, 4, 6, 8):
        source_charge = 4 * q_value
        sources = {
            (source_charge, horizontal_1, horizontal_2)
            for horizontal_1 in range(q_value)
            for horizontal_2 in range(q_value)
        }
        # Use q nonzero shifts so the daughter sheets are disjoint from
        # the source sheet as required by the abstract star.
        shifts = tuple(range(-q_value // 2, 0)) + tuple(
            range(1, q_value // 2 + 1)
        )
        assert set(shifts) == {-shift for shift in shifts}
        assert 0 not in shifts
        targets = {
            (charge + shift, horizontal_1, horizontal_2)
            for charge, horizontal_1, horizontal_2 in sources
            for shift in shifts
        }
        assert len(sources) == q_value**2
        assert len(shifts) == q_value
        assert len(targets) == q_value**3
        assert targets.isdisjoint(sources)
        negative_sources = {
            tuple(-coordinate for coordinate in source)
            for source in sources
        }
        negative_targets = {
            tuple(-coordinate for coordinate in target)
            for target in targets
        }
        assert targets.isdisjoint(negative_sources)
        assert targets.isdisjoint(negative_targets)
        assert sources.isdisjoint(negative_sources)

        # Every (source,shift) has a unique target label in this product
        # chart, so the support map has no first-step collisions.
        labels = {
            (charge + shift, horizontal_1, horizontal_2): (
                charge, horizontal_1, horizontal_2, shift
            )
            for charge, horizontal_1, horizontal_2 in sources
            for shift in shifts
        }
        assert len(labels) == q_value**3


def check_pure_normal_gate_darkness():
    # All gate wavevectors are parallel to N and all polarizations to a
    # fixed E perpendicular to N.  Hence every ordered gate-gate Euler
    # coefficient contains E dot (rN)=0.
    normal = (1, 1, 1)
    polarization = (1, -1, 0)
    assert dot(normal, polarization) == 0
    for left_charge in (-7, -2, 1, 5):
        for right_charge in (-5, -1, 3, 11):
            left_wavevector = tuple(left_charge * entry for entry in normal)
            right_wavevector = tuple(right_charge * entry for entry in normal)
            assert dot(polarization, left_wavevector) == 0
            assert dot(polarization, right_wavevector) == 0


def check_schedule_and_microgate_ledger():
    # A real q-shift gate chart uses a symmetric nonzero shift set, so take
    # the cofinal even-n subsequence.  This changes no power exponent.
    for n in (2, 4, 6):
        q_value = n**8
        b = F(1, n**2)
        pulses = n**2
        assert pulses == F(1, b)
        source_modes = q_value**2
        target_modes = q_value**3

        seed_l2 = F(1, n**28)
        gain = n**26
        final_l2 = gain * seed_l2
        assert final_l2 == b

        # Equal q^2-mode long-gain packet.
        seed_coefficient = seed_l2 / q_value
        source_coefficient = final_l2 / q_value
        source_l1 = source_modes * source_coefficient
        assert seed_coefficient == F(1, n**36)
        assert source_coefficient == F(1, n**10)
        assert source_l1 == n**6

        # Normalized q-way star endpoint.
        square_root_q = n**4
        assert square_root_q**2 == q_value
        daughter_coefficient = source_coefficient / square_root_q
        target_l2 = daughter_coefficient * n**12
        target_l1 = target_modes * daughter_coefficient
        assert daughter_coefficient == F(1, n**14)
        assert target_l2 == b
        assert target_l1 == n**10

        # J=b^-1 normalized-time factors, theta=b, so J*theta=1.
        # The exact factors are exp((pi/2J)G); the rational calculation
        # below strips the fixed pi/2 and checks only power exponents.
        theta = b
        assert pulses * theta == 1
        linearized_ledger_daughter = (
            theta * source_coefficient / square_root_q
        )
        assert linearized_ledger_daughter == F(1, n**16)
        assert pulses * linearized_ledger_daughter == daughter_coefficient
        one_factor_target_scale = n**12 * linearized_ledger_daughter
        assert one_factor_target_scale == b * theta == b**2
        assert pulses * one_factor_target_scale == b

        # C146 conditional residual powers at the collective star norm.
        per_pulse_wake = b * theta**2
        per_pulse_return = b * theta**3
        assert pulses * per_pulse_wake == b**2
        assert pulses * per_pulse_return == b**3

        # q sequential independent daughters do not fit J pulses, and
        # using unnormalized size theta per daughter overfills energy.
        assert q_value > pulses
        assert F(q_value, pulses) == n**6
        independent_energy_multiplier = q_value * theta**2
        assert independent_energy_multiplier == n**4 > 1


def main():
    check_normalized_star()
    check_charge_translation_cardinality()
    check_pure_normal_gate_darkness()
    check_schedule_and_microgate_ledger()
    print("PASS C161: normalized q-way star is an exact energy-preserving bright rotation")
    print("PASS C161: q shifts map q^2 sources to q^3 unique first-step targets")
    print("PASS C161: J=b^-1 normalized factors meet b, b^2, b^3 exponents")
    print("Coefficient l1 scales are not physical point-value lower bounds")
    print(
        "Uniform Leray star/coherence realization, depletion, collisions, "
        "collars, pressure, and BAFL remain OPEN"
    )


if __name__ == "__main__":
    main()
