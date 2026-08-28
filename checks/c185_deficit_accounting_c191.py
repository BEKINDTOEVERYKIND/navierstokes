#!/usr/bin/env python3
"""Dependency-free exact arithmetic for C191.

This checker verifies the C188 power ledger, the exact logarithmic action
required by the C185 operator-norm floor, the positive scalar
clock/collar/heat feasibility arithmetic, and the fixed-energy
normalization cancellation.  It does not prove the C185 spectral
inclusion or any finite-frequency, viscous, concentration, or nonlinear
bridge.
"""

from decimal import Decimal, getcontext
from fractions import Fraction as F


getcontext().prec = 60


def c188_stage_powers() -> None:
    # All entries are powers of n on the fixed q=n^8 C188 schedule.
    q = F(8)
    b = F(-5, 2)
    focus = F(3, 2) * q
    gain = b + focus
    deficit = F(3, 2) - F(9, 8)
    amplifier = deficit * q

    assert focus == 12
    assert gain == F(19, 2)
    assert deficit == F(3, 8)
    assert amplifier == 3

    preparation = F(51, 2)
    seed = F(-36)
    source_coefficient = preparation + seed
    assert source_coefficient == F(-21, 2)
    assert q + source_coefficient == b

    daughter_coefficient = source_coefficient - q / 2
    assert daughter_coefficient == F(-29, 2)
    assert 3 * q / 2 + daughter_coefficient == b
    assert 3 * q + daughter_coefficient == gain

    total_required = preparation + amplifier
    assert total_required == F(57, 2)


def certified_orbit_constants() -> None:
    t_lo = Decimal("3.0361377700939")
    t_hi = Decimal("3.0361377700941")
    beta_lo = Decimal("2.1108753438878")
    beta_hi = Decimal("2.1108753438881")

    assert t_lo > Decimal(3)
    assert t_hi < Decimal(76) / Decimal(25)
    assert beta_lo > Decimal(21) / Decimal(10)
    assert beta_hi > beta_lo

    # sqrt(3)>5/3 because 3>25/9.  Hence m=sqrt(3) beta>7/2.
    assert F(3) > F(25, 9)
    assert F(5, 3) * F(21, 10) == F(7, 2)


def required_window_thresholds() -> None:
    q_power = F(8)
    missing_q_power = F(3, 8)
    missing_n_power = q_power * missing_q_power
    terminal_periods_per_log_n = 5 * missing_n_power
    terminal_periods_per_log_q = terminal_periods_per_log_n / q_power

    assert missing_n_power == 3
    assert terminal_periods_per_log_n == 15
    assert terminal_periods_per_log_q == F(15, 8)

    # T>3 and T<76/25 give the exact lower/upper action coefficients.
    terminal_action_lower_per_log_q = 3 * terminal_periods_per_log_q
    terminal_action_upper_per_log_q = F(76, 25) * terminal_periods_per_log_q
    assert terminal_action_lower_per_log_q == F(45, 8)
    assert terminal_action_upper_per_log_q == F(57, 10)

    # If the same floor is assigned to H=n^(51/2) and q^(3/8)=n^3.
    total_n_power = F(51, 2) + 3
    total_periods_per_log_n = 5 * total_n_power
    total_periods_per_log_q = total_periods_per_log_n / q_power
    total_action_upper_per_log_q = F(76, 25) * total_periods_per_log_q

    assert total_n_power == F(57, 2)
    assert total_periods_per_log_n == F(285, 2)
    assert total_periods_per_log_q == F(285, 16)
    assert total_action_upper_per_log_q == F(1083, 20)


def scalar_feasibility() -> None:
    # Exact dyadic substitute for every n>=2:
    # (6/5)^116 > 2^29, hence for L=ceil(log_2 n),
    # (6/5)^(116L)>2^(29L)>=n^29>n^(57/2).
    assert F(6, 5) ** 4 > 2
    assert F(6, 5) ** 116 > 2**29
    assert F(29) > F(57, 2)

    total_periods_per_log_n = F(285, 2)

    # C188's log n <=6 n^(1/14) plus the ceiling:
    # R_* <= (6*285/2+1)n^(1/14)=856 n^(1/14).
    assert 6 * total_periods_per_log_n + 1 == 856
    collar_n_power = F(7, 2) * F(1, 14) - F(1, 2)
    assert collar_n_power == F(-1, 4)

    # T<76/25 and heat exponent 2 mu T R_*.
    assert 2 * F(76, 25) == F(152, 25)

    # If x=2 mu T R_*<=1/100, then 1/(1-x)<=100/99.
    assert 1 / (1 - F(1, 100)) == F(100, 99)

    # Actual C188 physical-time unit and C176 viscous-collar exponent.
    assert F(-35, 2) < 0
    viscous_collar_n_power = F(3, 2) * F(1, 14) + 8 + F(15, 2)
    assert viscous_collar_n_power == F(437, 28)


def normalization_fork() -> None:
    q = F(8)
    b = F(-5, 2)
    amplifier = F(3)  # q^(3/8)=n^3

    # Extra norm gain multiplies child energy by q^(3/4)=n^6.
    energy_overrun = 2 * amplifier
    scheduled_energy = 2 * b
    raw_energy = scheduled_energy + energy_overrun
    assert energy_overrun == 6
    assert scheduled_energy == -5
    assert raw_energy == 1
    assert raw_energy - scheduled_energy == 6

    # Per-mode seed arithmetic.
    preparation = F(51, 2)
    seed = F(-36)
    original_output = preparation + seed
    raw_output = amplifier + preparation + seed
    rescaled_seed = seed - amplifier
    normalized_output = amplifier + preparation + rescaled_seed
    assert original_output == F(-21, 2)
    assert raw_output == F(-15, 2)
    assert rescaled_seed == -39
    assert normalized_output == original_output

    # Packet seed and conditional C147 same-support writer after uniform
    # rescaling.  The point seed also loses A=n^3.
    packet_seed = q + seed
    packet_seed_rescaled = q + rescaled_seed
    point_seed_rescaled = F(-16) - amplifier
    full_log_gain = F(51, 2) + amplifier
    writer_l2 = point_seed_rescaled + q - 3 * q / 2
    writer_relative = writer_l2 - packet_seed_rescaled
    assert packet_seed == -28
    assert packet_seed_rescaled == -31
    assert point_seed_rescaled == -19
    assert full_log_gain == F(57, 2)
    assert writer_l2 == -23
    assert writer_relative == 8

    # The formal point exponent returns to 9/8 after final L2 normalization.
    assert F(9, 8) + F(3, 8) - F(3, 8) == F(9, 8)
    assert F(9, 8) + F(3, 8) == F(3, 2)


def operator_norm_nonimplication() -> None:
    # For G=A I on R^2, the squared concentration quotient
    # ||v||_infty^2/||v||_2^2 is unchanged although ||G||_2=A.
    vector = (F(3), F(-4))
    base_l2_sq = sum(entry * entry for entry in vector)
    base_linf_sq = max(abs(entry) for entry in vector) ** 2
    assert base_l2_sq == 25
    assert base_linf_sq == 16

    for amplitude in (F(2), F(17, 3), F(1000)):
        image = tuple(amplitude * entry for entry in vector)
        image_l2_sq = sum(entry * entry for entry in image)
        image_linf_sq = max(abs(entry) for entry in image) ** 2
        assert image_l2_sq == amplitude**2 * base_l2_sq
        assert image_linf_sq == amplitude**2 * base_linf_sq
        assert image_linf_sq * base_l2_sq == base_linf_sq * image_l2_sq

    # A retained fraction eta would add 5 log(eta^-1) full-return units.
    assert 5 * F(1, 5) == 1


def passive_two_d_three_c_identity() -> None:
    # With n=N/sqrt(3), (-sqrt(6) f)n=-sqrt(2) f N.
    assert F(6, 3) == 2

    # For N=(1,1,1), (N x grad f).grad f is the polynomial
    # (c-b)a+(a-c)b+(b-a)c.  Collect its quadratic monomials exactly.
    coefficients = {
        "ac": F(1) - F(1),
        "ab": -F(1) + F(1),
        "bc": -F(1) + F(1),
    }
    assert all(value == 0 for value in coefficients.values())


def main() -> None:
    c188_stage_powers()
    certified_orbit_constants()
    required_window_thresholds()
    scalar_feasibility()
    normalization_fork()
    operator_norm_nonimplication()
    passive_two_d_three_c_identity()
    print("C191 C185 deficit accounting: PASS")
    print("POSITIVE: certified floor has explicit logarithmic scalar capacity")
    print("FAILURES: no lower coverage/band bridge; common scalar gain cancels")
    print("CLASS SCOPE: universal secular lock is contradicted by C159/C185")


if __name__ == "__main__":
    main()
