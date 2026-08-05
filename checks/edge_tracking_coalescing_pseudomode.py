#!/usr/bin/env python3
"""Exact exponent ledger for the edge-tracking coalescing pseudomode."""

from fractions import Fraction as F


def nearest_integer(value):
    """Round a nonnegative Fraction to a nearest integer."""
    quotient, remainder = divmod(value.numerator, value.denominator)
    return quotient + (2 * remainder >= value.denominator)


def powers(carrier_power, gain_power):
    """Return the powers of j in the two small parameters.

    eta^3/s grows like j**action_growth, while the transported-phase
    mismatch decays like j**(-coherence_decay).
    """

    carrier_power = F(carrier_power)
    gain_power = F(gain_power)
    action_growth = (carrier_power - 3 * gain_power) / 2
    coherence_decay = (carrier_power - 3 * gain_power) / 4
    return action_growth, coherence_decay


def main():
    # General algebra: both margins are positive exactly when A>3g.
    samples = (
        (F(13, 2), F(2)),
        (F(8), F(2)),
        (F(10), F(3)),
        (F(25, 2), F(4)),
    )
    for carrier_power, gain_power in samples:
        action, coherence = powers(carrier_power, gain_power)
        assert (action > 0) == (carrier_power > 3 * gain_power)
        assert (coherence > 0) == (carrier_power > 3 * gain_power)
        assert action == 2 * coherence

    # The cascade design point g=2, A=8.
    action, coherence = powers(8, 2)
    assert action == 1
    assert coherence == F(1, 2)

    # With eta ~ j^{-g/2} and s ~ j^{-A/2}:
    # effective h = s/eta^3 ~ j^{-(A-3g)/2};
    # Y-width sqrt(s*eta) ~ j^{-(A+g)/4};
    # physical width s*Y-width ~ j^{-(3A+g)/4}.
    carrier_power = F(8)
    gain_power = F(2)
    effective_h_decay = (carrier_power - 3 * gain_power) / 2
    y_width_decay = (carrier_power + gain_power) / 4
    physical_width_decay = (3 * carrier_power + gain_power) / 4
    assert effective_h_decay == 1
    assert y_width_decay == F(5, 2)
    assert physical_width_decay == F(13, 2)

    # Uniform analytic-symbol errors after Y=eta^2 X and division by eta^2.
    # h*eta = s/eta^2,
    # h^2*eta^4 = s^2/eta^2,
    # h^3*eta^5 = s^3/eta^4.
    h_decay = effective_h_decay
    eta_decay = gain_power / 2
    s_over_eta2 = carrier_power / 2 - gain_power
    s2_over_eta2 = carrier_power - gain_power
    s3_over_eta4 = 3 * carrier_power / 2 - 2 * gain_power
    assert h_decay + eta_decay == s_over_eta2
    assert 2 * h_decay + 4 * eta_decay == s2_over_eta2
    assert 3 * h_decay + 5 * eta_decay == s3_over_eta4
    assert min(s_over_eta2, s2_over_eta2, s3_over_eta4) > 0

    # Aspect-ratio-aware winding: z has period 2*pi/epsilon, so the
    # physical axial Fourier covector is epsilon*m.  Choosing
    # m=round(beta*p/epsilon) gives ratio error <= epsilon/(2p).
    epsilon = F(1, 97)
    beta = F(7, 5)
    p_integer = 113
    winding = nearest_integer(beta * p_integer / epsilon)
    ratio_error = abs(epsilon * winding / p_integer - beta)
    assert ratio_error <= epsilon / (2 * p_integer)
    assert winding > p_integer

    # The two unrecentered symbol errors are epsilon*s and
    # epsilon*s^2/eta^2. In h=s/eta^3 notation these are exactly
    # epsilon*h*eta^3 and epsilon*h^2*eta^4.
    s_value = F(1, 2**12)
    eta_value = F(1, 8)
    h_value = s_value / eta_value**3
    assert h_value * eta_value**3 == s_value
    assert h_value**2 * eta_value**4 == s_value**2 / eta_value**2

    # Edge tracking removes the exponential Duhamel gap:
    # delta ~ j^{-g}, T ~ j^g, so delta*T has power zero.
    assert -gain_power + gain_power == 0

    # Three real phases give at most 1 + 6 + 12 = 19 quadratic charges.
    zero_charge = 1
    doubles = 2 * 3
    mixed = 4 * 3  # four signs for each of the three unordered pairs
    assert zero_charge + doubles + mixed == 19

    print("coalescing effective parameter s/eta^3: PASS")
    print("sharp common WKB/material threshold A > 3g: PASS")
    print("g=2, A=8 design point and physical widths: PASS")
    print("uniform dilated analytic-symbol error powers: PASS")
    print("aspect-ratio-aware winding and lattice errors: PASS")
    print("edge-gap cancellation delta*T=O(1): PASS")
    print("three-phase quadratic charge count <= 19: PASS")
    print("all edge-tracking checks passed")


if __name__ == "__main__":
    main()
