#!/usr/bin/env python3
"""Exact checks for C160: two-radius detuning of limiting DACR.

All symbolic reconstruction uses rational arithmetic in
Q(sqrt(2),sqrt(3)) and finite trigonometric polynomials.  Floating point
is used only for the final schedule-limit sanity check.
"""

from fractions import Fraction as F
from math import log, sqrt

import floquet_averaged_dacr_c155_c156 as base


def rational_circle(t):
    t = F(t)
    denominator = 1 + t * t
    return (
        (1 - t * t) / denominator,
        2 * t / denominator,
    )


SINE_TWO = base.trig_scale(
    base.trig_multiply(base.TRIG_SIN, base.TRIG_COS),
    base.q(2),
)
COSINE_TWO = base.trig_sub(
    base.trig_multiply(base.TRIG_COS, base.TRIG_COS),
    base.trig_multiply(base.TRIG_SIN, base.TRIG_SIN),
)


def second_harmonic_components(value):
    """Return g0,gc,gs for g0+gc*cos(2s)+gs*sin(2s)."""
    constant = base.trig_mean(value)
    cosine = base.qscale(
        base.trig_mean(base.trig_multiply(value, COSINE_TWO)), 2
    )
    sine = base.qscale(
        base.trig_mean(base.trig_multiply(value, SINE_TWO)), 2
    )
    reconstructed = base.trig_add(
        base.trig_constant(constant),
        base.trig_add(
            base.trig_scale(COSINE_TWO, cosine),
            base.trig_scale(SINE_TWO, sine),
        ),
    )
    assert reconstructed == value
    return constant, cosine, sine


def trig_kelvin_apply(wavevector, value):
    keys = set()
    for component in value:
        keys.update(component)
    output = [dict(), dict(), dict()]
    for key in keys:
        vector = tuple(
            component.get(key, base.ZERO) for component in value
        )
        image = base.kelvin_apply(wavevector, vector)
        for index, coefficient in enumerate(image):
            if coefficient != base.ZERO:
                output[index][key] = coefficient
    return tuple(base.trig_normalize(component) for component in output)


def trig_value_at_zero(value):
    """Evaluate a normalized trig polynomial at s=0 exactly."""
    output = base.ZERO
    for (sine_power, _cosine_power), coefficient in value.items():
        if sine_power == 0:
            output = base.qadd(output, coefficient)
    return output


def assert_transverse_trig(wavevector, value):
    assert base.trig_dot_constant(value, wavevector) == base.TRIG_ZERO


def assert_kelvin_plane_identity(wavevector, vector, omega_squared):
    assert base.vdot(wavevector, vector) == base.ZERO
    twice = base.kelvin_apply(
        wavevector, base.kelvin_apply(wavevector, vector)
    )
    assert twice == base.vscale(base.q(-omega_squared), vector)


def periodic_wake(wavevector, forcing):
    """Unique transverse periodic wake for 0<omega<2, omega!=1."""
    assert_transverse_trig(wavevector, forcing)
    components = [
        second_harmonic_components(component) for component in forcing
    ]
    forcing_zero = tuple(component[0] for component in components)
    forcing_cosine = tuple(component[1] for component in components)
    forcing_sine = tuple(component[2] for component in components)

    norm_squared = base.vdot(wavevector, wavevector)
    vertical_squared = base.qmul(wavevector[2], wavevector[2])
    assert norm_squared[1:] == (F(0), F(0), F(0))
    assert vertical_squared[1:] == (F(0), F(0), F(0))
    omega_squared = 4 * vertical_squared[0] / norm_squared[0]
    assert 0 < omega_squared < 4
    # In C160 the interior-angle identities sharpen this to
    # omega_- in (0,1), omega_+ in (1,2).  Excluding omega=1 is what
    # makes the 2*pi-periodic solution unique, not merely omega<2.
    assert omega_squared != 1

    for component in (forcing_zero, forcing_cosine, forcing_sine):
        assert_kelvin_plane_identity(
            wavevector, component, omega_squared
        )

    wake_zero = tuple(
        base.qscale(coefficient, F(1, omega_squared))
        for coefficient in base.kelvin_apply(
            wavevector, forcing_zero
        )
    )
    kelvin_cosine = base.kelvin_apply(
        wavevector, forcing_cosine
    )
    wake_cosine = tuple(
        base.qscale(
            base.qadd(kelvin_cosine[index],
                      base.qscale(forcing_sine[index], 2)),
            F(1, omega_squared - 4),
        )
        for index in range(3)
    )
    kelvin_wake_cosine = base.kelvin_apply(
        wavevector, wake_cosine
    )
    wake_sine = tuple(
        base.qscale(
            base.qadd(kelvin_wake_cosine[index],
                      forcing_cosine[index]),
            F(1, 2),
        )
        for index in range(3)
    )

    wake = tuple(
        base.trig_add(
            base.trig_constant(wake_zero[index]),
            base.trig_add(
                base.trig_scale(COSINE_TWO, wake_cosine[index]),
                base.trig_scale(SINE_TWO, wake_sine[index]),
            ),
        )
        for index in range(3)
    )
    derivative = tuple(
        base.trig_derivative(component) for component in wake
    )
    right_hand_side = tuple(
        base.trig_add(component, forcing[index])
        for index, component in enumerate(
            trig_kelvin_apply(wavevector, wake)
        )
    )
    assert derivative == right_hand_side
    assert_transverse_trig(wavevector, wake)

    wake_at_zero = tuple(trig_value_at_zero(component) for component in wake)
    assert_kelvin_plane_identity(
        wavevector, wake_at_zero, omega_squared
    )
    kelvin_at_zero = base.kelvin_apply(wavevector, wake_at_zero)
    # C_k is skew-adjoint on k^perp, hence has norm omega_k there.
    assert base.vdot(kelvin_at_zero, kelvin_at_zero) == base.qscale(
        base.vdot(wake_at_zero, wake_at_zero), omega_squared
    )
    return wake, omega_squared


def scaled_mode(radius, t):
    cosine, sine = rational_circle(t)
    wavevector = base.vscale(
        base.q(radius),
        (base.q(cosine), base.q(sine), base.INV_SQRT3),
    )
    parent = base.rotating_parent(cosine, sine)
    assert_transverse_trig(wavevector, parent)
    return wavevector, parent


def directed_return_rows(output_wavevector, partner_wavevector,
                         target_parent, partner_parent):
    """Rows pairing a general output wake back into the target mode."""
    rows = []
    allowed_harmonics = {(0, 0), (1, 1), (2, 0)}
    for index in range(3):
        basis_vector = tuple(
            base.trig_constant(base.ONE if index == j else base.ZERO)
            for j in range(3)
        )
        returned = base.trig_vector_neg(
            base.trig_symmetric_symbol(
                output_wavevector,
                basis_vector,
                partner_wavevector,
                partner_parent,
            )
        )
        row = base.trig_dot(target_parent, returned)
        # This harmonic support is load-bearing for the Cesaro and
        # small-divisor statements in the note.
        assert set(row).issubset(allowed_harmonics)
        rows.append(row)
    return tuple(rows)


def detuned_means(radius_a, radius_b, t):
    """Return exact sum/difference means and inertial frequencies."""
    radius_a = F(radius_a)
    radius_b = F(radius_b)
    t = F(t)
    assert radius_a > 0 and radius_b > 0 and radius_a != radius_b
    assert t != 0

    target_wavevector, target_parent = scaled_mode(radius_a, F(0))
    partner_wavevector, partner_parent = scaled_mode(radius_b, t)
    sum_wavevector = base.vadd(target_wavevector, partner_wavevector)
    difference_wavevector = base.vadd(
        target_wavevector, base.vneg(partner_wavevector)
    )

    directed_return_rows(
        sum_wavevector,
        base.vneg(partner_wavevector),
        target_parent,
        partner_parent,
    )
    directed_return_rows(
        difference_wavevector,
        partner_wavevector,
        target_parent,
        partner_parent,
    )

    sum_forcing = base.trig_symmetric_symbol(
        target_wavevector,
        target_parent,
        partner_wavevector,
        partner_parent,
    )
    difference_forcing = base.trig_symmetric_symbol(
        target_wavevector,
        target_parent,
        base.vneg(partner_wavevector),
        partner_parent,
    )
    sum_wake, omega_plus_squared = periodic_wake(
        sum_wavevector, sum_forcing
    )
    difference_wake, omega_minus_squared = periodic_wake(
        difference_wavevector, difference_forcing
    )

    sum_return = base.trig_vector_neg(
        base.trig_symmetric_symbol(
            sum_wavevector,
            sum_wake,
            base.vneg(partner_wavevector),
            partner_parent,
        )
    )
    difference_return = base.trig_vector_neg(
        base.trig_symmetric_symbol(
            difference_wavevector,
            difference_wake,
            partner_wavevector,
            partner_parent,
        )
    )
    sum_mean = base.trig_mean(
        base.trig_dot(target_parent, sum_return)
    )
    difference_mean = base.trig_mean(
        base.trig_dot(target_parent, difference_return)
    )
    assert sum_mean[1:] == (F(0), F(0), F(0))
    assert difference_mean[1:] == (F(0), F(0), F(0))
    return (
        sum_mean[0],
        difference_mean[0],
        omega_plus_squared,
        omega_minus_squared,
    )


def closed_mean(radius_a, radius_b, t):
    radius_a = F(radius_a)
    radius_b = F(radius_b)
    t = F(t)
    return (
        -radius_a
        * radius_b
        * t**3
        * (t * t + 2)
        * (t * t - t + 2)
        * (t * t + t + 2)
        / (2 * (1 + t * t) ** 5)
    )


def closed_frequencies(radius_a, radius_b, t):
    radius_a = F(radius_a)
    radius_b = F(radius_b)
    t = F(t)
    common = (1 + t * t) * (radius_a**2 + radius_b**2)
    denominator_plus = common + radius_a * radius_b * (2 - t * t)
    denominator_minus = common + radius_a * radius_b * (t * t - 2)
    return (
        (radius_a + radius_b) ** 2 * (1 + t * t)
        / denominator_plus,
        (radius_a - radius_b) ** 2 * (1 + t * t)
        / denominator_minus,
    )


def check_detuned_zero_mean():
    for radius_a, radius_b, t in (
        (F(2), F(1), F(1)),
        (F(3, 2), F(1), F(1, 2)),
        (F(11, 10), F(1), F(-3)),
        (F(101, 100), F(1), F(-1)),
        (F(2), F(3), F(1, 3)),
        (F(7, 5), F(9, 4), F(-2, 3)),
        (F(1001, 1000), F(1), F(7, 4)),
        (F(5), F(1, 4), F(-5)),
    ):
        sum_mean, difference_mean, omega_plus, omega_minus = (
            detuned_means(radius_a, radius_b, t)
        )
        expected = closed_mean(radius_a, radius_b, t)
        expected_plus, expected_minus = closed_frequencies(
            radius_a, radius_b, t
        )
        assert sum_mean == expected
        assert difference_mean == -expected
        assert sum_mean + difference_mean == 0
        assert omega_plus == expected_plus
        assert omega_minus == expected_minus
        common = (1 + t * t) * (radius_a**2 + radius_b**2)
        denominator_plus = common + radius_a * radius_b * (2 - t * t)
        denominator_minus = common + radius_a * radius_b * (t * t - 2)

        # Exact nonresonance factorization from (1.9)--(1.10).
        assert (omega_plus - 1) * denominator_plus == (
            3 * radius_a * radius_b * t * t
        )
        assert (1 - omega_minus) * denominator_minus == (
            3 * radius_a * radius_b * t * t
        )
        assert (4 - omega_plus) * denominator_plus == 3 * (
            (radius_a + radius_b) ** 2
            + t * t * (radius_a - radius_b) ** 2
        )
        assert 0 < omega_minus < 1 < omega_plus < 4

    # A!=B alone does not imply uniqueness of the periodic wake: at zero
    # angular separation both inertial frequencies are exactly one.
    assert closed_frequencies(F(2), F(1), F(0)) == (F(1), F(1))


def check_adjacent_layer_schedule():
    # At quarter separation and A=1+q^-1,B=1,
    # omega_-^2=2(A-B)^2/(2A^2-AB+2B^2), hence omega_-~q^-1.
    for n in (2, 3, 5, 8):
        q_value = n**8
        radius_a = F(q_value + 1, q_value)
        radius_b = F(1)
        _, omega_minus_squared = closed_frequencies(
            radius_a, radius_b, F(1)
        )
        expected = (
            2 * (radius_a - radius_b) ** 2
            / (2 * radius_a**2 - radius_a * radius_b
               + 2 * radius_b**2)
        )
        assert omega_minus_squared == expected
        omega_minus = sqrt(float(omega_minus_squared))
        stage_time = log(q_value)
        assert omega_minus * stage_time < 0.05
        assert (1 / omega_minus) / stage_time > 20


def rational_norm_squared(vector):
    value = base.vdot(vector, vector)
    assert value[1:] == (F(0), F(0), F(0))
    return value[0]


def check_causal_upper_scale_inputs():
    """Audit the uniform-amplitude inputs, without asserting a lower bound."""
    for t in (F(1, 2), F(1), F(2), F(-1, 2)):
        for q_value in (F(10), F(100), F(1000)):
            radius_a = 1 + 1 / q_value
            radius_b = F(1)
            target_wavevector, target_parent = scaled_mode(
                radius_a, F(0)
            )
            partner_wavevector, partner_parent = scaled_mode(
                radius_b, t
            )
            difference_wavevector = base.vadd(
                target_wavevector, base.vneg(partner_wavevector)
            )
            forcing = base.trig_symmetric_symbol(
                target_wavevector,
                target_parent,
                base.vneg(partner_wavevector),
                partner_parent,
            )
            components = [
                second_harmonic_components(component)
                for component in forcing
            ]
            forcing_zero = tuple(component[0] for component in components)
            wake, omega_squared = periodic_wake(
                difference_wavevector, forcing
            )
            wake_at_zero = tuple(
                trig_value_at_zero(component) for component in wake
            )

            # On this exact compact audit grid g_0=O(A-B) and the selected
            # periodic initial value stays uniformly bounded.  These are
            # inputs to the O(min(T,1/omega_-)) upper estimate only.
            gap_squared = (radius_a - radius_b) ** 2
            assert rational_norm_squared(forcing_zero) < 4 * gap_squared
            assert rational_norm_squared(wake_at_zero) < 16

            relative_gap_squared = (
                (radius_a - radius_b) / (radius_a + radius_b)
            ) ** 2
            frequency_ratio = omega_squared / relative_gap_squared
            assert F(1, 10) < frequency_ratio < F(100)

            directed_return_rows(
                difference_wavevector,
                partner_wavevector,
                target_parent,
                partner_parent,
            )


def main():
    check_detuned_zero_mean()
    check_adjacent_layer_schedule()
    check_causal_upper_scale_inputs()
    print("PASS C160: unequal-radius periodic sum/difference DACR means cancel exactly")
    print("PASS C160: omega_- is linear in the normalized relative gap")
    print("PASS C160: adjacent-layer detuning is too slow for O(log q)")
    print("PASS C160: causal estimate is certified only as an upper scale")
    print("Causal finite-time transient, finite epsilon, additive quartets, and localization remain OPEN")


if __name__ == "__main__":
    main()
