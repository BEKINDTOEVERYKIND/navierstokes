#!/usr/bin/env python3
"""Exact checks for C157--C158 on the limiting A2 resonant cone.

C157 reconstructs an explicit two-radius additive quartet in
Q(sqrt(2),sqrt(3)) and checks that its two projected quadratic wake
vectors are linearly independent whenever the radii differ.

C158 is only an exact reduced-half-lattice diagonal upper normalization.
It does not identify the coefficient l1 scale with physical point
amplitude without a coherence factor, and it does not assert that the
physical thick-packet kernel is coercive; C156's former
aggregate-coercivity claim has been retracted.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement

from a2_elliptic_gain_c149_c151 import (
    INV_SQRT2,
    INV_SQRT3,
    ONE,
    SQRT3,
    ZERO,
    q,
    qadd,
    qmul,
    qneg,
    qscale,
    symmetric_symbol,
    vadd,
    vcross,
    vdot,
    vneg,
    vscale,
)


def resonant_amplitude(cosine, sine):
    """The limiting selected unit line a(phi), exactly."""
    cosine = F(cosine)
    sine = F(sine)
    normalization = qscale(INV_SQRT2, F(1, 2))
    return (
        qmul(
            q(cosine * cosine + 2 * sine * sine - sine * cosine),
            normalization,
        ),
        qmul(
            q(sine * sine + 2 * cosine * cosine - sine * cosine),
            normalization,
        ),
        qmul(
            qneg(SQRT3),
            qmul(q(sine + cosine), normalization),
        ),
    )


def resonant_wavevector(radius, cosine, sine):
    return vscale(
        q(radius),
        (q(cosine), q(sine), INV_SQRT3),
    )


def determinant_formula(radius_a, radius_b):
    radius_a = F(radius_a)
    radius_b = F(radius_b)
    denominator = radius_a * radius_a + radius_b * radius_b
    quartic = (
        3 * radius_a**4
        - 4 * radius_a**3 * radius_b
        + 14 * radius_a * radius_a * radius_b * radius_b
        - 4 * radius_a * radius_b**3
        + 3 * radius_b**4
    )
    coefficient = (
        -radius_a
        * radius_b
        * (radius_a - radius_b)
        * (radius_a + radius_b) ** 2
        * quartic
        / (12 * denominator**3)
    )
    return qscale(SQRT3, coefficient)


def determinant_gap_coefficient_at_equal(common_radius):
    """Evaluate the determinant after cancelling A-B, then set A=B."""
    common_radius = F(common_radius)
    denominator = 2 * common_radius * common_radius
    quartic = 12 * common_radius**4
    coefficient = (
        -common_radius**2
        * (2 * common_radius) ** 2
        * quartic
        / (12 * denominator**3)
    )
    return qscale(SQRT3, coefficient)


def check_two_radius_additive_quartet():
    for radius_a, radius_b in (
        (F(1), F(2)),
        (F(2), F(3)),
        (F(3, 4), F(5, 6)),
        (F(11, 10), F(1)),
        (F(101, 100), F(1)),
    ):
        denominator = radius_a * radius_a + radius_b * radius_b

        # First decomposition: A at angle 0 and B at angle pi/2.
        p = resonant_wavevector(radius_a, F(1), F(0))
        q_wave = resonant_wavevector(radius_b, F(0), F(1))
        amplitude_p = resonant_amplitude(F(1), F(0))
        amplitude_q = resonant_amplitude(F(0), F(1))

        # The second decomposition is reflection across the horizontal
        # output direction (A,B).
        cosine_alpha = (
            radius_a * radius_a - radius_b * radius_b
        ) / denominator
        sine_alpha = 2 * radius_a * radius_b / denominator
        cosine_beta = 2 * radius_a * radius_b / denominator
        sine_beta = (
            radius_b * radius_b - radius_a * radius_a
        ) / denominator
        p_star = resonant_wavevector(
            radius_a, cosine_alpha, sine_alpha
        )
        q_star = resonant_wavevector(
            radius_b, cosine_beta, sine_beta
        )
        amplitude_p_star = resonant_amplitude(
            cosine_alpha, sine_alpha
        )
        amplitude_q_star = resonant_amplitude(
            cosine_beta, sine_beta
        )

        # The selected polarization is genuinely a unit transverse line.
        for radius, cosine, sine, wavevector, amplitude in (
            (radius_a, F(1), F(0), p, amplitude_p),
            (radius_b, F(0), F(1), q_wave, amplitude_q),
            (
                radius_a,
                cosine_alpha,
                sine_alpha,
                p_star,
                amplitude_p_star,
            ),
            (
                radius_b,
                cosine_beta,
                sine_beta,
                q_star,
                amplitude_q_star,
            ),
        ):
            assert cosine * cosine + sine * sine == 1
            assert wavevector == resonant_wavevector(radius, cosine, sine)
            assert vdot(wavevector, amplitude) == ZERO
            assert vdot(amplitude, amplitude) == ONE

        output = vadd(p, q_wave)
        assert output == vadd(p_star, q_star)
        assert p != p_star
        assert q_wave != q_star

        wake = symmetric_symbol(p, amplitude_p, q_wave, amplitude_q)
        wake_star = symmetric_symbol(
            p_star,
            amplitude_p_star,
            q_star,
            amplitude_q_star,
        )
        assert vdot(output, wake) == ZERO
        assert vdot(output, wake_star) == ZERO
        determinant = vdot(output, vcross(wake, wake_star))
        assert determinant == determinant_formula(radius_a, radius_b)
        assert determinant != ZERO

        # Exact positivity factorization of the quartic determinant
        # factor.  Put y=A/B+B/A >= 2.
        y = radius_a / radius_b + radius_b / radius_a
        quartic = (
            3 * radius_a**4
            - 4 * radius_a**3 * radius_b
            + 14 * radius_a * radius_a * radius_b * radius_b
            - 4 * radius_a * radius_b**3
            + 3 * radius_b**4
        )
        assert quartic == (
            radius_a * radius_a * radius_b * radius_b
            * (3 * y * y - 4 * y + 8)
        )
        assert y >= 2
        assert 3 * y * y - 4 * y + 8 > 0
        assert (
            3 * y * y - 4 * y + 8
            == 3 * (y - F(2, 3)) ** 2 + F(20, 3)
        )

        # Reality adds no third unordered decomposition of +output inside
        # this eight-mode closure.  Negative partners give -output.
        positive_modes = {
            "p": p,
            "q": q_wave,
            "p_star": p_star,
            "q_star": q_star,
        }
        reality_closed = dict(positive_modes)
        reality_closed.update(
            {
                f"minus_{name}": vneg(value)
                for name, value in positive_modes.items()
            }
        )
        output_pairs = {
            tuple(sorted((left_name, right_name)))
            for (left_name, left), (right_name, right) in combinations_with_replacement(
                reality_closed.items(), 2
            )
            if vadd(left, right) == output
        }
        assert output_pairs == {
            tuple(sorted(("p", "q"))),
            tuple(sorted(("p_star", "q_star"))),
        }

    # The determinant vanishes linearly, rather than uniformly, when the
    # radii coalesce.  Cancel A-B in the exact factorization and then set
    # A=B=a: the coefficient is exactly -(sqrt(3)/2)*a^2.
    for common_radius in (F(1), F(3, 2), F(7, 3)):
        limit = qscale(SQRT3, -common_radius * common_radius / 2)
        assert determinant_gap_coefficient_at_equal(common_radius) == limit


def check_diagonal_packet_normalization():
    # In one-positive-representative-per-reality-pair coordinates, equal
    # moduli c give A_l1=M*c and E_+=M*c^2=A_l1^2/M exactly.  The physical
    # real-field Fourier energy is 2E_+; point amplitude needs coherence.
    for mode_count in (1, 8, 27, 125):
        coefficient = F(3, 17 * mode_count)
        coefficient_l1 = mode_count * coefficient
        positive_energy = mode_count * coefficient * coefficient
        physical_fourier_energy = 2 * positive_energy
        assert positive_energy == coefficient_l1**2 / mode_count
        assert physical_fourier_energy == 2 * coefficient_l1**2 / mode_count

    # The selected cone line has a uniform positive projection.  Thus
    # aligned positive real coefficients and their conjugates do give a
    # physical point scale comparable, but not equal, to A_l1.
    coherence_direction = (INV_SQRT2, INV_SQRT2, ZERO)
    angles = (
        (F(1), F(0)),
        (F(0), F(1)),
        (F(3, 5), F(4, 5)),
        (F(-5, 13), F(12, 13)),
    )
    coefficient = F(7, 19)
    coefficient_l1 = len(angles) * coefficient
    physical_peak = (ZERO, ZERO, ZERO)
    for cosine, sine in angles:
        amplitude = resonant_amplitude(cosine, sine)
        projection = vdot(coherence_direction, amplitude)
        assert projection[1:] == (F(0), F(0), F(0))
        assert F(1, 2) <= projection[0] <= F(1)
        # A reality pair with positive real coefficient contributes 2ca.
        physical_peak = vadd(
            physical_peak,
            vscale(q(2 * coefficient), amplitude),
        )
    physical_projection = vdot(coherence_direction, physical_peak)
    assert physical_projection[1:] == (F(0), F(0), F(0))
    assert coefficient_l1 <= physical_projection[0] <= 2 * coefficient_l1

    # If |kappa_ij|<=kappa_star in a diagonal cubic map, its exact upper
    # parameter is Q^2*kappa_star*E_+.  This is an upper normalization, not
    # a lower/coercive physical claim.
    squared_amplitudes = (F(1, 7), F(2, 9), F(5, 11), F(3, 8))
    energy = sum(squared_amplitudes)
    carrier = F(13)
    kappa_star = F(5, 3)
    diagonal_squared_norm = F(0)
    input_squared_norm = energy
    for amplitude_squared in squared_amplitudes:
        row_mass = energy - amplitude_squared
        multiplier = carrier**2 * kappa_star * row_mass
        diagonal_squared_norm += amplitude_squared * multiplier**2
    upper_multiplier = carrier**2 * kappa_star * energy
    assert diagonal_squared_norm <= upper_multiplier**2 * input_squared_norm

    # C147 exponent ledger in A_l1 coordinates: q=n^8, M=q^3=n^24,
    # Q=q=n^8.  A uniform order-one coherence factor preserves these
    # exponents but is not asserted by this arithmetic alone.
    for n in (2, 3, 5):
        q_value = n**8
        mode_count = q_value**3
        carrier = q_value
        seed_l1 = F(1, n**16)
        final_l1 = F(n**10)
        unit_upper_parameter_l1 = F(n**4)
        seed_parameter = (
            carrier**2 * seed_l1**2 / mode_count
        )
        final_parameter = (
            carrier**2 * final_l1**2 / mode_count
        )
        assert seed_parameter == F(1, n**40)
        assert final_parameter == F(n**12)
        assert unit_upper_parameter_l1**2 == F(mode_count, carrier**2)
        assert unit_upper_parameter_l1 / seed_l1 == n**20
        assert final_l1 / seed_l1 == n**26
        assert final_l1 / unit_upper_parameter_l1 == n**6


def main():
    check_two_radius_additive_quartet()
    check_diagonal_packet_normalization()
    print("PASS C157: explicit two-radius additive quartet has rank two for A!=B")
    print("PASS C158: reduced-half-lattice diagonal upper normalization and exponents")
    print("No thick-packet coercivity, finite-epsilon persistence, or stage closure is claimed")


if __name__ == "__main__":
    main()
