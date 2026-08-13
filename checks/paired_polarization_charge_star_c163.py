#!/usr/bin/env python3
"""Exact arithmetic for C163's dual-helicity charge compression.

The Fourier-symbol path uses Q(sqrt(2),sqrt(3),i) arithmetic inherited
from C162.  It checks the selected return, signed-height positivity, reverse
leakage, and the equal-radius synchronization polynomial.  No invariant
full-symbol block, finite band, localization, or nonlinear stage is tested.
"""

from fractions import Fraction as F

import circular_normal_gate_star_no_go_c162 as base


def qreal(value):
    """Extract a rational real number from the exact complex field."""
    assert value[1] == base.ZERO
    assert value[0][1:] == (F(0), F(0), F(0))
    return value[0][0]


def closed_quantities(cosine, sine, y):
    cosine = F(cosine)
    sine = F(sine)
    y = F(y)
    h_squared = (cosine - sine) ** 2
    ell_squared = (cosine + sine) ** 2
    assert h_squared + ell_squared == 2
    denominator = y * y + 2 * y + 4
    forward_norm_squared = (
        ell_squared
        * h_squared
        * (4 - y * y) ** 2
        / (16 * denominator)
        + (2 * h_squared + y * ell_squared) ** 2 / 16
    )
    bracket = 1 + (
        ell_squared * y * y * (y - 2) / (4 * denominator)
    )
    directed_return = h_squared * bracket / 2
    return forward_norm_squared, directed_return, bracket


def exact_symbol_quantities(cosine, sine, y):
    cosine = F(cosine)
    sine = F(sine)
    y = F(y)
    source = base.real_vector(
        base.q(cosine), base.q(sine), base.INV_SQRT3
    )
    gate = base.real_vector(
        base.ZERO,
        base.ZERO,
        base.qscale(base.INV_SQRT3, y),
    )
    amplitude = base.selected_amplitude(cosine, sine)
    paired = base.real_vector(
        base.INV_SQRT2, base.qneg(base.INV_SQRT2), base.ZERO
    )
    forward = base.symmetric_symbol(source, amplitude, gate, paired)
    target = base.vadd(source, gate)
    negative_gate = tuple(base.cneg(component) for component in gate)
    reverse_of_forward = base.symmetric_symbol(
        target, forward, negative_gate, paired
    )
    norm_squared = qreal(base.hermitian_dot(forward, forward))
    directed_return = qreal(
        base.hermitian_dot(amplitude, reverse_of_forward)
    )
    assert base.bilinear_dot(target, forward) == base.CZERO
    assert base.bilinear_dot(source, reverse_of_forward) == base.CZERO
    return norm_squared, directed_return


def polynomial_multiply(left, right):
    output = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i + j] += a * b
    return output


def polynomial_scale(value, scalar):
    return [F(scalar) * coefficient for coefficient in value]


def check_closed_symbol():
    # Rational-circle angles are exact and exercise both signs of h and ell.
    for t in (F(-5), F(-2), F(-1), F(0), F(1, 3), F(1), F(2), F(7)):
        cosine, sine = base.rational_circle(t)
        for y in (
            F(-1),
            F(-3, 4),
            F(-1, 5),
            F(0),
            F(1, 5),
            F(1),
            F(4, 3),
            F(2),
            F(3),
            F(8),
        ):
            actual_norm, actual_return = exact_symbol_quantities(
                cosine, sine, y
            )
            closed_norm, closed_return, bracket = closed_quantities(
                cosine, sine, y
            )
            assert actual_norm == closed_norm
            assert actual_return == closed_return
            if cosine != sine and y >= -1:
                assert actual_norm > 0
                assert actual_return > 0
                assert bracket > 0

        # At the tuned height the forward norm is exactly one and the
        # directed return is beta^2=(c-s)^2/2.
        tuned_norm, tuned_return, _ = closed_quantities(
            cosine, sine, F(2)
        )
        assert tuned_norm == 1
        assert tuned_return == (cosine - sine) ** 2 / 2


def check_algebra_and_global_sign():
    # Cross-multiplying the raw return formula by 8D gives
    #   2 ell2(4-y2)+2(2-ell2)D+y ell2 D.
    # Cross-multiplying the closed formula gives
    #   4D+ell2 y2(y-2).
    # Compare coefficients in y separately at ell2=0 and ell2=1; both
    # sides are affine in ell2, so this is the exact polynomial identity.
    for ell_squared in (F(0), F(1)):
        h_squared = 2 - ell_squared
        for y in (F(-4), F(-1), F(0), F(2), F(5)):
            denominator = y * y + 2 * y + 4
            raw = (
                2 * ell_squared * (4 - y * y)
                + 2 * h_squared * denominator
                + y * ell_squared * denominator
            )
            closed = 4 * denominator + ell_squared * y * y * (y - 2)
            assert raw == closed

    # Exact factorization behind max_[0,2] y^2(2-y)=32/27:
    # 32/27-y^2(2-y)=(3y-4)^2(3y+2)/27.
    left = [F(32, 27), F(0), F(-2), F(1)]
    right = polynomial_scale(
        polynomial_multiply([F(16), F(-24), F(9)], [F(2), F(3)]),
        F(1, 27),
    )
    assert left == right

    # The note's global lower bound is then purely elementary:
    # ell^2<=2, D>=4, and y^2(2-y)<=32/27 on 0<=y<=2.
    assert 1 - F(2) * F(32, 27) / (4 * 4) == F(23, 27)
    for y in (F(0), F(1, 7), F(2, 3), F(4, 3), F(2), F(3), F(20)):
        for ell_squared in (F(0), F(1, 3), F(1), F(2)):
            denominator = y * y + 2 * y + 4
            bracket = 1 + ell_squared * y * y * (y - 2) / (
                4 * denominator
            )
            assert bracket >= F(23, 27)

    # On -1<=y<=0, y^2(2-y)<=3 and D=(y+1)^2+3>=3.
    assert 1 - F(2) * 3 / (4 * 3) == F(1, 2)
    for y in (F(-1), F(-6, 7), F(-1, 2), F(-1, 9), F(0)):
        assert y * y * (2 - y) <= 3
        denominator = y * y + 2 * y + 4
        assert denominator >= 3
        for ell_squared in (F(0), F(1, 3), F(1), F(2)):
            bracket = 1 + ell_squared * y * y * (y - 2) / (
                4 * denominator
            )
            assert bracket >= F(1, 2)


def source_complement(cosine, sine):
    """The exact unit source polarization orthogonal to selected a."""
    cosine = F(cosine)
    sine = F(sine)
    alpha = base.qscale(
        base.qmul(base.q(cosine + sine), base.INV_SQRT2), F(1, 2)
    )
    beta = base.qmul(base.q(cosine - sine), base.INV_SQRT2)
    e_s = base.real_vector(
        base.q(-cosine / 2),
        base.q(-sine / 2),
        base.qscale(base.SQRT3, F(1, 2)),
    )
    e_t = base.real_vector(
        base.q(-sine), base.q(cosine), base.ZERO
    )
    return base.vadd(
        base.vscale(base.c(beta), e_s),
        base.vscale(base.c(base.qscale(alpha, 2)), e_t),
    )


def check_reverse_leakage():
    # After division by alpha*beta, the raw source-complement coefficient is
    #   [(y+b(2-y))D-(y+2b)(4-y^2)]/D, b=beta^2.
    # The claimed numerator is y^2[(2-b)y+2(1+b)].  Both sides are affine
    # in b, so coefficient checks at b=0,1 prove the polynomial identity.
    for beta_squared in (F(0), F(1)):
        first = polynomial_multiply(
            [2 * beta_squared, 1 - beta_squared], [F(4), F(2), F(1)]
        )
        second = polynomial_multiply(
            [2 * beta_squared, F(1)], [F(4), F(0), F(-1)]
        )
        numerator = [left - right for left, right in zip(first, second)]
        claimed = [F(0), F(0), 2 * (1 + beta_squared), 2 - beta_squared]
        assert numerator == claimed

    paired = base.real_vector(
        base.INV_SQRT2, base.qneg(base.INV_SQRT2), base.ZERO
    )
    for t in (F(-3), F(-1), F(0), F(1, 3), F(2), F(4)):
        cosine, sine = base.rational_circle(t)
        source = base.real_vector(
            base.q(cosine), base.q(sine), base.INV_SQRT3
        )
        amplitude = base.selected_amplitude(cosine, sine)
        complement = source_complement(cosine, sine)
        assert base.hermitian_dot(amplitude, complement) == base.CZERO
        assert base.hermitian_dot(complement, complement) == base.CONE

        for y in (F(-1), F(-1, 2), F(0), F(1), F(2), F(5)):
            gate = base.real_vector(
                base.ZERO,
                base.ZERO,
                base.qscale(base.INV_SQRT3, y),
            )
            forward = base.symmetric_symbol(
                source, amplitude, gate, paired
            )
            reverse = base.symmetric_symbol(
                base.vadd(source, gate),
                forward,
                tuple(base.cneg(component) for component in gate),
                paired,
            )
            leakage = qreal(base.hermitian_dot(complement, reverse))
            alpha_beta = (cosine + sine) * (cosine - sine) / 4
            beta_squared = (cosine - sine) ** 2 / 2
            denominator = y * y + 2 * y + 4
            expected = (
                alpha_beta
                * y
                * y
                * (
                    (2 - beta_squared) * y
                    + 2 * (1 + beta_squared)
                )
                / denominator
            )
            assert leakage == expected

    # One exact generic witness: at phi=0,y=2, the reverse of the forward
    # daughter has coefficient 1/2 on the orthogonal source polarization.
    cosine, sine = F(1), F(0)
    source = base.real_vector(base.ONE, base.ZERO, base.INV_SQRT3)
    amplitude = base.selected_amplitude(cosine, sine)
    complement = source_complement(cosine, sine)
    gate = base.real_vector(
        base.ZERO,
        base.ZERO,
        base.qscale(base.INV_SQRT3, 2),
    )
    forward = base.symmetric_symbol(source, amplitude, gate, paired)
    reverse = base.symmetric_symbol(
        base.vadd(source, gate),
        forward,
        tuple(base.cneg(component) for component in gate),
        paired,
    )
    assert qreal(base.hermitian_dot(complement, reverse)) == F(1, 2)


def check_equal_radius_synchronization_obstruction():
    # For equal radius, every gate height y_j is independent of source
    # angle.  With x=sin(2phi), the compressed squared rate is
    # ((W+C)-W*x-C*x^2)/2.  Its rigid linear coefficient -W/2 rules out
    # constancy on an angular interval whenever the intensities are nonzero.
    gates = (
        (F(-1), F(1, 5)),
        (F(-1, 3), F(2, 7)),
        (F(1, 2), F(3, 11)),
        (F(3), F(5, 13)),
    )
    total_intensity = sum(weight for _, weight in gates)
    height_moment = sum(
        weight
        * y
        * y
        * (y - 2)
        / (4 * (y * y + 2 * y + 4))
        for y, weight in gates
    )
    assert total_intensity > 0
    assert -total_intensity / 2 < 0
    expanded = polynomial_multiply(
        [F(1), F(-1)],
        [total_intensity + height_moment, height_moment],
    )
    assert expanded == [
        total_intensity + height_moment,
        -total_intensity,
        -height_moment,
    ]

    for t in (F(-5), F(-2), F(-1, 2), F(0), F(1, 3), F(2), F(7)):
        cosine, sine = base.rational_circle(t)
        x = 2 * cosine * sine
        direct = sum(
            weight * closed_quantities(cosine, sine, y)[1]
            for y, weight in gates
        )
        polynomial = (
            total_intensity
            + height_moment
            - total_intensity * x
            - height_moment * x * x
        ) / 2
        assert direct == polynomial


def check_single_source_and_two_color_boundaries():
    plus = (
        base.c(base.INV_SQRT2),
        (base.ZERO, base.INV_SQRT2),
        base.CZERO,
    )
    minus = tuple(base.cconj(component) for component in plus)
    half_plus_i_half = (base.q(F(1, 2)), base.q(F(1, 2)))
    half_minus_i_half = base.cconj(half_plus_i_half)
    coherent_pair = base.vadd(
        base.vscale(half_plus_i_half, plus),
        base.vscale(half_minus_i_half, minus),
    )
    expected_pair = base.real_vector(
        base.INV_SQRT2, base.qneg(base.INV_SQRT2), base.ZERO
    )
    assert coherent_pair == expected_pair

    for t in (F(-3), F(-1), F(0), F(1, 3), F(1), F(4)):
        cosine, sine = base.rational_circle(t)
        source = base.real_vector(
            base.q(cosine), base.q(sine), base.INV_SQRT3
        )
        gate = base.real_vector(
            base.ZERO,
            base.ZERO,
            base.qscale(base.INV_SQRT3, F(2)),
        )
        negative_gate = tuple(base.cneg(component) for component in gate)
        amplitude = base.selected_amplitude(cosine, sine)
        tangent = base.real_vector(
            base.q(-sine), base.q(cosine), base.ZERO
        )

        # A source-specific radial linear polarization makes the selected
        # two-coordinate compression skew whenever beta is nonzero.  It
        # does not assert invariance against the other source polarization.
        radial = base.real_vector(
            base.q(cosine), base.q(sine), base.ZERO
        )
        radial_forward = base.symmetric_symbol(
            source, amplitude, gate, radial
        )
        forward_edge = base.hermitian_dot(tangent, radial_forward)
        radial_reverse = base.symmetric_symbol(
            base.vadd(source, gate), tangent, negative_gate, radial
        )
        reverse_edge = base.hermitian_dot(amplitude, radial_reverse)
        assert reverse_edge == base.cconj(forward_edge)

        # If the circular colors are orthogonalized and their products are
        # added with equal unit intensities, their imaginary parts cancel.
        directed_sum = base.CZERO
        norm_sum = base.CZERO
        for color in (plus, minus):
            forward = base.symmetric_symbol(source, amplitude, gate, color)
            scalar_forward = base.hermitian_dot(tangent, forward)
            reverse = base.symmetric_symbol(
                base.vadd(source, gate),
                tangent,
                negative_gate,
                tuple(base.cconj(component) for component in color),
            )
            scalar_reverse = base.hermitian_dot(amplitude, reverse)
            directed_sum = base.cadd(
                directed_sum,
                base.cmul(scalar_reverse, scalar_forward),
            )
            norm_sum = base.cadd(
                norm_sum,
                base.cmul(base.cconj(scalar_forward), scalar_forward),
            )
        beta_squared = (cosine - sine) ** 2 / 2
        assert directed_sum == base.c(base.q(beta_squared))
        assert norm_sum == base.CONE

    # At c=s=1/sqrt(2), beta=0.  At tuned height every horizontal gate maps
    # the selected source only into the tangent, while the reverse
    # coefficient from that tangent into the selected source is zero.
    # The reverse vector may occupy the other source polarization.
    cosine = base.INV_SQRT2
    sine = base.INV_SQRT2
    source = base.real_vector(cosine, sine, base.INV_SQRT3)
    amplitude = base.real_vector(
        base.qmul(base.q(F(1, 2)), base.INV_SQRT2),
        base.qmul(base.q(F(1, 2)), base.INV_SQRT2),
        base.qscale(base.SQRT3, F(-1, 2)),
    )
    assert base.bilinear_dot(source, amplitude) == base.CZERO
    assert base.hermitian_dot(amplitude, amplitude) == base.CONE
    tangent = base.real_vector(
        base.qneg(base.INV_SQRT2), base.INV_SQRT2, base.ZERO
    )
    gate = base.real_vector(
        base.ZERO,
        base.ZERO,
        base.qscale(base.INV_SQRT3, F(2)),
    )
    negative_gate = tuple(base.cneg(component) for component in gate)
    basis_gates = (
        base.real_vector(base.ONE, base.ZERO, base.ZERO),
        base.real_vector(base.ZERO, base.ONE, base.ZERO),
    )
    for color in basis_gates:
        forward = base.symmetric_symbol(source, amplitude, gate, color)
        # The tuned forward image has no component outside e_t.
        scalar = base.hermitian_dot(tangent, forward)
        assert forward == base.vscale(scalar, tangent)
        reverse = base.symmetric_symbol(
            base.vadd(source, gate), tangent, negative_gate, color
        )
        assert base.hermitian_dot(amplitude, reverse) == base.CZERO
        # The full reverse vector is not zero: it occupies the orthogonal
        # source polarization (which equals e_t at beta=0).  This is why
        # the triangular statement applies only after compression.
        assert base.hermitian_dot(tangent, reverse) != base.CZERO


def main():
    check_closed_symbol()
    check_algebra_and_global_sign()
    check_reverse_leakage()
    check_equal_radius_synchronization_obstruction()
    check_single_source_and_two_color_boundaries()
    print("PASS C163: dual-helicity E0 has selected P=-A^2 R<0 for y>=-1")
    print(
        "PASS C163: weighted-skew identity is a compression; "
        "reverse leakage is exact"
    )
    print(
        "PASS C163: equal-radius shared intensities cannot synchronize "
        "an angular interval"
    )
    print("PASS C163: balanced colors retain the selected beta=0 rank defect")


if __name__ == "__main__":
    main()
