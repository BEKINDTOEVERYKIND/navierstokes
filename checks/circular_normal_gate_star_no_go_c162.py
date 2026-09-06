#!/usr/bin/env python3
"""Exact Q(sqrt(2),sqrt(3),i) checks for C162.

The checker proves the explicit tuned one-polarization forward/reverse
product and checks the closed forward-norm formula.  No statement about
other polarizations/heights, paired-polarization, or time-dependent gate
bundles is made.
"""

from fractions import Fraction as F

from a2_elliptic_gain_c149_c151 import (
    INV_SQRT2,
    INV_SQRT3,
    ONE,
    SQRT2,
    SQRT3,
    ZERO,
    q,
    qadd,
    qmul,
    qneg,
    qscale,
    qsub,
)


# Complexification of Q(sqrt(2),sqrt(3)).
CZERO = (ZERO, ZERO)
CONE = (ONE, ZERO)
CI = (ZERO, ONE)
CMINUS_I = (ZERO, q(-1))


def c(value):
    return (value, ZERO)


def cadd(left, right):
    return (qadd(left[0], right[0]), qadd(left[1], right[1]))


def cneg(value):
    return (qneg(value[0]), qneg(value[1]))


def csub(left, right):
    return cadd(left, cneg(right))


def cmul(left, right):
    return (
        qsub(qmul(left[0], right[0]), qmul(left[1], right[1])),
        qadd(qmul(left[0], right[1]), qmul(left[1], right[0])),
    )


def cconj(value):
    return (value[0], qneg(value[1]))


def cscale(value, scalar):
    return (qmul(value[0], scalar), qmul(value[1], scalar))


def vadd(left, right):
    return tuple(cadd(a, b) for a, b in zip(left, right))


def vscale(scalar, value):
    return tuple(cmul(scalar, entry) for entry in value)


def bilinear_dot(left, right):
    output = CZERO
    for a, b in zip(left, right):
        output = cadd(output, cmul(a, b))
    return output


def hermitian_dot(left, right):
    output = CZERO
    for a, b in zip(left, right):
        output = cadd(output, cmul(cconj(a), b))
    return output


def real_vector(*entries):
    return tuple(c(entry) for entry in entries)


def leray_project(wavevector, value):
    norm_squared = bilinear_dot(wavevector, wavevector)
    assert norm_squared[1] == ZERO
    assert norm_squared[0][1:] == (F(0), F(0), F(0))
    coefficient = cscale(
        bilinear_dot(wavevector, value),
        q(F(1, norm_squared[0][0])),
    )
    return tuple(
        csub(component, cmul(coefficient, wave_component))
        for component, wave_component in zip(value, wavevector)
    )


def symmetric_symbol(p, amplitude_p, wavevector_q, amplitude_q):
    raw = vadd(
        vscale(bilinear_dot(amplitude_p, wavevector_q), amplitude_q),
        vscale(bilinear_dot(amplitude_q, p), amplitude_p),
    )
    return leray_project(vadd(p, wavevector_q), raw)


def rational_circle(t):
    t = F(t)
    denominator = 1 + t * t
    return (
        (1 - t * t) / denominator,
        2 * t / denominator,
    )


def polynomial_multiply(left, right):
    """Multiply coefficient lists in ascending powers."""
    output = [F(0)] * (len(left) + len(right) - 1)
    for left_power, left_coefficient in enumerate(left):
        for right_power, right_coefficient in enumerate(right):
            output[left_power + right_power] += (
                left_coefficient * right_coefficient
            )
    return output


def polynomial_add(left, right):
    output = [F(0)] * max(len(left), len(right))
    for power, coefficient in enumerate(left):
        output[power] += coefficient
    for power, coefficient in enumerate(right):
        output[power] += coefficient
    return output


def bivariate_add(left, right):
    output = dict(left)
    for monomial, coefficient in right.items():
        output[monomial] = output.get(monomial, F(0)) + coefficient
        if output[monomial] == 0:
            del output[monomial]
    return output


def bivariate_scale(value, scalar):
    return {
        monomial: F(scalar) * coefficient
        for monomial, coefficient in value.items()
        if F(scalar) * coefficient != 0
    }


def bivariate_multiply(left, right):
    output = {}
    for (c_left, s_left), left_coefficient in left.items():
        for (c_right, s_right), right_coefficient in right.items():
            monomial = (c_left + c_right, s_left + s_right)
            output[monomial] = output.get(monomial, F(0)) + (
                left_coefficient * right_coefficient
            )
            if output[monomial] == 0:
                del output[monomial]
    return output


def selected_amplitude(cosine, sine):
    cosine = F(cosine)
    sine = F(sine)
    normalization = qscale(INV_SQRT2, F(1, 2))
    return real_vector(
        qmul(
            q(cosine**2 + 2 * sine**2 - sine * cosine),
            normalization,
        ),
        qmul(
            q(sine**2 + 2 * cosine**2 - sine * cosine),
            normalization,
        ),
        qmul(
            qneg(SQRT3),
            qmul(q(sine + cosine), normalization),
        ),
    )


def forward_symbol(y, t):
    """A=1 and gate height x=y/sqrt(3)."""
    y = F(y)
    cosine, sine = rational_circle(t)
    source = real_vector(q(cosine), q(sine), INV_SQRT3)
    gate = real_vector(ZERO, ZERO, qscale(INV_SQRT3, y))
    amplitude = selected_amplitude(cosine, sine)
    circular = (
        c(INV_SQRT2),
        (ZERO, INV_SQRT2),
        CZERO,
    )
    return (
        symmetric_symbol(source, amplitude, gate, circular),
        source,
        gate,
        amplitude,
        circular,
    )


def check_forward_norm_formula():
    # The global lower bound uses the exact polynomial identity
    #
    # P(y)=y^4+y^3-3y^2-2y+4
    #     =(y-1)^2(y^2+3y+2)+(2-y),       0<=y<=2,
    #     =y^2(y^2-3)+y(y^2-2)+4,         y>=2.
    #
    # Each displayed summand is nonnegative on its stated range and at
    # least one is positive.  These coefficient checks keep that proof
    # independent of numerical sampling.
    polynomial = [F(4), F(-2), F(-3), F(1), F(1)]
    first_decomposition = polynomial_add(
        polynomial_multiply([F(1), F(-2), F(1)], [F(2), F(3), F(1)]),
        [F(2), F(-1)],
    )
    second_decomposition = polynomial_add(
        polynomial_multiply([F(0), F(0), F(1)], [F(-3), F(0), F(1)]),
        polynomial_add(
            polynomial_multiply([F(0), F(1)], [F(-2), F(0), F(1)]),
            [F(4)],
        ),
    )
    assert first_decomposition == polynomial
    assert second_decomposition == polynomial

    # For c^2+s^2=1, the selected line has
    # |a_z|^2=3(c+s)^2/8 <= 3/4.  The inequality follows from the exact
    # sum-of-squares identity 2(c^2+s^2)-(c+s)^2=(c-s)^2.
    cosine = {(1, 0): F(1)}
    sine = {(0, 1): F(1)}
    cosine_plus_sine = bivariate_add(cosine, sine)
    cosine_minus_sine = bivariate_add(cosine, bivariate_scale(sine, -1))
    unit_twice = bivariate_scale(
        bivariate_add(
            bivariate_multiply(cosine, cosine),
            bivariate_multiply(sine, sine),
        ),
        2,
    )
    assert bivariate_add(
        unit_twice,
        bivariate_scale(
            bivariate_multiply(cosine_plus_sine, cosine_plus_sine),
            -1,
        ),
    ) == bivariate_multiply(cosine_minus_sine, cosine_minus_sine)

    # y=0 is retained as an algebraic endpoint check only.  A circular
    # complex coefficient at Fourier mode zero is not a real physical gate.
    for y in (F(0), F(1, 3), F(1), F(2), F(3), F(7)):
        shape = y * (y - 2) * (y + 1) * (y + 2) / (
            3 * (y * y + 2 * y + 4)
        )
        for t in (F(-3), F(-1), F(0), F(1, 2), F(1), F(2)):
            forward, source, gate, amplitude, circular = forward_symbol(
                y, t
            )
            assert bilinear_dot(source, amplitude) == CZERO
            assert hermitian_dot(amplitude, amplitude) == CONE
            assert bilinear_dot(gate, circular) == CZERO
            assert hermitian_dot(circular, circular) == CONE
            norm_squared = hermitian_dot(forward, forward)
            assert norm_squared[1] == ZERO
            assert norm_squared[0][1:] == (F(0), F(0), F(0))
            vertical_squared = qmul(
                amplitude[2][0], amplitude[2][0]
            )
            assert vertical_squared[1:] == (F(0), F(0), F(0))
            expected = F(1, 2) + shape * vertical_squared[0]
            assert norm_squared[0][0] == expected

        # Exact elementary lower bound f>=-1/3.  With y>=0,
        # f+1/3=P(y)/(3(y^2+2y+4)), P>0 by the note's two-range proof.
        polynomial = y**4 + y**3 - 3 * y**2 - 2 * y + 4
        assert shape + F(1, 3) == polynomial / (
            3 * (y * y + 2 * y + 4)
        )
        assert polynomial > 0
        assert shape >= F(-1, 3)

    # At the tuned y=2 (x=2/sqrt(3)), the angular term vanishes.
    for t in (F(-5), F(-1), F(0), F(1, 3), F(2), F(7)):
        forward, _, _, _, _ = forward_symbol(F(2), t)
        assert hermitian_dot(forward, forward) == (q(F(1, 2)), ZERO)


def check_tuned_forward_reverse_product():
    forward, source, gate, source_amplitude, circular = forward_symbol(
        F(2), F(0)
    )
    norm_squared = hermitian_dot(forward, forward)
    assert norm_squared == (q(F(1, 2)), ZERO)
    tuned_component = (q(F(1, 2)), q(F(-1, 2)))
    assert forward == (CZERO, tuned_component, CZERO)
    # |F|=1/sqrt(2), so v=sqrt(2)F.
    target_amplitude = tuple(
        cscale(component, SQRT2) for component in forward
    )
    assert hermitian_dot(target_amplitude, target_amplitude) == CONE
    assert target_amplitude == (
        CZERO,
        (INV_SQRT2, qneg(INV_SQRT2)),
        CZERO,
    )

    negative_gate = tuple(cneg(component) for component in gate)
    conjugate_circular = tuple(cconj(component) for component in circular)
    reverse_vector = symmetric_symbol(
        vadd(source, gate),
        target_amplitude,
        negative_gate,
        conjugate_circular,
    )
    assert reverse_vector == forward

    forward_ode = cmul(CMINUS_I, c(INV_SQRT2))
    reverse_projection = hermitian_dot(source_amplitude, reverse_vector)
    reverse_leakage = tuple(
        csub(component, projected)
        for component, projected in zip(
            reverse_vector,
            vscale(reverse_projection, source_amplitude),
        )
    )
    assert hermitian_dot(reverse_leakage, reverse_leakage) != CZERO
    reverse_ode = cmul(CMINUS_I, reverse_projection)
    half_inv_sqrt2 = qscale(INV_SQRT2, F(1, 2))
    assert forward_ode == (ZERO, qneg(INV_SQRT2))
    assert reverse_projection == (
        half_inv_sqrt2,
        qneg(half_inv_sqrt2),
    )
    assert reverse_ode == (
        qneg(half_inv_sqrt2),
        qneg(half_inv_sqrt2),
    )
    product = cmul(forward_ode, reverse_ode)
    expected = (q(F(-1, 4)), q(F(1, 4)))
    assert product == expected

    # Source/daughter phase changes and reciprocal positive coordinate
    # scalings leave the off-diagonal product invariant.  A skew-Hermitian
    # 2x2 block would instead have product -|u|^2, a nonpositive real.
    assert product[1] != ZERO


def main():
    check_forward_norm_formula()
    check_tuned_forward_reverse_product()
    print("PASS C162: circular pure-N forward edge has a uniform nonzero bound")
    print("PASS C162: tuned forward/reverse product is exactly -1/4+i/4")
    print("PASS C162: the reverse vector also leaks off the selected source line")
    print("The specified tuned circular pure-N block is not a skew star")
    print("Other heights/polarizations and paired/time-dependent repairs remain OPEN")


if __name__ == "__main__":
    main()
