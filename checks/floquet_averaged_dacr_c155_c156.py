#!/usr/bin/env python3
"""Exact/algebraic checks for C155--C156.

The checker has three layers.

1. It reconstructs the fixed-ring pair symbols in Q(sqrt(2),sqrt(3))
   using the C149--C151 exact arithmetic.
2. It checks the reduced quarter-period wake equations and the exact
   full-period coefficient, including a rigorous elementary sign bound.
3. It checks the rational directed general-angle secular kernel, reconstructs
   physical forward/reverse returns for general rational pairs, validates an
   exact positive balanced triple, and checks the fixed-ring support logic.
   The reverse check retracts the formerly claimed symmetric coercive part.

Floating point is used only to print/check the displayed transcendental
decimal and to run an independent Simpson cross-check of the raw-derived
reduced integrand; nonvanishing follows from pi>3 and sqrt(10)<4.
"""

from fractions import Fraction as F
from math import cos, factorial, pi, sin, sqrt

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
    vdot,
    vneg,
    vscale,
)


def qratio(numerator, denominator):
    """Division in Q for an element known to be rational."""
    assert denominator[1:] == (F(0), F(0), F(0))
    return qscale(numerator, F(1, 1) / denominator[0])


def fixed_ring_mode(t):
    """Return eta_phi,a_phi for t=tan(phi/2), exactly."""
    t = F(t)
    denominator = 1 + t * t
    cosine = (1 - t * t) / denominator
    sine = 2 * t / denominator
    eta = (q(cosine), q(sine), INV_SQRT3)

    x_component = cosine * cosine + 2 * sine * sine - sine * cosine
    y_component = sine * sine + 2 * cosine * cosine - sine * cosine
    z_component = qmul(
        qneg(SQRT3),
        q(sine + cosine),
    )
    amplitude = (
        qmul(q(x_component), qscale(INV_SQRT2, F(1, 2))),
        qmul(q(y_component), qscale(INV_SQRT2, F(1, 2))),
        qmul(z_component, qscale(INV_SQRT2, F(1, 2))),
    )
    assert vdot(eta, amplitude) == ZERO
    assert vdot(amplitude, amplitude) == ONE
    return eta, amplitude


def check_pair_symbols():
    eta_zero = (ONE, ZERO, INV_SQRT3)
    amplitude_zero = (
        qscale(INV_SQRT2, F(1, 2)),
        INV_SQRT2,
        qneg(qmul(SQRT3, qscale(INV_SQRT2, F(1, 2)))),
    )

    for t in (F(1, 5), F(1, 2), F(1), F(3, 2), F(3)):
        eta, amplitude = fixed_ring_mode(t)
        plus = symmetric_symbol(
            eta_zero,
            amplitude_zero,
            eta,
            amplitude,
        )
        minus = symmetric_symbol(
            eta_zero,
            amplitude_zero,
            tuple(qneg(entry) for entry in eta),
            amplitude,
        )

        plus_scalar = (
            2 * t * t * (t * t + 2)
            / ((1 + t * t) ** 2 * (t * t + 4))
        )
        minus_scalar = (
            t * (t * t - 1) / (1 + t * t) ** 2
        )
        direction = (ONE, q(t), qneg(SQRT3))
        expected_plus = tuple(
            qmul(q(plus_scalar), entry) for entry in direction
        )
        expected_minus = tuple(
            qmul(q(minus_scalar), entry) for entry in direction
        )
        assert plus == expected_plus
        assert minus == expected_minus

    # The plus symbol is nonzero for every finite t>0.  At t=1 the
    # difference symbol vanishes, but the plus symbol remains nonzero.
    eta_quarter, amplitude_quarter = fixed_ring_mode(F(1))
    assert symmetric_symbol(
        eta_zero,
        amplitude_zero,
        eta_quarter,
        amplitude_quarter,
    ) != (ZERO, ZERO, ZERO)
    assert symmetric_symbol(
        eta_zero,
        amplitude_zero,
        tuple(qneg(entry) for entry in eta_quarter),
        amplitude_quarter,
    ) == (ZERO, ZERO, ZERO)

    # The leading degrees in (3.5)--(3.6) give the antipodal limits:
    # plus -> 0 and minus -> (0,1,0).  Exact rational samples verify the
    # stated convergence without pretending that an assigned tuple is a
    # derivation.
    previous_error = None
    for t in (F(10), F(100), F(1000)):
        plus_scalar = (
            2 * t * t * (t * t + 2)
            / ((1 + t * t) ** 2 * (t * t + 4))
        )
        minus_scalar = t * (t * t - 1) / (1 + t * t) ** 2
        plus_error = max(
            abs(plus_scalar),
            abs(t * plus_scalar),
            abs(plus_scalar),
        )
        minus_error = max(
            abs(minus_scalar),
            abs(t * minus_scalar - 1),
        )
        error = max(plus_error, minus_error)
        if previous_error is not None:
            assert error < previous_error
        previous_error = error


def matrix_multiply(left, right):
    return [
        [
            sum(
                (left[row][inner] * right[inner][column]
                 for inner in range(len(right))),
                F(0),
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def check_quarter_reduced_system():
    # Reconstruct the displayed bases, generator, forcing, and return rows
    # from the raw Leray-projected vector symbols.  This is deliberately
    # called before using the reduced coordinates below: a previous version
    # merely asserted the displayed forcing vector against itself.
    check_quarter_raw_reduction()

    a_sum = [[F(0), F(2)], [F(-4, 5), F(0)]]
    a_sum_squared = matrix_multiply(a_sum, a_sum)
    assert a_sum_squared == [
        [F(-8, 5), F(0)],
        [F(0), F(-8, 5)],
    ]

    # The sum wake does not reset after one parent period.  If both of its
    # endpoint coordinates vanished, sqrt(8/5) would be an integer; its
    # squared value is not even integral.  This is the exact non-Markov
    # boundary for the parent-only projected map.
    omega_squared = F(8, 5)
    assert omega_squared.denominator != 1

    # Difference-wake contribution.  Write S=sin(2s), C=cos(2s).
    # A polynomial is a dict (power(S),power(C))->coefficient.  Its
    # period mean uses the exact beta-integral moments.
    def poly_add(left, right):
        out = dict(left)
        for key, value in right.items():
            out[key] = out.get(key, F(0)) + value
            if out[key] == 0:
                del out[key]
        return out

    def poly_scale(value, scalar):
        return {
            key: coefficient * scalar
            for key, coefficient in value.items()
            if coefficient * scalar
        }

    def poly_multiply(left, right):
        out = {}
        for (s_left, c_left), a in left.items():
            for (s_right, c_right), b in right.items():
                key = (s_left + s_right, c_left + c_right)
                out[key] = out.get(key, F(0)) + a * b
        return {key: value for key, value in out.items() if value}

    def trig_mean(value):
        # Only degrees <=2 occur here.
        moments = {
            (0, 0): F(1),
            (1, 0): F(0),
            (0, 1): F(0),
            (2, 0): F(1, 2),
            (0, 2): F(1, 2),
            (1, 1): F(0),
        }
        return sum(
            coefficient * moments[key]
            for key, coefficient in value.items()
        )

    one = {(0, 0): F(1)}
    sine = {(1, 0): F(1)}
    cosine = {(0, 1): F(1)}

    w_delta_1 = poly_scale(poly_add(one, poly_scale(cosine, -1)), F(1, 4))
    # Remove sqrt(3) from both h_delta_2 and w_delta_2 and account for
    # their product by the factor 3.
    w_delta_2_over_sqrt3 = poly_scale(w_delta_1, -1)
    h_delta_1 = poly_add(
        poly_add(poly_scale(sine, F(-1, 2)),
                 poly_scale(cosine, F(-3, 4))),
        {(0, 0): F(-1, 2)},
    )
    h_delta_2_over_sqrt3 = poly_add(
        poly_add(poly_scale(sine, F(1, 4)),
                 poly_scale(cosine, F(-1, 4))),
        {(0, 0): F(1, 12)},
    )
    delta_mean = trig_mean(
        poly_add(
            poly_multiply(h_delta_1, w_delta_1),
            poly_scale(
                poly_multiply(
                    h_delta_2_over_sqrt3,
                    w_delta_2_over_sqrt3,
                ),
                3,
            ),
        )
    )
    assert delta_mean == F(-3, 16)
    # Multiplication by the period 2*pi gives I_delta=-3*pi/8.

    # The physical forward/reverse means are reconstructed later directly
    # from the raw vector symbols.  In units of the period factor pi they
    # split as sum plus difference contributions (-3/4,-3/8) and
    # (+3/4,+3/8), respectively.  This replaces hard-coded sum means and a
    # vacuous scalar "X-X=0" check from the earlier audit draft.
    direct_components = physical_directed_kernel_components(F(0), F(1))
    reverse_components = physical_directed_kernel_components(F(1), F(0))
    assert direct_components == (F(-3, 4), F(-3, 8))
    assert reverse_components == (F(3, 4), F(3, 8))
    assert sum(direct_components) == F(-9, 8)
    assert sum(reverse_components) == F(9, 8)
    assert sum(direct_components) + sum(reverse_components) == 0


def first_period_coefficient():
    theta = 4 * sqrt(10) * pi / 5
    i_sum = (
        3
        * (-16 * pi - 3 + 3 * cos(theta) + 4 * sqrt(10) * sin(theta))
        / 64
    )
    i_delta = -3 * pi / 8
    return i_sum + i_delta


def check_first_period_coefficient():
    coefficient = first_period_coefficient()
    assert abs(coefficient - (-3.097745763580765)) < 2e-14

    # Independent numerical integration of the raw-derived rows and the
    # causal wake solutions.  This does not supply the rigorous sign (the
    # rational enclosure below does); it catches a missing path, sign, or
    # factor in the closed transcendental expression.
    omega = sqrt(F(8, 5))

    def reduced_integrand(s):
        h_sum_1 = cos(2 * s) / 4 - 1
        h_sum_2 = -5 * sin(2 * s) / 4 - F(3, 4)
        wake_sum_1 = F(3, 8) * (1 - cos(omega * s))
        wake_sum_2 = 3 * omega * sin(omega * s) / 16
        h_difference_1 = (
            -sin(2 * s) / 2 - 3 * cos(2 * s) / 4 - F(1, 2)
        )
        h_difference_2 = sqrt(3) * (
            F(1, 12) + sin(2 * s) / 4 - cos(2 * s) / 4
        )
        wake_difference_1 = (1 - cos(2 * s)) / 4
        wake_difference_2 = -sqrt(3) * (1 - cos(2 * s)) / 4
        return float(
            h_sum_1 * wake_sum_1
            + h_sum_2 * wake_sum_2
            + h_difference_1 * wake_difference_1
            + h_difference_2 * wake_difference_2
        )

    panels = 16384
    step = 2 * pi / panels
    simpson = reduced_integrand(0) + reduced_integrand(2 * pi)
    simpson += 4 * sum(
        reduced_integrand(index * step)
        for index in range(1, panels, 2)
    )
    simpson += 2 * sum(
        reduced_integrand(index * step)
        for index in range(2, panels, 2)
    )
    simpson *= step / 3
    assert abs(simpson - coefficient) < 2e-12

    # Rigorous elementary enclosure:
    # cos(theta)<=1, sin(theta)<=1, pi>3, sqrt(10)<4.
    # The bracket in (2.9) is therefore <-56.
    bracket_upper = F(-24) * 3 - 3 + 3 + F(4) * 4
    assert bracket_upper == F(-56)
    coefficient_upper = F(3, 64) * bracket_upper
    assert coefficient_upper == F(-21, 8)
    assert coefficient < float(coefficient_upper)

    # The monomial has phases z1*z2*conj(z2)=z1*|z2|^2.
    # Test this identity on exact Gaussian-integer scalar phases.
    for z1 in (complex(1, 0), complex(2, -3), complex(-4, 5)):
        for z2 in (complex(1, 1), complex(-2, 3), complex(4, 0)):
            assert abs(z1 * z2 * z2.conjugate() - z1 * abs(z2) ** 2) < 1e-12


# Trigonometric polynomials in sin(s),cos(s), with coefficients in
# Q(sqrt(2),sqrt(3)).  Reduction modulo sin(s)^2+cos(s)^2=1 leaves
# cos-power zero or one and makes the raw general-angle check exact.
def trig_normalize(value):
    pending = list(value.items())
    out = {}
    while pending:
        (sin_power, cos_power), coefficient = pending.pop()
        if coefficient == ZERO:
            continue
        if cos_power >= 2:
            pending.append(((sin_power, cos_power - 2), coefficient))
            pending.append(((sin_power + 2, cos_power - 2), qneg(coefficient)))
            continue
        key = (sin_power, cos_power)
        out[key] = qadd(out.get(key, ZERO), coefficient)
        if out[key] == ZERO:
            del out[key]
    return out


def trig_constant(value):
    return {} if value == ZERO else {(0, 0): value}


TRIG_ZERO = {}
TRIG_ONE = trig_constant(ONE)
TRIG_SIN = {(1, 0): ONE}
TRIG_COS = {(0, 1): ONE}


def trig_add(left, right):
    out = dict(left)
    for key, coefficient in right.items():
        out[key] = qadd(out.get(key, ZERO), coefficient)
        if out[key] == ZERO:
            del out[key]
    return trig_normalize(out)


def trig_neg(value):
    return {key: qneg(coefficient) for key, coefficient in value.items()}


def trig_sub(left, right):
    return trig_add(left, trig_neg(right))


def trig_scale(value, scalar):
    return trig_normalize(
        {
            key: qmul(scalar, coefficient)
            for key, coefficient in value.items()
            if qmul(scalar, coefficient) != ZERO
        }
    )


def trig_multiply(left, right):
    out = {}
    for (sin_left, cos_left), coefficient_left in left.items():
        for (sin_right, cos_right), coefficient_right in right.items():
            key = (sin_left + sin_right, cos_left + cos_right)
            coefficient = qmul(coefficient_left, coefficient_right)
            out[key] = qadd(out.get(key, ZERO), coefficient)
    return trig_normalize(out)


def trig_derivative(value):
    out = {}
    for (sin_power, cos_power), coefficient in value.items():
        if sin_power:
            key = (sin_power - 1, cos_power + 1)
            term = qscale(coefficient, sin_power)
            out[key] = qadd(out.get(key, ZERO), term)
        if cos_power:
            key = (sin_power + 1, cos_power - 1)
            term = qscale(coefficient, -cos_power)
            out[key] = qadd(out.get(key, ZERO), term)
    return trig_normalize(out)


def trig_mean(value):
    out = ZERO
    for (sin_power, cos_power), coefficient in value.items():
        if sin_power % 2 or cos_power % 2:
            continue
        a = sin_power // 2
        b = cos_power // 2
        moment = F(
            factorial(2 * a) * factorial(2 * b),
            4 ** (a + b)
            * factorial(a)
            * factorial(b)
            * factorial(a + b),
        )
        out = qadd(out, qscale(coefficient, moment))
    return out


def trig_vector_constant(value):
    return tuple(trig_constant(entry) for entry in value)


def trig_vector_add(left, right):
    return tuple(trig_add(a, b) for a, b in zip(left, right))


def trig_vector_neg(value):
    return tuple(trig_neg(entry) for entry in value)


def trig_vector_scale(scalar, value):
    return tuple(trig_multiply(scalar, entry) for entry in value)


def trig_dot(left, right):
    out = TRIG_ZERO
    for a, b in zip(left, right):
        out = trig_add(out, trig_multiply(a, b))
    return out


def trig_dot_constant(left, right):
    out = TRIG_ZERO
    for polynomial, constant in zip(left, right):
        out = trig_add(out, trig_scale(polynomial, constant))
    return out


def trig_divide_rational(value, denominator):
    return {
        key: qratio(coefficient, denominator)
        for key, coefficient in value.items()
    }


def trig_leray_project(wavevector, value):
    norm_squared = vdot(wavevector, wavevector)
    coefficient = trig_divide_rational(
        trig_dot_constant(value, wavevector),
        norm_squared,
    )
    return tuple(
        trig_sub(component, trig_scale(coefficient, wave_component))
        for component, wave_component in zip(value, wavevector)
    )


def trig_symmetric_symbol(p, amplitude_p, wavevector_q, amplitude_q):
    output = trig_vector_add(
        trig_vector_scale(
            trig_dot_constant(amplitude_p, wavevector_q),
            amplitude_q,
        ),
        trig_vector_scale(
            trig_dot_constant(amplitude_q, p),
            amplitude_p,
        ),
    )
    return trig_leray_project(vadd(p, wavevector_q), output)


def kelvin_apply(wavevector, value):
    j_value = (qneg(value[1]), value[0], ZERO)
    norm_squared = vdot(wavevector, wavevector)
    projection_scalar = qratio(vdot(wavevector, j_value), norm_squared)
    return vadd(
        vscale(q(-2), j_value),
        vscale(qmul(q(2), projection_scalar), wavevector),
    )


def rotating_parent(cosine_phi, sine_phi):
    """Exact b_phi(s)=R(-s)a(s+phi) as a trig-polynomial vector."""
    cosine_phi = q(cosine_phi)
    sine_phi = q(sine_phi)
    sine_theta = trig_add(
        trig_scale(TRIG_SIN, cosine_phi),
        trig_scale(TRIG_COS, sine_phi),
    )
    cosine_theta = trig_sub(
        trig_scale(TRIG_COS, cosine_phi),
        trig_scale(TRIG_SIN, sine_phi),
    )
    sine_squared = trig_multiply(sine_theta, sine_theta)
    cosine_squared = trig_multiply(cosine_theta, cosine_theta)
    sine_cosine = trig_multiply(sine_theta, cosine_theta)
    normalization = qscale(INV_SQRT2, F(1, 2))
    amplitude_x = trig_scale(
        trig_sub(
            trig_add(cosine_squared, trig_scale(sine_squared, q(2))),
            sine_cosine,
        ),
        normalization,
    )
    amplitude_y = trig_scale(
        trig_sub(
            trig_add(sine_squared, trig_scale(cosine_squared, q(2))),
            sine_cosine,
        ),
        normalization,
    )
    amplitude_z = trig_scale(
        trig_add(sine_theta, cosine_theta),
        qneg(qmul(SQRT3, normalization)),
    )
    return (
        trig_add(
            trig_multiply(TRIG_COS, amplitude_x),
            trig_multiply(TRIG_SIN, amplitude_y),
        ),
        trig_add(
            trig_neg(trig_multiply(TRIG_SIN, amplitude_x)),
            trig_multiply(TRIG_COS, amplitude_y),
        ),
        amplitude_z,
    )


def check_quarter_raw_reduction():
    """Derive (2.2)--(2.6) from the un-reduced vector symbols."""
    eta_zero = (ONE, ZERO, INV_SQRT3)
    eta_quarter = (ZERO, ONE, INV_SQRT3)
    eta_sum = vadd(eta_zero, eta_quarter)
    eta_difference = vadd(eta_zero, vneg(eta_quarter))
    parent_zero = rotating_parent(F(1), F(0))
    parent_quarter = rotating_parent(F(0), F(1))

    sum_basis = (
        (ONE, q(-1), ZERO),
        (ONE, ONE, qneg(SQRT3)),
    )
    difference_basis = (
        (ONE, ONE, ZERO),
        (ZERO, ZERO, ONE),
    )

    # Generator columns in the displayed transverse bases.
    assert kelvin_apply(eta_sum, sum_basis[0]) == vscale(
        q(F(-4, 5)), sum_basis[1]
    )
    assert kelvin_apply(eta_sum, sum_basis[1]) == vscale(
        q(2), sum_basis[0]
    )
    assert kelvin_apply(
        eta_difference, difference_basis[0]
    ) == (ZERO, ZERO, ZERO)
    assert kelvin_apply(
        eta_difference, difference_basis[1]
    ) == (ZERO, ZERO, ZERO)

    sum_forcing = trig_symmetric_symbol(
        eta_zero, parent_zero, eta_quarter, parent_quarter
    )
    expected_sum_forcing = trig_vector_constant(
        vscale(q(F(3, 10)), sum_basis[1])
    )
    assert sum_forcing == expected_sum_forcing

    difference_forcing = trig_symmetric_symbol(
        eta_zero, parent_zero, vneg(eta_quarter), parent_quarter
    )
    sine_cosine = trig_multiply(TRIG_SIN, TRIG_COS)
    expected_difference_forcing = trig_vector_scale(
        sine_cosine,
        trig_vector_constant((ONE, ONE, qneg(SQRT3))),
    )
    assert difference_forcing == expected_difference_forcing

    def return_row(wake_wavevector, basis, partner_wavevector):
        return tuple(
            trig_dot(
                parent_zero,
                trig_vector_neg(
                    trig_symmetric_symbol(
                        wake_wavevector,
                        trig_vector_constant(vector),
                        partner_wavevector,
                        parent_quarter,
                    )
                ),
            )
            for vector in basis
        )

    sum_row = return_row(eta_sum, sum_basis, vneg(eta_quarter))
    difference_row = return_row(
        eta_difference, difference_basis, eta_quarter
    )
    sine_squared = trig_multiply(TRIG_SIN, TRIG_SIN)

    expected_sum_row = (
        trig_add(
            trig_constant(q(F(-3, 4))),
            trig_scale(sine_squared, q(F(-1, 2))),
        ),
        trig_add(
            trig_constant(q(F(-3, 4))),
            trig_scale(sine_cosine, q(F(-5, 2))),
        ),
    )
    expected_difference_row = (
        trig_add(
            trig_add(
                trig_constant(q(F(-5, 4))),
                trig_scale(sine_cosine, q(-1)),
            ),
            trig_scale(sine_squared, q(F(3, 2))),
        ),
        trig_add(
            trig_add(
                trig_constant(qmul(SQRT3, q(F(-1, 6)))),
                trig_scale(
                    sine_cosine, qmul(SQRT3, q(F(1, 2)))
                ),
            ),
            trig_scale(
                sine_squared, qmul(SQRT3, q(F(1, 2)))
            ),
        ),
    )
    assert sum_row == expected_sum_row
    assert difference_row == expected_difference_row


def raw_kernel_ratio(t):
    """Derive K(t)/pi from the raw vector symbols and exact time means."""
    t = F(t)
    denominator = 1 + t * t
    cosine_phi = (1 - t * t) / denominator
    sine_phi = 2 * t / denominator
    eta_zero = (ONE, ZERO, INV_SQRT3)
    eta_phi = (q(cosine_phi), q(sine_phi), INV_SQRT3)
    eta_sum = vadd(eta_zero, eta_phi)
    eta_difference = vadd(eta_zero, vneg(eta_phi))
    parent_zero = rotating_parent(F(1), F(0))
    parent_phi = rotating_parent(cosine_phi, sine_phi)

    forcing_sum = trig_symmetric_symbol(
        eta_zero,
        parent_zero,
        eta_phi,
        parent_phi,
    )
    assert all(set(component).issubset({(0, 0)}) for component in forcing_sum)
    forcing_sum_constant = tuple(
        component.get((0, 0), ZERO) for component in forcing_sum
    )
    equilibrium = (
        q(t**3 * (t * t + 2) / denominator**3),
        q(-t * t * (t * t + 2) / denominator**3),
        ZERO,
    )
    assert vadd(
        kelvin_apply(eta_sum, equilibrium),
        forcing_sum_constant,
    ) == (ZERO, ZERO, ZERO)

    sine_two = trig_scale(
        trig_multiply(TRIG_SIN, TRIG_COS),
        q(2),
    )
    cosine_two = trig_sub(
        trig_multiply(TRIG_COS, TRIG_COS),
        trig_multiply(TRIG_SIN, TRIG_SIN),
    )
    f_scalar = trig_scale(
        trig_add(
            trig_add(
                trig_scale(sine_two, q(t * t - 1)),
                trig_scale(cosine_two, q(-2 * t)),
            ),
            trig_constant(q(2 * t)),
        ),
        q(F(1, 2) / denominator**2),
    )
    difference_direction = (q(t), q(t * t), qmul(q(-t), SQRT3))
    wake_difference = tuple(
        trig_scale(f_scalar, component)
        for component in difference_direction
    )
    assert kelvin_apply(
        eta_difference,
        difference_direction,
    ) == (ZERO, ZERO, ZERO)
    forcing_difference = trig_symmetric_symbol(
        eta_zero,
        parent_zero,
        vneg(eta_phi),
        parent_phi,
    )
    assert tuple(
        trig_derivative(component) for component in wake_difference
    ) == forcing_difference

    return_sum = trig_symmetric_symbol(
        eta_sum,
        trig_vector_constant(equilibrium),
        vneg(eta_phi),
        parent_phi,
    )
    return_difference = trig_symmetric_symbol(
        eta_difference,
        wake_difference,
        eta_phi,
        parent_phi,
    )
    cubic_forcing = trig_vector_neg(
        trig_vector_add(return_sum, return_difference)
    )
    mean = trig_mean(trig_dot(parent_zero, cubic_forcing))
    ratio = qscale(mean, 2)  # one 2*pi period, divided by pi
    assert ratio[1:] == (F(0), F(0), F(0))
    return ratio[0]


def rational_circle(t):
    """Return (cos(phi),sin(phi)) for t=tan(phi/2)."""
    t = F(t)
    denominator = 1 + t * t
    return (
        (1 - t * t) / denominator,
        2 * t / denominator,
    )


def periodic_second_harmonic_primitive(value):
    """Zero-at-zero primitive of a mean-zero 0/2-harmonic polynomial."""
    sine_two = trig_scale(
        trig_multiply(TRIG_SIN, TRIG_COS),
        q(2),
    )
    cosine_two = trig_sub(
        trig_multiply(TRIG_COS, TRIG_COS),
        trig_multiply(TRIG_SIN, TRIG_SIN),
    )
    assert trig_mean(value) == ZERO
    sine_coefficient = qscale(
        trig_mean(trig_multiply(value, sine_two)), 2
    )
    cosine_coefficient = qscale(
        trig_mean(trig_multiply(value, cosine_two)), 2
    )
    primitive = trig_add(
        trig_scale(
            trig_sub(TRIG_ONE, cosine_two),
            qscale(sine_coefficient, F(1, 2)),
        ),
        trig_scale(
            sine_two,
            qscale(cosine_coefficient, F(1, 2)),
        ),
    )
    assert trig_derivative(primitive) == value
    return primitive


def physical_directed_kernel_components(target_t, partner_t):
    """Return the physical (sum,difference) coefficients divided by pi.

    Unlike raw_kernel_ratio(-t), this routine changes both the target mode
    and its selected left line.  It derives the persistent sum equilibrium
    and zero-initial periodic difference wake from the vector symbols.
    """
    target_cosine, target_sine = rational_circle(target_t)
    partner_cosine, partner_sine = rational_circle(partner_t)
    target_wavevector = (
        q(target_cosine), q(target_sine), INV_SQRT3,
    )
    partner_wavevector = (
        q(partner_cosine), q(partner_sine), INV_SQRT3,
    )
    target_parent = rotating_parent(target_cosine, target_sine)
    partner_parent = rotating_parent(partner_cosine, partner_sine)
    sum_wavevector = vadd(target_wavevector, partner_wavevector)
    difference_wavevector = vadd(
        target_wavevector, vneg(partner_wavevector)
    )

    sum_forcing = trig_symmetric_symbol(
        target_wavevector,
        target_parent,
        partner_wavevector,
        partner_parent,
    )
    assert all(
        set(component).issubset({(0, 0)})
        for component in sum_forcing
    )
    sum_forcing_constant = tuple(
        component.get((0, 0), ZERO) for component in sum_forcing
    )
    sum_norm_squared = vdot(sum_wavevector, sum_wavevector)
    sum_vertical_squared = qmul(sum_wavevector[2], sum_wavevector[2])
    assert sum_norm_squared[1:] == (F(0), F(0), F(0))
    assert sum_vertical_squared[1:] == (F(0), F(0), F(0))
    inertial_frequency_squared = (
        4 * sum_vertical_squared[0] / sum_norm_squared[0]
    )
    sum_equilibrium = tuple(
        qscale(component, F(1, inertial_frequency_squared))
        for component in kelvin_apply(
            sum_wavevector, sum_forcing_constant
        )
    )
    assert vadd(
        kelvin_apply(sum_wavevector, sum_equilibrium),
        sum_forcing_constant,
    ) == (ZERO, ZERO, ZERO)

    difference_forcing = trig_symmetric_symbol(
        target_wavevector,
        target_parent,
        vneg(partner_wavevector),
        partner_parent,
    )
    difference_wake = tuple(
        periodic_second_harmonic_primitive(component)
        for component in difference_forcing
    )

    sum_return = trig_symmetric_symbol(
        sum_wavevector,
        trig_vector_constant(sum_equilibrium),
        vneg(partner_wavevector),
        partner_parent,
    )
    difference_return = trig_symmetric_symbol(
        difference_wavevector,
        difference_wake,
        partner_wavevector,
        partner_parent,
    )
    ratios = []
    for return_term in (sum_return, difference_return):
        cubic_forcing = trig_vector_neg(return_term)
        mean = trig_mean(trig_dot(target_parent, cubic_forcing))
        ratio = qscale(mean, 2)
        assert ratio[1:] == (F(0), F(0), F(0))
        ratios.append(ratio[0])
    return tuple(ratios)


def physical_directed_kernel_ratio(target_t, partner_t):
    """Physical directed secular coefficient divided by pi."""
    return sum(physical_directed_kernel_components(target_t, partner_t))


def kernel_ratio(t):
    """Return K(t)/pi as an exact rational number."""
    t = F(t)
    return (
        -t**3
        * (t * t + 2)
        * (t**4 + 2 * t**3 + 3 * t * t + 2 * t + 4)
        / (1 + t * t) ** 5
    )


def check_general_secular_kernel():
    assert kernel_ratio(F(1)) == F(-9, 8)

    for t in (F(1, 10), F(1, 3), F(1, 2), F(1), F(2), F(7, 2)):
        assert raw_kernel_ratio(t) == kernel_ratio(t)
        assert raw_kernel_ratio(-t) == kernel_ratio(-t)
        assert physical_directed_kernel_ratio(F(0), t) == kernel_ratio(t)
        assert physical_directed_kernel_ratio(t, F(0)) == -kernel_ratio(t)
        omega_squared = 4 * (1 + t * t) / (t * t + 4)
        assert 1 < omega_squared < 4
        # The distance from the second-harmonic resonance tends to zero
        # only at the excluded antipodal endpoint t=infinity.
        assert 4 - omega_squared == F(12, 1) / (t * t + 4)
        assert kernel_ratio(t) < 0

    # raw_kernel_ratio(-t) still targets the phase-zero mode and is not
    # the physical reverse return.  At quarter separation the two numbers
    # differ exactly; check_quarter_reduced_system reconstructs the latter.
    assert kernel_ratio(F(-1)) == F(3, 8)
    physical_reverse_quarter_ratio = F(9, 8)
    assert kernel_ratio(F(-1)) != physical_reverse_quarter_ratio


def check_physical_antisymmetry_and_balanced_triple():
    # For the isolated two-parent shared-wake normal form, the order-four
    # energy identity forces physical K_ij=-K_ji.  Reconstruct several
    # pairs directly, including the quarter pair.
    for left, right in (
        (F(0), F(1)),
        (F(0), F(-1)),
        (F(1, 2), F(2)),
        (F(-3), F(-1)),
        (F(-1), F(1, 2)),
    ):
        forward = physical_directed_kernel_ratio(left, right)
        reverse = physical_directed_kernel_ratio(right, left)
        assert forward + reverse == 0

    assert physical_directed_kernel_ratio(F(0), F(1)) == F(-9, 8)
    assert physical_directed_kernel_ratio(F(1), F(0)) == F(9, 8)

    # An exact positive balanced three-ray state.  The phases are encoded
    # by half-angle parameters (-3,-1,1/2).  With K_ij measured in units
    # of pi, the positive vector below lies in the kernel of K.
    phases = (F(-3), F(-1), F(1, 2))
    k12 = physical_directed_kernel_ratio(phases[0], phases[1])
    k13 = physical_directed_kernel_ratio(phases[0], phases[2])
    k23 = physical_directed_kernel_ratio(phases[1], phases[2])
    assert k12 == F(-1746, 3125)
    assert k13 == F(227409, 78125000)
    assert k23 == F(-3861, 25000)
    weights = (k23, -k13, k12)
    if all(weight < 0 for weight in weights):
        weights = tuple(-weight for weight in weights)
    assert weights == (
        F(3861, 25000),
        F(227409, 78125000),
        F(1746, 3125),
    )
    assert all(weight > 0 for weight in weights)
    matrix = (
        (F(0), k12, k13),
        (-k12, F(0), k23),
        (-k13, -k23, F(0)),
    )
    for row in matrix:
        assert sum(entry * weight for entry, weight in zip(row, weights)) == 0

def check_fixed_ring_support_logic():
    # A nonzero sum of two unit vectors determines the unordered pair.
    # This exact rational-circle sample checks the reconstruction:
    # u,v are the two roots around midpoint s/2, with their difference
    # perpendicular to s and of fixed magnitude.
    points = []
    for t in (F(0), F(1, 3), F(1, 2), F(1), F(2), F(3)):
        denominator = 1 + t * t
        points.append(((1 - t * t) / denominator, 2 * t / denominator))

    sums = {}
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            output = (
                points[i][0] + points[j][0],
                points[i][1] + points[j][1],
            )
            if output != (F(0), F(0)):
                assert output not in sums
                sums[output] = (i, j)

    # No cubic aggregate-sign assertion is made here.  The exact
    # quarter-pair reverse calculation disproves the former one.


def main():
    check_pair_symbols()
    check_quarter_reduced_system()
    check_first_period_coefficient()
    check_general_secular_kernel()
    check_physical_antisymmetry_and_balanced_triple()
    check_fixed_ring_support_logic()
    print("PASS C155: full-period quarter-pair DACR coefficient is exact and < -21/8")
    print("PASS C156: fixed-ring quadratic null and corrected directed-cubic boundary")
    print("PASS AUDIT: quarter physical returns are -9*pi/8 and +9*pi/8; aggregate zero")
    print("PASS REPAIR: physical K is antisymmetric and an exact positive 3-ray balance exists")
    print("Thick multi-charge packet, finite-epsilon normal form, and wake graph remain OPEN")


if __name__ == "__main__":
    main()
