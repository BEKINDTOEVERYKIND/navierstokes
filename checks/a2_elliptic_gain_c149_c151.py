#!/usr/bin/env python3
"""Exact checks for C149--C151 on the existing A2 pump.

The Floquet first variation is evaluated by integrating a finite
trigonometric polynomial over one period.  Coefficients live in
Q(sqrt(2),sqrt(3)); no external symbolic package is used.  Floating point
appears only in the limiting lattice-approximation sanity check.  The
Q^-1 ledger is explicitly conditional and is not a time-coherence check.
"""

from fractions import Fraction as F
from math import factorial, sqrt


# Elements a+b*sqrt(2)+c*sqrt(3)+d*sqrt(6).
ZERO = (F(0), F(0), F(0), F(0))
ONE = (F(1), F(0), F(0), F(0))
SQRT2 = (F(0), F(1), F(0), F(0))
SQRT3 = (F(0), F(0), F(1), F(0))
SQRT6 = (F(0), F(0), F(0), F(1))


def q(value=0):
    return (F(value), F(0), F(0), F(0))


def qadd(left, right):
    return tuple(a + b for a, b in zip(left, right))


def qneg(value):
    return tuple(-entry for entry in value)


def qsub(left, right):
    return qadd(left, qneg(right))


def qmul(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e + 2 * b * f + 3 * c * g + 6 * d * h,
        a * f + b * e + 3 * c * h + 3 * d * g,
        a * g + c * e + 2 * b * h + 2 * d * f,
        a * h + d * e + b * g + c * f,
    )


def qscale(value, scalar):
    return tuple(F(scalar) * entry for entry in value)


def qis_zero(value):
    return value == ZERO


INV_SQRT2 = qscale(SQRT2, F(1, 2))
INV_SQRT3 = qscale(SQRT3, F(1, 3))
INV_SQRT6 = qscale(SQRT6, F(1, 6))


def vadd(left, right):
    return tuple(qadd(a, b) for a, b in zip(left, right))


def vneg(value):
    return tuple(qneg(entry) for entry in value)


def vscale(scalar, value):
    return tuple(qmul(scalar, entry) for entry in value)


def vdot(left, right):
    out = ZERO
    for a, b in zip(left, right):
        out = qadd(out, qmul(a, b))
    return out


def vcross(left, right):
    return (
        qsub(qmul(left[1], right[2]), qmul(left[2], right[1])),
        qsub(qmul(left[2], right[0]), qmul(left[0], right[2])),
        qsub(qmul(left[0], right[1]), qmul(left[1], right[0])),
    )


def rational_div(value, denominator):
    return qscale(value, F(1, denominator))


def leray_project(wavevector, value):
    norm_sq = vdot(wavevector, wavevector)
    assert norm_sq[1:] == (F(0), F(0), F(0))
    coefficient = rational_div(vdot(wavevector, value), norm_sq[0])
    return vadd(value, vscale(qneg(coefficient), wavevector))


def symmetric_symbol(p, amplitude_p, wavevector_q, amplitude_q):
    output = vadd(
        vscale(vdot(amplitude_p, wavevector_q), amplitude_q),
        vscale(vdot(amplitude_q, p), amplitude_p),
    )
    return leray_project(vadd(p, wavevector_q), output)


# Polynomial in sin(s), cos(s), with Q(sqrt2,sqrt3) coefficients.
def pconstant(value):
    return {} if qis_zero(value) else {(0, 0): value}


PZERO = {}
PONE = pconstant(ONE)
PSIN = {(1, 0): ONE}
PCOS = {(0, 1): ONE}


def padd(left, right):
    out = dict(left)
    for key, value in right.items():
        out[key] = qadd(out.get(key, ZERO), value)
        if qis_zero(out[key]):
            del out[key]
    return out


def pneg(value):
    return {key: qneg(coefficient) for key, coefficient in value.items()}


def psub(left, right):
    return padd(left, pneg(right))


def pmul(left, right):
    out = {}
    for (sin_left, cos_left), coefficient_left in left.items():
        for (sin_right, cos_right), coefficient_right in right.items():
            key = (sin_left + sin_right, cos_left + cos_right)
            value = qmul(coefficient_left, coefficient_right)
            out[key] = qadd(out.get(key, ZERO), value)
            if qis_zero(out[key]):
                del out[key]
    return out


def pscale(value, scalar):
    return {
        key: qmul(scalar, coefficient)
        for key, coefficient in value.items()
        if not qis_zero(qmul(scalar, coefficient))
    }


def psum(values):
    out = {}
    for value in values:
        out = padd(out, value)
    return out


def matrix_add(left, right):
    return [
        [padd(left[row][column], right[row][column])
         for column in range(len(left[0]))]
        for row in range(len(left))
    ]


def matrix_neg(value):
    return [[pneg(entry) for entry in row] for row in value]


def matrix_scale(value, scalar):
    return [[pscale(entry, scalar) for entry in row] for row in value]


def matrix_multiply(left, right):
    return [
        [
            psum(
                pmul(left[row][inner], right[inner][column])
                for inner in range(len(right))
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def transpose(value):
    return [list(row) for row in zip(*value)]


def outer(left, right):
    return [
        [pmul(left[row], right[column]) for column in range(len(right))]
        for row in range(len(left))
    ]


def polynomial_matrix_from_constants(value):
    return [[pconstant(entry) for entry in row] for row in value]


def trig_moment_over_pi(sin_power, cos_power):
    """Return pi^{-1} integral_0^{2pi} sin^i cos^j."""
    if sin_power % 2 or cos_power % 2:
        return F(0)
    a = sin_power // 2
    b = cos_power // 2
    return F(
        2 * factorial(2 * a) * factorial(2 * b),
        4 ** (a + b) * factorial(a) * factorial(b) * factorial(a + b),
    )


def integrate_over_period_div_pi(value):
    out = ZERO
    for (sin_power, cos_power), coefficient in value.items():
        out = qadd(
            out,
            qscale(
                coefficient,
                trig_moment_over_pi(sin_power, cos_power),
            ),
        )
    return out


def floquet_first_variation_div_pi(detuning):
    """Exact coefficient M'(0)/pi for r=1/sqrt(3)+detuning*eps."""
    sin_2s = pscale(pmul(PSIN, PCOS), q(2))
    cos_2s = psub(pmul(PCOS, PCOS), pmul(PSIN, PSIN))

    zero = {}
    j_matrix = [
        [zero, pneg(PONE), zero],
        [PONE, zero, zero],
        [zero, zero, zero],
    ]
    s_rotated = [
        [pneg(sin_2s), pneg(cos_2s), zero],
        [pneg(cos_2s), sin_2s, zero],
        [zero, zero, zero],
    ]

    eta_zero = [PONE, zero, pconstant(INV_SQRT3)]
    eta_one = [
        pscale(psub(PONE, cos_2s), q(F(1, 2))),
        pmul(PSIN, PCOS),
        pconstant(detuning),
    ]
    eta_dot = psum(
        pmul(eta_zero[index], eta_one[index]) for index in range(3)
    )

    # |eta_0|^2=4/3.  Differentiate P=eta eta^T/|eta|^2.
    p_zero = matrix_scale(outer(eta_zero, eta_zero), q(F(3, 4)))
    d_projection_first = matrix_scale(
        matrix_add(outer(eta_one, eta_zero), outer(eta_zero, eta_one)),
        q(F(3, 4)),
    )
    d_projection_second = [
        [
            pmul(
                entry,
                pscale(eta_dot, q(F(9, 8))),
            )
            for entry in row
        ]
        for row in outer(eta_zero, eta_zero)
    ]
    d_projection = matrix_add(
        d_projection_first,
        matrix_neg(d_projection_second),
    )

    c_star = matrix_add(
        matrix_scale(j_matrix, q(-2)),
        matrix_scale(matrix_multiply(p_zero, j_matrix), q(2)),
    )
    c_one = matrix_add(
        matrix_scale(matrix_multiply(d_projection, j_matrix), q(2)),
        matrix_add(
            matrix_neg(s_rotated),
            matrix_scale(matrix_multiply(p_zero, s_rotated), q(2)),
        ),
    )

    identity = [
        [PONE if row == column else zero for column in range(3)]
        for row in range(3)
    ]
    c_star_sq = matrix_multiply(c_star, c_star)
    one_minus_cosine = psub(PONE, PCOS)
    exponential_plus = matrix_add(
        identity,
        matrix_add(
            [[pmul(PSIN, entry) for entry in row] for row in c_star],
            [[pmul(one_minus_cosine, entry) for entry in row]
             for row in c_star_sq],
        ),
    )
    exponential_minus = matrix_add(
        identity,
        matrix_add(
            [[pneg(pmul(PSIN, entry)) for entry in row] for row in c_star],
            [[pmul(one_minus_cosine, entry) for entry in row]
             for row in c_star_sq],
        ),
    )

    # E1=(0,1,0), E2=(-1/2,0,sqrt(3)/2) at epsilon=0.
    frame = polynomial_matrix_from_constants(
        [
            [ZERO, q(F(-1, 2))],
            [ONE, ZERO],
            [ZERO, qscale(SQRT3, F(1, 2))],
        ]
    )
    integrand = matrix_multiply(
        matrix_multiply(exponential_minus, c_one),
        exponential_plus,
    )
    transverse = matrix_multiply(
        matrix_multiply(transpose(frame), integrand),
        frame,
    )
    return [
        [integrate_over_period_div_pi(entry) for entry in row]
        for row in transverse
    ]


def integer_dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def check_c149_elliptic_cocycle():
    n_vector = (1, 1, 1)
    r1 = (1, -1, 0)
    r2 = (0, 1, -1)
    r3 = (-1, 0, 1)
    roots = (r1, r2, r3)
    assert tuple(sum(root[index] for root in roots) for index in range(3)) == (0, 0, 0)
    assert all(integer_dot(root, n_vector) == 0 for root in roots)
    assert all(integer_dot(root, root) == 2 for root in roots)

    # In the ex,ey,ez frame, N cross Hess(f_delta) is exactly
    # sqrt(3)(2+delta) A_epsilon.  The two nonzero entries simplify to
    # -3sqrt(3) and sqrt(3)(1+2delta).
    ex = vscale(INV_SQRT2, tuple(q(a + b) for a, b in zip(r1, r2)))
    ey = vscale(INV_SQRT6, tuple(q(a - b) for a, b in zip(r1, r2)))
    ez = vscale(INV_SQRT3, tuple(q(entry) for entry in n_vector))
    frame = (ex, ey, ez)
    assert vcross(ex, ey) == vneg(ez)
    n_q = tuple(q(entry) for entry in n_vector)
    root_q = [tuple(q(entry) for entry in root) for root in roots]

    for delta in (F(1, 5), F(1, 2), F(4, 5)):
        weights = (F(1), F(1), delta)
        columns = []
        for basis_vector in frame:
            hessian_action = (ZERO, ZERO, ZERO)
            for weight, root in zip(weights, root_q):
                hessian_action = vadd(
                    hessian_action,
                    vscale(q(-weight), vscale(vdot(root, basis_vector), root)),
                )
            columns.append(vcross(n_q, hessian_action))
        represented = [
            [vdot(frame[row], columns[column]) for column in range(3)]
            for row in range(3)
        ]
        expected = [
            [ZERO, qscale(SQRT3, -3), ZERO],
            [qscale(SQRT3, 1 + 2 * delta), ZERO, ZERO],
            [ZERO, ZERO, ZERO],
        ]
        assert represented == expected

        epsilon = (1 - delta) / (2 + delta)
        assert 0 < epsilon < F(1, 2)
        assert (1 - 2 * epsilon) / (1 + epsilon) == delta

    detunings = (
        ZERO,
        qscale(SQRT3, F(1, 6)),  # 1/(2sqrt(3))
        qscale(SQRT3, F(-1, 12)),
        qscale(SQRT3, F(1, 4)),
    )
    for detuning in detunings:
        actual = floquet_first_variation_div_pi(detuning)
        root_term = qmul(qscale(SQRT3, F(3, 2)), detuning)
        expected = [
            [ZERO, qadd(q(F(-15, 8)), root_term)],
            [qneg(qadd(q(F(3, 8)), root_term)), ZERO],
        ]
        assert actual == expected

    c_zero = qscale(SQRT3, F(1, 6))
    split = floquet_first_variation_div_pi(c_zero)
    assert split == [
        [ZERO, q(F(-9, 8))],
        [q(F(-9, 8)), ZERO],
    ]

    # The open first-order tongue is exactly the intersection where both
    # factors in (1.19) are positive.
    lower = -1 / (4 * sqrt(3))
    upper = 5 / (4 * sqrt(3))
    c_numeric = 1 / (2 * sqrt(3))
    assert lower < c_numeric < upper
    # Divide the 9*pi/8 multiplier splitting by the 2*pi period.
    predicted_dimensionless_exponent = F(9, 8) / 2
    assert predicted_dimensionless_exponent == F(9, 16)


def check_c150_charge_darkness():
    for delta in range(0, 12):
        for margin in range(1, 8):
            m_zero = 3 * delta + margin
            retained_positive = set(range(m_zero - delta, m_zero + delta + 1))
            retained = retained_positive | {-value for value in retained_positive}
            quadratic = {
                left + right
                for left in retained
                for right in retained
            }
            assert retained.isdisjoint(quadratic)

            # The zero-charge pump preserves every charge exactly.
            pump_shifts = {charge + 0 for charge in retained}
            assert pump_shifts == retained

            # Cubic return is support-allowed: m0+m0-m0=m0.
            assert m_zero + m_zero - m_zero in retained

    # The concrete C144 packet has normal charges 4q,...,5q-1.
    for q_value in range(1, 20):
        band = set(range(4 * q_value, 5 * q_value))
        real_band = band | {-charge for charge in band}
        quadratic = {left + right for left in real_band for right in real_band}
        assert real_band.isdisjoint(quadratic)


def float_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def float_scale(scalar, value):
    return tuple(scalar * entry for entry in value)


def float_dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def float_project(wavevector, value):
    coefficient = float_dot(wavevector, value) / float_dot(wavevector, wavevector)
    return float_add(value, float_scale(-coefficient, wavevector))


def float_normalize(value):
    return float_scale(1 / sqrt(float_dot(value, value)), value)


def float_symmetric_symbol(p, amplitude_p, wavevector_q, amplitude_q):
    output = float_add(
        float_scale(float_dot(amplitude_p, wavevector_q), amplitude_q),
        float_scale(float_dot(amplitude_q, p), amplitude_p),
    )
    return float_project(float_add(p, wavevector_q), output)


def float_norm(value):
    return sqrt(float_dot(value, value))


def check_c151_derivative_amplified_return():
    r = INV_SQRT3
    k = (ONE, ZERO, r)
    ell = (ZERO, ONE, r)
    amplitude_a = (
        qscale(INV_SQRT2, F(1, 2)),
        INV_SQRT2,
        qneg(qmul(SQRT3, qscale(INV_SQRT2, F(1, 2)))),
    )
    amplitude_b = (
        INV_SQRT2,
        qscale(INV_SQRT2, F(1, 2)),
        qneg(qmul(SQRT3, qscale(INV_SQRT2, F(1, 2)))),
    )
    assert vdot(k, amplitude_a) == ZERO
    assert vdot(ell, amplitude_b) == ZERO
    assert vdot(amplitude_a, amplitude_a) == ONE
    assert vdot(amplitude_b, amplitude_b) == ONE

    wake = symmetric_symbol(k, amplitude_a, ell, amplitude_b)
    expected_wake = (
        q(F(3, 10)),
        q(F(3, 10)),
        qscale(SQRT3, F(-3, 10)),
    )
    assert wake == expected_wake

    # Self, opposite, and difference channels vanish, so the real four-mode
    # first derivative has only the sum wake and its reality partner.
    zero_vector = (ZERO, ZERO, ZERO)
    assert symmetric_symbol(k, amplitude_a, k, amplitude_a) == zero_vector
    assert symmetric_symbol(ell, amplitude_b, ell, amplitude_b) == zero_vector
    assert symmetric_symbol(k, amplitude_a, vneg(ell), amplitude_b) == zero_vector
    opposite_k = vadd(
        vscale(vdot(amplitude_a, vneg(k)), amplitude_a),
        vscale(vdot(amplitude_a, k), amplitude_a),
    )
    opposite_ell = vadd(
        vscale(vdot(amplitude_b, vneg(ell)), amplitude_b),
        vscale(vdot(amplitude_b, ell), amplitude_b),
    )
    assert opposite_k == zero_vector
    assert opposite_ell == zero_vector
    assert symmetric_symbol(
        vneg(k), amplitude_a, vneg(ell), amplitude_b
    ) == vneg(wake)

    returned = symmetric_symbol(
        vadd(k, ell),
        wake,
        vneg(ell),
        amplitude_b,
    )
    expected_return = vscale(qscale(INV_SQRT2, F(1, 2)), wake)
    assert returned == expected_return
    assert vdot(amplitude_a, vneg(returned)) == q(F(-9, 40))
    negative_return_symbol = symmetric_symbol(
        vneg(vadd(k, ell)),
        wake,
        ell,
        amplitude_b,
    )
    # The negative wake Fourier coefficient is +i*wake, so the extra
    # phase i together with the Euler factor -i yields the same real
    # returned coefficient as the positive-frequency path.
    assert negative_return_symbol == vneg(returned)

    # Actual A2-lattice sequences with equal normal charge approach both
    # epsilon=0 limiting rays and preserve the nonzero limiting symbol.
    # This does not certify exact finite-epsilon Floquet eigenpackets.
    # Coordinates here are in the ex,ey,ez frame; the integer coefficients
    # are B_m and D_m.
    target_k = (1.0, 0.0, 1 / sqrt(3))
    target_ell = (0.0, 1.0, 1 / sqrt(3))
    target_a = (
        1 / (2 * sqrt(2)),
        1 / sqrt(2),
        -sqrt(3) / (2 * sqrt(2)),
    )
    target_b = (
        1 / sqrt(2),
        1 / (2 * sqrt(2)),
        -sqrt(3) / (2 * sqrt(2)),
    )
    target_wake = (0.3, 0.3, -0.3 * sqrt(3))
    errors = []
    for m in (1000, 10000, 100000):
        b_index = round(3 * m / sqrt(2))
        d_index = round(3 * m / sqrt(6))
        p_m = (
            b_index * sqrt(2) / (3 * m),
            0.0,
            1 / sqrt(3),
        )
        q_m = (
            0.0,
            d_index * sqrt(6) / (3 * m),
            1 / sqrt(3),
        )
        a_m = float_normalize(float_project(p_m, target_a))
        b_m = float_normalize(float_project(q_m, target_b))
        wake_m = float_symmetric_symbol(p_m, a_m, q_m, b_m)
        return_m = float_symmetric_symbol(
            float_add(p_m, q_m),
            wake_m,
            float_scale(-1, q_m),
            b_m,
        )
        errors.append(float_norm(float_add(wake_m, float_scale(-1, target_wake))))
        assert float_dot(a_m, float_scale(-1, return_m)) < -0.2
        assert abs(float_dot(p_m, a_m)) < 1e-12
        assert abs(float_dot(q_m, b_m)) < 1e-12
    assert errors[-1] < 1e-5

    # Explicit fixed-epsilon retuning of the two same-charge directions.
    # This checks only the lattice-direction asymptotic, not a PDE packet.
    for epsilon in (0.01, 0.05, 0.1):
        r_epsilon = 1 / sqrt(3) + epsilon / (2 * sqrt(3))
        alpha_epsilon = sqrt((1 + epsilon) / (1 - epsilon))
        direction_errors = []
        for m in (1000, 10000, 100000):
            q_scale = sqrt(3) * m / r_epsilon
            b_index = round(sqrt(3) * m / (sqrt(2) * r_epsilon))
            d_index = round(alpha_epsilon * m / (sqrt(2) * r_epsilon))
            p_m = (
                sqrt(2) * b_index / q_scale,
                0.0,
                sqrt(3) * m / q_scale,
            )
            q_m = (
                0.0,
                sqrt(6) * d_index / q_scale,
                sqrt(3) * m / q_scale,
            )
            direction_errors.append(
                max(
                    float_norm(float_add(p_m, (-1.0, 0.0, -r_epsilon))),
                    float_norm(
                        float_add(
                            q_m,
                            (0.0, -alpha_epsilon, -r_epsilon),
                        )
                    ),
                )
            )
        assert direction_errors[-1] < direction_errors[0]
        assert direction_errors[-1] < 1e-5

    # Conditional coherent scalar ledger: an uncancelled, signed
    # Q^2 A^3 channel reaches the raw threshold A=Q^-1.  This arithmetic
    # does not prove temporal/Floquet coherence of the instantaneous symbol.
    for n in (2, 5, 11):
        q_value = n**8
        seed = F(1, n**16)
        threshold = F(1, q_value)
        assert seed == F(1, q_value**2)
        assert threshold / seed == q_value

        final_amplitude = n**10
        relative_cubic_scale = q_value**2 * final_amplitude**2
        assert relative_cubic_scale == n**36


def main():
    check_c149_elliptic_cocycle()
    check_c150_charge_darkness()
    check_c151_derivative_amplified_return()
    print("PASS C149: exact elliptic A2 Kelvin splitting and positive Floquet derivative")
    print("PASS C150: normal-charge quadratic darkness and cubic support grading")
    print("PASS C151: nonzero instantaneous Q^2 A^3 return; conditional Q^-1 ledger")
    print("Local principal cocycle only; nonlinear localized broad-band stage remains open")


if __name__ == "__main__":
    main()
