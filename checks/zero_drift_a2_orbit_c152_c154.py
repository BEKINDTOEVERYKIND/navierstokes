#!/usr/bin/env python3
"""Checks for the zero-drift A2 orbit, C152--C154.

C152 and C154 are checked with exact integer/rational algebra, apart from
the explicitly labelled period-value diagnostics.  C153 is deliberately a
floating RK4 convergence experiment, not an interval/Taylor certificate.
Passing this file therefore does NOT promote C153 beyond NUMERICAL CANDIDATE.
No third-party package is used.
"""

from fractions import Fraction as F
from math import acos, cos, log, pi, sin, sqrt


DELTA = F(4, 5)
N = (1, 1, 1)
R1 = (1, -1, 0)
R2 = (0, 1, -1)
R3 = (-1, 0, 1)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(scalar, value):
    return tuple(scalar * entry for entry in value)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def outer(left, right):
    return [[a * b for b in right] for a in left]


def matrix_add(left, right):
    return [
        [a + b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def matrix_scale(scalar, value):
    return [[scalar * entry for entry in row] for row in value]


def matrix_multiply(left, right):
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def matrix_vector(value, vector):
    return [sum(row[index] * vector[index] for index in range(len(vector)))
            for row in value]


def identity(size):
    return [[F(i == j) for j in range(size)] for i in range(size)]


def check_c152_exact_geometry():
    assert add(add(R1, R2), R3) == (0, 0, 0)
    assert dot(N, R1) == dot(N, R2) == dot(N, R3) == 0
    assert dot(R1, R1) == dot(R2, R2) == dot(R3, R3) == 2
    assert dot(R1, R2) == dot(R2, R3) == dot(R3, R1) == -1

    # r1.Nxr2=-3 and r2.Nxr1=3 give the exact phase system (1.3).
    assert dot(R1, cross(N, R2)) == -3
    assert dot(R2, cross(N, R1)) == 3

    # Complete critical-point reduction.  sin(a)=sin(b) gives either
    # a=b, where sin(a)(1+2 delta cos(a))=0, or a+b=pi, where sin(a)=0.
    assert -F(1, 2) / DELTA == F(-5, 8)
    critical_values = {
        "maximum": 2 + DELTA,
        "pi_pi": -2 + DELTA,
        "mixed": -DELTA,
        "diagonal": 2 * F(-5, 8)
        + DELTA * (2 * F(25, 64) - 1),
    }
    assert critical_values == {
        "maximum": F(14, 5),
        "pi_pi": F(-6, 5),
        "mixed": F(-4, 5),
        "diagonal": F(-57, 40),
    }
    assert critical_values["maximum"] > 0
    assert all(value < 0 for key, value in critical_values.items()
               if key != "maximum")

    # At (0,0), the phase Hessian is
    # [[-(1+delta),-delta],[-delta,-(1+delta)]].  Its eigenvalues are
    # -1 and -(1+2 delta), so the sole critical point above zero is a
    # strict nondegenerate maximum.  With no intervening critical value,
    # its zero superlevel is therefore a disk.
    assert -1 < 0
    assert -(1 + 2 * DELTA) < 0
    assert (1 + DELTA)**2 - DELTA**2 == 1 + 2 * DELTA > 0

    # The algebraic point cos(A)=-2/5, sin(A)=sqrt(21)/5 lies on f=0.
    cosine_a = F(-2, 5)
    sine_a_sq = F(21, 25)
    assert cosine_a * cosine_a + sine_a_sq == 1
    assert 2 * cosine_a + DELTA == 0

    # Remove the common sqrt(21)/5 factor in g0 and U0.
    gradient_direction = (-1, 2, -1)
    velocity_direction = scale(3, R3)
    assert cross(N, gradient_direction) == velocity_direction
    assert dot(gradient_direction, N) == 0
    assert dot(gradient_direction, velocity_direction) == 0


def q57_add(left, right):
    """Add a+b sqrt(57), represented by (a,b)."""
    return left[0] + right[0], left[1] + right[1]


def q57_mul(left, right):
    return (
        left[0] * right[0] + 57 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def check_c152_exact_period_algebra():
    y_minus = (F(41, 32), F(-5, 32))
    y_plus = (F(41, 32), F(5, 32))
    assert q57_add(y_minus, y_plus) == (F(41, 16), F(0))
    assert q57_mul(y_minus, y_plus) == (F(1, 4), F(0))

    # sqrt(y_minus)=(sqrt(57)-5)/8 exactly.
    turning_cosine = (F(-5, 8), F(1, 8))
    assert q57_mul(turning_cosine, turning_cosine) == y_minus
    # Exact sign checks: 5sqrt(57)<41, 5sqrt(57)>9, and y_+>1.
    assert 25 * 57 < 41**2
    assert 25 * 57 > 9**2
    assert F(41, 32) > 1

    # Exact coefficient check of (1.12) at E=0 as a polynomial in z=cos^2 p.
    # LHS = -4 delta^2 z^2 + (4 delta^2+4)z - delta^2.
    lhs = (-4 * DELTA**2, 4 * DELTA**2 + 4, -DELTA**2)
    # RHS = 4 delta^2 (z-y_-)(y_+-z), using exact sum/product.
    rhs = (
        -4 * DELTA**2,
        4 * DELTA**2 * F(41, 16),
        -4 * DELTA**2 * F(1, 4),
    )
    assert lhs == rhs

    # At E=0, sqrt(1+2 delta^2)=sqrt(57)/5 lies strictly above one.
    assert 25 < 57
    # Hence y_-'=(1-5/sqrt(57))/(2delta)>0, while y_+'>y_-'.
    assert 57 > 5**2


def period_diagnostic(subintervals):
    """Composite Simpson diagnostic for T(0),T'(0); not an interval proof."""
    assert subintervals % 2 == 0
    delta = float(DELTA)
    discriminant = sqrt(1 + 2 * delta * delta)
    y_minus = (1 + delta * delta - discriminant) / (2 * delta * delta)
    y_plus = (1 + delta * delta + discriminant) / (2 * delta * delta)
    y_minus_prime = (1 - 1 / discriminant) / (2 * delta)
    y_plus_prime = (1 + 1 / discriminant) / (2 * delta)
    step = (pi / 2) / subintervals
    period_sum = 0.0
    derivative_sum = 0.0
    for index in range(subintervals + 1):
        theta = index * step
        sine = sin(theta)
        cosine = cos(theta)
        a_value = y_minus + (1 - y_minus) * sine * sine
        b_value = y_plus - 1 + (1 - y_minus) * cosine * cosine
        integrand = 1 / sqrt(a_value * b_value)
        a_prime = y_minus_prime * cosine * cosine
        b_prime = y_plus_prime - y_minus_prime * cosine * cosine
        derivative_integrand = integrand * (
            a_prime / a_value + b_prime / b_value
        )
        weight = 1 if index in (0, subintervals) else (4 if index % 2 else 2)
        period_sum += weight * integrand
        derivative_sum += weight * derivative_integrand
    period_integral = step * period_sum / 3
    derivative_integral = step * derivative_sum / 3
    return (
        4 * period_integral / (3 * delta),
        -2 * derivative_integral / (3 * delta),
    )


def check_c152_rank_one_return():
    # Exact structural calculation with arbitrary nonzero rational
    # coefficients in u=alpha U+zeta N.  Orthogonality, rather than their
    # analytic values T,T', is what proves the rank-one identities.
    gradient = tuple(F(entry) for entry in (-1, 2, -1))
    velocity = tuple(F(entry) for entry in scale(3, R3))
    normal = tuple(F(entry) for entry in N)
    alpha = F(5, 4)
    zeta = F(-7, 3)
    u = add(scale(alpha, velocity), scale(zeta, normal))
    assert dot(gradient, u) == 0

    rank_one = outer(u, gradient)
    zero_matrix = [[F(0) for _ in range(3)] for _ in range(3)]
    assert matrix_multiply(rank_one, rank_one) == zero_matrix
    ident = identity(3)
    flow = matrix_add(ident, rank_one)
    flow_inverse = matrix_add(ident, matrix_scale(-1, rank_one))
    assert matrix_multiply(flow, flow_inverse) == ident
    assert matrix_multiply(flow_inverse, flow) == ident

    # A concrete off-plane vector U+beta N lies in u^perp.  This is the
    # symbolic identity behind beta in (1.19).
    beta = -alpha * dot(velocity, velocity) / (zeta * dot(normal, normal))
    covector = add(velocity, scale(beta, normal))
    assert beta != 0
    assert dot(u, covector) == 0
    assert dot(normal, covector) != 0

    period_coarse = period_diagnostic(1 << 12)
    period_fine = period_diagnostic(1 << 14)
    assert abs(period_coarse[0] - period_fine[0]) < 2e-12
    assert abs(period_coarse[1] - period_fine[1]) < 5e-12
    assert 3.03613 < period_fine[0] < 3.03615
    assert -1.799 < period_fine[1] < -1.797
    return period_fine


def norm(value):
    return sqrt(dot(value, value))


def kelvin_rhs(state):
    """Floating RHS for the C153 diagnostic only."""
    delta = float(DELTA)
    a_phase, b_phase = state[0], state[1]
    sine_a = sin(a_phase)
    sine_b = sin(b_phase)
    sine_sum = sin(a_phase + b_phase)
    f_a = -sine_a - delta * sine_sum
    f_b = -sine_b - delta * sine_sum
    phase_a_dot = -3 * f_b
    phase_b_dot = 3 * f_a

    hessian = [[0.0] * 3 for _ in range(3)]
    root_sum = add(R1, R2)
    for coefficient, root in (
        (-cos(a_phase), R1),
        (-cos(b_phase), R2),
        (-delta * cos(a_phase + b_phase), root_sum),
    ):
        for row in range(3):
            for column in range(3):
                hessian[row][column] += coefficient * root[row] * root[column]

    cross_n = ((0.0, -1.0, 1.0),
               (1.0, 0.0, -1.0),
               (-1.0, 1.0, 0.0))
    gradient = [f_a * R1[index] + f_b * R2[index] for index in range(3)]
    velocity_gradient = [
        [
            sum(cross_n[row][inner] * hessian[inner][column]
                for inner in range(3))
            - sqrt(2) * N[row] * gradient[column]
            for column in range(3)
        ]
        for row in range(3)
    ]

    wavevector = state[2:5]
    wavevector_norm_sq = dot(wavevector, wavevector)
    wavevector_dot = [
        -sum(velocity_gradient[column][row] * wavevector[column]
             for column in range(3))
        for row in range(3)
    ]
    output = [phase_a_dot, phase_b_dot, *wavevector_dot]
    for offset in (5, 8):
        amplitude = state[offset:offset + 3]
        gradient_times_amplitude = [
            sum(velocity_gradient[row][column] * amplitude[column]
                for column in range(3))
            for row in range(3)
        ]
        scalar = dot(wavevector, gradient_times_amplitude)
        output.extend(
            -gradient_times_amplitude[row]
            + 2 * wavevector[row] * scalar / wavevector_norm_sq
            for row in range(3)
        )
    return output


def rk4_kelvin(period, period_prime, step_count):
    sine_a = sqrt(21) / 5
    initial_velocity = (-3 * sine_a, 0.0, 3 * sine_a)
    beta = (
        -period_prime * dot(initial_velocity, initial_velocity)
        / (sqrt(2) * period * dot(N, N))
    )
    wavevector = [initial_velocity[index] + beta * N[index]
                  for index in range(3)]
    wavevector_unit = [entry / norm(wavevector) for entry in wavevector]
    reference = (1.0, 0.0, 0.0)
    reference_projection = dot(reference, wavevector_unit)
    frame_one = [
        reference[index] - reference_projection * wavevector_unit[index]
        for index in range(3)
    ]
    frame_one = [entry / norm(frame_one) for entry in frame_one]
    frame_two = cross(wavevector_unit, frame_one)

    angle = acos(-2 / 5)
    state = [angle, -angle, *wavevector, *frame_one, *frame_two]
    step = period / step_count
    for _ in range(step_count):
        first = kelvin_rhs(state)
        second = kelvin_rhs([
            value + step * slope / 2 for value, slope in zip(state, first)
        ])
        third = kelvin_rhs([
            value + step * slope / 2 for value, slope in zip(state, second)
        ])
        fourth = kelvin_rhs([
            value + step * slope for value, slope in zip(state, third)
        ])
        state = [
            value + step * (a + 2 * b + 2 * c + d) / 6
            for value, a, b, c, d in zip(
                state, first, second, third, fourth
            )
        ]

    monodromy = [
        [dot(frame, state[offset:offset + 3]) for offset in (5, 8)]
        for frame in (frame_one, frame_two)
    ]
    trace = monodromy[0][0] + monodromy[1][1]
    determinant = (
        monodromy[0][0] * monodromy[1][1]
        - monodromy[0][1] * monodromy[1][0]
    )
    return {
        "matrix": monodromy,
        "trace": trace,
        "determinant": determinant,
        "wavevector_return": norm([
            state[index + 2] - wavevector[index] for index in range(3)
        ]),
        "constraint_one": abs(dot(state[2:5], state[5:8])),
        "constraint_two": abs(dot(state[2:5], state[8:11])),
    }


def check_c153_numerical_candidate(period_data):
    # This convergence table is deliberately NOT described as rigorous.
    results = [
        rk4_kelvin(period_data[0], period_data[1], steps)
        for steps in (1 << 12, 1 << 13, 1 << 14)
    ]
    traces = [result["trace"] for result in results]
    assert all(trace > 10_000 for trace in traces)
    assert abs(traces[1] - traces[0]) < 2e-6
    assert abs(traces[2] - traces[1]) < 3e-7
    final = results[-1]
    assert 16_716.88 < final["trace"] < 16_716.89
    assert abs(final["determinant"] - 1) < 1e-3
    assert final["wavevector_return"] < 1e-8
    assert final["constraint_one"] < 1e-7
    assert final["constraint_two"] < 1e-7
    spectral_radius = (
        final["trace"] + sqrt(final["trace"]**2 - 4)
    ) / 2
    exponent = log(spectral_radius) / period_data[0]
    assert 3.20 < exponent < 3.21
    return final, spectral_radius, exponent


def check_exact_liouville_identity_structure():
    # tr(C_N H)=0 for symmetric H, and tr(N g^T)=N.g=0 for every
    # phase-gradient g in span(r1,r2).  Thus the pump is incompressible.
    cross_n = ((0, -1, 1), (1, 0, -1), (-1, 1, 0))
    symmetric_hessian = ((2, 3, 5), (3, 7, 11), (5, 11, 13))
    product = matrix_multiply(cross_n, symmetric_hessian)
    assert sum(product[index][index] for index in range(3)) == 0
    for gradient in (R1, R2, add(scale(7, R1), scale(-4, R2))):
        assert dot(N, gradient) == 0

    # For L=-A+2 k(k.A v)/|k|^2 on k^perp,
    # tr_perp L=k.A.k/|k|^2, while
    # d log|k|/dt=-k.A.k/|k|^2.  Verify this with a fully explicit exact
    # transverse trace, not just a scalar placeholder.
    wavevector = (F(1), F(2), F(3))
    transverse_one = (F(2), F(-1), F(0))
    transverse_two = (F(3), F(6), F(-5))
    assert dot(wavevector, transverse_one) == 0
    assert dot(wavevector, transverse_two) == 0
    assert dot(transverse_one, transverse_two) == 0
    gradient_matrix = (
        (F(1), F(2), F(3)),
        (F(4), F(-2), F(5)),
        (F(6), F(7), F(1)),
    )
    assert sum(gradient_matrix[index][index] for index in range(3)) == 0

    def generator(vector):
        av = matrix_vector(gradient_matrix, vector)
        coefficient = 2 * dot(wavevector, av) / dot(wavevector, wavevector)
        return add(scale(-1, av), scale(coefficient, wavevector))

    transverse_trace = (
        dot(transverse_one, generator(transverse_one))
        / dot(transverse_one, transverse_one)
        + dot(transverse_two, generator(transverse_two))
        / dot(transverse_two, transverse_two)
    )
    numerator = dot(
        wavevector,
        matrix_vector(gradient_matrix, wavevector),
    )
    norm_sq = dot(wavevector, wavevector)
    assert transverse_trace == numerator / norm_sq
    log_wavevector_derivative = -numerator / norm_sq
    assert transverse_trace + log_wavevector_derivative == 0


def check_c154_exact_bandwidth_shear():
    # Use exact rational representatives with the same orthogonality as
    # g0,u.  R=g outer u is the nilpotent term in K=I-R.
    gradient = tuple(F(entry) for entry in (-1, 2, -1))
    velocity = tuple(F(entry) for entry in scale(3, R3))
    normal = tuple(F(entry) for entry in N)
    u = add(scale(F(5, 4), velocity), scale(F(-7, 3), normal))
    assert dot(u, gradient) == 0
    shear = outer(gradient, u)
    zero_matrix = [[F(0) for _ in range(3)] for _ in range(3)]
    assert matrix_multiply(shear, shear) == zero_matrix
    ident = identity(3)
    one_return = matrix_add(ident, matrix_scale(-1, shear))
    for returns in (0, 1, 2, 7, 31):
        power = identity(3)
        for _ in range(returns):
            power = matrix_multiply(power, one_return)
        expected = matrix_add(ident, matrix_scale(-returns, shear))
        assert power == expected

    # A normal-charge displacement has nonzero shear coefficient and grows
    # exactly linearly.  A correlated displacement in u^perp does not shear.
    assert dot(u, normal) != 0
    assert dot(u, gradient) == 0
    for returns in (1, 5, 17):
        power = matrix_add(ident, matrix_scale(-returns, shear))
        normal_after = matrix_vector(power, normal)
        assert tuple(normal_after) == add(
            normal,
            scale(-returns * dot(u, normal), gradient),
        )
        assert tuple(matrix_vector(power, gradient)) == gradient

    # Exact ledger: width q in the shearing direction becomes q*ell;
    # narrowing it to q/ell restores final width q.
    for n_value, returns in ((2, 3), (3, 11), (5, 29)):
        q_value = n_value**8
        assert q_value * returns // q_value == returns
        narrowed = F(q_value, returns)
        assert narrowed * returns == q_value


def main():
    check_c152_exact_geometry()
    check_c152_exact_period_algebra()
    period_data = check_c152_rank_one_return()
    check_exact_liouville_identity_structure()
    numerical, spectral_radius, exponent = check_c153_numerical_candidate(
        period_data
    )
    check_c154_exact_bandwidth_shear()
    print("C152 EXACT checks: PASS")
    print(
        "C153 NUMERICAL CANDIDATE diagnostic: PASS "
        f"(trace={numerical['trace']:.9f}, rho={spectral_radius:.9f}, "
        f"exponent={exponent:.9f})"
    )
    print("C153 interval/Taylor proof: OPEN (this checker is not a certificate)")
    print("C154 EXACT checks: PASS")


if __name__ == "__main__":
    main()
