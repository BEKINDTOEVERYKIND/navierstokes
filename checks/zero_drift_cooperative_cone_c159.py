#!/usr/bin/env python3
"""Outward-rational cooperative-cone certificate for C159.

No floating-point arithmetic or third-party package is used in the proof
path.  Decimal endpoints are finite decimals (hence rational numbers), and
every interval operation is performed in a directed-rounding Context.
Taylor coefficients are merely reference data generated deterministically;
the checker independently encloses their ODE residuals.
"""

from decimal import (
    Context, Decimal as D, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN,
    localcontext,
)
from fractions import Fraction as F
from math import isqrt


PREC = 70
NEAR = Context(prec=PREC, rounding=ROUND_HALF_EVEN)
DOWN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)


class I:
    """Closed decimal-rational interval with directed operations."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo if isinstance(lo, D) else D(lo)
        self.hi = self.lo if hi is None else (hi if isinstance(hi, D) else D(hi))
        assert self.lo <= self.hi

    def __add__(self, other):
        other = other if isinstance(other, I) else I(other)
        return I(DOWN.add(self.lo, other.lo), UP.add(self.hi, other.hi))

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return I(other) + (-self)

    def __mul__(self, other):
        other = other if isinstance(other, I) else I(other)
        lows = [
            DOWN.multiply(a, b)
            for a in (self.lo, self.hi) for b in (other.lo, other.hi)
        ]
        highs = [
            UP.multiply(a, b)
            for a in (self.lo, self.hi) for b in (other.lo, other.hi)
        ]
        return I(min(lows), max(highs))

    __rmul__ = __mul__

    def reciprocal(self):
        assert not (self.lo <= 0 <= self.hi)
        return I(DOWN.divide(D(1), self.hi), UP.divide(D(1), self.lo))

    def __truediv__(self, other):
        other = other if isinstance(other, I) else I(other)
        return self * other.reciprocal()

    def __rtruediv__(self, other):
        return I(other) / self

    def square(self):
        if self.lo <= 0 <= self.hi:
            upper = max(UP.multiply(self.lo, self.lo), UP.multiply(self.hi, self.hi))
            return I(0, upper)
        lows = [DOWN.multiply(self.lo, self.lo), DOWN.multiply(self.hi, self.hi)]
        highs = [UP.multiply(self.lo, self.lo), UP.multiply(self.hi, self.hi)]
        return I(min(lows), max(highs))

    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def widen(self, radius):
        radius = radius if isinstance(radius, D) else D(radius)
        return I(DOWN.subtract(self.lo, radius), UP.add(self.hi, radius))


class Ball:
    """Decimal center-radius ball with an explicit correctly-rounded ulp charge."""

    __slots__ = ("center", "radius")

    def __init__(self, center, radius=D(0)):
        self.center = center if isinstance(center, D) else D(center)
        self.radius = radius if isinstance(radius, D) else D(radius)
        assert self.radius >= 0

    @staticmethod
    def ulp(value):
        if not value:
            return D(1).scaleb(-PREC - 10)
        return D(1).scaleb(value.adjusted() - PREC + 1)

    def __add__(self, other):
        other = other if isinstance(other, Ball) else Ball(other)
        center = NEAR.add(self.center, other.center)
        radius = UP.add(UP.add(self.radius, other.radius), self.ulp(center))
        return Ball(center, radius)

    __radd__ = __add__

    def __neg__(self):
        return Ball(-self.center, self.radius)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return Ball(other) + (-self)

    def __mul__(self, other):
        other = other if isinstance(other, Ball) else Ball(other)
        center = NEAR.multiply(self.center, other.center)
        radius = UP.add(
            UP.add(
                UP.multiply(abs(self.center), other.radius),
                UP.multiply(abs(other.center), self.radius),
            ),
            UP.multiply(self.radius, other.radius),
        )
        radius = UP.add(radius, self.ulp(center))
        return Ball(center, radius)

    __rmul__ = __mul__

    def abs_upper(self):
        return UP.add(abs(self.center), self.radius)


def as_ball(value):
    """Enclose a point or interval by a center-radius ball.

    For an interval we deliberately use its lower endpoint as center.  This
    avoids introducing an unaccounted midpoint rounding; the directed width
    is then an immediate valid radius.
    """
    if isinstance(value, Ball):
        return value
    if isinstance(value, I):
        return Ball(value.lo, UP.subtract(value.hi, value.lo))
    return Ball(value)


def decimal_fraction(value):
    numerator, denominator = value.as_integer_ratio()
    return F(numerator, denominator)


def fraction_interval(value):
    numerator = D(value.numerator)
    denominator = D(value.denominator)
    return I(DOWN.divide(numerator, denominator), UP.divide(numerator, denominator))


def sqrt_fraction_bound(value, digits=PREC):
    """Exact floor/ceiling decimal bounds for sqrt(Fraction(value))."""
    value = F(value)
    scale = 10 ** digits
    scaled_num = value.numerator * scale * scale
    floor_arg = scaled_num // value.denominator
    ceil_arg = -(-scaled_num // value.denominator)
    lower_int = isqrt(floor_arg)
    upper_int = isqrt(ceil_arg)
    if upper_int * upper_int < ceil_arg:
        upper_int += 1
    scale_decimal = D(scale)
    lower = DOWN.divide(D(lower_int), scale_decimal)
    upper = UP.divide(D(upper_int), scale_decimal)
    assert decimal_fraction(lower) ** 2 <= value <= decimal_fraction(upper) ** 2
    return I(lower, upper)


def interval_sqrt(value):
    assert value.lo >= 0
    lo = sqrt_fraction_bound(decimal_fraction(value.lo)).lo
    hi = sqrt_fraction_bound(decimal_fraction(value.hi)).hi
    return I(lo, hi)


def atan_alternating_bounds(x, terms=220):
    x = F(x)
    total = F(0)
    power = x
    sign = 1
    for index in range(terms):
        total += sign * power / (2 * index + 1)
        power *= x * x
        sign *= -1
    next_term = sign * power / (2 * terms + 1)
    return min(total, total + next_term), max(total, total + next_term)


def parameter_intervals():
    # Machin: pi=16 atan(1/5)-4 atan(1/239).
    a_lo, a_hi = atan_alternating_bounds(F(1, 5))
    b_lo, b_hi = atan_alternating_bounds(F(1, 239))
    pi_interval = I(
        fraction_interval(16 * a_lo - 4 * b_hi).lo,
        fraction_interval(16 * a_hi - 4 * b_lo).hi,
    )
    assert D("3.141592653589793238462643383279") < pi_interval.hi
    assert pi_interval.lo < D("3.141592653589793238462643383280")

    delta = I(D(4) / 5)
    root57 = sqrt_fraction_bound(F(57))
    discriminant = root57 / 5
    y_minus = (1 + delta * delta - discriminant) / (2 * delta * delta)
    y_plus = (1 + delta * delta + discriminant) / (2 * delta * delta)
    c_star = y_plus - y_minus
    d_star = y_plus - 1
    parameter = 1 - y_minus * d_star / c_star

    # Derivatives with respect to the energy E at E=0.
    y_minus_prime = (1 - 1 / discriminant) / (2 * delta)
    y_plus_prime = (1 + 1 / discriminant) / (2 * delta)
    c_prime = y_plus_prime - y_minus_prime
    d_prime = y_plus_prime
    parameter_prime = -(
        y_minus_prime * d_star + y_minus * d_prime
    ) / c_star + y_minus * d_star * c_prime / (c_star * c_star)

    # AGM plus the E/K sum.  c_0^2=m, hence the initial sum is m/2.
    agm_a = I(1)
    agm_b = interval_sqrt(1 - parameter)
    first_c = (agm_a - agm_b) / 2
    assert agm_b.lo > D("0.25")
    assert first_c.hi < D("0.5")
    e_sum = parameter / 2
    power_two = D(1)
    last_c = None
    for _ in range(12):
        c_value = (agm_a - agm_b) / 2
        e_sum = e_sum + power_two * c_value * c_value
        new_a = (agm_a + agm_b) / 2
        new_b = interval_sqrt(agm_a * agm_b)
        agm_a, agm_b = new_a, new_b
        power_two *= 2
        last_c = c_value
    # For positive AGM iterates,
    # c_(j+1)=c_j^2/(sqrt(a_j)+sqrt(b_j))^2.  Since b_j is increasing,
    # b_0>0 and c_1<1/2, the omitted weighted terms decrease at least
    # geometrically with ratio 1/2.  The factor 4 is deliberately twice the
    # resulting factor-2 tail bound.
    next_c = (agm_a - agm_b) / 2
    tail = I(0, UP.multiply(D(4) * power_two, next_c.square().hi))
    assert tail.hi < D("1e-35"), tail.hi
    e_over_k = 1 - (e_sum + tail)
    agm_limit = I(agm_b.lo, agm_a.hi)
    elliptic_k = pi_interval / (2 * agm_limit)
    period = (4 / (3 * delta)) * elliptic_k / interval_sqrt(c_star)

    k_log_derivative = (
        e_over_k / (2 * parameter * (1 - parameter))
        - 1 / (2 * parameter)
    )
    t_log_derivative = (
        k_log_derivative * parameter_prime - c_prime / (2 * c_star)
    )
    sqrt2 = sqrt_fraction_bound(F(2))
    carrier = I(F(378, 25).numerator) / F(378, 25).denominator
    beta = -t_log_derivative * carrier / (3 * sqrt2)

    assert period.lo > D("3.0361377700939")
    assert period.hi < D("3.0361377700941")
    assert beta.lo > D("2.1108753438878")
    assert beta.hi < D("2.1108753438881")
    assert UP.subtract(period.hi, period.lo) < D("1e-60")
    assert UP.subtract(beta.hi, beta.lo) < D("1e-60")
    assert UP.subtract(sqrt2.hi, sqrt2.lo) < D("1e-60")
    return period, beta, sqrt2


# Small polynomial helpers.  Coefficients can be Decimal points (reference
# generation) or intervals (certificate verification).
def p_add(left, right, zero):
    out = [zero for _ in range(max(len(left), len(right)))]
    for index, value in enumerate(left):
        out[index] = out[index] + value
    for index, value in enumerate(right):
        out[index] = out[index] + value
    return out


def p_scale(value, scalar):
    return [scalar * entry for entry in value]


def p_mul(left, right, zero, limit=None):
    size = len(left) + len(right) - 1
    if limit is not None:
        size = min(size, limit)
    out = [zero for _ in range(size)]
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            if i + j >= size:
                break
            out[i + j] = out[i + j] + a_value * b_value
    return out


def p_pad(value, size, zero):
    return value + [zero for _ in range(size - len(value))]


DELTA = D(4) / 5
R1 = (D(1), D(-1), D(0))
R2 = (D(0), D(1), D(-1))
RS = tuple(a + b for a, b in zip(R1, R2))
CN = ((D(0), D(-1), D(1)), (D(1), D(0), D(-1)), (D(-1), D(1), D(0)))
CARRIER = D(378) / 25


def point_rhs(coefficients, period, sqrt2, beta, limit):
    """Polynomial augmented phase,r=1/h,gamma RHS at point decimals."""
    zero = D(0)
    ca, sa, cb, sb, inverse_h, gamma = coefficients
    sab = p_add(p_mul(sa, cb, zero, limit), p_mul(ca, sb, zero, limit), zero)
    cab = p_add(p_mul(ca, cb, zero, limit), p_scale(p_mul(sa, sb, zero, limit), -1), zero)
    fa = p_add(p_scale(sa, -1), p_scale(sab, -DELTA), zero)
    fb = p_add(p_scale(sb, -1), p_scale(sab, -DELTA), zero)
    adot = p_scale(p_add(sb, p_scale(sab, DELTA), zero), 3)
    bdot = p_scale(p_add(sa, p_scale(sab, DELTA), zero), -3)
    hessian = [[
        p_add(
            p_add(p_scale(ca, -R1[i] * R1[j]), p_scale(cb, -R2[i] * R2[j]), zero),
            p_scale(cab, -DELTA * RS[i] * RS[j]), zero,
        ) for j in range(3)
    ] for i in range(3)]
    gradient = [p_add(p_scale(fa, R1[i]), p_scale(fb, R2[i]), zero) for i in range(3)]
    velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_scale(gradient[j], CN[i][j]), zero)
        velocity.append(entry)
    h_velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_scale(p_mul(hessian[i][j], velocity[j], zero, limit), D(1)), zero)
        h_velocity.append(entry)
    g_h_u = [zero]
    u_h_u = [zero]
    for i in range(3):
        g_h_u = p_add(g_h_u, p_mul(gradient[i], h_velocity[i], zero, limit), zero)
        u_h_u = p_add(u_h_u, p_mul(velocity[i], h_velocity[i], zero, limit), zero)
    inverse_h_sq = p_mul(inverse_h, inverse_h, zero, limit)
    return [
        p_scale(p_mul(sa, adot, zero, limit), -period),
        p_scale(p_mul(ca, adot, zero, limit), period),
        p_scale(p_mul(sb, bdot, zero, limit), -period),
        p_scale(p_mul(cb, bdot, zero, limit), period),
        p_scale(p_mul(inverse_h_sq, g_h_u, zero, limit), -2 * period),
        p_add(
            [period * 3 * sqrt2 * beta],
            p_scale(p_mul(u_h_u, inverse_h_sq, zero, limit), period * (2 * CARRIER / 3)),
            zero,
        ),
    ]


def generate_reference(period, beta, sqrt2, panels=32, order=40):
    """Generate exact-decimal reference polynomials; no claim uses generation accuracy."""
    with localcontext(NEAR):
        period_point = +(period.lo + period.hi) / 2
        beta_point = +(beta.lo + beta.hi) / 2
        sqrt2_point = +(sqrt2.lo + sqrt2.hi) / 2
        sqrt21 = +(D(21).sqrt())
        state = [D(-2) / 5, sqrt21 / 5, D(-2) / 5, -sqrt21 / 5, D(25) / 126, D(0)]
        step = D(1) / panels
        records = []
        for _ in range(panels):
            coefficients = [[D(0)] * (order + 1) for _ in range(6)]
            for index in range(6):
                coefficients[index][0] = state[index]
            for degree in range(order):
                short = [row[:degree + 1] for row in coefficients]
                rhs = point_rhs(short, period_point, sqrt2_point, beta_point, degree + 1)
                for index in range(6):
                    coefficients[index][degree + 1] = rhs[index][degree] / D(degree + 1)
            records.append(coefficients)
            powers = [D(1)]
            for _degree in range(order):
                powers.append(powers[-1] * step)
            state = [sum((a * b for a, b in zip(row, powers)), D(0)) for row in coefficients]
    return records, period_point, beta_point, sqrt2_point


def interval_polynomial(value):
    return [I(entry) for entry in value]


def interval_phase_rhs(coefficients, period_point):
    zero = I(0)
    ca, sa, cb, sb = coefficients
    sab = p_add(p_mul(sa, cb, zero), p_mul(ca, sb, zero), zero)
    adot = p_scale(p_add(sb, p_scale(sab, I(DELTA)), zero), I(3))
    bdot = p_scale(p_add(sa, p_scale(sab, I(DELTA)), zero), I(-3))
    return [
        p_scale(p_mul(sa, adot, zero), I(-period_point)),
        p_scale(p_mul(ca, adot, zero), I(period_point)),
        p_scale(p_mul(sb, bdot, zero), I(-period_point)),
        p_scale(p_mul(cb, bdot, zero), I(period_point)),
    ]


def interval_gamma_numerator(coefficients, period_point, sqrt2_point, beta_point):
    """Return h^2*(Pgamma'-const)-T*(2c/3)U.H.U as an interval polynomial."""
    zero = I(0)
    ca, sa, cb, sb = coefficients[:4]
    gamma = coefficients[5]
    sab = p_add(p_mul(sa, cb, zero), p_mul(ca, sb, zero), zero)
    cab = p_add(p_mul(ca, cb, zero), p_scale(p_mul(sa, sb, zero), I(-1)), zero)
    fa = p_add(p_scale(sa, I(-1)), p_scale(sab, I(-DELTA)), zero)
    fb = p_add(p_scale(sb, I(-1)), p_scale(sab, I(-DELTA)), zero)
    hessian = [[
        p_add(
            p_add(p_scale(ca, I(-R1[i] * R1[j])), p_scale(cb, I(-R2[i] * R2[j])), zero),
            p_scale(cab, I(-DELTA * RS[i] * RS[j])), zero,
        ) for j in range(3)
    ] for i in range(3)]
    gradient = [p_add(p_scale(fa, I(R1[i])), p_scale(fb, I(R2[i])), zero) for i in range(3)]
    velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_scale(gradient[j], I(CN[i][j])), zero)
        velocity.append(entry)
    h_value = [zero]
    for entry in gradient:
        h_value = p_add(h_value, p_mul(entry, entry, zero), zero)
    h_velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_mul(hessian[i][j], velocity[j], zero), zero)
        h_velocity.append(entry)
    u_h_u = [zero]
    for i in range(3):
        u_h_u = p_add(u_h_u, p_mul(velocity[i], h_velocity[i], zero), zero)
    gamma_derivative = [I(index + 1) * gamma[index + 1] for index in range(len(gamma) - 1)]
    gamma_derivative[0] = gamma_derivative[0] - I(period_point * 3 * sqrt2_point * beta_point)
    return p_add(
        p_mul(p_mul(h_value, h_value, zero), gamma_derivative, zero),
        p_scale(u_h_u, I(-period_point * (2 * CARRIER / 3))),
        zero,
    ), h_value


def ball_phase_rhs(coefficients, period_value):
    zero = Ball(0)
    period_ball = as_ball(period_value)
    ca, sa, cb, sb = coefficients
    sab = p_add(p_mul(sa, cb, zero), p_mul(ca, sb, zero), zero)
    adot = p_scale(p_add(sb, p_scale(sab, Ball(DELTA)), zero), Ball(3))
    bdot = p_scale(p_add(sa, p_scale(sab, Ball(DELTA)), zero), Ball(-3))
    return [
        p_scale(p_mul(sa, adot, zero), -period_ball),
        p_scale(p_mul(ca, adot, zero), period_ball),
        p_scale(p_mul(sb, bdot, zero), -period_ball),
        p_scale(p_mul(cb, bdot, zero), period_ball),
    ]


def ball_gamma_numerator(coefficients, period_value, sqrt2_value, beta_value):
    zero = Ball(0)
    period_ball = as_ball(period_value)
    sqrt2_ball = as_ball(sqrt2_value)
    beta_ball = as_ball(beta_value)
    ca, sa, cb, sb = coefficients[:4]
    gamma = coefficients[5]
    sab = p_add(p_mul(sa, cb, zero), p_mul(ca, sb, zero), zero)
    cab = p_add(p_mul(ca, cb, zero), p_scale(p_mul(sa, sb, zero), Ball(-1)), zero)
    fa = p_add(p_scale(sa, Ball(-1)), p_scale(sab, Ball(-DELTA)), zero)
    fb = p_add(p_scale(sb, Ball(-1)), p_scale(sab, Ball(-DELTA)), zero)
    hessian = [[
        p_add(
            p_add(p_scale(ca, Ball(-R1[i] * R1[j])), p_scale(cb, Ball(-R2[i] * R2[j])), zero),
            p_scale(cab, Ball(-DELTA * RS[i] * RS[j])), zero,
        ) for j in range(3)
    ] for i in range(3)]
    gradient = [p_add(p_scale(fa, Ball(R1[i])), p_scale(fb, Ball(R2[i])), zero) for i in range(3)]
    velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_scale(gradient[j], Ball(CN[i][j])), zero)
        velocity.append(entry)
    h_value = [zero]
    for entry in gradient:
        h_value = p_add(h_value, p_mul(entry, entry, zero), zero)
    h_velocity = []
    for i in range(3):
        entry = [zero]
        for j in range(3):
            entry = p_add(entry, p_mul(hessian[i][j], velocity[j], zero), zero)
        h_velocity.append(entry)
    u_h_u = [zero]
    for i in range(3):
        u_h_u = p_add(u_h_u, p_mul(velocity[i], h_velocity[i], zero), zero)
    gamma_derivative = [Ball(index + 1) * gamma[index + 1] for index in range(len(gamma) - 1)]
    gamma_derivative[0] = gamma_derivative[0] - period_ball * Ball(3) * sqrt2_ball * beta_ball
    return p_add(
        p_mul(p_mul(h_value, h_value, zero), gamma_derivative, zero),
        p_scale(u_h_u, -(period_ball * Ball(2 * CARRIER / 3))),
        zero,
    )


def ball_polynomial_abs_bound(coefficients, step):
    power = I(1)
    total = I(0)
    for coefficient in coefficients:
        total = total + I(0, coefficient.abs_upper()) * power
        power = power * I(step)
    return total.hi


def polynomial_abs_bound(coefficients, step):
    power = I(1)
    total = I(0)
    for coefficient in coefficients:
        total = total + I(0, coefficient.abs_upper()) * power
        power = power * I(step)
    return total.hi


def evaluate_polynomial(coefficients, argument):
    output = I(0)
    for coefficient in reversed(coefficients):
        output = output * argument + coefficient
    return output


class AD:
    __slots__ = ("value", "gradient")

    def __init__(self, value, gradient=None):
        self.value = value if isinstance(value, I) else I(value)
        self.gradient = [I(0)] * 4 if gradient is None else gradient

    def __add__(self, other):
        other = other if isinstance(other, AD) else AD(other)
        return AD(self.value + other.value, [a + b for a, b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return AD(-self.value, [-entry for entry in self.gradient])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return AD(other) + (-self)

    def __mul__(self, other):
        other = other if isinstance(other, AD) else AD(other)
        return AD(
            self.value * other.value,
            [a * other.value + self.value * b for a, b in zip(self.gradient, other.gradient)],
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, AD) else AD(other)
        denominator = other.value * other.value
        return AD(
            self.value / other.value,
            [(a * other.value - self.value * b) / denominator for a, b in zip(self.gradient, other.gradient)],
        )

    def __rtruediv__(self, other):
        return AD(other) / self


def vector_dot(left, right):
    output = left[0] * right[0]
    for a_value, b_value in zip(left[1:], right[1:]):
        output = output + a_value * b_value
    return output


def vector_cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def matrix_vector(matrix, vector):
    return [vector_dot(row, vector) for row in matrix]


def square_norm(vector):
    output = vector[0].square()
    for entry in vector[1:]:
        output = output + entry.square()
    return output


def phase_and_gamma_derivative_bounds(phase, period, beta, sqrt2):
    variables = []
    for index, value in enumerate(phase):
        gradient = [I(0)] * 4
        gradient[index] = I(1)
        variables.append(AD(value, gradient))
    ca, sa, cb, sb = variables
    sab = sa * cb + ca * sb
    adot = 3 * (sb + D(DELTA) * sab)
    bdot = -3 * (sa + D(DELTA) * sab)
    period_ad = AD(period)
    phase_rhs = [period_ad * (-sa * adot), period_ad * (ca * adot), period_ad * (-sb * bdot), period_ad * (cb * bdot)]
    logarithmic_norm = D("-1e100")
    for row, value in enumerate(phase_rhs):
        bound = value.gradient[row].hi
        for column, derivative in enumerate(value.gradient):
            if column != row:
                bound = UP.add(bound, derivative.abs_upper())
        logarithmic_norm = max(logarithmic_norm, bound)

    fa = -sa - D(DELTA) * sab
    fb = -sb - D(DELTA) * sab
    cab = ca * cb - sa * sb
    hessian = [[
        -ca * R1[i] * R1[j] - cb * R2[i] * R2[j] - D(DELTA) * cab * RS[i] * RS[j]
        for j in range(3)
    ] for i in range(3)]
    gradient = [fa * R1[i] + fb * R2[i] for i in range(3)]
    velocity = vector_cross((AD(1), AD(1), AD(1)), gradient)
    h_value = vector_dot(gradient, gradient)
    u_h_u = vector_dot(velocity, matrix_vector(hessian, velocity))
    gamma_rhs = period_ad * (AD(3 * sqrt2 * beta) + (2 * CARRIER / 3) * u_h_u / (h_value * h_value))
    gamma_lipschitz = D(0)
    for entry in gamma_rhs.gradient:
        gamma_lipschitz = UP.add(gamma_lipschitz, entry.abs_upper())
    return logarithmic_norm, gamma_lipschitz


def cone_coefficients(phase, gamma, period, beta, sqrt2, sqrt3):
    ca, sa, cb, sb = phase
    sab = sa * cb + ca * sb
    cab = ca * cb - sa * sb
    fa = -sa - I(DELTA) * sab
    fb = -sb - I(DELTA) * sab
    hessian = [[
        -ca * R1[i] * R1[j] - cb * R2[i] * R2[j] - I(DELTA) * cab * RS[i] * RS[j]
        for j in range(3)
    ] for i in range(3)]
    gradient = [fa * R1[i] + fb * R2[i] for i in range(3)]
    velocity = vector_cross((I(1), I(1), I(1)), gradient)
    h_value = square_norm(gradient)
    p_vector = [I(CARRIER / 3) * velocity[i] / h_value + gamma * gradient[i] for i in range(3)]
    normal = tuple(I(1) / sqrt3 for _ in range(3))
    m_value = sqrt3 * beta
    d_value = square_norm(p_vector)
    q_value = d_value + m_value * m_value
    tangent = vector_cross(p_vector, normal)
    velocity_gradient = [[
        sum((I(CN[i][inner]) * hessian[inner][j] for inner in range(3)), I(0))
        - sqrt2 * gradient[j]
        for j in range(3)
    ] for i in range(3)]
    p_s_t = vector_dot(p_vector, matrix_vector(velocity_gradient, tangent))
    l_p = vector_dot(normal, matrix_vector(velocity_gradient, p_vector))
    p_s_p = vector_dot(p_vector, matrix_vector(velocity_gradient, p_vector))
    b11 = m_value * l_p / d_value
    b22 = (2 * p_s_p + m_value * l_p) / d_value
    b21 = m_value * m_value * sqrt2 * I(CARRIER) / (q_value * d_value)
    b12 = (2 * m_value * p_s_t + sqrt2 * I(CARRIER) * (d_value - m_value * m_value)) / d_value
    b11, b12, b21, b22 = [period * value for value in (b11, b12, b21, b22)]
    return b12, b21, b11 + I(F(3, 20).numerator) / F(3, 20).denominator * b12, I(F(20, 3).numerator) / F(20, 3).denominator * b21 + b22


def exp_upper_42():
    # Exact rational Taylor sum plus a geometric tail (ratio <1/4 after k=168).
    total = F(1)
    term = F(1)
    for index in range(1, 169):
        term *= F(42, index)
        total += term
    tail = term * F(42, 169) * F(4, 3)
    return total + tail


def certify_path(records, period, beta, sqrt2):
    panels = len(records)
    order = len(records[0][0]) - 1
    step = D(1) / panels
    phase_residual_integral = D(0)
    gamma_residual_integral = D(0)
    jump_phase = D(0)
    jump_gamma = D(0)

    for panel, coefficients in enumerate(records):
        interval_coefficients = [interval_polynomial(row) for row in coefficients]
        ball_coefficients = [[Ball(entry) for entry in row] for row in coefficients]
        # Insert the full outward parameter intervals into the residual.
        # Thus no unproved point-parameter allowance enters the certificate.
        phase_rhs = ball_phase_rhs(ball_coefficients[:4], period)
        panel_phase_residual = D(0)
        for index in range(4):
            derivative = [Ball(degree + 1) * ball_coefficients[index][degree + 1] for degree in range(order)]
            residual = p_add(phase_rhs[index], p_scale(derivative, Ball(-1)), Ball(0))
            panel_phase_residual = max(panel_phase_residual, ball_polynomial_abs_bound(residual, step))
        phase_residual_integral = UP.add(
            phase_residual_integral,
            UP.multiply(panel_phase_residual, step),
        )

        numerator = ball_gamma_numerator(ball_coefficients, period, sqrt2, beta)
        h_lower = D("1e100")
        for subcell in range(16):
            sub_argument = I(step * D(subcell) / 16, step * D(subcell + 1) / 16)
            phase_panel = [evaluate_polynomial(interval_coefficients[index], sub_argument) for index in range(4)]
            ca, sa, cb, sb = phase_panel
            sab = sa * cb + ca * sb
            fa = -sa - I(DELTA) * sab
            fb = -sb - I(DELTA) * sab
            gradient = [fa * R1[index] + fb * R2[index] for index in range(3)]
            h_lower = min(h_lower, square_norm(gradient).lo)
        assert h_lower > D(3)
        gamma_residual = UP.divide(
            ball_polynomial_abs_bound(numerator, step),
            DOWN.multiply(h_lower, h_lower),
        )
        gamma_residual_integral = UP.add(gamma_residual_integral, UP.multiply(gamma_residual, step))

        if panel + 1 < panels:
            for index in range(6):
                endpoint = evaluate_polynomial(interval_coefficients[index], I(step))
                next_value = records[panel + 1][index][0]
                distance = (endpoint - I(next_value)).abs_upper()
                if index < 4:
                    jump_phase = UP.add(jump_phase, distance)
                elif index == 5:
                    jump_gamma = UP.add(jump_gamma, distance)

    # Enclose the algebraic initial sine explicitly.  The other relevant
    # initial coordinates (-2/5,-2/5 and gamma=0) are exact decimals.
    initial_sine = sqrt_fraction_bound(F(21)) / 5
    exact_initial_phase = [I(-2) / 5, initial_sine, I(-2) / 5, -initial_sine]
    initial_phase_error = max(
        (I(records[0][index][0]) - exact_initial_phase[index]).abs_upper()
        for index in range(4)
    )
    assert records[0][5][0] == 0

    phase_total = UP.add(
        UP.add(phase_residual_integral, jump_phase),
        initial_phase_error,
    )
    gamma_total = UP.add(gamma_residual_integral, jump_gamma)
    assert phase_total < D("7e-25"), (phase_total, phase_residual_integral, jump_phase)
    assert gamma_total < D("4e-17"), gamma_total
    assert exp_upper_42() < F(2 * 10**18)
    assert UP.multiply(D(2 * 10**18), phase_total) < D("2e-6")

    phase_radius = D("2e-6")
    gamma_radius = D("8e-4")
    evaluation_subdivisions = 32
    min_values = [D("1e100")] * 4
    maximum_mu = D("-1e100")
    maximum_gamma_lipschitz = D(0)
    sqrt3 = sqrt_fraction_bound(F(3))

    for coefficients in records:
        interval_coefficients = [interval_polynomial(row) for row in coefficients]
        for subcell in range(evaluation_subdivisions):
            left = step * D(subcell) / evaluation_subdivisions
            right = step * D(subcell + 1) / evaluation_subdivisions
            argument = I(left, right)
            phase = [evaluate_polynomial(interval_coefficients[index], argument).widen(phase_radius) for index in range(4)]
            gamma = evaluate_polynomial(interval_coefficients[5], argument).widen(gamma_radius)
            mu_value, gamma_lipschitz = phase_and_gamma_derivative_bounds(phase, period, beta, sqrt2)
            maximum_mu = max(maximum_mu, mu_value)
            maximum_gamma_lipschitz = max(maximum_gamma_lipschitz, gamma_lipschitz)
            values = cone_coefficients(phase, gamma, period, beta, sqrt2, sqrt3)
            for index, value in enumerate(values):
                min_values[index] = min(min_values[index], value.lo)

    assert maximum_mu < D(42)
    assert maximum_gamma_lipschitz < D(300)
    gamma_tube_bound = UP.add(
        gamma_total,
        UP.multiply(maximum_gamma_lipschitz, phase_radius),
    )
    assert gamma_tube_bound < gamma_radius
    assert min_values[0] > D(32)
    assert min_values[1] > D("0.9")
    assert min_values[2] > D("0.7"), min_values
    assert min_values[3] > D("0.2"), min_values
    return phase_total, gamma_total, maximum_mu, maximum_gamma_lipschitz, min_values


def exact_structural_checks():
    # The key Beltrami/axial identities behind the scalar reduction.
    n = (F(1), F(1), F(1))
    r1 = tuple(F(value) for value in R1)
    r2 = tuple(F(value) for value in R2)
    assert sum(n[i] * r1[i] for i in range(3)) == 0
    assert sum(n[i] * r2[i] for i in range(3)) == 0
    assert F(3) * F(126, 25) == F(378, 25)
    # Cooperative comparison uses w=(1,3/20): the two normalized rows are
    # B11+(3/20)B12 and (20/3)B21+B22.
    assert F(3, 20) * F(20, 3) == 1


def main():
    exact_structural_checks()
    period, beta, sqrt2 = parameter_intervals()
    print("C159 parameters: enclosed", flush=True)
    records, _period_point, _beta_point, _sqrt2_point = generate_reference(period, beta, sqrt2)
    print("C159 reference: generated", flush=True)
    result = certify_path(records, period, beta, sqrt2)
    phase_total, gamma_total, maximum_mu, gamma_lipschitz, margins = result
    print("C159 outward cooperative-cone certificate: PASS")
    print(f"phase total defect budget < {phase_total}")
    print(f"gamma total defect budget < {gamma_total}")
    print(f"mu_infinity < {maximum_mu}; gamma phase-Lipschitz < {gamma_lipschitz}")
    print("outward lower margins:", ", ".join(str(value) for value in margins))
    print("Conclusion: rho(M)>1, det(M)=1, and tr(M)>2 (no floating trace premise).")


if __name__ == "__main__":
    main()
