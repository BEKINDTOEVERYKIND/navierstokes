#!/usr/bin/env python3
"""Outward certificate for quantitative off-ray C159 dominated cones.

This checker does four things, without floating-point arithmetic or third
party packages.

1. It reruns the complete C159 phase/covector certificate.
2. It fattens every one of the 2048 C192 cells by explicit phase, gamma,
   vertical-charge, and horizontal-charge boxes.
3. It proves strict finite-horizon cone invariance and >3000 coefficient
   gain for both the forward cone and the sign-reflected inverse cone.
4. It checks that the C194 action-angle bounds put q^(-1/4), and also the
   proposed q^(-1/12)/q^(-1/3) two-width packets, inside those fat boxes at
   explicit stage thresholds.

The result is a principal phase-space dominated-cone certificate.  It does
not construct a canonical invariant bundle, one closed periodic phase, a
real finite-frequency packet, a viscous solution, or a nonlinear return.
"""

from decimal import Decimal as D
from fractions import Fraction as F

import c159_stable_line_c193 as stable
import strong_zero_drift_gain_c192 as gain
import zero_drift_cooperative_cone_c159 as base


# Extra physical tube around C159's already-certified reference tube.
EXTRA_PHASE = D("1e-7")
EXTRA_GAMMA = D("1e-5")
EXTRA_BETA = D("1e-7")
EXTRA_CHARGE = D("1e-5")
EXTRA_ENERGY = D("1e-11")

CENTRAL_PHASE_RADIUS = D("2e-6")
CENTRAL_GAMMA_RADIUS = D("8e-4")

CONE_LO = D("0.137")
CONE_HI = D("0.2")
IMAGE_LO = D("0.138")
IMAGE_HI = D("0.194")


def exact_general_block_formula_checks():
    """Check C159 (2.3) without the zero-level symmetry of S.

    Away from f=0 the horizontal block S need not be symmetric, but it is
    still trace free.  The load-bearing formula only needs tr S=0.  Compute
    the moving-basis Kelvin generator directly, including E1' and E2'.  The
    zero-level B21 row acquires the necessary skew correction

        m {t.S.p-p.S.t}/(qD),  t=p x n,

    while the other three C159 rows are unchanged.  Compare all four rows
    for several exact rational inputs with b!=c.
    """

    def dot(left, right):
        return sum((a * b for a, b in zip(left, right)), F(0))

    def mv(matrix, vector):
        return tuple(dot(row, vector) for row in matrix)

    def add(left, right):
        return tuple(a + b for a, b in zip(left, right))

    def scale(scalar, vector):
        return tuple(scalar * value for value in vector)

    cases = (
        (F(2), F(-1), F(3), F(1), F(4), F(-2), F(5)),
        (F(-3, 2), F(5, 3), F(7, 4), F(-2), F(1, 3), F(5), F(-4)),
        (F(4, 3), F(7, 5), F(-5, 2), F(3, 2), F(-7, 3), F(2, 5), F(9, 4)),
    )
    for p1, p2, m, a, b, c_value, ell1 in cases:
        # S=[[a,b],[c,-a]] and l=(ell1,ell2) are deliberately nonsymmetric.
        ell2 = F(11, 7) - ell1
        p = (p1, p2)
        tangent = (p2, -p1)
        s_matrix = ((a, b), (c_value, -a))
        ell = (ell1, ell2)
        d_value = dot(p, p)
        q_value = d_value + m * m
        p_dot = add(scale(-1, mv(tuple(zip(*s_matrix)), p)), scale(-m, ell))
        d_dot = 2 * dot(p, p_dot)
        q_dot = d_dot

        k = (p1, p2, m)
        e1 = (-m * p1 / q_value, -m * p2 / q_value, d_value / q_value)
        e2 = (p2, -p1, F(0))
        e1_dot = (
            -m * (p_dot[0] * q_value - p1 * q_dot) / q_value**2,
            -m * (p_dot[1] * q_value - p2 * q_dot) / q_value**2,
            (d_dot * q_value - d_value * q_dot) / q_value**2,
        )
        e2_dot = (p_dot[1], -p_dot[0], F(0))
        full_a = (
            (a, b, F(0)),
            (c_value, -a, F(0)),
            (ell1, ell2, F(0)),
        )

        def kelvin(vector):
            av = mv(full_a, vector)
            return add(scale(-1, av), scale(2 * dot(k, av) / q_value, k))

        direct = []
        for vector, derivative in ((e1, e1_dot), (e2, e2_dot)):
            residual = add(kelvin(vector), scale(-1, derivative))
            direct.append((
                dot(e1, residual) / dot(e1, e1),
                dot(e2, residual) / dot(e2, e2),
            ))
        # direct[column][row] -> standard B[row][column].
        direct_entries = (
            direct[0][0], direct[1][0], direct[0][1], direct[1][1]
        )
        lp = dot(ell, p)
        psp = dot(p, mv(s_matrix, p))
        pst = dot(p, mv(s_matrix, tangent))
        tsp = dot(tangent, mv(s_matrix, p))
        minus_l_t = -dot(ell, tangent)
        claimed_entries = (
            m * lp / d_value,
            (2 * m * pst + minus_l_t * (d_value - m * m)) / d_value,
            (
                m * (tsp - pst) + m * m * minus_l_t
            ) / (q_value * d_value),
            (2 * psp + m * lp) / d_value,
        )
        assert direct_entries == claimed_entries


def coefficientwise_general_block_formula_check():
    """Prove all four general moving-frame rows coefficient by coefficient.

    A tiny exact multivariate polynomial/rational-function ring independently
    repeats the direct E1/E2 differentiation with symbolic
    p1,p2,m,a,b,c,l1,l2 and S=[[a,b],[c,-a]].  Equality is tested only after
    cross multiplication and exact coefficient collection; numerical
    specialization is not a premise.
    """

    variable_count = 8
    zero_monomial = (0,) * variable_count

    def clean(value):
        return {key: coefficient for key, coefficient in value.items() if coefficient}

    def poly_constant(value):
        value = F(value)
        return {} if value == 0 else {zero_monomial: value}

    def poly_variable(index):
        power = [0] * variable_count
        power[index] = 1
        return {tuple(power): F(1)}

    def poly_add(left, right):
        out = dict(left)
        for key, coefficient in right.items():
            out[key] = out.get(key, F(0)) + coefficient
        return clean(out)

    def poly_neg(value):
        return {key: -coefficient for key, coefficient in value.items()}

    def poly_mul(left, right):
        out = {}
        for left_key, left_coefficient in left.items():
            for right_key, right_coefficient in right.items():
                key = tuple(a + b for a, b in zip(left_key, right_key))
                out[key] = out.get(key, F(0)) + left_coefficient * right_coefficient
        return clean(out)

    class Rat:
        __slots__ = ("numerator", "denominator")

        def __init__(self, numerator, denominator=None):
            if isinstance(numerator, Rat):
                self.numerator = numerator.numerator
                self.denominator = numerator.denominator
                return
            self.numerator = (
                numerator if isinstance(numerator, dict) else poly_constant(numerator)
            )
            self.denominator = (
                poly_constant(1) if denominator is None else denominator
            )

        def __add__(self, other):
            other = Rat(other)
            return Rat(
                poly_add(
                    poly_mul(self.numerator, other.denominator),
                    poly_mul(other.numerator, self.denominator),
                ),
                poly_mul(self.denominator, other.denominator),
            )

        __radd__ = __add__

        def __neg__(self):
            return Rat(poly_neg(self.numerator), self.denominator)

        def __sub__(self, other):
            return self + (-Rat(other))

        def __rsub__(self, other):
            return Rat(other) - self

        def __mul__(self, other):
            other = Rat(other)
            return Rat(
                poly_mul(self.numerator, other.numerator),
                poly_mul(self.denominator, other.denominator),
            )

        __rmul__ = __mul__

        def __truediv__(self, other):
            other = Rat(other)
            return Rat(
                poly_mul(self.numerator, other.denominator),
                poly_mul(self.denominator, other.numerator),
            )

        def __rtruediv__(self, other):
            return Rat(other) / self

        def __pow__(self, power):
            assert power >= 0
            out = Rat(1)
            for _ in range(power):
                out = out * self
            return out

        def coefficientwise_equal(self, other):
            other = Rat(other)
            return poly_mul(self.numerator, other.denominator) == poly_mul(
                other.numerator, self.denominator
            )

    def variable(index):
        return Rat(poly_variable(index))

    def dot(left, right):
        return sum((a * b for a, b in zip(left, right)), Rat(0))

    def mv(matrix, vector):
        return tuple(dot(row, vector) for row in matrix)

    def add(left, right):
        return tuple(a + b for a, b in zip(left, right))

    def scale(scalar, vector):
        return tuple(scalar * value for value in vector)

    p1, p2, m, a, b, c_value, ell1, ell2 = (
        variable(index) for index in range(variable_count)
    )
    p = (p1, p2)
    tangent = (p2, -p1)
    s_matrix = ((a, b), (c_value, -a))
    ell = (ell1, ell2)
    d_value = dot(p, p)
    q_value = d_value + m * m
    p_dot = add(scale(-1, mv(tuple(zip(*s_matrix)), p)), scale(-m, ell))
    d_dot = 2 * dot(p, p_dot)
    q_dot = d_dot

    k = (p1, p2, m)
    e1 = (-m * p1 / q_value, -m * p2 / q_value, d_value / q_value)
    e2 = (p2, -p1, Rat(0))
    e1_dot = (
        -m * (p_dot[0] * q_value - p1 * q_dot) / q_value**2,
        -m * (p_dot[1] * q_value - p2 * q_dot) / q_value**2,
        (d_dot * q_value - d_value * q_dot) / q_value**2,
    )
    e2_dot = (p_dot[1], -p_dot[0], Rat(0))
    full_a = (
        (a, b, Rat(0)),
        (c_value, -a, Rat(0)),
        (ell1, ell2, Rat(0)),
    )

    def kelvin(vector):
        av = mv(full_a, vector)
        return add(scale(-1, av), scale(2 * dot(k, av) / q_value, k))

    direct = []
    for vector, derivative in ((e1, e1_dot), (e2, e2_dot)):
        residual = add(kelvin(vector), scale(-1, derivative))
        direct.append((
            dot(e1, residual) / dot(e1, e1),
            dot(e2, residual) / dot(e2, e2),
        ))
    direct_entries = (
        direct[0][0], direct[1][0], direct[0][1], direct[1][1]
    )
    lp = dot(ell, p)
    psp = dot(p, mv(s_matrix, p))
    pst = dot(p, mv(s_matrix, tangent))
    tsp = dot(tangent, mv(s_matrix, p))
    minus_l_t = -dot(ell, tangent)
    claimed_entries = (
        m * lp / d_value,
        (2 * m * pst + minus_l_t * (d_value - m * m)) / d_value,
        (m * (tsp - pst) + m * m * minus_l_t) / (q_value * d_value),
        (2 * psp + m * lp) / d_value,
    )
    assert all(
        actual.coefficientwise_equal(expected)
        for actual, expected in zip(direct_entries, claimed_entries)
    )


def exact_a2_skew_identity_checks():
    """Verify the A2 reduction of the off-level B21 skew correction.

    For t=p x n and S=C_N H,

        t.S.p-p.S.t = 2 sqrt(3) f |p|^2.

    Multiplying by sqrt(3) and using t=(p x N)/sqrt(3) gives the rational
    identity checked below.
    """

    def dot(left, right):
        return sum((a * b for a, b in zip(left, right)), F(0))

    def cross(left, right):
        return (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )

    def mv(matrix, vector):
        return tuple(dot(row, vector) for row in matrix)

    normal = (F(1), F(1), F(1))
    r1 = tuple(F(value) for value in (1, -1, 0))
    r2 = tuple(F(value) for value in (0, 1, -1))
    root_sum = tuple(a + b for a, b in zip(r1, r2))
    cross_n = tuple(tuple(F(value) for value in row) for row in base.CN)
    for ca, cb, cab in (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    ):
        hessian = [[
            -ca * r1[i] * r1[j]
            - cb * r2[i] * r2[j]
            - F(4, 5) * cab * root_sum[i] * root_sum[j]
            for j in range(3)
        ] for i in range(3)]
        s_matrix = [[
            sum((cross_n[i][ell] * hessian[ell][j] for ell in range(3)), F(0))
            for j in range(3)
        ] for i in range(3)]
        energy = ca + cb + F(4, 5) * cab
        for p in (
            tuple(F(value) for value in (1, -1, 0)),
            tuple(F(value) for value in (1, 1, -2)),
            tuple(F(value) for value in (2, -3, 1)),
        ):
            tangent_n = cross(p, normal)
            lhs = dot(tangent_n, mv(s_matrix, p)) - dot(
                p, mv(s_matrix, tangent_n)
            )
            assert lhs == 6 * energy * dot(p, p)


def generalized_raw_coefficients(
    phase, gamma, period, beta, sqrt2, sqrt3, horizontal_charge, energy
):
    """C159 (2.3), with beta and p.V allowed to range independently.

    On a nearby invariant level E=f, write V=N x grad f and

        p = (c_h/(3 |grad f|^2)) V + gamma grad f,
        c_h = p.V.

    The C159 derivation uses c_h, not the full three-dimensional invariant
    k.U=c_h-sqrt(6)mE.  Off the zero level the exact A2 skew identity above
    reduces the additional B21 row to 2 sqrt(3) m E/q.  E is bounded as an
    invariant, rather than through the much wider nonphysical independent
    phase-validation box.
    """

    ca, sa, cb, sb = phase
    sab = sa * cb + ca * sb
    cab = ca * cb - sa * sb
    fa = -sa - base.I(base.DELTA) * sab
    fb = -sb - base.I(base.DELTA) * sab
    hessian = [[
        -ca * base.R1[i] * base.R1[j]
        - cb * base.R2[i] * base.R2[j]
        - base.I(base.DELTA) * cab * base.RS[i] * base.RS[j]
        for j in range(3)
    ] for i in range(3)]
    gradient = [fa * base.R1[i] + fb * base.R2[i] for i in range(3)]
    velocity = base.vector_cross((base.I(1), base.I(1), base.I(1)), gradient)
    h_value = base.square_norm(gradient)
    p_vector = [
        horizontal_charge * velocity[i] / (3 * h_value) + gamma * gradient[i]
        for i in range(3)
    ]
    normal = tuple(base.I(1) / sqrt3 for _ in range(3))
    m_value = sqrt3 * beta
    d_value = base.square_norm(p_vector)
    q_value = d_value + m_value * m_value
    tangent = base.vector_cross(p_vector, normal)
    velocity_gradient = [[
        sum(
            (base.I(base.CN[i][inner]) * hessian[inner][j] for inner in range(3)),
            base.I(0),
        ) - sqrt2 * gradient[j]
        for j in range(3)
    ] for i in range(3)]
    p_s_t = base.vector_dot(
        p_vector, base.matrix_vector(velocity_gradient, tangent)
    )
    l_p = base.vector_dot(
        normal, base.matrix_vector(velocity_gradient, p_vector)
    )
    p_s_p = base.vector_dot(
        p_vector, base.matrix_vector(velocity_gradient, p_vector)
    )
    b11 = m_value * l_p / d_value
    b22 = (2 * p_s_p + m_value * l_p) / d_value
    b21 = (
        m_value * m_value * sqrt2 * horizontal_charge / (q_value * d_value)
        + 2 * sqrt3 * m_value * energy / q_value
    )
    b12 = (
        2 * m_value * p_s_t
        + sqrt2 * horizontal_charge * (d_value - m_value * m_value)
    ) / d_value
    return tuple(period * value for value in (b11, b12, b21, b22))


def fattened_cells():
    """Return all outward B boxes after rerunning the C159 proof."""

    exact_general_block_formula_checks()
    coefficientwise_general_block_formula_check()
    exact_a2_skew_identity_checks()
    base.exact_structural_checks()
    period, beta, sqrt2 = base.parameter_intervals()
    records, *_ = base.generate_reference(period, beta, sqrt2)
    base.certify_path(records, period, beta, sqrt2)

    sqrt3 = base.sqrt_fraction_bound(base.F(3))
    beta_box = beta.widen(EXTRA_BETA)
    charge_box = base.I(base.CARRIER).widen(EXTRA_CHARGE)
    energy_box = base.I(-EXTRA_ENERGY, EXTRA_ENERGY)
    panel_step = D(1) / len(records)
    subdivisions = 64
    entries = []
    minima = [D("1e100")] * 4
    for coefficients in records:
        polynomials = [base.interval_polynomial(row) for row in coefficients]
        for subcell in range(subdivisions):
            left = panel_step * D(subcell) / subdivisions
            right = panel_step * D(subcell + 1) / subdivisions
            argument = base.I(left, right)
            phase = [
                base.evaluate_polynomial(polynomials[index], argument).widen(
                    CENTRAL_PHASE_RADIUS + EXTRA_PHASE
                )
                for index in range(4)
            ]
            gamma = base.evaluate_polynomial(polynomials[5], argument).widen(
                CENTRAL_GAMMA_RADIUS + EXTRA_GAMMA
            )
            assert gamma.abs_upper() < D(1)
            cell = generalized_raw_coefficients(
                phase,
                gamma,
                period,
                beta_box,
                sqrt2,
                sqrt3,
                charge_box,
                energy_box,
            )
            assert cell[1].lo > 0 and cell[2].lo > 0
            for index, value in enumerate(cell):
                minima[index] = min(minima[index], value.lo)
            entries.append(cell)
    assert len(entries) == 2048
    return entries, minima, period, beta


def reflected_inverse_cell(cell):
    """Generator for R Phi^{-1} R, R=diag(1,-1), in reverse time."""

    b11, b12, b21, b22 = cell
    return -b11, b12, b21, -b22


def lower_cone_gain(cells):
    """Lower z1 gain at the worst (smallest) positive cone slope."""

    dt_lo = base.DOWN.divide(D(1), D(len(cells)))
    dt_hi = base.UP.divide(D(1), D(len(cells)))
    z1, z2 = D(1), CONE_LO
    maximum_alpha = D(0)
    for cell in cells:
        z1, z2, alpha = gain.lower_step(z1, z2, cell, dt_lo, dt_hi)
        maximum_alpha = max(maximum_alpha, alpha)
    return z1, z2, maximum_alpha


def riccati_cone_image(cells):
    """Outward image of the complete input slope interval."""

    dt_lo = base.DOWN.divide(D(1), D(len(cells)))
    dt_hi = base.UP.divide(D(1), D(len(cells)))
    lower = base.I(CONE_LO)
    upper = base.I(CONE_HI)
    for b11, b12, b21, b22 in cells:
        # For r>=0 and B in an interval box,
        # l21+(l22-u11)r-u12 r^2 <= r' <=
        # u21+(u22-l11)r-l12 r^2.
        lower = stable.interval_taylor_step(
            lower,
            b21.lo,
            base.DOWN.subtract(b22.lo, b11.hi),
            b12.hi,
            dt_lo,
            dt_hi,
        )
        upper = stable.interval_taylor_step(
            upper,
            b21.hi,
            base.UP.subtract(b22.hi, b11.lo),
            b12.lo,
            dt_lo,
            dt_hi,
        )
    return lower.lo, upper.hi


def exponential_and_cone_geometry_checks():
    """Exact gain and separation arithmetic for arbitrary cone-line pairs."""

    # Rational Taylor upper bound for exp(8), as in C192.
    term = F(1)
    exponential_upper = F(1)
    for order in range(1, 33):
        term *= F(8, order)
        exponential_upper += term
    exponential_upper += term * F(8, 33) * F(4, 3)
    assert exponential_upper < 3000

    # At comparison-section times the fat tube has 5<|k|<6. Thus a coefficient
    # slope in [137/1000,1/5] has physical orthonormal slope in the safer
    # interval [13/20,6/5].  For y in this interval, y+1/y is maximized at
    # 13/20 and equals 569/260.  Hence for positive y_u,y_s,
    #
    # sin(angle)^2 >= 4/[(y_u+1/y_u)(y_s+1/y_s)],
    #
    # and sin(angle)>=520/569>9/10.
    lower_y = F(13, 20)
    upper_y = F(6, 5)
    lower_sum = lower_y + 1 / lower_y
    upper_sum = upper_y + 1 / upper_y
    assert lower_sum == F(569, 260)
    assert upper_sum == F(61, 30)
    assert lower_sum > upper_sum
    angle_sine = F(520, 569)
    assert angle_sine > F(9, 10)
    projector_norm = 1 / angle_sine
    assert projector_norm == F(569, 520)
    assert projector_norm < F(10, 9)

    # The coordinate determinant of a Kelvin block is D_start/D_end,
    # where D=|P_n k|^2.  The section tube has 3<|P_n k|<4, so det<16/9<2.
    # Since the first coordinate expands by >3000 throughout the cone, the
    # associated projective Mobius derivative is below 2/3000^2.
    graph_contraction = F(2, 3000**2)
    assert F(16, 9) < 2
    assert graph_contraction < F(1, 4_000_000)

    # Between section times the physical Kelvin generator has norm at most
    # 3||DU||<=18.  Since T<76/25 and e<3, one whole-block propagator and
    # its inverse have physical operator norm below the explicit fixed
    # factor 3^55. Thus O(log q) blocks introduce no within-block power
    # loss beyond the already-certified comparison-section multipliers.
    assert 3 * 6 == 18
    assert F(18 * 76, 25) < 55
    continuous_block_factor = 3**55
    return (
        exponential_upper,
        angle_sine,
        projector_norm,
        graph_contraction,
        continuous_block_factor,
    )


def action_angle_and_stage_threshold_checks(period, beta):
    """Verify the explicit block-closeness criterion and q thresholds."""

    # On one common Euclidean torus lift, r is by definition
    # max(|X0-X0_*|_2, |k0-k0_*|_2).  Every row below is imported from the
    # displayed C193/C194 physical bounds or derived here from that r.

    # C194: J1<=216 L, J2<=286992 L^2, L=1+t.  K1 is an Eulerian
    # derivative.  Comparing rays launched from distinct initial base
    # points therefore pays |dX|<=J1 r before using K1<=7J2.  Hence
    #
    # |dk| <=216 L r +216*286992*7 L^3 r <434000000 L^3 r.
    assert 216 * 286_992 * 7 == 433_931_904
    assert 216 + 433_931_904 < 434_000_000

    # Phase variables cost sqrt(2).  For gamma=(p.g)/|g|^2, use the safe
    # physical rows |g|<4, |H|<6, |p_*|<12, |gamma_*|<1, and h>=3/2.
    # Direct comparison with the central ray gives
    # |dgamma|<=2|dk|+100|dX|.
    assert 306**2 > 2 * 216**2
    assert 2 * 434_000_000 + 100 * 216 == 868_021_600
    assert 868_021_600 < 869_000_000

    # The bound |p_*|<12 follows from c_h<76/5, sqrt(3h)>2,
    # |gamma_*|<1, and |g|<4.
    assert F(76, 5) / 2 + 4 == F(58, 5) < 12

    closeness_denominator = 87_000_000_000_000  # 8.7e13
    assert F(869_000_000, closeness_denominator) < F(1, 100_000)
    assert F(306, closeness_denominator) < F(1, 10_000_000)
    assert F(1, closeness_denominator) < F(1, 10_000_000)
    # At the initial comparison section, |p0_*|^2=CARRIER<16. Hence the
    # invariant horizontal charge obeys
    # |dc_h|<(|V|+|p0_*||DV|)r<(7+4*10)r=47r<200r.
    assert base.CARRIER < D(16)
    assert 7 + 4 * 10 == 47 < 200
    assert F(200, closeness_denominator) < F(1, 100_000)
    # The invariant energy obeys |E|<=4r.
    assert F(4, closeness_denominator) < F(1, 100_000_000_000)

    # The same criterion keeps the base in |f|<1/10 from f=0.  The local
    # global physical row |grad f|<4 gives |E|<4r, and f is then invariant
    # along the perturbed base trajectory.
    assert F(4, closeness_denominator) < F(1, 10)

    # At a comparison section the exact central k0 has 5.3<|k0|<5.4. The action-angle
    # criterion gives |dk|<434000000/(8.7e13)<1e-5, hence 5<|k|<6 and
    # 3<|P_n k|<4.  Use the actual outward beta interval here rather than
    # C193's deliberately coarse (21/10,9/4) display bounds.
    k_squared = base.I(base.CARRIER) + 3 * beta * beta
    assert k_squared.lo > D("28.09")       # (53/10)^2
    assert k_squared.hi < D("29.16")       # (27/5)^2
    assert F(434_000_000, closeness_denominator) < F(1, 100_000)
    assert base.CARRIER > D(9)              # central |p|>3
    assert base.CARRIER < D(16)             # central |p|<4

    # Current C193 clock:
    # R_filt=ceil((3/8)log n)+1, T<76/25, q=n^8.
    # Thus L=1+T R_filt <=177/25+(57/50)log n.
    # C188 gives log n<=6 n^(1/14), whence L<14 n^(1/14).
    assert period.hi < D(76) / 25
    assert F(1) + F(76, 25) * 2 == F(177, 25)
    assert F(76, 25) * F(3, 8) == F(57, 50)
    assert F(57, 50) * 6 == F(171, 25)
    assert F(177 + 171, 25) == F(348, 25) < 14

    # It is enough that 2744*n^power <=1/(8.7e13).  Raise to the common
    # denominator so every threshold is checked with integers only.
    threshold = 2744 * closeness_denominator

    # r=q^(-1/4)=n^-2: r L^3 <=2744 n^(-25/14).
    n_quarter = 10**10
    assert threshold**14 <= n_quarter**25

    # r_s=q^(-1/12)=n^(-2/3): <=2744 n^(-19/42).
    n_stable = 10**39
    assert threshold**42 <= n_stable**19

    # r_u=q^(-1/3)=n^(-8/3): <=2744 n^(-103/42).
    n_unstable = 10**8
    assert threshold**42 <= n_unstable**103
    assert n_stable >= n_unstable
    return closeness_denominator, n_quarter, n_stable, n_unstable


def main():
    cells, minima, period, beta = fattened_cells()

    forward_gain = lower_cone_gain(cells)
    inverse_cells = [reflected_inverse_cell(cell) for cell in reversed(cells)]
    inverse_gain = lower_cone_gain(inverse_cells)
    forward_image = riccati_cone_image(cells)
    inverse_image = riccati_cone_image(inverse_cells)

    assert forward_gain[0] > 3000
    assert inverse_gain[0] > 3000
    assert forward_image[0] > IMAGE_LO
    assert forward_image[1] < IMAGE_HI
    assert inverse_image[0] > IMAGE_LO
    assert inverse_image[1] < IMAGE_HI

    exponential_upper, angle_sine, projector_norm, contraction, block_factor = (
        exponential_and_cone_geometry_checks()
    )
    denominator, n_quarter, n_stable, n_unstable = (
        action_angle_and_stage_threshold_checks(period, beta)
    )

    print("A2 off-ray finite-horizon dominated cones: PASS")
    print("fat-box raw B lower endpoints:", *(str(value) for value in minima))
    print("forward lower endpoint:", forward_gain[0], forward_gain[1])
    print("inverse-reflected lower endpoint:", inverse_gain[0], inverse_gain[1])
    print("forward cone image:", *forward_image)
    print("inverse-reflected cone image:", *inverse_image)
    print("exact exp(8) upper < 3000:", exponential_upper)
    print(
        "any forward/backward cone-line angle sine >=",
        angle_sine,
        "; two-line oblique-projector norm <=",
        projector_norm,
    )
    print("each block projective slope contraction <", contraction)
    print("continuous one-block physical operator norm <", block_factor)
    print("block-closeness criterion: r*(1+T*R_filt)^3 <= 1/", denominator)
    print("q=n^8 threshold for r=q^(-1/4): n >=", n_quarter)
    print(
        "two-width thresholds: r_s=q^(-1/12), n >=",
        n_stable,
        "; r_u=q^(-1/3), n >=",
        n_unstable,
    )
    print(
        "BOUNDARY: finite-horizon principal dominated cones only; "
        "no invariant bundle/closed phase/packet/viscosity"
    )


if __name__ == "__main__":
    main()
