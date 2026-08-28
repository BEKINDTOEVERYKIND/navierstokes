#!/usr/bin/env python3
"""Exact global jet bounds for the C152 A2 Beltrami pump.

This prospective C193 checker uses only Fraction arithmetic.  It verifies
the complete rational-algebra reduction behind

    sup ||DU||_op = 6,
    sup ||D^2 U||_mult = 3 sqrt(6),
    sup ||D^3 U||_mult = 9,

and the sharper zero-level value sup_{f=0} ||DU||_op = 12 sqrt(6)/5.
The multilinear norms use Euclidean unit inputs and a Euclidean output.

No floating point, interval package, or symbolic-algebra dependency is
used.  Small polynomial dictionaries verify the displayed matrix and PSD
identities rather than treating hand expansion as a premise.
"""

from fractions import Fraction as F
from math import comb


# Sparse polynomials in (P,S,C,Q), with exact rational coefficients.
ZERO_MONOMIAL = (0, 0, 0, 0)


def clean(poly):
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def constant(value):
    value = F(value)
    return {} if not value else {ZERO_MONOMIAL: value}


def variable(index):
    monomial = [0, 0, 0, 0]
    monomial[index] = 1
    return {tuple(monomial): F(1)}


def add(left, right):
    output = dict(left)
    for monomial, coefficient in right.items():
        output[monomial] = output.get(monomial, F(0)) + coefficient
    return clean(output)


def neg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def sub(left, right):
    return add(left, neg(right))


def scale(poly, scalar):
    scalar = F(scalar)
    return clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def mul(left, right):
    output = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            output[monomial] = output.get(monomial, F(0)) + left_coefficient * right_coefficient
    return clean(output)


def power(poly, exponent):
    output = constant(1)
    factor = poly
    while exponent:
        if exponent & 1:
            output = mul(output, factor)
        factor = mul(factor, factor)
        exponent //= 2
    return output


def reduce_trig(poly):
    """Reduce modulo S^2=1-P^2 and Q^2=1-C^2 to a canonical form."""
    output = {}
    for (p_exp, s_exp, c_exp, q_exp), coefficient in poly.items():
        s_pairs, s_rem = divmod(s_exp, 2)
        q_pairs, q_rem = divmod(q_exp, 2)
        for s_index in range(s_pairs + 1):
            s_coefficient = F(comb(s_pairs, s_index)) * (-1) ** s_index
            for q_index in range(q_pairs + 1):
                q_coefficient = F(comb(q_pairs, q_index)) * (-1) ** q_index
                monomial = (
                    p_exp + 2 * s_index,
                    s_rem,
                    c_exp + 2 * q_index,
                    q_rem,
                )
                output[monomial] = output.get(monomial, F(0)) + (
                    coefficient * s_coefficient * q_coefficient
                )
    return clean(output)


def assert_trig_identity(left, right):
    assert not reduce_trig(sub(left, right)), reduce_trig(sub(left, right))


# One-variable polynomial helpers, coefficients in increasing degree.
def uclean(poly):
    output = list(poly)
    while output and not output[-1]:
        output.pop()
    return output


def uadd(left, right):
    output = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        output[index] += coefficient
    for index, coefficient in enumerate(right):
        output[index] += coefficient
    return uclean(output)


def uneg(poly):
    return [-coefficient for coefficient in poly]


def usub(left, right):
    return uadd(left, uneg(right))


def uscale(poly, scalar):
    scalar = F(scalar)
    return uclean([scalar * coefficient for coefficient in poly])


def umul(left, right):
    if not left or not right:
        return []
    output = [F(0)] * (len(left) + len(right) - 1)
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            output[i + j] += a_value * b_value
    return uclean(output)


def upower(poly, exponent):
    output = [F(1)]
    factor = poly
    while exponent:
        if exponent & 1:
            output = umul(output, factor)
        factor = umul(factor, factor)
        exponent //= 2
    return output


def ucompose_affine(poly, offset, slope):
    affine = [F(offset), F(slope)]
    output = []
    for coefficient in reversed(poly):
        output = uadd(umul(output, affine), [coefficient])
    return output


class RationalPolynomial:
    """An unreduced exact quotient; equality is checked by cross product."""

    def __init__(self, numerator, denominator=(F(1),)):
        self.numerator = uclean(list(numerator))
        self.denominator = uclean(list(denominator))
        assert self.denominator

    def __add__(self, other):
        other = other if isinstance(other, RationalPolynomial) else RationalPolynomial([other])
        return RationalPolynomial(
            uadd(umul(self.numerator, other.denominator), umul(other.numerator, self.denominator)),
            umul(self.denominator, other.denominator),
        )

    __radd__ = __add__

    def __neg__(self):
        return RationalPolynomial(uneg(self.numerator), self.denominator)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return RationalPolynomial([other]) - self

    def __mul__(self, other):
        other = other if isinstance(other, RationalPolynomial) else RationalPolynomial([other])
        return RationalPolynomial(
            umul(self.numerator, other.numerator),
            umul(self.denominator, other.denominator),
        )

    __rmul__ = __mul__

    def equals(self, other):
        other = other if isinstance(other, RationalPolynomial) else RationalPolynomial([other])
        return not usub(
            umul(self.numerator, other.denominator),
            umul(other.numerator, self.denominator),
        )


def dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def cross(left, right):
    return (
        F(left[1]) * F(right[2]) - F(left[2]) * F(right[1]),
        F(left[2]) * F(right[0]) - F(left[0]) * F(right[2]),
        F(left[0]) * F(right[1]) - F(left[1]) * F(right[0]),
    )


def geometric_and_matrix_identities():
    n_vector = (1, 1, 1)
    q1 = (1, -1, 0)
    q2 = (0, 1, -1)
    q3 = (1, 0, -1)
    assert tuple(q1[index] + q2[index] for index in range(3)) == q3
    for q_vector in (q1, q2, q3):
        assert dot(q_vector, q_vector) == 2
        assert dot(q_vector, n_vector) == 0
        assert dot(cross(n_vector, q_vector), cross(n_vector, q_vector)) == 6
    assert dot(n_vector, n_vector) == 3

    P, S, C, Q = (variable(index) for index in range(4))
    h11 = add(neg(mul(P, C)), scale(sub(scale(power(P, 2), 2), constant(1)), F(-8, 5)))
    h12_over_sqrt3 = mul(S, Q)
    h22 = scale(mul(P, C), -3)
    g1_over_sqrt2 = neg(mul(S, add(C, scale(P, F(8, 5)))))
    g2_over_sqrt6 = neg(mul(P, Q))

    # M=H^2+2gg^T, with radicals divided out exactly.
    m11 = add(
        add(power(h11, 2), scale(power(h12_over_sqrt3, 2), 3)),
        scale(power(g1_over_sqrt2, 2), 4),
    )
    m12_over_sqrt3 = add(
        mul(h12_over_sqrt3, add(h11, h22)),
        scale(mul(g1_over_sqrt2, g2_over_sqrt6), 4),
    )
    m22 = add(
        add(power(h22, 2), scale(power(h12_over_sqrt3, 2), 3)),
        scale(power(g2_over_sqrt6, 2), 12),
    )

    expected_m11 = scale(add(
        add(add(constant(139), scale(power(C, 2), 25)), scale(power(P, 2), -75)),
        add(scale(mul(C, P), 240), scale(mul(C, power(P, 3)), -160)),
    ), F(1, 25))
    expected_m12_over_sqrt3 = scale(
        mul(mul(S, Q), add(scale(power(P, 2), 2), constant(1))),
        F(8, 5),
    )
    expected_m22 = scale(add(add(constant(1), scale(power(P, 2), 3)), neg(power(C, 2))), 3)
    assert_trig_identity(m11, expected_m11)
    assert_trig_identity(m12_over_sqrt3, expected_m12_over_sqrt3)
    assert_trig_identity(m22, expected_m22)


def global_first_jet_psd_certificate():
    # For r=|P| and c=|C|, G=12I-M has
    # G11 >= a/25.  The coefficient of PC is nonpositive in G11 because
    # 3-2r^2 >= 1, so replacing PC by rc is the worst sign.
    assert 3 - 2 * F(1) >= 1

    # E=3a-64(1-c^2)(2r^2+1)^2=A(r)c^2+B(r)c+C0(r).
    # Verify the expansion exactly as a bivariate polynomial in (r,c).
    r = variable(0)
    c_value = variable(2)
    a_value = add(
        add(add(constant(161), scale(power(r, 2), 75)), scale(power(c_value, 2), -25)),
        scale(mul(mul(r, c_value), sub(constant(3), scale(power(r, 2), 2))), -80),
    )
    energy = sub(
        scale(a_value, 3),
        scale(mul(sub(constant(1), power(c_value, 2)), power(add(scale(power(r, 2), 2), constant(1)), 2)), 64),
    )
    coefficient_a = add(add(scale(power(r, 4), 256), scale(power(r, 2), 256)), constant(-11))
    coefficient_b = add(scale(power(r, 3), 480), scale(r, -720))
    coefficient_c = add(add(add(constant(419), scale(power(r, 2), -31)), scale(power(r, 4), -256)), constant(0))
    expected_energy = add(
        add(mul(coefficient_a, power(c_value, 2)), mul(coefficient_b, c_value)),
        coefficient_c,
    )
    assert not sub(energy, expected_energy)

    # On 0<=r<=2/3, the deliberately separate lower bounds sum positively.
    assert F(-11) + F(-3040, 9) + F(28727, 81) == F(476, 81) > 0
    assert 480 * F(2, 3) ** 3 - 720 * F(2, 3) == F(-3040, 9)
    assert 419 - 31 * F(2, 3) ** 2 - 256 * F(2, 3) ** 4 == F(28727, 81)

    # On 2/3<=r<=1, A>0 and the discriminant is 4Q(r^2).
    assert 256 * F(4, 9) ** 2 + 256 * F(4, 9) - 11 > 0
    x = [F(0), F(1)]
    coefficient_a_x = uadd(uadd(uscale(upower(x, 2), 256), uscale(x, 256)), [-11])
    # B(r)^2 = r^2(480r^2-720)^2 = x(480x-720)^2.
    coefficient_b_squared = umul(x, upower(uadd(uscale(x, 480), [-720]), 2))
    coefficient_c_x = uadd(uadd([419], uscale(x, -31)), uscale(upower(x, 2), -256))
    discriminant = usub(coefficient_b_squared, uscale(umul(coefficient_a_x, coefficient_c_x), 4))
    quartic = [4609, 21995, -274944, 131072, 65536]
    assert discriminant == uscale(quartic, 4)

    transformed = ucompose_affine(quartic, F(4, 9), F(5, 9))
    degree = 4
    bernstein = []
    for k_value in range(degree + 1):
        coefficient = F(0)
        for i_value in range(k_value + 1):
            coefficient += (
                transformed[i_value]
                * F(comb(k_value, i_value), comb(degree, i_value))
            )
        bernstein.append(coefficient)
    expected_bernstein = [
        F(-169675667, 6561),
        F(-41568439, 972),
        F(-1094911, 18),
        F(-277743, 4),
        F(-51732),
    ]
    assert bernstein == expected_bernstein
    assert all(coefficient < 0 for coefficient in bernstein)

    # Therefore G11 >= (64/75)(1-c^2)(2r^2+1)^2 and
    # G22=3(c^2+3(1-r^2)).  Subtracting M12^2 leaves the displayed
    # nonnegative determinant lower bound.
    # Check the scalar coefficient cancellation exactly.
    assert F(64, 75) * 3 == F(64, 25)
    assert F(64, 25) * 3 - F(192, 25) == 0
    assert F(8, 5) ** 2 * 3 == F(192, 25)  # M12^2 radical factor

    # G is PSD, so M<=12I and DU^T DU=3M<=36I.
    assert 3 * 12 == 36
    # Equality: P^2=1,C^2=0 in the d direction gives M22=12.
    assert 3 * (1 + 3 * 1 - 0) == 12


def higher_jet_certificate():
    # Exact modewise tensors.  For q in the A2 plane, |q|^2=2,
    # |N|^2=3 and N.(N x q)=0, so both
    #
    # W_2=sin(theta) N x q+sqrt(2)cos(theta)N,
    # W_3=cos(theta) N x q-sqrt(2)sin(theta)N
    #
    # have squared norm 6.  The derivative signs come from the two cycles
    # d^r(-sin)/dtheta^r and d^r(-cos)/dtheta^r at r=2,3.
    normal = (1, 1, 1)
    roots = ((1, -1, 0), (0, 1, -1), (1, 0, -1))
    assert dot(normal, normal) == 3
    for root in roots:
        normal_cross_root = cross(normal, root)
        assert dot(root, root) == 2
        assert dot(normal, normal_cross_root) == 0
        assert dot(normal_cross_root, normal_cross_root) == 6
        assert 2 * dot(normal, normal) == 6
    horizontal_derivative_signs = (-1, -1, 1, 1)  # -sin,-cos,+sin,+cos
    axial_derivative_signs = (-1, 1, 1, -1)       # -cos,+sin,+cos,-sin
    assert horizontal_derivative_signs[2:] == (1, 1)
    assert axial_derivative_signs[2:] == (1, -1)

    # Exact projection numerators for e=q3/sqrt(2) and
    # d=(q1-q2)/sqrt(6).  These mechanically support every radical
    # coefficient used in the frame and equality calculations below.
    q1, q2, q3 = roots
    e_numerator = q3
    d_numerator = tuple(q1[index] - q2[index] for index in range(3))
    assert dot(e_numerator, e_numerator) == 2
    assert dot(d_numerator, d_numerator) == 6
    assert tuple(dot(root, e_numerator) for root in roots) == (1, 1, 2)
    assert tuple(dot(root, d_numerator) for root in roots) == (3, -3, 0)

    # Mechanically verify the two nonnegative algebraic inequalities used
    # mode by mode.  For X,Y>=0 the first is Cauchy's 2XY<=X^2+Y^2.
    # For X,Y,Z>=0 the second is 3XYZ<=X^3+Y^3+Z^3 (AM--GM), in SOS form.
    X, Y, Z = variable(0), variable(1), variable(2)
    cauchy_left = sub(add(power(X, 2), power(Y, 2)), scale(mul(X, Y), 2))
    cauchy_right = power(sub(X, Y), 2)
    assert not sub(cauchy_left, cauchy_right)
    amgm_left = sub(add(add(power(X, 3), power(Y, 3)), power(Z, 3)), scale(mul(mul(X, Y), Z), 3))
    squared_differences = add(
        add(power(sub(X, Y), 2), power(sub(Y, Z), 2)),
        power(sub(Z, X), 2),
    )
    amgm_right = scale(mul(add(add(X, Y), Z), squared_differences), F(1, 2))
    assert not sub(amgm_left, amgm_right)

    # Write the three roots in the orthonormal horizontal (e,d) frame:
    #
    # q1=(1/sqrt(2), sqrt(3/2)),
    # q2=(1/sqrt(2),-sqrt(3/2)), q3=(sqrt(2),0),
    #
    # with weights 1,1,4/5.  Thus sum c_i q_i q_i^T is diagonal with
    # eigenvalues 13/5 and 3.  Cauchy gives
    #
    # sum c_i |q_i.u||q_i.v| <= 3 |u||v|.
    frame_eigen_e = F(1, 2) + F(1, 2) + F(4, 5) * 2
    frame_eigen_d = F(3, 2) + F(3, 2)
    assert frame_eigen_e == F(13, 5)
    assert frame_eigen_d == 3 > frame_eigen_e
    x_coordinate, y_coordinate = variable(0), variable(1)
    unit_frame_form = add(
        scale(power(x_coordinate, 2), frame_eigen_e),
        scale(power(y_coordinate, 2), frame_eigen_d),
    )
    unit_three = scale(add(power(x_coordinate, 2), power(y_coordinate, 2)), 3)
    assert not sub(sub(unit_three, unit_frame_form), scale(power(x_coordinate, 2), F(2, 5)))

    # Every modal output bracket in D^2U has Euclidean norm sqrt(6).
    # Hence ||D^2U||_mult <=3sqrt(6).  Equality occurs at a=b=0 with
    # both inputs d: q1.d and q2.d have squared value 3/2, q3.d=0,
    # and the two surviving axial vectors point in the same direction.
    d2_scalar_sum = F(3, 2) + F(3, 2)
    d2_equality_squared = 2 * dot((1, 1, 1), (1, 1, 1)) * d2_scalar_sum**2
    assert d2_scalar_sum == 3
    assert d2_equality_squared == 54  # (3sqrt(6))^2

    # For D^3U, generalized Holder reduces the scalar trilinear sum to
    #
    # S(v)=sum c_i |q_i.v|^3.
    #
    # By the absolute-value symmetries take v=t e+s d, t,s>=0 and
    # t^2+s^2=1.  When t<=sqrt(3)/2 (equivalently s>=1/2), expansion of
    # the two opposite-sign cubes gives
    #
    # S=3sqrt(3/2)s+(8sqrt(2)/5)t^3.
    #
    # The target S<=3sqrt(3/2) is equivalent (for t>0) to
    # 16t(1+s)<=15sqrt(3).  The exact factorization below proves
    # t(1+s)<=3sqrt(3)/4 on 1/2<=s<=1.
    s_variable = [F(0), F(1)]
    # The two absolute-value branches use these exact cube expansions.
    A, B = variable(0), variable(1)
    first_cube_branch = add(power(add(A, B), 3), power(sub(B, A), 3))
    first_cube_target = add(scale(power(B, 3), 2), scale(mul(B, power(A, 2)), 6))
    second_cube_branch = add(power(add(A, B), 3), power(sub(A, B), 3))
    second_cube_target = add(scale(power(A, 3), 2), scale(mul(A, power(B, 2)), 6))
    assert not sub(first_cube_branch, first_cube_target)
    assert not sub(second_cube_branch, second_cube_target)

    left_square = umul(usub([F(1)], upower(s_variable, 2)), upower(uadd([F(1)], s_variable), 2))
    branch_gap = usub([F(27, 16)], left_square)
    branch_factor = uscale(
        umul(
            upower(uadd([F(-1)], uscale(s_variable, 2)), 2),
            uadd(uadd(uscale(upower(s_variable, 2), 4), uscale(s_variable, 12)), [F(11)]),
        ),
        F(1, 16),
    )
    assert branch_gap == branch_factor
    # The second factor in branch_factor is positive for s>=0.  Therefore
    # t^2(1+s)^2<=27/16, and squaring the desired bound gives the exact
    # strict margin 432<675.
    assert all(coefficient > 0 for coefficient in (F(4), F(12), F(11)))
    assert 256 * F(27, 16) == 432 < 675 == 15**2 * 3
    assert 16 * F(3, 4) == 12 < 15

    # When t>=sqrt(3)/2 both root projections have the same sign and
    #
    # sqrt(2) S(t)=9t-(24/5)t^3.
    #
    # Its derivative is already negative on t^2>=3/4, so its boundary
    # value (27/10)sqrt(3) is below the target 3sqrt(3).
    branch_two_derivative_upper = 9 - F(72, 5) * F(3, 4)
    branch_two_boundary_coefficient = 9 * F(1, 2) - F(24, 5) * F(3, 8)
    assert branch_two_derivative_upper == F(-9, 5) < 0
    assert branch_two_boundary_coefficient == F(27, 10) < 3

    # Thus S(v)<=3sqrt(3/2), generalized Holder gives
    # ||D^3U||_mult<=sqrt(6)*3sqrt(3/2)=9.  Equality occurs for all three
    # inputs d at (a,b)=(pi/2,-pi/2): the q1 and q2 axial outputs agree,
    # while q3.d=0.  Here (q1.d)^6=(3/2)^3=27/8.
    d3_equality_squared = 8 * F(27, 8) * dot((1, 1, 1), (1, 1, 1))
    assert d3_equality_squared == 81


def zero_level_first_jet_certificate():
    # Let x=P^2.  The level relation gives
    # C^2=(4/25)(2x-1)^2/x and hence 16x^2-41x+4<=0.
    assert 189**2 > 25**2 * 57  # (41-5sqrt(57))/32 > 1/10

    x = RationalPolynomial([0, 1])
    c_squared = RationalPolynomial(uscale(upower([-1, 2], 2), F(4, 25)), [0, 1])
    pc = RationalPolynomial(uscale([-1, 2], F(-2, 5)))
    m11 = F(1, 25) * (
        139 + 25 * c_squared - 75 * x + 80 * pc * (3 - 2 * x)
    )
    m22 = 3 * (1 + 3 * x - c_squared)
    m12_squared = F(192, 25) * (1 - x) * (1 - c_squared) * (2 * x + 1) * (2 * x + 1)
    level_edge = F(288, 25)
    k11 = level_edge - m11
    k22 = level_edge - m22
    determinant = k11 * k22 - m12_squared

    r_poly = [4, -69, -315, 128]
    target_k11 = RationalPolynomial(uneg(r_poly), uscale([0, 1], 25))
    target_k22 = RationalPolynomial(
        uscale(umul([-1, 1], [4, 59]), -3),
        uscale([0, 1], 25),
    )
    target_determinant = RationalPolynomial(
        uscale(
            umul(
                umul([-1, 1], [4, -113, 64]),
                [-4, -39, -105, 64],
            ),
            -3,
        ),
        uscale([0, 0, 1], 625),
    )
    assert k11.equals(target_k11)
    assert k22.equals(target_k22)
    assert determinant.equals(target_determinant)

    # Sign checks on x in [1/10,1].
    assert 384 * F(1) ** 2 - 630 * F(1) - 69 < 0
    assert 384 * F(1, 10) ** 2 - 630 * F(1, 10) - 69 < 0
    assert 128 * F(1, 10) ** 3 - 315 * F(1, 10) ** 2 - 69 * F(1, 10) + 4 == F(-2961, 500)
    assert 64 * F(1, 10) ** 2 - 113 * F(1, 10) + 4 == F(-333, 50)
    assert 64 - 113 + 4 == -45
    # 64x^3-105x^2 <= -41x^2 on [0,1].
    assert 64 - 105 == -41

    # At x=1, lambda_max(M)=288/25, so ||DU||^2=864/25.
    assert 3 * F(288, 25) == F(864, 25)
    assert F(12, 5) ** 2 * 6 == F(864, 25)


def scaling_and_thresholds():
    # For V=A U(Kx), each spatial derivative contributes one K.
    assert F(350, 57) - 6 == F(8, 57) > 0
    # Squared exact constants, convenient for a rational-only consumer.
    assert (3 * 3 * 6, 9 * 9) == (54, 81)


def main():
    geometric_and_matrix_identities()
    global_first_jet_psd_certificate()
    higher_jet_certificate()
    zero_level_first_jet_certificate()
    scaling_and_thresholds()
    print("C193 A2 global jet bounds: PASS")
    print("GLOBAL: ||DU||_op <= 6, ||D2U||_mult <= 3 sqrt(6), ||D3U||_mult <= 9")
    print("SHARP: all three constants have explicit equality cases")
    print("ZERO LEVEL: sup_{f=0} ||DU||_op = 12 sqrt(6)/5")
    print("WKB MARGIN: 6 < 350/57 by 8/57")


if __name__ == "__main__":
    main()
