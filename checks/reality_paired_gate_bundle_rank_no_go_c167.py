#!/usr/bin/env python3
"""Dependency-free exact checks for C167's static gate-bundle no-go.

The checker verifies with rational polynomial arithmetic:

* the full complex-polarization round-trip coefficients;
* the non-tangential eigenline/common-metric slope;
* the reality-pair overdetermination and circular-polarization identity;
* the sign-definite diagonal gap of every positive pair mixture; and
* strict radial detuning throughout the oscillatory band |y|<2.

No time-dependent control, finite-band, collar, wake, viscosity, or
Navier--Stokes stage claim is tested.
"""

from fractions import Fraction as Q


def trim(p):
    p = [Q(x) for x in p]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, q):
    out = [Q(0)] * max(len(p), len(q))
    for j, value in enumerate(p):
        out[j] += value
    for j, value in enumerate(q):
        out[j] += value
    return trim(out)


def neg(p):
    return trim([-value for value in p])


def sub(p, q):
    return add(p, neg(q))


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for j, left in enumerate(p):
        for k, right in enumerate(q):
            out[j + k] += left * right
    return trim(out)


def scale(p, value):
    return trim([Q(value) * coefficient for coefficient in p])


def derivative(p):
    if len(p) == 1:
        return (Q(0),)
    return trim([j * p[j] for j in range(1, len(p))])


ONE = (Q(1),)
TWO = (Q(2),)
FOUR = (Q(4),)
Y = (Q(0), Q(1))
Y2 = mul(Y, Y)
Y4 = mul(Y2, Y2)
D_PLUS = add(add(FOUR, scale(Y, 2)), Y2)
D_MINUS = add(add(FOUR, scale(Y, -2)), Y2)
P_EVEN = add(add((Q(16),), scale(Y2, 4)), Y4)


# Complex rational scalars are stored as (real, imaginary).  The three-vector
# coordinates below use (e_r,e_t,e_z/sqrt(3)); their Euclidean metric is
# diag(1,1,1/3).  This makes a direct evaluation of the projected Euler
# symbol rational and independent of the closed formulas being checked.
CZ = (Q(0), Q(0))


def cadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def cneg(value):
    return -value[0], -value[1]


def csub(left, right):
    return cadd(left, cneg(right))


def cmul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def cscale(value, coefficient):
    coefficient = Q(coefficient)
    return coefficient * value[0], coefficient * value[1]


def cconj(value):
    return value[0], -value[1]


def cdiv_real(value, denominator):
    return cscale(value, Q(1) / Q(denominator))


def vadd(left, right):
    return tuple(cadd(a, b) for a, b in zip(left, right))


def vcscale(value, coefficient):
    return tuple(cmul(coefficient, item) for item in value)


def real_vector(*entries):
    return tuple((Q(entry), Q(0)) for entry in entries)


def bilinear_dot(left, right):
    weights = (Q(1), Q(1), Q(1, 3))
    out = CZ
    for weight, a, b in zip(weights, left, right):
        out = cadd(out, cscale(cmul(a, b), weight))
    return out


def symmetric_symbol_unprojected(k, a, q, b):
    """Return (a.q)b+(b.k)a in the scaled rational coordinates."""
    return vadd(
        vcscale(b, bilinear_dot(a, q)),
        vcscale(a, bilinear_dot(b, k)),
    )


def check_projected_complex_blocks_directly():
    # Dotting the unprojected symbol with a divergence-free output vector is
    # exactly the corresponding coefficient after Leray projection.  Use the
    # unnormalised daughter vector sqrt(D)e_perp so no square root occurs.
    # Several signed rational heights and genuinely complex u,v make the
    # conjugation and orientation tests non-tautological.
    u = (Q(2, 3), Q(-1, 5))
    v = (Q(-3, 7), Q(4, 9))
    ubar, vbar = cconj(u), cconj(v)
    et = real_vector(0, 1, 0)
    e_sigma = real_vector(Q(1, 2), 0, Q(-3, 2))

    for y in (Q(-7, 4), Q(-2, 3), Q(1, 5), Q(7, 4)):
        p = real_vector(1, 0, 1)
        g = real_vector(0, 0, y)
        minus_g = real_vector(0, 0, -y)
        daughter_wave = real_vector(1, 0, 1 + y)
        e_perp_tilde = real_vector(-(1 + y), 0, 3)
        gate = (u, v, CZ)
        reverse_gate = (ubar, vbar, CZ)

        forward_sigma = symmetric_symbol_unprojected(
            p, e_sigma, g, gate
        )
        forward_tangent = symmetric_symbol_unprojected(p, et, g, gate)
        f11_tilde = bilinear_dot(e_perp_tilde, forward_sigma)
        f21 = bilinear_dot(et, forward_sigma)
        f12_tilde = bilinear_dot(e_perp_tilde, forward_tangent)
        f22 = bilinear_dot(et, forward_tangent)
        assert f11_tilde == cscale(u, -(4 - y * y) / 2)
        assert f21 == cscale(v, -y / 2)
        assert f12_tilde == CZ
        assert f22 == u

        # Feeding sqrt(D)e_perp scales the first reverse column by sqrt(D).
        reverse_perp_tilde = symmetric_symbol_unprojected(
            daughter_wave, e_perp_tilde, minus_g, reverse_gate
        )
        reverse_tangent = symmetric_symbol_unprojected(
            daughter_wave, et, minus_g, reverse_gate
        )
        r11_tilde = bilinear_dot(e_sigma, reverse_perp_tilde)
        r21_tilde = bilinear_dot(et, reverse_perp_tilde)
        r12 = bilinear_dot(e_sigma, reverse_tangent)
        r22 = bilinear_dot(et, reverse_tangent)
        assert r11_tilde == cscale(ubar, -(y + 2))
        assert r21_tilde == cscale(vbar, -y)
        assert r12 == CZ
        assert r22 == ubar

        # Multiply the directly obtained blocks.  The two tilde entries each
        # carry one sqrt(D), so their product is divided by D.
        d = y * y + 2 * y + 4
        m11 = cdiv_real(cmul(r11_tilde, f11_tilde), d)
        m21 = cadd(cdiv_real(cmul(r21_tilde, f11_tilde), d), cmul(r22, f21))
        m22 = cmul(r22, f22)
        abs_u2 = cmul(ubar, u)
        expected_r = (4 - y * y) * (y + 2) / (2 * d)
        expected_q = cadd(
            cscale(cmul(u, vbar), y * (4 - y * y) / (2 * d)),
            cscale(cmul(ubar, v), -y / 2),
        )
        assert m11 == cscale(abs_u2, expected_r)
        assert m21 == expected_q
        assert m22 == abs_u2


def check_round_trip_and_slope():
    # Strip the common |u|^2.  The upper return eigenvalue is
    # r(y)=(4-y^2)(y+2)/(2D), and the lower entry is
    # y[(4-y^2) conjugate(rho)-D rho]/(2D).
    four_minus_y2 = sub(FOUR, Y2)
    y_plus_2 = add(Y, TWO)
    y_plus_4 = add(Y, FOUR)
    r_num = mul(four_minus_y2, y_plus_2)

    # 2D(r-1)=-y^2(y+4).
    assert sub(r_num, scale(D_PLUS, 2)) == neg(mul(Y2, y_plus_4))

    # With rho=a+i b,
    # (4-y^2)conj(rho)-D rho
    # =-2y(y+1)a - 2i(y+4)b.
    coefficient_a = sub(four_minus_y2, D_PLUS)
    coefficient_ib = neg(add(four_minus_y2, D_PLUS))
    assert coefficient_a == scale(mul(Y, add(Y, ONE)), -2)
    assert coefficient_ib == scale(y_plus_4, -2)

    # Hence m=2(y+1)/(y+4) Re(rho) + 2i/y Im(rho).
    # The real-coefficient difference between the y and -y reality
    # members has numerator 12y and denominator 16-y^2.
    numerator = scale(
        sub(
            mul(add(Y, ONE), sub(FOUR, Y)),
            mul(sub(ONE, Y), add(Y, FOUR)),
        ),
        2,
    )
    assert numerator == scale(Y, 12)
    assert mul(add(Y, FOUR), sub(FOUR, Y)) == sub((Q(16),), Y2)

    # The common Hermitian source-metric equation for
    # M=[[lambda,0],[q,mu]] is
    # (mu-lambda) z = conjugate(q) q22.
    # Substitution z/q22=-conjugate(m), m=q/(lambda-mu), is exact.
    lam = Q(7, 5)
    mu = Q(11, 6)
    qr, qi = Q(3, 7), Q(-5, 9)
    # Complex pairs (real,imag).
    m = (qr / (lam - mu), qi / (lam - mu))
    z = (-m[0], m[1])  # -conjugate(m)
    left = ((mu - lam) * z[0], (mu - lam) * z[1])
    right = (qr, -qi)
    assert left == right

    # Check the entire Hermitian relation QM=M*Q, not just its (1,2)
    # scalar consequence.  A large q11 makes Q positive; its value cancels
    # from the only nontrivial compatibility equation.
    q11, q22 = Q(17), Q(3, 2)
    z_pair = (q22 * z[0], q22 * z[1])
    q_pair = (qr, qi)

    def cmatmul(left, right):
        return tuple(
            tuple(
                cadd(cmul(left[row][0], right[0][column]),
                     cmul(left[row][1], right[1][column]))
                for column in range(2)
            )
            for row in range(2)
        )

    def cadjoint(matrix):
        return tuple(
            tuple(cconj(matrix[column][row]) for column in range(2))
            for row in range(2)
        )

    matrix_q = (
        ((q11, Q(0)), z_pair),
        (cconj(z_pair), (q22, Q(0))),
    )
    matrix_m = (
        ((lam, Q(0)), CZ),
        (q_pair, (mu, Q(0))),
    )
    assert q11 * q22 > z_pair[0] ** 2 + z_pair[1] ** 2
    assert cmatmul(matrix_q, matrix_m) == cmatmul(
        cadjoint(matrix_m), matrix_q
    )


def madd(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, Q(0)) + coefficient
        if out[monomial] == 0:
            del out[monomial]
    return out


def mneg(value):
    return {monomial: -coefficient for monomial, coefficient in value.items()}


def mscale(value, coefficient):
    coefficient = Q(coefficient)
    return {
        monomial: coefficient * item
        for monomial, item in value.items()
        if coefficient * item
    }


def mmul(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            monomial = tuple(j + k for j, k in zip(a, b))
            out[monomial] = out.get(monomial, Q(0)) + x * y
            if out[monomial] == 0:
                del out[monomial]
    return out


def cs_add(left, right):
    return madd(left[0], right[0]), madd(left[1], right[1])


def cs_scale(value, coefficient):
    return mscale(value[0], coefficient), mscale(value[1], coefficient)


def cs_mul(left, right):
    return (
        madd(mmul(left[0], right[0]), mneg(mmul(left[1], right[1]))),
        madd(mmul(left[0], right[1]), mmul(left[1], right[0])),
    )


def cs_conj(value):
    return value[0], mneg(value[1])


def check_angular_circular_rigidity():
    # Variables are c,s,Xr,Xi,Yr,Yi.
    def var(index):
        powers = [0] * 6
        powers[index] = 1
        return {tuple(powers): Q(1)}

    zero = {}
    c, s = var(0), var(1)
    xr, xi, yr, yi = var(2), var(3), var(4), var(5)
    X = (xr, xi)
    Yc = (yr, yi)
    c_scalar = (c, zero)
    s_scalar = (s, zero)

    u = cs_add(cs_mul(c_scalar, X), cs_mul(s_scalar, Yc))
    v = cs_add(
        cs_scale(cs_mul(s_scalar, X), -1),
        cs_mul(c_scalar, Yc),
    )
    actual = cs_mul(v, cs_conj(u))[0]

    abs_x2 = cs_mul(X, cs_conj(X))[0]
    abs_y2 = cs_mul(Yc, cs_conj(Yc))[0]
    re_y_bar_x = cs_mul(Yc, cs_conj(X))[0]
    expected = madd(
        mmul(mmul(c, s), madd(abs_y2, mneg(abs_x2))),
        mmul(madd(mmul(c, c), mneg(mmul(s, s))), re_y_bar_x),
    )
    assert actual == expected

    # If Y=sigma*iX, verify v=sigma*i*u for both helicities as a formal
    # polynomial identity, without choosing an angle or phase.
    for sigma in (-1, 1):
        i_x = (mneg(xi), xr)
        y_circular = cs_scale(i_x, sigma)
        u_circular = cs_add(
            cs_mul(c_scalar, X), cs_mul(s_scalar, y_circular)
        )
        v_circular = cs_add(
            cs_scale(cs_mul(s_scalar, X), -1),
            cs_mul(c_scalar, y_circular),
        )
        i_u = (mneg(u_circular[1]), u_circular[0])
        assert cs_add(v_circular, cs_scale(i_u, -sigma)) == ({}, {})

    # After dividing by nonzero X, the two open-arc coefficient conditions
    # say Re(Y/X)=0 and |Y/X|=1.  With Y/X=i*b this leaves b^2-1=0; record
    # its exact two-helicity factorization rather than merely sampling the
    # two advertised solutions.
    b_minus_1 = (Q(-1), Q(1))
    b_plus_1 = (Q(1), Q(1))
    assert mul(b_minus_1, b_plus_1) == (Q(-1), Q(0), Q(1))


def check_distinct_height_incompatibility():
    # For circular polarizations, equality of the source slopes at positive
    # heights g_a,g_b is sigma_a*g_b-sigma_b*g_a=0.  Same helicity reduces
    # to g_b-g_a=0; opposite helicity reduces to +/-(g_a+g_b)=0, impossible
    # for positive heights.  Encode all four sign cases exactly.
    # Linear forms use coefficients (constant, g_a, g_b).
    for sigma_a in (-1, 1):
        for sigma_b in (-1, 1):
            relation = (Q(0), Q(-sigma_b), Q(sigma_a))
            if sigma_a == sigma_b:
                assert relation == (
                    Q(0), Q(-sigma_a), Q(sigma_a)
                )
            else:
                assert relation == (
                    Q(0), Q(sigma_a), Q(sigma_a)
                )

    # A nontrivial exact palette confirms that no accidental equality was
    # introduced by reciprocal or sign conventions.
    positive_heights = (Q(1, 3), Q(2, 5), Q(7, 9))
    slopes = {
        (sigma, height): Q(sigma) / height
        for sigma in (-1, 1)
        for height in positive_heights
    }
    for (sigma_a, height_a), slope_a in slopes.items():
        for (sigma_b, height_b), slope_b in slopes.items():
            if slope_a == slope_b:
                assert sigma_a == sigma_b and height_a == height_b


def check_pair_gap_and_radial_monotonicity():
    four_minus_y2 = sub(FOUR, Y2)
    r_plus_num = mul(four_minus_y2, add(Y, TWO))
    r_minus_num = mul(four_minus_y2, sub(TWO, Y))

    # D(y)D(-y)=y^4+4y^2+16.
    assert mul(D_PLUS, D_MINUS) == P_EVEN

    # r(y)+r(-y)=8(4-y^2)/P_EVEN.  Before the common factor 1/2,
    # the cross numerator is 16(4-y^2).
    cross_num = add(mul(r_plus_num, D_MINUS), mul(r_minus_num, D_PLUS))
    assert cross_num == scale(four_minus_y2, 16)

    # 2-s(y)=2y^2(y^2+8)/P_EVEN, strictly positive for y != 0.
    gap_num = sub(scale(P_EVEN, 2), scale(four_minus_y2, 8))
    assert gap_num == scale(mul(Y2, add(Y2, (Q(8),))), 2)

    # For x=y^2, h(x)=8(4-x)/[x(x^2+4x+16)].  Verify the exact
    # derivative numerator 16(x^3-4x^2-16x-32).
    X = Y  # Reuse the univariate polynomial variable under the name x.
    X2 = mul(X, X)
    denominator = mul(X, add(add(X2, scale(X, 4)), (Q(16),)))
    numerator = scale(sub((Q(4),), X), 8)
    derivative_num = sub(
        mul(derivative(numerator), denominator),
        mul(numerator, derivative(denominator)),
    )
    sign_poly = (Q(-32), Q(-16), Q(-4), Q(1))
    assert derivative_num == scale(sign_poly, 16)

    # p(x)=x^2(x-4)-16x-32 is negative term-by-term on 0<x<4.
    assert sign_poly == add(
        mul(X2, sub(X, (Q(4),))),
        (Q(-32), Q(-16)),
    )


def check_normalized_positive_gap_example():
    # An exact q=6 normalized reality palette: three positive weights with
    # 2 sum w=1.  The aggregate gap remains a positive rational number for
    # any nonzero heights and radial projections.
    weights = [Q(1, 6)] * 3
    assert 2 * sum(weights, Q(0)) == 1
    y_values = [Q(1, 5), Q(1, 3), Q(1, 2)]
    u2_values = [Q(2, 7), Q(3, 8), Q(5, 9)]

    gap = Q(0)
    for weight, y, u2 in zip(weights, y_values, u2_values):
        denominator = y**4 + 4 * y**2 + 16
        gap += weight * u2 * 2 * y**2 * (y**2 + 8) / denominator
    assert gap > 0


def main():
    check_projected_complex_blocks_directly()
    check_round_trip_and_slope()
    check_angular_circular_rigidity()
    check_distinct_height_incompatibility()
    check_pair_gap_and_radial_monotonicity()
    check_normalized_positive_gap_example()
    print("C167 exact reality-paired full-block checks passed.")
    print("  one reality pair on an open angular arc forces circular E")
    print("  distinct charge magnitudes have incompatible source lines/metrics")
    print("  every nondegenerate normalized positive palette retains a strict gap")
    print("  both oscillatory branches detune on every open radial interval")
    print("  static pure-normal terminal star: ruled out in the stated class")


if __name__ == "__main__":
    main()
