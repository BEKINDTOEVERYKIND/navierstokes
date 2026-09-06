#!/usr/bin/env python3
"""Exact arithmetic for the pressure-resolving A2 curl/WKB bridge.

This checker is dependency-free.  It verifies the algebraic and numerical
constants used by the C194 bridge claim; it does not integrate a
ray, a polarization column, or a PDE.

Norm conventions in the accompanying proof are fixed as follows.

* Vectors use the Euclidean norm.
* Matrices use the induced Euclidean operator norm.
* ``D^m F`` uses the corresponding multilinear operator norm.
* Spatial gradients of vector fields in the curl estimates use the
  Frobenius norm.  The proof charges ``|curl c| <= sqrt(2)|Dc|`` and
  ``|A|_F <= sqrt(2)|A|_op`` because ``A N=0`` has rank at most two.
* Function norms are ordinary Lebesgue L2 norms on R^3.  The periodic A2
  field is used as a bounded coefficient field; no torus normalization is
  invoked.

The exact C193 premises consumed here are

    |DU|_op <= 6,  |D2 U|_mult <= 3 sqrt(6),
    |D3 U|_mult <= 9.

The checker proves the invariant-annulus lower bound, the adapted-frame
flow-jet constants J1/J2, a safe J3 ledger for a second corrector, the
fully expanded first-corrector residual constants, and the C192 power
margins.  The structural curl and pressure identities are also checked at
the scalar/sign level.  A checker cannot replace their displayed
product-rule proof.
"""

from fractions import Fraction as F


# Sparse bivariate polynomials in (X,E), with exact rational coefficients.
def clean(poly):
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def const(value):
    value = F(value)
    return {} if not value else {(0, 0): value}


def var(index):
    exponent = [0, 0]
    exponent[index] = 1
    return {tuple(exponent): F(1)}


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
    for (i, j), left_coefficient in left.items():
        for (k, ell), right_coefficient in right.items():
            monomial = (i + k, j + ell)
            output[monomial] = output.get(monomial, F(0)) + left_coefficient * right_coefficient
    return clean(output)


def power(poly, exponent):
    output = const(1)
    factor = poly
    while exponent:
        if exponent & 1:
            output = mul(output, factor)
        factor = mul(factor, factor)
        exponent //= 2
    return output


# Univariate polynomials in t, coefficients in increasing order.
def padd(left, right):
    output = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        output[index] += F(coefficient)
    for index, coefficient in enumerate(right):
        output[index] += F(coefficient)
    while output and not output[-1]:
        output.pop()
    return output


def pscale(poly, scalar):
    return [F(scalar) * F(coefficient) for coefficient in poly]


def pmul(left, right):
    output = [F(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            output[i + j] += F(left_coefficient) * F(right_coefficient)
    return output


def pintegral(poly):
    return [F(0)] + [F(coefficient, index + 1) for index, coefficient in enumerate(poly)]


def coefficientwise_leq(left, right):
    length = max(len(left), len(right))
    left = list(left) + [F(0)] * (length - len(left))
    right = list(right) + [F(0)] * (length - len(right))
    return all(a <= b for a, b in zip(left, right))


# Generic sparse polynomials and unreduced rational functions.  These are
# used below to verify the action--angle frame identities symbolically,
# without importing a computer-algebra package.  A polynomial monomial is
# an exponent tuple of fixed length; equality of rational functions is
# checked by exact cross multiplication.
def gclean(poly):
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def gconst(value, dimension=5):
    value = F(value)
    return {} if not value else {(0,) * dimension: value}


def gvar(index, dimension=5):
    exponent = [0] * dimension
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def gadd(left, right):
    output = dict(left)
    for monomial, coefficient in right.items():
        output[monomial] = output.get(monomial, F(0)) + coefficient
    return gclean(output)


def gneg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def gsub(left, right):
    return gadd(left, gneg(right))


def gscale(poly, scalar):
    scalar = F(scalar)
    return gclean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def gmul(left, right):
    output = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            output[monomial] = output.get(monomial, F(0)) + left_coefficient * right_coefficient
    return gclean(output)


def gpower(poly, exponent):
    dimension = len(next(iter(poly), (0,) * 5))
    output = gconst(1, dimension)
    factor = poly
    while exponent:
        if exponent & 1:
            output = gmul(output, factor)
        factor = gmul(factor, factor)
        exponent //= 2
    return output


def rat(numerator, denominator=None):
    dimension = len(next(iter(numerator), (0,) * 5))
    return numerator, denominator if denominator is not None else gconst(1, dimension)


def radd(left, right):
    left_numerator, left_denominator = left
    right_numerator, right_denominator = right
    return rat(
        gadd(gmul(left_numerator, right_denominator), gmul(right_numerator, left_denominator)),
        gmul(left_denominator, right_denominator),
    )


def rneg(value):
    numerator, denominator = value
    return rat(gneg(numerator), denominator)


def rsub(left, right):
    return radd(left, rneg(right))


def rmul(left, right):
    return rat(gmul(left[0], right[0]), gmul(left[1], right[1]))


def rscale(value, scalar):
    return rat(gscale(value[0], scalar), value[1])


def requal(left, right):
    return not gsub(gmul(left[0], right[1]), gmul(right[0], left[1]))


def rzero(dimension=5):
    return rat(gconst(0, dimension))


def rdot(left, right):
    output = rzero()
    for left_entry, right_entry in zip(left, right):
        output = radd(output, rmul(left_entry, right_entry))
    return output


def rmatvec(matrix, vector):
    return [rdot(row, vector) for row in matrix]


def rmatmul(left, right):
    columns = list(zip(*right))
    return [[rdot(row, column) for column in columns] for row in left]


def assert_rvector_equal(left, right):
    assert len(left) == len(right)
    assert all(requal(a, b) for a, b in zip(left, right))


def assert_rmatrix_equal(left, right):
    assert len(left) == len(right)
    for left_row, right_row in zip(left, right):
        assert_rvector_equal(left_row, right_row)


def fraction_matmul(left, right):
    columns = list(zip(*right))
    return [
        [sum((F(a) * F(b) for a, b in zip(row, column)), F(0)) for column in columns]
        for row in left
    ]


def fraction_transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def fraction_dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def fraction_cross(left, right):
    return [
        F(left[1]) * F(right[2]) - F(left[2]) * F(right[1]),
        F(left[2]) * F(right[0]) - F(left[0]) * F(right[2]),
        F(left[0]) * F(right[1]) - F(left[1]) * F(right[0]),
    ]


def action_angle_frame_identities():
    """Verify the exact A2 first-integral frame and its shear equation."""
    # Variables are (g1,g2,H11,H12,H22).  Here g=grad_(a,b) f and H=Dg.
    g1, g2, h11, h12, h22 = (gvar(index) for index in range(5))
    one = rat(gconst(1))
    zero = rzero()
    g_vector = [rat(g1), rat(g2)]
    h_value = gadd(gpower(g1, 2), gpower(g2, 2))
    h_rat = rat(h_value)
    h_squared = gpower(h_value, 2)
    jg = [rat(gneg(g2)), rat(g1)]
    velocity = [rscale(jg[0], 3), rscale(jg[1], 3)]  # V=3Jg.
    hessian = [
        [rat(h11), rat(h12)],
        [rat(h12), rat(h22)],
    ]

    # The phase equations are (a',b')=3J grad f, hence E'=g.V=0.
    assert requal(rdot(g_vector, velocity), zero)

    # The physical-to-(a,b,z) map has exact singular values
    # sqrt(3), 1, 1/sqrt(3).  This mechanically records both coordinate
    # factors used when passing phase-coordinate jets back to physical x.
    coordinate = [
        [F(1), F(-1), F(0)],
        [F(0), F(1), F(-1)],
        [F(1, 3), F(1, 3), F(1, 3)],
    ]
    gram = fraction_matmul(coordinate, fraction_transpose(coordinate))
    assert gram == [
        [F(2), F(-1), F(0)],
        [F(-1), F(2), F(0)],
        [F(0), F(0), F(1, 3)],
    ]
    eigenvectors = ([F(1), F(-1), F(0)], [F(1), F(1), F(0)], [F(0), F(0), F(1)])
    eigenvalues = (F(3), F(1), F(1, 3))
    for vector, eigenvalue in zip(eigenvectors, eigenvalues):
        image = [sum((row[index] * vector[index] for index in range(3)), F(0)) for row in gram]
        assert image == [eigenvalue * entry for entry in vector]

    # The axial equation is z'=-sqrt(2)E.  The integer geometry checks
    # N.q_j=0 and |N|^2=3, so N.(N cross grad f)/3=0 while
    # N.(-sqrt(2) f N)/3=-sqrt(2)f.
    n_vector = (F(1), F(1), F(1))
    roots = ((F(1), F(-1), F(0)), (F(0), F(1), F(-1)), (F(1), F(0), F(-1)))
    assert fraction_dot(n_vector, n_vector) == 3
    for root in roots:
        assert fraction_dot(n_vector, root) == 0
    # These scalar triples mechanically produce a'=-3 f_b and b'=3 f_a.
    assert fraction_dot(roots[0], fraction_cross(n_vector, roots[0])) == 0
    assert fraction_dot(roots[1], fraction_cross(n_vector, roots[1])) == 0
    assert fraction_dot(roots[0], fraction_cross(n_vector, roots[1])) == -3
    assert fraction_dot(roots[1], fraction_cross(n_vector, roots[0])) == 3

    # Q=[(V,0),(e,0),e_z], e=g/h.  Check Q^{-1}Q=I exactly and check
    # the orthogonal column lengths |V|^2=9h, |e|^2=h^{-1}.
    e_vector = [rat(g1, h_value), rat(g2, h_value)]
    q_matrix = [
        [velocity[0], e_vector[0], zero],
        [velocity[1], e_vector[1], zero],
        [zero, zero, one],
    ]
    q_inverse = [
        [rat(gneg(g2), gscale(h_value, 3)), rat(g1, gscale(h_value, 3)), zero],
        [rat(g1), rat(g2), zero],
        [zero, zero, one],
    ]
    identity = [[one if i == j else zero for j in range(3)] for i in range(3)]
    assert_rmatrix_equal(rmatmul(q_inverse, q_matrix), identity)
    assert_rmatrix_equal(rmatmul(q_matrix, q_inverse), identity)
    assert requal(rdot(velocity, e_vector), zero)
    assert requal(rdot(g_vector, e_vector), one)
    assert requal(rdot(velocity, velocity), rscale(h_rat, 9))
    assert requal(rdot(e_vector, e_vector), rat(gconst(1), h_value))

    # Along a trajectory, D_t g=H V and DV=3JH.  Direct differentiation
    # of e=g/h verifies
    #
    #   (DV)e-D_t e = alpha' V,
    #   alpha'=[g.Hg-(Jg).H(Jg)]/h^2.
    #
    # This is the exact shear identity behind F=Q_t S Q_0^{-1}.
    dg_dt = rmatvec(hessian, velocity)
    dh_dt = rscale(rdot(g_vector, dg_dt), 2)
    de_dt = [
        rsub(rmul(dg_dt[index], rat(gconst(1), h_value)),
             rmul(g_vector[index], rmul(dh_dt, rat(gconst(1), h_squared))))
        for index in range(2)
    ]
    dv_matrix = [
        [rscale(hessian[1][0], -3), rscale(hessian[1][1], -3)],
        [rscale(hessian[0][0], 3), rscale(hessian[0][1], 3)],
    ]
    dv_e = rmatvec(dv_matrix, e_vector)
    hg = rmatvec(hessian, g_vector)
    hjg = rmatvec(hessian, jg)
    alpha_numerator = rsub(rdot(g_vector, hg), rdot(jg, hjg))
    alpha_prime = rmul(alpha_numerator, rat(gconst(1), h_squared))
    assert_rvector_equal(
        [rsub(dv_e[index], de_dt[index]) for index in range(2)],
        [rmul(alpha_prime, velocity[index]) for index in range(2)],
    )

    # The tangent column also solves the variational equation exactly:
    # D_t V=(DV)V.  This is an independent symbolic check of the first
    # column of the shear representation.
    dt_velocity_from_g = [rscale(dg_dt[1], -3), rscale(dg_dt[0], 3)]
    assert_rvector_equal(dt_velocity_from_g, rmatvec(dv_matrix, velocity))

    # For W=(V,-lambda E), the preceding identities also verify the axial
    # shear in the second column.  Namely DE(e)=g.e=1, so DW(q_2) has
    # axial component -lambda, while DE(V)=g.V=0.  Therefore
    #
    #   q_2(t)+alpha(t)q_1(t)-lambda*t*q_3
    #
    # solves the variational equation when alpha'=the expression above.
    # This is exactly the second column of
    # S=[[1,alpha,0],[0,1,0],[0,-lambda*t,1]].
    assert requal(rdot(g_vector, e_vector), one)
    assert requal(rdot(g_vector, velocity), zero)


def annulus_identity_and_lower_bound():
    """Verify |grad_(a,b) f|^2 > 3/2 on |f| <= 1/10."""
    X = var(0)
    E = var(1)

    # With p=(a+b)/2, X=cos^2 p, y=cos((a-b)/2), and E=f,
    # Z=x*y=(E+4/5-(8/5)X)/2.  The phase-gradient identity is
    # h=2[(1-X)(Z^2/X+(16/5)Z+(64/25)X)+X-Z^2].
    Z = scale(add(add(E, const(F(4, 5))), scale(X, F(-8, 5))), F(1, 2))
    one_minus_X = sub(const(1), X)
    # Form 50 X h without introducing a rational-function class.
    fifty_X_h = scale(
        add(
            add(
                mul(one_minus_X, power(Z, 2)),
                mul(mul(X, one_minus_X), add(scale(Z, F(16, 5)), scale(X, F(64, 25)))),
            ),
            sub(power(X, 2), mul(X, power(Z, 2))),
        ),
        100,
    )
    fifty_X_h_minus_three_halves = sub(fifty_X_h, scale(X, 75))
    expected = add(
        add(
            add(scale(power(X, 3), -128), scale(power(X, 2), 164)),
            add(scale(X, -43), const(16)),
        ),
        add(mul(sub(const(25), scale(X, 50)), power(E, 2)), scale(E, 40)),
    )
    assert fifty_X_h_minus_three_halves == expected

    # The X-polynomial is at least -7 on [0,1].
    # P(X)+7=(1-X)(128X^2-36X+7), and the quadratic is positive.
    polynomial_plus_seven = add(
        add(add(scale(power(X, 3), -128), scale(power(X, 2), 164)), scale(X, -43)),
        const(7),
    )
    factorized_polynomial = mul(
        sub(const(1), X),
        add(add(scale(power(X, 2), 128), scale(X, -36)), const(7)),
    )
    assert polynomial_plus_seven == factorized_polynomial
    sign_factorized = add(
        add(
            const(9),
            mul(
                sub(const(1), X),
                add(add(scale(power(X, 2), 128), scale(X, -36)), const(7)),
            ),
        ),
        add(mul(sub(const(25), scale(X, 50)), power(E, 2)), scale(E, 40)),
    )
    # This is the exact sign-ready form of (1.1), not merely a numerical
    # sample of it:
    # 50X(h-3/2)=9+(1-X)(128X^2-36X+7)+(25-50X)E^2+40E.
    assert expected == sign_factorized
    discriminant = 36**2 - 4 * 128 * 7
    assert discriminant == -2288 < 0
    assert 128 > 0

    # For |E|<=1/10, (25-50X)E^2 >= -1/4 and 40E>=-4.
    numerator_floor = -7 + 16 - F(1, 4) - 4
    assert numerator_floor == F(19, 4) > 0
    # Since 0<X<=1, division by 50X leaves at least 19/200.
    assert F(3, 2) + numerator_floor / 50 == F(319, 200) > F(3, 2)

    # X=0 is actually disjoint from the annulus: cos p=0 gives
    # E=(4/5)cos(2p)=-4/5, independently of the other phase.
    assert abs(F(-4, 5)) > F(1, 10)


def curl_pressure_signs_and_first_jet_constants():
    # (1/i)*i=1 and i*(2i)=-2: these are the two load-bearing signs.
    # Store Gaussian integers as pairs (real,imaginary) to avoid floats.
    def gaussian_mul(left, right):
        a, b = left
        c, d = right
        return (a * c - b * d, a * d + b * c)

    inv_i = (0, -1)
    i_value = (0, 1)
    assert gaussian_mul(inv_i, i_value) == (1, 0)
    assert gaussian_mul(i_value, (0, 2)) == (-2, 0)

    # Exact global coefficient premises and elementary comparisons.
    assert 6 < F(350, 57)
    assert F(350, 57) - 6 == F(8, 57)
    assert F(550, 57) - 6 == F(208, 57)

    # The C1/H0/H1 phase-coordinate chain.  For
    # f=cos(a)+cos(b)+(4/5)cos(a+b), each Hessian row has one diagonal
    # contribution <=9/5 and one off-diagonal contribution <=4/5.
    # Hence H0<=13/5 in induced l2 operator norm by the symmetric row-sum
    # bound.  The third derivative splits into the two unit coordinate
    # cubes and (4/5)(e_a+e_b)^tensor3, so
    # H1<=2+(8/5)sqrt(2)<22/5.
    diagonal_hessian = F(1) + F(4, 5)
    off_diagonal_hessian = F(4, 5)
    hessian = diagonal_hessian + off_diagonal_hessian
    assert diagonal_hessian == F(9, 5)
    assert hessian == F(13, 5)
    # sqrt(2)<3/2 is checked after squaring positive quantities.
    assert 2 < F(3, 2) ** 2
    third_f_radical_majorant = F(2) + F(8, 5) * F(3, 2)
    third_f = F(22, 5)
    assert third_f_radical_majorant == third_f

    # C1=sqrt(6) is a safe upper bound for |g|.  Put
    # p=(a+b)/2.  The only cross term in |g|^2 obeys
    # |sin(a+b)(sin(a)+sin(b))|<=8/(3sqrt(3))<17/10:
    # max_{0<=X<=1} X(1-X)^2=4/27.  Consequently
    # |g|^2<=2+(8/5)(17/10)+32/25=6.
    assert F(64, 27) < F(17, 10) ** 2
    gradient_squared_upper = F(2) + F(8, 5) * F(17, 10) + F(32, 25)
    assert gradient_squared_upper == 6
    c1_squared = F(6)
    assert c1_squared == gradient_squared_upper

    # Phase-coordinate derivative bounds on the invariant annulus.
    h_floor = F(3, 2)
    # |alpha'|<=2H0/h.
    assert 2 * hessian / h_floor == F(52, 15) < F(7, 2)
    # sqrt(h_floor)>6/5 gives h_floor^(3/2)>9/5.
    assert h_floor > F(6, 5) ** 2
    # If n=g.Hg-(Jg).H(Jg), then
    # |n|<=2H0 h, |Dn|<=4H0^2 sqrt(h)+2H1 h, and |Dh|<=2H0 sqrt(h).
    # Differentiating alpha'=n/h^2 therefore gives
    # |D alpha'|<=12H0^2/h^(3/2)+2H1/h.
    derivative_a = 12 * hessian**2 / F(9, 5) + 2 * third_f / h_floor
    assert derivative_a == F(764, 15) < 51

    # Q=[(3Jg,0),(g/|g|^2,0),e_z].
    assert 3 * hessian == F(39, 5) < 8
    assert 3 * hessian / h_floor == F(26, 5) < 6
    assert 8 + 6 == 14
    assert 6 * 14 == 84  # |D Q^{-1}| <= |Q^{-1}|^2 |DQ|.
    # |Q|<=3sqrt(6), |Q^{-1}|<=sqrt(6), so their product is 18.
    assert 3 * 6 == 18
    # sqrt((7/2)^2+2)=sqrt(57)/2<4.
    assert F(7, 2) ** 2 + 2 == F(57, 4) < 16

    # The coordinate map L:x->(a,b,z) has squared singular values
    # 3,1,1/3 (verified symbolically above).  Thus a first derivative
    # pays ||L^-1||||L||=3 and a second derivative pays
    # ||L^-1||||L||^2=3sqrt(3)<6.
    assert 3**2 * 3 < 6**2


def explicit_j1_j2_and_residual_arithmetic():
    one_plus_four_t = [F(1), F(4)]
    t_plus_two_t2 = [F(0), F(1), F(2)]

    # In phase coordinates F_y=Q_t S Q_0^{-1}.  The preceding function
    # proves the exact formula, while the C1/H0/H1 chain gives
    # ||Q||<=3sqrt(6), ||Q^-1||<=sqrt(6), and
    # ||S||<=1+sqrt(alpha^2+2t^2)<=1+4t.  Thus J1_y<=18(1+4t).
    assert 3 * 6 == 18
    j1_phase = pscale(one_plus_four_t, 18)
    # Conjugating by x->(a,b,z) pays the exact factor
    # ||L^-1||||L||=3.
    j1_physical = pscale(j1_phase, 3)
    assert j1_physical == [54, 216]

    # D_{y0} alpha <= integral_0^t 51 J1_y(s) ds
    # =918(t+2t^2).
    d_alpha = pscale(t_plus_two_t2, 918)
    assert 51 * 18 == 918

    # Differentiate Q_t S Q_0^-1.  Its three product-rule terms are
    #
    #  ||DQ|| J1_y ||S|| ||Q^-1||
    #      <=14*18*sqrt(6)(1+4t)^2
    #      <=756(1+4t)^2,
    #  ||Q|| ||D S|| ||Q^-1||
    #      <=18*918(t+2t^2)=16524(t+2t^2),
    #  ||Q|| ||S|| ||DQ^-1||
    #      <=252sqrt(6)(1+4t)<=672(1+4t).
    #
    # The first term uses sqrt(6)<3 and the third sqrt(6)<8/3; both
    # comparisons are verified by squaring positive rationals.
    assert 6 < 3**2
    assert 6 < F(8, 3) ** 2
    assert 14 * 18 * 3 == 756
    assert 18 * 918 == 16_524
    assert 3 * 84 * F(8, 3) == 672

    # Phase-coordinate J2, followed by the safe physical factor 6>3sqrt(3).
    j2_phase = padd(
        padd(pscale(pmul(one_plus_four_t, one_plus_four_t), 756), pscale(t_plus_two_t2, 16524)),
        pscale(one_plus_four_t, 672),
    )
    j2_physical = pscale(j2_phase, 6)
    assert j2_physical == [8568, 151560, 270864]

    one_plus_t = [F(1), F(1)]
    j1_simple = pscale(one_plus_t, 216)
    j2_simple = pscale(pmul(one_plus_t, one_plus_t), 286_992)
    assert coefficientwise_leq(j1_physical, j1_simple)
    assert coefficientwise_leq(j2_physical, j2_simple)

    # Audited coefficient chain for (2.4).  With c=-k cross b/|k|^2,
    #
    # |Dc| <= B1/K_- + 3 K1 B0/K_-^2,
    # |D(D_t c)| <= 4D B1/K_- + 4A1 B0/K_-
    #                  +20D K1 B0/K_-^2.
    #
    # The curl/material commutator together with A curl(c) costs another
    # 2D|Dc|.  Therefore D_t d+Ad has coefficients
    # sqrt(2)*(6D,4A1,26DK1).  Direct differentiation of
    # pi=2i(k.Ab)/|k|^2 contributes (2D,2A1,6DK1).  We retain the exact
    # B1 coefficient and use the audited safe A1/K1 row
    #
    # (6sqrt(2)+2)D B1/K_- +[10A1/K_-+58DK1/K_-^2]B0.       (*)
    # Dc has one Db term; differentiating |k|^-2 supplies two of its
    # three Dk terms.
    dc_b1 = F(1)
    dc_k1 = F(1) + F(2)
    assert (dc_b1, dc_k1) == (1, 3)

    # D_t c is the sum
    # (A^T k cross b)/r +(k cross Ab)/r
    # -2(k cross b)(k.A^T k)/r^2.  Upon one spatial derivative the
    # (B1,A1,K1) coefficient triples are, term by term,
    # (1,1,3), (1,1,3), and (2,2,14).
    dt_c_term_1 = (F(1), F(1), F(1) + F(2))
    dt_c_term_2 = (F(1), F(1), F(1) + F(2))
    # In the last term K1 differentiates k cross b (2), the quadratic
    # scalar k.A^T k (4), and r^-2 (8).
    dt_c_term_3 = (F(2), F(2), F(2) + F(4) + F(8))
    dt_dc_b1, dt_dc_a1, dt_dc_k1 = (
        sum((term[index] for term in (dt_c_term_1, dt_c_term_2, dt_c_term_3)), F(0))
        for index in range(3)
    )
    assert (dt_dc_b1, dt_dc_a1, dt_dc_k1) == (4, 4, 20)
    material_curl_b1 = dt_dc_b1 + 2 * dc_b1
    material_curl_a1 = dt_dc_a1
    material_curl_k1 = dt_dc_k1 + 2 * dc_k1
    assert (material_curl_b1, material_curl_a1, material_curl_k1) == (6, 4, 26)

    # For theta=(k.Ab)/r, Dtheta has (B1,A1,K1)=(1,1,3), with the
    # denominator again contributing two of the K1 units.  pi=2i theta.
    pressure_theta = (F(1), F(1), F(1) + F(2))
    pressure_pi = tuple(2 * coefficient for coefficient in pressure_theta)
    assert pressure_pi == (2, 2, 6)
    def rational_plus_sqrt2_lt(rational_part, radical_part, target):
        """Check a+b*sqrt(2)<target exactly for nonnegative b."""
        rational_part = F(rational_part)
        radical_part = F(radical_part)
        target = F(target)
        assert radical_part >= 0 and target > rational_part
        return 2 * radical_part**2 < (target - rational_part) ** 2

    # After the curl's sqrt(2) and the axial rank-two Frobenius conversion,
    # the DA channel costs 4*2=8, while its scalar pressure gradient costs
    # 2 with no Frobenius conversion.  The Dk channel similarly costs
    # 20*2=40 from curl(D_t c), 2*3*2=12 from the material-curl/A d
    # terms, and 6 from the scalar pressure gradient.
    assert 4 * 2 + pressure_pi[1] == 10
    assert 20 * 2 + 2 * 3 * 2 + pressure_pi[2] == 58

    # Algebraic sentinel: controlling only P_q D beta loses a real
    # longitudinal derivative.  At q=e1, beta=e2, take dq=e2 and
    # d beta=-e1.  Then d(q.beta)=0 and P_q d beta=0, although d beta is
    # nonzero.  Thus q.D beta=-(Dq).beta must be charged separately.
    q_sentinel = (F(1), F(0), F(0))
    beta_sentinel = (F(0), F(1), F(0))
    dq_sentinel = (F(0), F(1), F(0))
    dbeta_sentinel = (F(-1), F(0), F(0))
    assert fraction_dot(q_sentinel, beta_sentinel) == 0
    assert (
        fraction_dot(dq_sentinel, beta_sentinel)
        + fraction_dot(q_sentinel, dbeta_sentinel)
    ) == 0
    projected_dbeta = tuple(
        dbeta_sentinel[index]
        - q_sentinel[index] * fraction_dot(q_sentinel, dbeta_sentinel)
        for index in range(3)
    )
    assert projected_dbeta == (0, 0, 0)
    assert fraction_dot(dbeta_sentinel, dbeta_sentinel) == 1

    # The Kelvin generator is L=(-I+2nn^T)A.  The reflection has norm one,
    # while its derivative costs 4|Dn|.  With |Dn|<=J1 J2 this gives
    # |D_a L|<=A1 J1+4D J1 J2
    #            =3sqrt(6)J1+24J1J2.
    assert 4 * 6 == 24

    # D_b = integral (3sqrt(6) J1 + 24 J1 J2).
    # Use sqrt(6)<5/2 and the simple polynomial envelopes.
    assert 3 * F(5, 2) * 216 == 1620
    cubic_coefficient = 24 * 216 * 286_992
    assert cubic_coefficient == 1_487_766_528
    d_b_constant = F(1620, 2) + F(cubic_coefficient, 4)
    assert d_b_constant == 371_942_442

    # The transverse beta derivative is controlled by D_b.  Differentiating
    # q.beta=0 supplies the missing longitudinal term J2, and the declared
    # Frobenius conversion gives the honest full-gradient estimate
    #
    # B1<=e^(6t)|b_*|[
    #       J1 eps^-1||Dchi||+sqrt(2)J1(D_b+J2)||chi||].
    #
    # Before this Frobenius factor, (*) has B1 coefficient
    # 6(6sqrt(2)+2)=36sqrt(2)+12<90.  After multiplication by sqrt(2)
    # it is 72+12sqrt(2)<90.  The direct K1 row is 58*6=348; adding the
    # longitudinal J2 contribution (<90) gives <438<462.  Thus the final
    # displayed 90/30/462 majorant is unchanged, but now includes the
    # load-bearing longitudinal derivative.
    d_global = 6
    a1_global_without_sqrt6 = 3
    assert rational_plus_sqrt2_lt(12, 36, 90)
    assert rational_plus_sqrt2_lt(72, 12, 90)
    a1_substituted = 10 * a1_global_without_sqrt6
    direct_k1_substituted = 58 * d_global
    assert a1_substituted == 30
    assert direct_k1_substituted == 348
    # 348+(72+12sqrt(2))=420+12sqrt(2)<438<462.
    assert rational_plus_sqrt2_lt(direct_k1_substituted + 72, 12, 438)
    assert 438 < 462
    assert j1_physical[0] >= 1

    # First-corrector residual before time integration:
    # 90 J1^2 eps^{-1} ||D chi||
    # +(90 J1^2 D_b + 30sqrt(6)J1 + 462J1^3J2)||chi||.
    gradient_constant = 90 * 216**2
    assert gradient_constant == 4_199_040
    amplitude_terms = (
        gradient_constant * 371_942_442,
        30 * F(5, 2) * 216,
        462 * 216**3 * 286_992,
    )
    assert amplitude_terms == (
        1_561_801_191_655_680,
        16_200,
        1_336_204_776_259_584,
    )
    amplitude_sum = sum(amplitude_terms, F(0))
    amplitude_constant = 2_898_006_000_000_000
    assert amplitude_sum == 2_898_005_967_931_464 < amplitude_constant

    # One time integration adds one polynomial power but no exponential:
    # e^{6(t-s)}e^{6s}=e^{6t} exactly.
    assert 6 - 6 == 0


def second_corrector_j3_ledger():
    """Safe explicit constants needed by, but not a proof of, corrector two."""
    hessian = F(13, 5)
    third_f = F(22, 5)
    fourth_f = F(26, 5)
    g_upper = F(5, 2)  # sqrt(6)<5/2.
    h_floor = F(3, 2)

    # Derivatives of h=|g|^2 and r=h^{-1}.
    dh = 2 * g_upper * hessian
    d2h = 2 * (hessian**2 + g_upper * third_f)
    assert dh == 13
    assert d2h == F(888, 25) < 36
    dr = F(4, 9) * dh
    d2r = 2 * F(8, 27) * dh**2 + F(4, 9) * 36
    assert dr == F(52, 9) < 6
    assert d2r < 117

    # e=g/h and Q second derivatives.
    d2e = F(2, 3) * third_f + 2 * 6 * hessian + g_upper * 117
    assert d2e < 327
    d2v = 3 * third_f
    assert d2v == F(66, 5) < 14
    d2q = 341
    assert 327 + 14 == d2q
    d2q_inverse = 2 * 15 * 14**2 + 6 * d2q
    assert d2q_inverse == 7926

    # Second derivative of a=(gHg-(Jg)H(Jg))/h^2.
    d2numerator = (
        12 * hessian * third_f * g_upper
        + 2 * fourth_f * g_upper**2
        + 4 * hessian**3
    )
    assert d2numerator < 479
    dnumerator = 4 * hessian**2 * g_upper + 2 * third_f * g_upper**2
    assert dnumerator < 123
    numerator = 2 * hessian * g_upper**2
    assert numerator == F(65, 2)
    d_h_minus_2 = 2 * F(8, 27) * dh
    assert d_h_minus_2 < 8
    d2_h_minus_2 = 6 * F(16, 81) * dh**2 + 2 * F(8, 27) * 36
    assert d2_h_minus_2 < 222
    d2a = F(4, 9) * 479 + 2 * 8 * 123 + F(65, 2) * 222
    assert d2a < 9396

    # An explicit, unexpanded phase-coordinate J3 formula.  Here j1,j2
    # are the phase-coordinate bounds already proved above.
    one_plus_four_t = [F(1), F(4)]
    t_plus_two_t2 = [F(0), F(1), F(2)]
    j1 = pscale(one_plus_four_t, 18)
    j2 = padd(
        padd(pscale(pmul(one_plus_four_t, one_plus_four_t), 756), pscale(t_plus_two_t2, 16524)),
        pscale(one_plus_four_t, 672),
    )
    d_alpha = pscale(t_plus_two_t2, 918)
    d2_alpha = pintegral(padd(pscale(pmul(j1, j1), 9396), pscale(j2, 51)))

    # D^2(Q_t S Q_0^{-1}), term by term; physical J3 is nine times this.
    a2 = padd(pscale(pmul(j1, j1), 341), pscale(j2, 14))
    j3_phase = padd(
        padd(
            padd(pscale(pmul(a2, one_plus_four_t), 3), pscale(d2_alpha, 24)),
            pscale(one_plus_four_t, 63_408),
        ),
        padd(
            padd(pscale(pmul(j1, d_alpha), 84), pscale(pmul(j1, one_plus_four_t), 2352)),
            pscale(d_alpha, 1344),
        ),
    )
    j3_physical = pscale(j3_phase, 9)
    assert all(coefficient >= 0 for coefficient in j3_physical)
    assert len(j3_physical) <= 4  # at most cubic growth.
    j3_simple = pscale(pmul(pmul([F(1), F(1)], [F(1), F(1)]), [F(1), F(1)]), 4_031_918_208)
    assert coefficientwise_leq(j3_physical, j3_simple)


def clock_margins():
    # At t <= (57/400)log Q+76/25 and Gamma=6:
    exponent_charge = F(57, 400) * 6
    assert exponent_charge == F(171, 200)
    # sqrt(h)=Q^{-1/2}, divided by the Q^{3/8} signal.
    assert F(1, 2) + F(3, 8) - exponent_charge == F(1, 50)
    # h=Q^{-1}, divided by the same signal.
    assert 1 + F(3, 8) - exponent_charge == F(13, 25)
    # The C193 concentration width eps=Q^-1/4 gives h/eps=Q^-3/4.
    assert F(3, 4) + F(3, 8) - exponent_charge == F(27, 100)


def main():
    action_angle_frame_identities()
    annulus_identity_and_lower_bound()
    curl_pressure_signs_and_first_jet_constants()
    explicit_j1_j2_and_residual_arithmetic()
    second_corrector_j3_ledger()
    clock_margins()
    print("C194 A2 curl/WKB bridge arithmetic: PASS")
    print("ANNULUS: |f|<=1/10 gives |grad_(a,b) f|^2 >= 319/200")
    print("FRAME: action-angle equations, Q/Q^-1, alpha shear checked exactly")
    print("JETS: J1,J2 explicit; J3 <= 4031918208*(1+t)^3 auxiliary ledger checked")
    print("FIRST CORRECTOR: hbar*(1+t)^7*exp(6t), Gamma=6")
    print("CLOCK: eps=Q^-1/4 margin Q^-27/100; eps=Q^-1/2 margin Q^-1/50")


if __name__ == "__main__":
    main()
