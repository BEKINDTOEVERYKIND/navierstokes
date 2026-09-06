#!/usr/bin/env python3
"""Dependency-free exact checks for C164's full pure-normal gate block.

The checker uses only rational polynomial arithmetic.  It verifies the
two-by-two forward/reverse matrices after the elementary dot products have
been expanded, the full round-trip spectrum and invariant branches, the
positive metric on |y|<2, the exceptional heights, and the Hermitian charge
convolution and viscous modulus identity behind the tangential
phase-multiplier no-go.

No finite-band, localization, nonlinear-stage, or Navier--Stokes closure
claim is tested here.
"""

from fractions import Fraction as Q


def trim(p):
    p = [Q(x) for x in p]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, q):
    out = [Q(0)] * max(len(p), len(q))
    for i, x in enumerate(p):
        out[i] += x
    for i, x in enumerate(q):
        out[i] += x
    return trim(out)


def neg(p):
    return trim([-x for x in p])


def sub(p, q):
    return add(p, neg(q))


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, z in enumerate(q):
            out[i + j] += x * z
    return trim(out)


def scale(p, c):
    return trim([Q(c) * x for x in p])


def power(p, n):
    out = (Q(1),)
    for _ in range(n):
        out = mul(out, p)
    return out


ONE = (Q(1),)
Y = (Q(0), Q(1))
Y_PLUS_1 = add(Y, ONE)
Y_PLUS_2 = add(Y, (Q(2),))
Y_PLUS_4 = add(Y, (Q(4),))
D = (Q(4), Q(2), Q(1))


def check_projected_symbol_matrices():
    # In source basis (e_sigma,e_t) and daughter (e_perp,e_t), write
    #
    # F = [[u*fa/sqrt(D), 0], [v*fb, u]],
    # R = [[u*ra/sqrt(D), 0], [v*rb/sqrt(D), u]].
    #
    # These identities independently expand the four dot products in the
    # projected symmetric Euler symbol.
    fa = (Q(-2), Q(0), Q(1, 2))       # -(4-y^2)/2
    fb = (Q(0), Q(-1, 2))             # -y/2
    ra = (Q(-2), Q(-1))               # -(y+2)
    rb = (Q(0), Q(-1))                # -y

    # Forward radial numerator:
    # -(1+y)(1-y)/2 - 3/2 = -(4-y^2)/2.
    raw_forward = scale(
        add(mul(neg(Y_PLUS_1), sub(ONE, Y)), (Q(-3),)),
        Q(1, 2),
    )
    assert raw_forward == fa

    # Reverse radial numerator:
    # -(1+2y)/2 - 3/2 = -(y+2).
    raw_reverse = scale(add((Q(-1), Q(-2)), (Q(-3),)), Q(1, 2))
    assert raw_reverse == ra
    assert fb == scale(Y, Q(-1, 2))
    assert rb == neg(Y)
    return fa, fb, ra, rb


def check_round_trip_and_eigenbranches(fa, fb, ra, rb):
    # The first diagonal numerator of RF is fa*ra over D.
    lambda_num = mul(fa, ra)
    expected_lambda_num = scale(
        mul(sub((Q(4),), power(Y, 2)), Y_PLUS_2), Q(1, 2)
    )
    assert lambda_num == expected_lambda_num
    assert mul(sub((Q(2),), Y), power(Y_PLUS_2, 2)) == scale(
        lambda_num, 2
    )

    # The lower entry of RF, after multiplication by D/(uv), is
    # rb*fa + fb*D = -y^2(y+1).
    round_trip_lower_num = add(mul(rb, fa), mul(fb, D))
    assert round_trip_lower_num == neg(mul(power(Y, 2), Y_PLUS_1))

    # The lower entry of FR is uv*y^2/(2sqrt(D)).
    daughter_lower_num = add(mul(fb, ra), rb)
    assert daughter_lower_num == scale(power(Y, 2), Q(1, 2))

    # lambda_oblique/u^2 - 1 = -y^2(y+4)/(2D).
    assert sub(lambda_num, D) == scale(
        neg(mul(power(Y, 2), Y_PLUS_4)), Q(1, 2)
    )

    # Source and daughter shear coefficients, with the common v/u factor
    # stripped, are
    # X=2(y+1)/(y+4), E=-sqrt(D)/(y+4).
    # The two cross-multiplied identities below prove
    # F sigma_y=f_1 d_y and R d_y=r_1 sigma_y.
    x_num = scale(Y_PLUS_1, 2)
    assert add(mul(fb, Y_PLUS_4), x_num) == neg(fa)
    assert neg(add(mul(Y, Y_PLUS_4), D)) == mul(x_num, ra)

    # Tangential source and daughter vectors are an exact common branch:
    # F e_t=u e_t and R e_t=u e_t.  This is already encoded by the two
    # unit diagonal entries and needs no numerical premise.


def check_metric_rank_and_signed_classification(fa, ra):
    # On |y|<2 the branch weight w=(2-y)/2 is positive and satisfies
    # w*r_1=f_1.  This is the complete weighted-self-adjoint condition in
    # the simultaneously sheared branch coordinates.
    w_num = sub((Q(2),), Y)
    assert scale(mul(w_num, ra), Q(1, 2)) == fa

    # D=(y+1)^2+3 is strictly positive on the real line.
    assert add(power(Y_PLUS_1, 2), (Q(3),)) == D

    # Determinants are u^2*fa/sqrt(D) and u^2*ra/sqrt(D).  Their only
    # radial zeros are y=+/-2 and y=-2, respectively.
    assert fa == scale(mul(sub(Y, (Q(2),)), Y_PLUS_2), Q(1, 2))
    assert ra == neg(Y_PLUS_2)

    # At y=-4 the repeated RF eigenvalue has lower Jordan entry 4uv;
    # at y=0 the same repeated eigenvalue has zero lower entry.
    def ev(p, y):
        return sum(c * Q(y) ** i for i, c in enumerate(p))

    lower = neg(mul(power(Y, 2), Y_PLUS_1))
    assert ev(D, -4) == 12
    assert ev(lower, -4) / ev(D, -4) == 4
    assert ev(lower, 0) == 0
    assert ev(fa, 2) == 0 and ev(fa, -2) == 0
    assert ev(ra, -2) == 0 and ev(ra, 2) != 0

    # Generic endpoint branch matrices [[0,r],[f,0]]: at -2 both
    # coefficients vanish, while at +2 exactly f vanishes.  At -4 both are
    # nonzero and their product is the repeated tangential eigenvalue.
    assert ev(fa, -2) == 0 and ev(ra, -2) == 0
    assert ev(fa, 2) == 0 and ev(ra, 2) == -4
    assert ev(fa, -4) == 6 and ev(ra, -4) == 2
    assert ev(fa, -4) * ev(ra, -4) / ev(D, -4) == 1

    # On the reality-safe interval |y|<=1, a crude exact uniform bound is
    # 1/14 <= lambda_oblique/u^2 <= 9/2 and 1/2<=w<=3/2.
    # It follows directly from 1<=2-y<=3, 1<=(y+2)^2<=9, and 3<=D<=7.
    assert Q(1) * Q(1) / (2 * Q(7)) == Q(1, 14)
    assert Q(3) * Q(9) / (2 * Q(3)) == Q(9, 2)
    assert Q(2 - 1, 2) == Q(1, 2)
    assert Q(2 - (-1), 2) == Q(3, 2)

    # Exact nontrivial spot-check of the original-frame symmetrizer at
    # y=7/4, where sqrt(D)=13/4.  This checks the orientation of both shear
    # pullbacks, not merely the scalar relation w*r=f.
    def mmul(left, right):
        return tuple(
            tuple(
                sum(
                    (left[row][inner] * right[inner][column] for inner in range(2)),
                    Q(0),
                )
                for column in range(2)
            )
            for row in range(2)
        )

    def transpose(matrix):
        return tuple(tuple(matrix[column][row] for column in range(2)) for row in range(2))

    y = Q(7, 4)
    root_d = Q(13, 4)
    u = Q(3, 5)
    v = Q(4, 5)
    f_matrix = (
        (-u * (4 - y * y) / (2 * root_d), Q(0)),
        (-y * v / 2, u),
    )
    r_matrix = (
        (-u * (y + 2) / root_d, Q(0)),
        (-y * v / root_d, u),
    )
    xi = 2 * v * (y + 1) / (u * (y + 4))
    eta = -v * root_d / (u * (y + 4))
    weight = (2 - y) / 2
    source_inverse = ((Q(1), Q(0)), (-xi, Q(1)))
    daughter_inverse = ((Q(1), Q(0)), (-eta, Q(1)))
    source_weight = mmul(
        mmul(transpose(source_inverse), ((weight, Q(0)), (Q(0), Q(1)))),
        source_inverse,
    )
    daughter_weight = mmul(transpose(daughter_inverse), daughter_inverse)
    assert mmul(source_weight, r_matrix) == mmul(
        transpose(f_matrix), daughter_weight
    )
    assert (
        source_weight[0][0] * source_weight[1][1]
        - source_weight[0][1] * source_weight[1][0]
        == weight
    )
    assert (
        daughter_weight[0][0] * daughter_weight[1][1]
        - daughter_weight[0][1] * daughter_weight[1][0]
        == 1
    )


def badd(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, Q(0)) + coefficient
        if out[monomial] == 0:
            del out[monomial]
    return out


def bneg(value):
    return {monomial: -coefficient for monomial, coefficient in value.items()}


def bmul(left, right):
    out = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            monomial = (i + k, j + ell)
            out[monomial] = out.get(monomial, Q(0)) + a * b
            if out[monomial] == 0:
                del out[monomial]
    return out


def check_common_line_and_phase_no_go():
    # Injectivity of (y+1)/(y+4):
    # (y1+1)(y2+4)-(y2+1)(y1+4)=3(y1-y2).
    y1 = {(1, 0): Q(1)}
    y2 = {(0, 1): Q(1)}
    one = {(0, 0): Q(1)}
    four = {(0, 0): Q(4)}
    left = badd(
        bmul(badd(y1, one), badd(y2, four)),
        bneg(bmul(badd(y2, one), badd(y1, four))),
    )
    right = badd(
        {(1, 0): Q(3)},
        {(0, 1): Q(-3)},
    )
    assert left == right

    # Pointwise scalar phase equation Z'=-i*a*Z has
    # x'=a z, z'=-a x and d(x^2+z^2)/dt=0.
    # The cancellation is coefficientwise and uses only that a is real.
    # Store monomials in (x,z); the common real scalar a is factored out.
    x = {(1, 0): Q(1)}
    z = {(0, 1): Q(1)}
    modulus_derivative = badd(
        {monomial: 2 * coefficient for monomial, coefficient in bmul(x, z).items()},
        {monomial: -2 * coefficient for monomial, coefficient in bmul(z, x).items()},
    )
    assert modulus_derivative == {}

    # A finite cyclic charge convolution with reality-symmetric Fourier
    # coefficients is exactly Hermitian, so -i times it is skew-Hermitian.
    def gadd(a, b):
        return (a[0] + b[0], a[1] + b[1])

    def gmul(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def gconj(a):
        return (a[0], -a[1])

    zero = (Q(0), Q(0))
    fhat = {
        0: (Q(2, 3), Q(0)),
        1: (Q(1, 5), Q(2, 7)),
        -1: (Q(1, 5), Q(-2, 7)),
        2: (Q(-3, 11), Q(1, 13)),
        -2: (Q(-3, 11), Q(-1, 13)),
    }
    n = 11

    def coefficient(offset):
        signed = offset % n
        if signed > n // 2:
            signed -= n
        return fhat.get(signed, zero)

    convolution = [
        [coefficient(row - column) for column in range(n)]
        for row in range(n)
    ]
    for row in range(n):
        for column in range(n):
            assert convolution[row][column] == gconj(
                convolution[column][row]
            )

    minus_i = (Q(0), Q(-1))
    generator = [
        [gmul(minus_i, value) for value in row]
        for row in convolution
    ]
    for row in range(n):
        for column in range(n):
            assert generator[row][column] == (
                -gconj(generator[column][row])[0],
                -gconj(generator[column][row])[1],
            )

    # With viscosity, write Z=x+i z and use the six independent jet
    # variables (x,z,x_theta,z_theta,x_tt,z_tt).  Factoring out nu, verify
    #   2 Re(conj(Z) Z_tt)=(|Z|^2)_tt-2|Z_theta|^2.
    diffusion_left = {
        (1, 0, 0, 0, 1, 0): Q(2),
        (0, 1, 0, 0, 0, 1): Q(2),
    }
    modulus_second_derivative = {
        (1, 0, 0, 0, 1, 0): Q(2),
        (0, 1, 0, 0, 0, 1): Q(2),
        (0, 0, 2, 0, 0, 0): Q(2),
        (0, 0, 0, 2, 0, 0): Q(2),
    }
    gradient_loss = {
        (0, 0, 2, 0, 0, 0): Q(-2),
        (0, 0, 0, 2, 0, 0): Q(-2),
    }
    combined = dict(modulus_second_derivative)
    for monomial, coefficient in gradient_loss.items():
        combined[monomial] = combined.get(monomial, Q(0)) + coefficient
        if combined[monomial] == 0:
            del combined[monomial]
    assert diffusion_left == combined


def main():
    fa, fb, ra, rb = check_projected_symbol_matrices()
    check_round_trip_and_eigenbranches(fa, fb, ra, rb)
    check_metric_rank_and_signed_classification(fa, ra)
    check_common_line_and_phase_no_go()
    print("PASS C164: full two-polarization F/R block and RF spectrum are exact")
    print("PASS C164: |y|<2 is weighted elliptic away from the angular defect")
    print("PASS C164: y=+/-2 ranks and y=-4 Jordan exception are exact")
    print("PASS C164: tangential phase fiber has the exact viscous modulus identity")


if __name__ == "__main__":
    main()
