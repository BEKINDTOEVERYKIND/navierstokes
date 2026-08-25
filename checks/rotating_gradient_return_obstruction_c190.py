#!/usr/bin/env python3
"""Exact checks for C190's orbit-specific PPRG obstruction.

The computation is closed form.  It verifies the moving-frame Kelvin
generator on C186's rotating-gradient orbit, the return-charge defect, and
the two-block verdict with formal polynomial arithmetic in T=2*pi.  It
does not claim a finite-frequency or complete-class PPRG theorem.
"""

from fractions import Fraction as F
from typing import TypeVar


Scalar = TypeVar("Scalar")
Matrix = tuple[tuple[Scalar, ...], ...]


def transpose(a: Matrix[Scalar]) -> Matrix[Scalar]:
    return tuple(tuple(a[i][j] for i in range(len(a))) for j in range(len(a[0])))


def matmul(a: Matrix[Scalar], b: Matrix[Scalar]) -> Matrix[Scalar]:
    assert len(a[0]) == len(b)
    return tuple(
        tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), start=0) for j in range(len(b[0])))
        for i in range(len(a))
    )


def matadd(a: Matrix[Scalar], b: Matrix[Scalar]) -> Matrix[Scalar]:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
        for i in range(len(a))
    )


def matscale(c: Scalar, a: Matrix[Scalar]) -> Matrix[Scalar]:
    return tuple(tuple(c * entry for entry in row) for row in a)


def eye(n: int, one: Scalar, zero: Scalar) -> Matrix[Scalar]:
    return tuple(
        tuple(one if i == j else zero for j in range(n)) for i in range(n)
    )


def trace(a: Matrix[Scalar]) -> Scalar:
    return sum((a[i][i] for i in range(len(a))), start=0)


def det2(a: Matrix[Scalar]) -> Scalar:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def exact_moving_frame_generator() -> None:
    # At t=0 for p0=-e2.  Columns of E0 are n and k0 x n=-e1.
    a0: Matrix[F] = (
        (F(0), F(-1), F(0)),
        (F(1), F(0), F(0)),
        (F(1), F(0), F(0)),
    )
    k0: Matrix[F] = ((F(0),), (F(-1),), (F(0),))
    e0: Matrix[F] = (
        (F(0), F(-1)),
        (F(0), F(0)),
        (F(1), F(0)),
    )
    omega: Matrix[F] = (
        (F(0), F(-1), F(0)),
        (F(1), F(0), F(0)),
        (F(0), F(0), F(0)),
    )

    kta = matmul(transpose(k0), a0)
    assert kta == ((F(-1), F(0), F(0)),)
    k0_sq = matmul(transpose(k0), k0)[0][0]
    assert k0_sq == 1
    kelvin0 = matadd(matscale(F(-1), a0), matscale(F(2), matmul(k0, kta)))
    assert kelvin0 == (
        (F(0), F(1), F(0)),
        (F(1), F(0), F(0)),
        (F(-1), F(0), F(0)),
    )

    connection = matadd(
        matmul(matmul(transpose(e0), kelvin0), e0),
        matscale(F(-1), matmul(matmul(transpose(e0), omega), e0)),
    )
    assert connection == ((F(0), F(1)), (F(0), F(0)))
    assert matmul(connection, connection) == ((F(0), F(0)), (F(0), F(0)))

    # Check the orthogonal rotational covariance at two exact rational
    # points of the unit circle.  This is a check on signs and frame order,
    # not a numerical approximation to the orbit.
    i3 = eye(3, F(1), F(0))
    for cosine, sine in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13))):
        q: Matrix[F] = (
            (cosine, -sine, F(0)),
            (sine, cosine, F(0)),
            (F(0), F(0), F(1)),
        )
        assert matmul(transpose(q), q) == i3
        at = matmul(matmul(q, a0), transpose(q))
        kt = matmul(q, k0)
        et = matmul(q, e0)
        kelvint = matmul(matmul(q, kelvin0), transpose(q))
        assert matmul(transpose(kt), kt) == ((F(1),),)
        assert matmul(matmul(transpose(et), kelvint), et) == matmul(
            matmul(transpose(e0), kelvin0), e0
        )
        assert at[0][2] == at[1][2] == at[2][2] == 0

    # Keep the passive gradient fixed at g0=e1 and vary the returning
    # horizontal covector.  This checks (2.9), rather than appealing to a
    # covariance that rotates the scalar gradient together with p0.  For a
    # unit h0=(q1,q2), the canonical frame has columns n and h0 x n, and
    # the exact connection is beta E12 with beta=e1.Jh0=-q2.
    for q1, q2 in (
        (F(1), F(0)),
        (F(0), F(1)),
        (F(0), F(-1)),
        (F(3, 5), F(4, 5)),
        (F(-5, 13), F(12, 13)),
    ):
        assert q1**2 + q2**2 == 1
        k: Matrix[F] = ((q1,), (q2,), (F(0),))
        frame: Matrix[F] = (
            (F(0), q2),
            (F(0), -q1),
            (F(1), F(0)),
        )
        assert matmul(transpose(frame), frame) == eye(2, F(1), F(0))
        kta = matmul(transpose(k), a0)
        kelvin = matadd(
            matscale(F(-1), a0),
            matscale(F(2), matmul(k, kta)),
        )
        connection = matadd(
            matmul(matmul(transpose(frame), kelvin), frame),
            matscale(
                F(-1),
                matmul(matmul(transpose(frame), omega), frame),
            ),
        )
        beta = -q2
        assert connection == ((F(0), beta), (F(0), F(0)))
        assert beta == -q2  # e1 dot J h0
        assert beta**2 <= 1


def return_charge_obstruction() -> None:
    # At a full coefficient return R_T=I, C183 gives
    # p(T)=p0-m*T*e1.  Positive T makes return equivalent to m=0.
    def endpoint(p0: tuple[F, F], m: F, period: F) -> tuple[F, F]:
        return (p0[0] - m * period, p0[1])

    for p0 in ((F(0), F(-1)), (F(3, 5), F(4, 5))):
        for period in (F(1), F(7, 3), F(44, 7)):
            assert endpoint(p0, F(0), period) == p0
            for charge in (F(1), F(-1), F(5, 9)):
                assert endpoint(p0, charge, period) != p0

    # Even projective return cannot help when m is nonzero: conservation of
    # the vertical component forces the proportionality factor to be one.
    for charge in (F(1), F(-2, 3)):
        for factor in (F(-1), F(1, 2), F(1), F(3)):
            if charge == factor * charge:
                assert factor == 1


def consecutive_quarter_determinant_gate() -> None:
    # C183 gives det Phi(t1,t0)=|k(t0)|/|k(t1)|.  Hence two
    # determinant-one h-episodes require equality of these three squared
    # norms.  The exact differences factor as displayed below.
    def norm_sq(a: F, b: F, m: F, t: F) -> F:
        return (a - m * t) ** 2 + b**2 + m**2

    # First verify the load-bearing identities universally, by comparing
    # their expanded coefficient dictionaries in the variables (a,b,m,h).
    # If both norm differences vanish, their difference gives
    # 2*m^2*h^2=0; since an episode has h>0, this forces m=0.
    Monomial4 = tuple[int, int, int, int]
    MVPoly = dict[Monomial4, F]

    def mvclean(value: MVPoly) -> MVPoly:
        return {monomial: coefficient for monomial, coefficient in value.items()
                if coefficient}

    def mvvar(index: int) -> MVPoly:
        powers = [0, 0, 0, 0]
        powers[index] = 1
        return {tuple(powers): F(1)}  # type: ignore[dict-item]

    def mvadd(left: MVPoly, right: MVPoly) -> MVPoly:
        result = dict(left)
        for monomial, coefficient in right.items():
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return mvclean(result)

    def mvscale(coefficient: F | int, value: MVPoly) -> MVPoly:
        return mvclean({monomial: F(coefficient) * entry
                        for monomial, entry in value.items()})

    def mvmul(left: MVPoly, right: MVPoly) -> MVPoly:
        result: MVPoly = {}
        for left_monomial, left_coefficient in left.items():
            for right_monomial, right_coefficient in right.items():
                monomial = tuple(
                    left_monomial[index] + right_monomial[index]
                    for index in range(4)
                )
                result[monomial] = result.get(monomial, F(0)) \
                    + left_coefficient * right_coefficient
        return mvclean(result)

    def mvsquare(value: MVPoly) -> MVPoly:
        return mvmul(value, value)

    a_poly, b_poly, m_poly, h_poly = (mvvar(index) for index in range(4))
    mh_poly = mvmul(m_poly, h_poly)

    def symbolic_norm(time_multiple: int) -> MVPoly:
        shifted_a = mvadd(a_poly, mvscale(-time_multiple, mh_poly))
        return mvadd(
            mvadd(mvsquare(shifted_a), mvsquare(b_poly)),
            mvsquare(m_poly),
        )

    norm_zero = symbolic_norm(0)
    norm_one = symbolic_norm(1)
    norm_two = symbolic_norm(2)
    first_formal = mvadd(norm_one, mvscale(-1, norm_zero))
    second_formal = mvadd(norm_two, mvscale(-1, norm_one))
    expected_first = mvmul(
        mh_poly,
        mvadd(mh_poly, mvscale(-2, a_poly)),
    )
    expected_second = mvmul(
        mh_poly,
        mvadd(mvscale(3, mh_poly), mvscale(-2, a_poly)),
    )
    contradiction = mvadd(second_formal, mvscale(-1, first_formal))
    expected_contradiction = mvscale(
        2,
        mvmul(mvsquare(m_poly), mvsquare(h_poly)),
    )
    assert first_formal == expected_first
    assert second_formal == expected_second
    assert contradiction == expected_contradiction
    assert contradiction == {(0, 0, 2, 2): F(2)}

    for a, b, m, h in (
        (F(2, 3), F(5, 7), F(4, 9), F(3, 5)),
        (F(-7, 4), F(1, 6), F(-5, 8), F(11, 9)),
    ):
        first = norm_sq(a, b, m, h) - norm_sq(a, b, m, F(0))
        second = norm_sq(a, b, m, 2 * h) - norm_sq(a, b, m, h)
        assert first == m * h * (m * h - 2 * a)
        assert second == m * h * (3 * m * h - 2 * a)

        # As a polynomial in t, |k(t)|^2 has these exact coefficients.
        coefficients = (a**2 + b**2 + m**2, -2 * a * m, m**2)
        for t in (F(-2), F(0), F(3, 7), F(5)):
            evaluated = sum(coefficients[j] * t**j for j in range(3))
            assert evaluated == norm_sq(a, b, m, t)
        assert coefficients[2] == m**2 > 0

    # If m*h is nonzero, the first equality forces a=m*h/2 while the
    # second forces a=3*m*h/2.  Their difference is m*h, so both cannot
    # hold.  Check the exact generic solution formulas on several
    # nonzero rational parameter pairs.
    for m, h in ((F(1), F(1)), (F(-2, 3), F(7, 5)), (F(11, 13), F(5, 8))):
        first_solution = m * h / 2
        second_solution = 3 * m * h / 2
        assert first_solution != second_solution
        assert second_solution - first_solution == m * h
        assert m * h * (m * h - 2 * first_solution) == 0
        assert m * h * (3 * m * h - 2 * second_solution) == 0


# Dependency-free formal polynomials in one indeterminate T.  A tuple is
# the coefficient list, so (a,b,c) denotes a+b*T+c*T^2.
Poly = tuple[F, ...]


def poly(*coefficients: F | int) -> Poly:
    values = [F(value) for value in coefficients]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [F(0)])


def padd(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return poly(
        *((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n))
    )


def pmul(a: Poly, b: Poly) -> Poly:
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return poly(*out)


def pscale(c: F, a: Poly) -> Poly:
    return poly(*(c * value for value in a))


def pmatadd(a: Matrix[Poly], b: Matrix[Poly]) -> Matrix[Poly]:
    return tuple(
        tuple(padd(a[i][j], b[i][j]) for j in range(len(a[0])))
        for i in range(len(a))
    )


def pmatmul(a: Matrix[Poly], b: Matrix[Poly]) -> Matrix[Poly]:
    assert len(a[0]) == len(b)
    return tuple(
        tuple(
            _poly_sum(pmul(a[i][k], b[k][j]) for k in range(len(b)))
            for j in range(len(b[0]))
        )
        for i in range(len(a))
    )


def _poly_sum(values: object) -> Poly:
    total = poly(0)
    for value in values:  # type: ignore[union-attr]
        total = padd(total, value)
    return total


def pmatscale(c: Poly, a: Matrix[Poly]) -> Matrix[Poly]:
    return tuple(tuple(pmul(c, entry) for entry in row) for row in a)


def ptrace(a: Matrix[Poly]) -> Poly:
    return _poly_sum(a[i][i] for i in range(len(a)))


def pdet2(a: Matrix[Poly]) -> Poly:
    return padd(pmul(a[0][0], a[1][1]), pscale(F(-1), pmul(a[0][1], a[1][0])))


def formal_two_block_verdict() -> None:
    zero, one, period = poly(0), poly(1), poly(0, 1)
    i2: Matrix[Poly] = ((one, zero), (zero, one))
    e12: Matrix[Poly] = ((zero, one), (zero, zero))
    phi = pmatadd(i2, pmatscale(period, e12))
    nilpotent = pmatadd(phi, pmatscale(poly(-1), i2))

    assert pdet2(phi) == one
    assert pmatmul(nilpotent, nilpotent) == ((zero, zero), (zero, zero))

    phi1 = phi
    phi2 = phi
    assert phi1 == phi2
    product = pmatmul(phi2, phi1)
    assert product == ((one, poly(0, 2)), (zero, one))
    assert ptrace(product) == poly(2)

    power = i2
    for n in range(65):
        assert power == ((one, poly(0, n)), (zero, one))
        power = pmatmul(phi, power)

    # Quarter periods are T/4 in the same moving frame.  Their maps are
    # again identical and their two-block trace is exactly two.
    quarter = pmatadd(i2, pmatscale(pscale(F(1, 4), period), e12))
    assert pdet2(quarter) == one
    assert ptrace(pmatmul(quarter, quarter)) == poly(2)


def fixed_success_criteria_fail() -> None:
    epsilon = F(1, 100)

    # The U and V boxes cannot contain the same common-frame matrix: their
    # (1,2) center entries differ by one, larger than twice the radius.
    assert F(1) > 2 * epsilon

    # In the displayed frame Phi_21=0 while V_21=1.
    assert abs(F(0) - F(1)) > epsilon

    # From ||I+N*T|| <= 1+|T| and 2*pi < 44/7, every N-return power obeys
    # the explicit majorant in C190.  The exact matrix power was checked
    # above; here the rational constant ledger is checked separately.
    period_upper = F(44, 7)
    for n in range(129):
        majorant = F(1) + period_upper * n
        assert majorant == F(1) + F(44 * n, 7)
        assert majorant >= 1


def main() -> None:
    exact_moving_frame_generator()
    return_charge_obstruction()
    consecutive_quarter_determinant_gate()
    formal_two_block_verdict()
    fixed_success_criteria_fail()
    print("C190 rotating-gradient return obstruction checks passed")
    print("VERDICT: pre-registered outcome (b) for the chosen C186 orbit")
    print("EXACT: Phi1=Phi2, square-zero gate, and trace(Phi2 Phi1)=2")
    print("BOUNDARY: no complete-class PPRG, finite-frequency, viscous, or UVSR claim")


if __name__ == "__main__":
    main()
