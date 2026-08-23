#!/usr/bin/env python3
"""Exact arithmetic checks for C186's conditional unipotent PPRG lemma.

The script certifies the matrix witness, rational Lyapunov lower bounds,
entrywise robustness box, passive 2D3C coefficient identities, and the
conditional common-conjugacy formulas.  It does not prove that C183's
Kelvin maps are unipotent or realize the two matrix boxes in one PDE orbit.
"""

from fractions import Fraction as F


Matrix = tuple[tuple[F, F], tuple[F, F]]
Vector = tuple[F, F]


I: Matrix = ((F(1), F(0)), (F(0), F(1)))
ZERO: Matrix = ((F(0), F(0)), (F(0), F(0)))


def add(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def sub(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(a[i][j] - b[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def mul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def mv(a: Matrix, v: Vector) -> Vector:
    return tuple(sum(a[i][j] * v[j] for j in range(2)) for i in range(2))  # type: ignore[return-value]


def det(a: Matrix) -> F:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def trace(a: Matrix) -> F:
    return a[0][0] + a[1][1]


def scale(c: F, a: Matrix) -> Matrix:
    return tuple(tuple(c * a[i][j] for j in range(2)) for i in range(2))  # type: ignore[return-value]


def transpose(a: Matrix) -> Matrix:
    return ((a[0][0], a[1][0]), (a[0][1], a[1][1]))


def witness_and_flags() -> None:
    u: Matrix = ((F(1), F(1)), (F(0), F(1)))
    v: Matrix = ((F(1), F(0)), (F(1), F(1)))
    nu = sub(u, I)
    nv = sub(v, I)

    assert mul(nu, nu) == ZERO
    assert mul(nv, nv) == ZERO
    assert det(u) == det(v) == 1

    # The Gram matrices have characteristic polynomial x^2-3x+1, so the
    # squared operator norm is (3+sqrt(5))/2 as used in (1.4).
    for block in (u, v):
        gram = mul(transpose(block), block)
        assert trace(gram) == 3
        assert det(gram) == 1

    e1: Vector = (F(1), F(0))
    e2: Vector = (F(0), F(1))
    assert mv(nu, e1) == (0, 0)
    assert mv(nu, e2) != (0, 0)
    assert mv(nv, e2) == (0, 0)
    assert mv(nv, e1) != (0, 0)

    # tau is the exact common-flag decision scalar; here tau=1>0.
    tau = trace(mul(nu, nv))
    assert tau == 1
    assert (2 + tau) ** 2 - 4 == tau * (tau + 4)

    p = mul(u, v)
    assert p == ((F(2), F(1)), (F(1), F(1)))
    assert det(p) == 1
    assert trace(p) == 3

    w: Vector = (F(13, 8), F(1))
    pw = mv(p, w)
    lower = (F(34, 13) * w[0], F(34, 13) * w[1])
    assert pw == (F(17, 4), F(21, 8))
    assert pw[0] >= lower[0] and pw[1] >= lower[1]

    # ||U||_2^2=(3+sqrt(5))/2 < (13/8)^2.  The last inequality follows
    # from sqrt(5)<73/32, certified by squaring positive rationals.
    assert F(5) < F(73, 32) ** 2
    assert (F(3) + F(73, 32)) / 2 == F(13, 8) ** 2


def logarithmic_growth_bounds() -> None:
    # log(x)=2(z+z^3/3+z^5/5+...), z=(x-1)/(x+1).
    z = F(21, 47)  # x=34/13
    lower = 2 * (z + z**3 / 3 + z**5 / 5)
    assert lower - F(24, 25) == F(1185042, 5733625175)
    assert lower > F(24, 25)
    assert lower / 2 > F(12, 25)

    # Exact schedule coefficients for q^(3/8) and q^(1/2).
    assert F(24, 25) * F(25, 64) == F(3, 8)
    assert F(24, 25) * F(25, 48) == F(1, 2)


def robustness_box() -> None:
    epsilon = F(1, 100)
    trace_error = 6 * epsilon + 4 * epsilon**2
    assert trace_error == F(151, 2500)
    trace_lower = 3 - trace_error
    assert trace_lower == F(7349, 2500)
    assert trace_lower > F(29, 10)

    # For det=1, trace>29/10 implies lambda_+>5/2 because
    # lambda + lambda^-1 is increasing for lambda>1.
    assert F(5, 2) + F(2, 5) == F(29, 10)

    z = F(3, 7)  # x=5/2
    log_lower = 2 * (z + z**3 / 3)
    assert log_lower == F(312, 343)
    assert log_lower > F(9, 10)
    assert log_lower / 2 > F(9, 20)


# Tiny exact commutative polynomial helper for the steady-Euler identities.
Monomial = tuple[int, int, int, int]  # powers of sx,cx,sy,cy
Polynomial = dict[Monomial, F]


def poly_term(index: int, coefficient: F = F(1)) -> Polynomial:
    powers = [0, 0, 0, 0]
    powers[index] = 1
    return {tuple(powers): coefficient}  # type: ignore[dict-item]


def poly_add(a: Polynomial, b: Polynomial) -> Polynomial:
    out = dict(a)
    for monomial, coefficient in b.items():
        out[monomial] = out.get(monomial, F(0)) + coefficient
        if out[monomial] == 0:
            del out[monomial]
    return out


def poly_mul(a: Polynomial, b: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            monomial = tuple(ma[i] + mb[i] for i in range(4))
            out[monomial] = out.get(monomial, F(0)) + ca * cb
    return {m: c for m, c in out.items() if c}


def poly_scale(c: F, a: Polynomial) -> Polynomial:
    return {m: c * value for m, value in a.items() if c * value}


def passive_pde_example() -> None:
    sx, cx, sy, cy = (poly_term(i) for i in range(4))
    v1 = poly_scale(F(-1), sy)
    v2 = sx
    omega_x = poly_scale(F(-1), sx)
    omega_y = poly_scale(F(-1), sy)
    advected_omega = poly_add(poly_mul(v1, omega_x), poly_mul(v2, omega_y))
    assert advected_omega == {}

    # (v.grad)v=(-sx*cy,-sy*cx); grad(-cx*cy)=(sx*cy,cx*sy).
    convective_1 = poly_scale(F(-1), poly_mul(sx, cy))
    convective_2 = poly_scale(F(-1), poly_mul(sy, cx))
    pressure_x = poly_mul(sx, cy)
    pressure_y = poly_mul(cx, sy)
    assert poly_add(convective_1, pressure_x) == {}
    assert poly_add(convective_2, pressure_y) == {}

    j: Matrix = ((F(0), F(-1)), (F(1), F(0)))
    assert mul(j, j) == scale(F(-1), I)
    e1: Vector = (F(1), F(0))
    r_zero = I
    r_quarter = j
    assert mv(r_zero, e1) == (F(1), F(0))
    assert mv(r_quarter, e1) == (F(0), F(1))


def common_conjugacy_boundary() -> None:
    # Verify C^{-j} M C^j formula on several exact rational instances.
    e12: Matrix = ((F(0), F(1)), (F(0), F(0)))
    for s, j_index, a, b, c in (
        (F(2, 3), 7, F(3, 5), F(-4, 7), F(5, 11)),
        (F(-5, 4), 8, F(-2, 9), F(7, 6), F(-3, 8)),
    ):
        m: Matrix = ((a, b), (c, -a))
        c_power = add(I, scale(F(j_index) * s, e12))
        c_inverse_power = add(I, scale(-F(j_index) * s, e12))
        actual = mul(mul(c_inverse_power, m), c_power)
        expected: Matrix = (
            (
                a - F(j_index) * s * c,
                b
                + 2 * F(j_index) * s * a
                - F(j_index) ** 2 * s**2 * c,
            ),
            (c, -a + F(j_index) * s * c),
        )
        assert actual == expected

    # Exact finite-difference coefficients behind |c| and |a| bounds.
    s = F(7, 5)
    a, b, c = F(2, 3), F(-4, 9), F(5, 8)

    def q(j_index: F) -> F:
        return b + 2 * j_index * s * a - j_index**2 * s**2 * c

    j_total = F(10)
    second_difference = q(j_total) - 2 * q(j_total / 2) + q(F(0))
    assert second_difference == -j_total**2 * s**2 * c / 2
    first_difference = q(j_total) - q(F(0))
    assert first_difference == 2 * j_total * s * a - j_total**2 * s**2 * c


def main() -> None:
    witness_and_flags()
    logarithmic_growth_bounds()
    robustness_box()
    passive_pde_example()
    common_conjugacy_boundary()
    print("C186 conditional unipotent PPRG checks passed")
    print("BOUNDARY: C183 does not imply unipotence of Kelvin blocks")
    print("BOUNDARY: the two robust matrix boxes are not realized by one PDE orbit")
    print("BOUNDARY: no finite-frequency, viscous, UVSR, or singularity theorem")


if __name__ == "__main__":
    main()
