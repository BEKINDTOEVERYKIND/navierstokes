#!/usr/bin/env python3
"""Dependency-free checks for C175's resonant invariant-graph obstruction.

The script checks exact A2 label/polarization arithmetic, material
resonance under rational covector maps, the graph Riccati/Melnikov matrix
identity over the rationals, the full-domain/restricted-domain distinction,
and the scalar Leray--heat graph formulas.

It does not prove the actual broad-packet Melnikov estimate, localization,
MCKC, BAFL, an unforced stage, or blow-up.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import exp, factorial, isclose


Vec = tuple[F, ...]
Mat = tuple[tuple[F, ...], ...]


def vadd(a: Vec, b: Vec) -> Vec:
    return tuple(x + y for x, y in zip(a, b))


def vsub(a: Vec, b: Vec) -> Vec:
    return tuple(x - y for x, y in zip(a, b))


def vscale(c: F, a: Vec) -> Vec:
    return tuple(c * x for x in a)


def dot(a: Vec, b: Vec) -> F:
    return sum((x * y for x, y in zip(a, b)), F(0))


def norm_sq(a: Vec) -> F:
    return dot(a, a)


def transpose(a: Mat) -> Mat:
    return tuple(tuple(a[j][i] for j in range(len(a)))
                 for i in range(len(a[0])))


def det3(a: Mat) -> F:
    assert len(a) == 3 and all(len(row) == 3 for row in a)
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def mmul(a: Mat, b: Mat) -> Mat:
    bt = transpose(b)
    return tuple(tuple(dot(row, col) for col in bt) for row in a)


def madd(a: Mat, b: Mat) -> Mat:
    return tuple(tuple(x + y for x, y in zip(ar, br))
                 for ar, br in zip(a, b))


def msub(a: Mat, b: Mat) -> Mat:
    return tuple(tuple(x - y for x, y in zip(ar, br))
                 for ar, br in zip(a, b))


def mscale(c: F, a: Mat) -> Mat:
    return tuple(tuple(c * x for x in row) for row in a)


def mvec(a: Mat, x: Vec) -> Vec:
    return tuple(dot(row, x) for row in a)


def row_mat(row: Vec, a: Mat) -> Vec:
    return tuple(dot(row, col) for col in transpose(a))


def outer(a: Vec, b: Vec) -> Mat:
    return tuple(tuple(x * y for y in b) for x in a)


def identity(n: int) -> Mat:
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def leray(k: Vec, value: Vec) -> Vec:
    return vsub(value, vscale(dot(k, value) / norm_sq(k), k))


def check_a2_resonance_and_leray() -> None:
    n = (F(1), F(1), F(1))
    k1 = (F(1), F(-1), F(0))
    kc = (F(1), F(0), F(-1))
    e1 = (F(2), F(-1), F(-1))
    a = (F(1), F(1), F(-2))

    assert vsub(e1, k1) == kc
    assert dot(a, k1) == 0
    assert dot(a, e1) == 3
    assert dot(n, e1) == dot(n, kc) == 0
    assert norm_sq(k1) == norm_sq(kc) == 2
    assert norm_sq(e1) == 6
    assert leray(kc, n) == n
    assert leray(vadd(e1, k1), n) == n

    # V.grad W has child sine coefficient -3/2.  Moving the linearized
    # convection term to the right of the evolution gives signed coupling
    # +3/2, while the reverse ordering is zero because N.k1=0.
    assert dot(a, e1) / 2 == F(3, 2)
    assert dot(n, k1) == 0

    # Two nonorthogonal rational material covector maps.  Linearity of
    # F^{-T} preserves the resonance exactly, not approximately.
    maps: tuple[Mat, ...] = (
        (
            (F(1), F(2), F(0)),
            (F(0), F(1), F(1)),
            (F(0), F(0), F(1)),
        ),
        (
            (F(2), F(1), F(0)),
            (F(1), F(1), F(1)),
            (F(0), F(0), F(1)),
        ),
    )
    for covector_map in maps:
        assert det3(covector_map) == 1
        assert vsub(mvec(covector_map, e1), mvec(covector_map, k1)) == mvec(
            covector_map, kc
        )

    # C161 translation labels obey the same relation for every source and
    # charge shift in this exact finite sample.
    r1 = k1
    r2 = (F(0), F(1), F(-1))
    for b in range(4):
        for c in range(4):
            source = vadd(vscale(F(9), n), vadd(vscale(F(b), r1), vscale(F(c), r2)))
            for charge in (-2, -1, 1, 2):
                gate = vscale(F(charge), n)
                daughter = vadd(source, gate)
                for covector_map in maps:
                    assert mvec(covector_map, daughter) == vadd(
                        mvec(covector_map, source), mvec(covector_map, gate)
                    )


def check_graph_and_melnikov_matrix_identity() -> None:
    # Rectangular active<-wake graph: active dimension 2, wake dimension 3.
    a: Mat = ((F(1), F(2)), (F(-1), F(3)))
    b: Mat = (
        (F(2), F(-1), F(4)),
        (F(0), F(3), F(1)),
    )
    c: Mat = (
        (F(1), F(0)),
        (F(-2), F(1)),
        (F(3), F(-1)),
    )
    d: Mat = (
        (F(0), F(2), F(-1)),
        (F(1), F(-1), F(0)),
        (F(2), F(3), F(1)),
    )
    k: Mat = (
        (F(1, 2), F(-1, 3), F(2, 5)),
        (F(-2, 7), F(3, 4), F(1, 6)),
    )

    # K' = A K + B - K D - K C K.
    kdot = msub(madd(mmul(a, k), b), madd(mmul(k, d), mmul(mmul(k, c), k)))

    ell: Vec = (F(2, 3), F(-5, 4))
    r: Vec = (F(3, 5), F(-2, 7), F(4, 9))
    elldot = vscale(F(-1), row_mat(ell, a))
    rdot = mvec(d, r)

    # Direct product-rule derivative of ell K r.
    direct = (
        dot(elldot, mvec(k, r))
        + dot(ell, mvec(kdot, r))
        + dot(ell, mvec(k, rdot))
    )
    melnikov = dot(ell, mvec(b, r)) - dot(
        ell, mvec(mmul(mmul(k, c), k), r)
    )
    assert direct == melnikov

    # If C=0, the graph identity is exactly first order.
    zero_c: Mat = tuple(tuple(F(0) for _ in range(2)) for _ in range(3))
    assert mmul(mmul(k, zero_c), k) == tuple(
        tuple(F(0) for _ in range(3)) for _ in range(2)
    )


def check_restricted_dark_subspace() -> None:
    # A=C=D=0 and rank-one B=(1,0).  The full wake-to-active block is
    # nonzero, but W_adm=ker B=span{(0,1)} is exactly dark.  This certifies
    # why a finite-codimension graph only requires a restricted Melnikov
    # bound plus membership of the physical wake.
    b: Mat = ((F(1), F(0)),)
    bright = (F(1), F(0))
    dark = (F(0), F(1))
    assert mvec(b, bright) == (F(1),)
    assert mvec(b, dark) == (F(0),)
    assert b != ((F(0), F(0)),)


def heat_integral(nu: float, t: float) -> float:
    if nu == 0.0:
        return t
    return -__import__("math").expm1(-6.0 * nu * t) / (6.0 * nu)


def terminal_graph(nu: float, t: float, horizon: float, b0: float = 1.5) -> float:
    if nu == 0.0:
        return -b0 * (horizon - t)
    return -b0 * (
        exp(-2.0 * nu * t) - exp(4.0 * nu * t - 6.0 * nu * horizon)
    ) / (6.0 * nu)


def check_heat_graph_formula() -> None:
    b0 = 1.5
    horizon = 0.7
    for nu in (0.0, 1.0e-8, 0.03, 0.2):
        k0 = terminal_graph(nu, 0.0, horizon, b0)
        expected = -b0 * heat_integral(nu, horizon)
        assert isclose(k0, expected, rel_tol=2e-9, abs_tol=2e-9)
        assert isclose(terminal_graph(nu, horizon, horizon, b0), 0.0,
                       rel_tol=0.0, abs_tol=2e-12)

        # Check K' = 4 nu K + B0 exp(-2 nu t) by a centered difference.
        # The closed form has a removable cancellation at nu=0; avoid
        # differentiating that cancellation numerically at 1e-8.
        if nu > 1.0e-6:
            for t in (0.1, 0.3, 0.6):
                step = 1.0e-6
                derivative = (
                    terminal_graph(nu, t + step, horizon, b0)
                    - terminal_graph(nu, t - step, horizon, b0)
                ) / (2.0 * step)
                rhs = 4.0 * nu * terminal_graph(nu, t, horizon, b0) + b0 * exp(
                    -2.0 * nu * t
                )
                assert isclose(derivative, rhs, rel_tol=2e-7, abs_tol=2e-8)

        born = b0 * exp(-2.0 * nu * horizon) * heat_integral(nu, horizon)
        assert born > 0.0

    # Inviscid limits of the graph tilt and Born response.
    for nu in (1.0e-3, 1.0e-5, 1.0e-7):
        assert isclose(
            terminal_graph(nu, 0.0, horizon, b0),
            -b0 * horizon,
            rel_tol=3.1 * nu,
            abs_tol=2e-9,
        )

    # Elementary lower bound 1-exp(-x) >= x/2 for x in [0,1].
    for j in range(1, 101):
        x = j / 100.0
        assert 1.0 - exp(-x) >= x / 2.0


def check_factorial_schedule_not_small_graph() -> None:
    # q=n^8, b=n^-2.  Factorial viscosity makes the fixed-time Melnikov
    # coefficient approach 3T/2, so its ratio to b diverges.
    physical_nu = 1.0
    horizon = 0.5
    previous_ratio = 0.0
    for n in range(4, 13):
        # C127 uses n=j+1 and mu_j=nu/(j!)^2=nu/((n-1)!)^2.
        mu = physical_nu / (factorial(n - 1) ** 2)
        b = n ** -2
        coefficient = 1.5 * heat_integral(mu, horizon)
        assert 6.0 * mu * horizon <= 1.0
        assert coefficient >= 0.75 * horizon
        ratio = coefficient / b
        assert ratio > previous_ratio
        previous_ratio = ratio
    assert previous_ratio > 100.0


def main() -> None:
    check_a2_resonance_and_leray()
    check_graph_and_melnikov_matrix_identity()
    check_restricted_dark_subspace()
    check_heat_graph_formula()
    check_factorial_schedule_not_small_graph()
    print("C175 Lagrangian resonant graph checks passed")


if __name__ == "__main__":
    main()
