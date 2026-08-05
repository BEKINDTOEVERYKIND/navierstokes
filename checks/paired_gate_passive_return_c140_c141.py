#!/usr/bin/env python3
"""Exact checks for C140--C141: 2D3C closure and cubic wake return.

The base helical arithmetic is imported from the C114--C117 checker.  All
coefficients remain in Q(sqrt(2)) + i Q(sqrt(2)); no floating point or
external symbolic package is used.
"""

from fractions import Fraction as F

from terminal_triad_hexagon_c114_c117 import (
    C2,
    I,
    INV_SQRT2,
    K1,
    K2,
    N,
    SQRT2,
    ZERO,
    convolution,
    dot,
    iadd,
    idot,
    ineg,
    mode,
    vadd,
    vconj,
    vscale,
)


ZV = (ZERO, ZERO, ZERO)
KC = iadd(K1, K2)
D12 = iadd(K1, ineg(K2))
E1 = iadd(K1, KC)
E2 = iadd(K2, KC)
ACTIVE = {K1, K2, KC, ineg(K1), ineg(K2), ineg(KC)}
NAMED_WAKE = {E1, E2, ineg(E1), ineg(E2)}


def clean(values):
    return {k: v for k, v in values.items() if v != ZV}


def add_dicts(*values):
    out = {}
    for value in values:
        for k, v in value.items():
            out[k] = vadd(out.get(k, ZV), v)
    return clean(out)


def scale_dict(value, scalar):
    return clean({k: vscale(scalar, v) for k, v in value.items()})


def subtract_dicts(left, right):
    return add_dicts(left, scale_dict(right, -1))


def restrict(value, support):
    return {k: v for k, v in value.items() if k in support and v != ZV}


def laplacian(value, viscosity):
    return {
        k: vscale(-viscosity * idot(k, k), v)
        for k, v in value.items()
        if v != ZV
    }


def explicit_gate():
    a1 = mode(K1, I, C2(1))
    a2 = mode(K2, -I, C2(1))
    return {
        K1: a1,
        K2: a2,
        ineg(K1): vconj(a1),
        ineg(K2): vconj(a2),
    }


def taylor_series(initial, viscosity, projected, last_order=3):
    """Coefficients c_n in u(t)=sum c_n t^n for NS or its P truncation."""
    series = [initial]
    for n in range(last_order):
        nonlinear = add_dicts(
            *(convolution(series[a], series[n - a]) for a in range(n + 1))
        )
        rhs = add_dicts(laplacian(series[n], viscosity), nonlinear)
        coefficient = scale_dict(rhs, F(1, n + 1))
        series.append(restrict(coefficient, ACTIVE) if projected else coefficient)
    return series


def scalar_times_n(scalar):
    return vscale(scalar, N)


def abs_sq(scalar):
    product = scalar.conj() * scalar
    assert product.im == ZERO.re
    return product.re


def check_c140_2d3c_and_viscous_return():
    initial = explicit_gate()

    # Orthogonal splitting u=v+phi*N on the plane k.N=0.  The explicit
    # planar component is monochromatic and has zero projected self-action.
    planar = {}
    vertical = {}
    for k, value in initial.items():
        assert idot(k, (1, 1, 1)) == 0
        phi = dot(N, value) / 3
        vertical[k] = scalar_times_n(phi)
        planar[k] = vadd(value, vscale(-1, vertical[k]))
        assert dot(N, planar[k]) == ZERO
    assert clean(convolution(planar, planar)) == {}
    assert clean(convolution(vertical, planar)) == {}
    assert clean(convolution(vertical, vertical)) == {}
    planar_advects_scalar = clean(convolution(planar, vertical))
    assert planar_advects_scalar
    assert all(vscale(1, value) == scalar_times_n(dot(N, value) / 3)
               for value in planar_advects_scalar.values())
    assert all(idot(k, k) == 2 for k in planar)

    expected_wake = {
        E1: scalar_times_n(C2(-9, -9)),
        E2: scalar_times_n(C2(-9, 9)),
        ineg(E1): scalar_times_n(C2(-9, 9)),
        ineg(E2): scalar_times_n(C2(-9, -9)),
    }
    expected_return = {
        KC: scalar_times_n(C2(0, 18 * SQRT2)),
        ineg(KC): scalar_times_n(C2(0, -18 * SQRT2)),
    }

    # Three exact viscosities check the coefficient recurrence, while the
    # proof uses only that Delta is diagonal and the active roots share one
    # shell.  Both leakage coefficients are in fact viscosity-independent.
    for viscosity in (F(0), F(1, 7), F(5, 3)):
        full = taylor_series(initial, viscosity, projected=False)
        retained = taylor_series(initial, viscosity, projected=True)

        assert restrict(full[0], set(full[0]) - ACTIVE) == {}
        assert restrict(full[1], set(full[1]) - ACTIVE) == {}
        assert restrict(full[2], set(full[2]) - ACTIVE) == expected_wake
        for order in range(3):
            assert restrict(subtract_dicts(full[order], retained[order]), ACTIVE) == {}
        assert restrict(subtract_dicts(full[3], retained[3]), ACTIVE) == expected_return

        # The cubic active difference is exactly one further interaction of
        # the quadratic wake with the original parents.
        return_formula = scale_dict(
            restrict(
                add_dicts(
                    convolution(initial, expected_wake),
                    convolution(expected_wake, initial),
                ),
                ACTIVE,
            ),
            F(1, 3),
        )
        assert return_formula == expected_return


def cancelled_family(a, c, y1=C2(1), y2=C2(1)):
    r1 = C2(a, c)
    r2 = C2(a, -c)
    x1 = y1 * r1
    x2 = y2 * r2
    a1 = mode(K1, x1, y1)
    a2 = mode(K2, x2, y2)
    initial = {
        K1: a1,
        K2: a2,
        ineg(K1): vconj(a1),
        ineg(K2): vconj(a2),
    }
    return initial, x1, x2


def nonlinear_derivatives(initial):
    """Return D1, D2 and the wake-fed portion of D3 for Euler."""
    d1 = clean(convolution(initial, initial))
    d2 = add_dicts(convolution(initial, d1), convolution(d1, initial))
    wake2 = restrict(d2, NAMED_WAKE)
    d3_wake = add_dicts(
        convolution(initial, wake2),
        convolution(wake2, initial),
    )
    return d1, d2, wake2, d3_wake


def check_c141_general_cubic_return():
    # With y1=y2=1, every component of D3_wake is a polynomial of degree at
    # most four in each of a,c.  Equality on a 5-by-5 rational grid is
    # therefore an exact interpolation certificate for the displayed
    # bivariate identity, not a floating-point sample.
    nodes = (-2, -1, 0, 1, 2)
    for a in nodes:
        for c in nodes:
            initial, _, _ = cancelled_family(a, c)
            d1, _, wake2, d3_wake = nonlinear_derivatives(initial)
            expected_child = scalar_times_n(C2(0, -6 * c * SQRT2))
            expected_return = scalar_times_n(
                C2(0, 54 * c * ((a - 1) ** 2 + c * c) * SQRT2)
            )
            assert d1.get(KC, ZV) == expected_child
            assert d3_wake.get(KC, ZV) == expected_return
            if c != 0:
                assert set(d1) == {KC, ineg(KC)}
                assert set(wake2) == NAMED_WAKE
                assert set(restrict(d3_wake, ACTIVE)) == {KC, ineg(KC)}
                assert expected_child != ZV
                assert expected_return != ZV
            else:
                assert d1 == {}
                assert wake2 == {}

    # Independently check the general complex-y formula and its two named
    # intermediate paths at exact rational-complex samples.
    samples = (
        (F(0), F(1), C2(1), C2(1)),
        (F(2), F(-3), C2(2, 1), C2(-1, 2)),
        (F(-1, 2), F(5, 3), C2(F(3, 2), F(-2, 3)), C2(F(4, 5), F(7, 4))),
    )
    for a, c, y1, y2 in samples:
        initial, x1, x2 = cancelled_family(a, c, y1, y2)
        d1, _, wake2, d3_wake = nonlinear_derivatives(initial)
        terminal = y1 * x2 - x1 * y2
        difference = x1 * y2.conj() - y1 * x2.conj()
        assert difference == ZERO
        assert terminal != ZERO

        d_1 = x1 - y1
        d_2 = x2 - y2
        assert d1[KC] == scalar_times_n(terminal * (3 * SQRT2))
        assert wake2[E1] == scalar_times_n(-9 * d_1 * terminal)
        assert wake2[E2] == scalar_times_n(9 * d_2 * terminal)
        expected = scalar_times_n(
            terminal * (-27 * INV_SQRT2 * (abs_sq(d_1) + abs_sq(d_2)))
        )
        assert d3_wake[KC] == expected
        assert expected != ZV
        assert set(d1) == {KC, ineg(KC)}
        assert set(wake2) == NAMED_WAKE
        assert set(restrict(d3_wake, ACTIVE)) == {KC, ineg(KC)}

    # Unparameterized finite-family audit.  Unlike cancelled_family(), this
    # loop starts from four independent Gaussian-rational amplitudes and
    # filters by the cancellation and terminal conditions.  It therefore
    # checks that the ratio parameterization has not hidden a degenerate
    # branch and ties the coordinate-free formula directly to convolution.
    gaussian_values = (
        ZERO,
        C2(1),
        C2(-1),
        I,
        -I,
        C2(1, 1),
    )
    witnesses = 0
    for x1 in gaussian_values:
        for y1 in gaussian_values:
            for x2 in gaussian_values:
                for y2 in gaussian_values:
                    difference = x1 * y2.conj() - y1 * x2.conj()
                    terminal = y1 * x2 - x1 * y2
                    if difference != ZERO or terminal == ZERO:
                        continue

                    # T != 0 together with D = 0 excludes y_i = 0.
                    assert y1 != ZERO and y2 != ZERO
                    a1 = mode(K1, x1, y1)
                    a2 = mode(K2, x2, y2)
                    initial = {
                        K1: a1,
                        K2: a2,
                        ineg(K1): vconj(a1),
                        ineg(K2): vconj(a2),
                    }
                    d1, _, wake2, d3_wake = nonlinear_derivatives(initial)
                    assert d1.get(D12, ZV) == ZV
                    assert d1.get(ineg(D12), ZV) == ZV
                    assert set(d1) == {KC, ineg(KC)}
                    assert set(wake2) == NAMED_WAKE
                    assert set(restrict(d3_wake, ACTIVE)) == {
                        KC,
                        ineg(KC),
                    }

                    delta1 = x1 - y1
                    delta2 = x2 - y2
                    expected = scalar_times_n(
                        terminal
                        * (
                            -27
                            * INV_SQRT2
                            * (abs_sq(delta1) + abs_sq(delta2))
                        )
                    )
                    assert expected != ZV
                    assert d3_wake[KC] == expected
                    witnesses += 1
    assert witnesses >= 40


def main():
    check_c140_2d3c_and_viscous_return()
    check_c141_general_cubic_return()
    print("PASS C140: exact 2D3C gate and viscosity-independent cubic return jet")
    print("PASS C141: every nondegenerate cancelled four-mode gate has cubic wake return")
    print("Bare fixed-plane gate only; localized conversion/focus BAFL remains open")


if __name__ == "__main__":
    main()
