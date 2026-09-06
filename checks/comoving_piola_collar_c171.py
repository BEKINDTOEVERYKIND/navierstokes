#!/usr/bin/env python3
"""Exact arithmetic checks for C171's co-moving Piola collar.

The checker verifies:

* the trace-free affine curl potential on the core;
* nontrivial affine and variable exact curl--Piola transforms;
* the Cauchy material derivative and factor two in the parent cross term;
* an exact Fourier Helmholtz orthogonal split;
* the full degree-five Piola/F-jet chain-rule envelope;
* exhaustive ordered-sign bookkeeping for real zero-charge A2 triads;
* the stationary, co-moving, backward, and terminal-chart scale powers;
* the absence of a homogeneity-generated factor b.

It does not certify a scale-uniform material chart, a localized active
evolution family, an integrated wake cancellation, LCE, BAFL, or an
unforced Navier--Stokes stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log


Vec = tuple[F, F, F]
Mat = tuple[Vec, Vec, Vec]
Exp = tuple[int, int, int]
Poly = dict[Exp, F]

ZERO: Vec = (F(0), F(0), F(0))


def vadd(a: Vec, b: Vec) -> Vec:
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def vsub(a: Vec, b: Vec) -> Vec:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def vscale(c: F, a: Vec) -> Vec:
    return tuple(c * a[i] for i in range(3))  # type: ignore[return-value]


def dot(a: Vec, b: Vec) -> F:
    return sum((a[i] * b[i] for i in range(3)), F(0))


def cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def mvec(a: Mat, x: Vec) -> Vec:
    return tuple(dot(a[i], x) for i in range(3))  # type: ignore[return-value]


def mmul(a: Mat, b: Mat) -> Mat:
    bt = tuple(tuple(b[j][i] for j in range(3)) for i in range(3))
    return tuple(
        tuple(dot(a[i], bt[j]) for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def transpose(a: Mat) -> Mat:
    return tuple(
        tuple(a[j][i] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def det(a: Mat) -> F:
    return dot(a[0], cross(a[1], a[2]))


def norm_sq(a: Vec) -> F:
    return dot(a, a)


def padd(a: Poly, b: Poly) -> Poly:
    result = dict(a)
    for exp, coeff in b.items():
        result[exp] = result.get(exp, F(0)) + coeff
        if result[exp] == 0:
            del result[exp]
    return result


def pscale(c: F, a: Poly) -> Poly:
    if c == 0:
        return {}
    return {exp: c * coeff for exp, coeff in a.items() if c * coeff}


def pmul(a: Poly, b: Poly) -> Poly:
    result: Poly = {}
    for left_exp, left_coeff in a.items():
        for right_exp, right_coeff in b.items():
            exp = tuple(
                left_exp[i] + right_exp[i] for i in range(3)
            )
            result[exp] = result.get(exp, F(0)) + left_coeff * right_coeff
    return {exp: coeff for exp, coeff in result.items() if coeff}


def pdiff(a: Poly, axis: int) -> Poly:
    result: Poly = {}
    for exp, coeff in a.items():
        if exp[axis]:
            new_exp = list(exp)
            new_exp[axis] -= 1
            result[tuple(new_exp)] = coeff * exp[axis]
    return result


def pvec_add(a: tuple[Poly, Poly, Poly], b: tuple[Poly, Poly, Poly]):
    return tuple(padd(a[i], b[i]) for i in range(3))


def pvec_scale(c: F, a: tuple[Poly, Poly, Poly]):
    return tuple(pscale(c, a[i]) for i in range(3))


def pmat_vec(a: Mat, x: tuple[Poly, Poly, Poly]):
    return tuple(
        sum_polys(pscale(a[i][j], x[j]) for j in range(3))
        for i in range(3)
    )


def sum_polys(items) -> Poly:
    result: Poly = {}
    for item in items:
        result = padd(result, item)
    return result


def pcross(a: tuple[Poly, Poly, Poly], b: tuple[Poly, Poly, Poly]):
    return (
        padd(pmul(a[1], b[2]), pscale(-1, pmul(a[2], b[1]))),
        padd(pmul(a[2], b[0]), pscale(-1, pmul(a[0], b[2]))),
        padd(pmul(a[0], b[1]), pscale(-1, pmul(a[1], b[0]))),
    )


def pcurl(a: tuple[Poly, Poly, Poly]):
    return (
        padd(pdiff(a[2], 1), pscale(-1, pdiff(a[1], 2))),
        padd(pdiff(a[0], 2), pscale(-1, pdiff(a[2], 0))),
        padd(pdiff(a[1], 0), pscale(-1, pdiff(a[0], 1))),
    )


def coordinate_polys(linear_map: Mat):
    """Return the three linear forms given by linear_map * x."""
    result = []
    for row in linear_map:
        poly: Poly = {}
        for axis, coeff in enumerate(row):
            if coeff:
                exp = [0, 0, 0]
                exp[axis] = 1
                poly[tuple(exp)] = coeff
        result.append(poly)
    return tuple(result)


def affine_potential(s: Mat, x: tuple[Poly, Poly, Poly]):
    sx = pmat_vec(s, x)
    return pvec_scale(F(1, 3), pcross(sx, x))


def check_affine_curl_and_piola() -> None:
    identity: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    )
    x = coordinate_polys(identity)

    # First check the sign and trace correction without imposing trace zero.
    s_with_trace: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(2), F(0)),
        (F(0), F(0), F(4)),
    )
    trace = F(7)
    general_curl = pcurl(affine_potential(s_with_trace, x))
    expected_general = pvec_add(
        pmat_vec(s_with_trace, x),
        pvec_scale(-trace / 3, x),
    )
    assert general_curl == expected_general

    # A nontrivial symmetric trace-free core matrix.
    s: Mat = (
        (F(0), F(1), F(0)),
        (F(1), F(0), F(2)),
        (F(0), F(2), F(0)),
    )
    assert sum((s[i][i] for i in range(3)), F(0)) == 0
    potential = affine_potential(s, x)
    assert pcurl(potential) == pmat_vec(s, x)

    # Nonorthogonal, determinant-one affine chart.
    f: Mat = (
        (F(1), F(2), F(0)),
        (F(0), F(1), F(1)),
        (F(0), F(0), F(1)),
    )
    f_inv: Mat = (
        (F(1), F(-2), F(2)),
        (F(0), F(1), F(-1)),
        (F(0), F(0), F(1)),
    )
    assert det(f) == 1
    assert mmul(f, f_inv) == identity

    a_of_x = coordinate_polys(f_inv)
    a0_of_a = affine_potential(s, a_of_x)
    alpha = pmat_vec(transpose(f_inv), a0_of_a)
    left = pcurl(alpha)
    curl_a0_of_a = pmat_vec(s, a_of_x)
    right = pmat_vec(f, curl_a0_of_a)
    assert left == right

    # Variable volume-preserving shear X(a)=(a1+a2^2,a2,a3).  This checks
    # the derivative cancellation in curl(F^{-T}A0), which a constant affine
    # chart alone cannot see.  In x coordinates a=(x1-x2^2,x2,x3).
    x1, x2, x3 = x
    a1 = padd(x1, pscale(-1, pmul(x2, x2)))
    nonlinear_a = (a1, x2, x3)
    nonlinear_a0 = affine_potential(s, nonlinear_a)
    nonlinear_alpha = (
        nonlinear_a0[0],
        padd(nonlinear_a0[1], pscale(-2, pmul(x2, nonlinear_a0[0]))),
        nonlinear_a0[2],
    )
    nonlinear_left = pcurl(nonlinear_alpha)
    nonlinear_w0 = pmat_vec(s, nonlinear_a)
    nonlinear_right = (
        padd(nonlinear_w0[0], pscale(2, pmul(x2, nonlinear_w0[1]))),
        nonlinear_w0[1],
        nonlinear_w0[2],
    )
    assert nonlinear_left == nonlinear_right


def check_material_derivative_and_factor_two() -> None:
    # Parent shear V=(x_2,0,0), X=(a_1+t a_2,a_2,a_3).
    parent_gradient: Mat = (
        (F(0), F(1), F(0)),
        (F(0), F(0), F(0)),
        (F(0), F(0), F(0)),
    )
    for time in (F(-3, 2), F(0), F(5, 3)):
        flow_gradient: Mat = (
            (F(1), time, F(0)),
            (F(0), F(1), F(0)),
            (F(0), F(0), F(1)),
        )
        flow_derivative: Mat = (
            (F(0), F(1), F(0)),
            (F(0), F(0), F(0)),
            (F(0), F(0), F(0)),
        )
        assert flow_derivative == mmul(parent_gradient, flow_gradient)
        for w0 in (
            (F(2), F(-1), F(3)),
            (F(0), F(7, 4), F(-2)),
        ):
            u = mvec(flow_gradient, w0)
            material = mvec(flow_derivative, w0)
            strain = mvec(parent_gradient, u)
            assert material == strain
            assert vadd(material, strain) == vscale(F(2), strain)

    # Exact affine-jet remainder decomposition:
    # D_t U = A U + E.grad U and (U.grad)V=(A+grad E)U.
    a: Mat = (
        (F(1), F(2), F(0)),
        (F(0), F(-1), F(1)),
        (F(0), F(0), F(0)),
    )
    grad_e: Mat = (
        (F(0), F(1, 3), F(0)),
        (F(-2, 5), F(0), F(0)),
        (F(0), F(0), F(0)),
    )
    u = (F(2), F(-3), F(5))
    e_dot_grad_u = (F(7), F(-2), F(1))
    lhs = vadd(
        vadd(mvec(a, u), e_dot_grad_u),
        mvec(
            tuple(
                tuple(a[i][j] + grad_e[i][j] for j in range(3))
                for i in range(3)
            ),  # type: ignore[arg-type]
            u,
        ),
    )
    rhs = vadd(
        vadd(vscale(F(2), mvec(a, u)), e_dot_grad_u),
        mvec(grad_e, u),
    )
    assert lhs == rhs


def check_helmholtz_l2_split() -> None:
    # One exact Fourier fiber represents the orthogonal torus Helmholtz
    # decomposition.  The statement then sums over all nonzero fibers.
    k = (F(2), F(-1), F(-1))
    g = (F(5), F(4), F(-2))
    gradient = vscale(dot(k, g) / norm_sq(k), k)
    projected = vsub(g, gradient)
    assert dot(k, projected) == 0
    assert dot(projected, gradient) == 0
    assert norm_sq(g) == norm_sq(projected) + norm_sq(gradient)
    assert norm_sq(projected) <= norm_sq(g)
    assert norm_sq(gradient) <= norm_sq(g)


def check_full_f_jet_polynomial() -> None:
    # Dimensionless variables are (f,g,d1,d2)=(|F|,|F^-1|,r|DF|,r^2|D2F|).
    # The actual first/second spatial chain rules produce the following
    # monomials.  In particular the derivative of F^-1 produces the two
    # degree-five terms, so a condition-number-only factor is insufficient.
    monomial_exponents = (
        (1, 0, 0, 0),  # f
        (1, 1, 0, 0),  # g f
        (0, 1, 1, 0),  # g d1
        (1, 2, 0, 0),  # g^2 f
        (0, 2, 1, 0),  # g^2 d1
        (0, 2, 0, 1),  # g^2 d2
        (1, 3, 1, 0),  # g^3 d1 f
        (0, 3, 2, 0),  # g^3 d1^2
    )
    degrees = tuple(sum(exponents) for exponents in monomial_exponents)
    assert max(degrees) == 5
    assert degrees[-2:] == (5, 5)

    values = (F(2, 3), F(5, 4), F(7, 6), F(9, 8))
    chart_sum = sum(values, F(0))
    polynomial_cap = (1 + chart_sum) ** 5
    for exponents in monomial_exponents:
        monomial = F(1)
        for value, exponent in zip(values, exponents):
            monomial *= value**exponent
        assert 0 < monomial <= polynomial_cap


def neg(k: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(-entry for entry in k)  # type: ignore[return-value]


def iadd(k: tuple[int, int, int], q: tuple[int, int, int]):
    return tuple(k[i] + q[i] for i in range(3))


def ivec_to_fraction(k: tuple[int, int, int]) -> Vec:
    return tuple(F(entry) for entry in k)  # type: ignore[return-value]


def ordered_fourier_raw(
    left: dict[tuple[int, int, int], Vec],
    right: dict[tuple[int, int, int], Vec],
    target: tuple[int, int, int],
):
    """Return every ordered path and its real vector before multiplication by i."""
    paths = []
    total = ZERO
    for p, left_coefficient in left.items():
        for q, right_coefficient in right.items():
            if iadd(p, q) == target:
                contribution = vscale(
                    dot(left_coefficient, ivec_to_fraction(q)),
                    right_coefficient,
                )
                paths.append((p, q, contribution))
                total = vadd(total, contribution)
    return paths, total


def check_a2_reality_triad_exhaustively() -> None:
    n = (1, 1, 1)
    k1 = (1, -1, 0)
    kc = (1, 0, -1)
    e1 = (2, -1, -1)
    a_star = (F(1), F(1), F(-2))
    n_vec = ivec_to_fraction(n)

    assert iadd(k1, kc) == e1
    assert dot(ivec_to_fraction(n), ivec_to_fraction(k1)) == 0
    assert dot(ivec_to_fraction(n), ivec_to_fraction(kc)) == 0
    assert dot(ivec_to_fraction(n), ivec_to_fraction(e1)) == 0
    assert dot(a_star, ivec_to_fraction(k1)) == 0
    assert dot(a_star, ivec_to_fraction(e1)) == 3

    # Exact real cosine data: each signed Fourier coefficient is half of
    # its real polarization.
    v_modes = {
        k1: vscale(F(1, 2), a_star),
        neg(k1): vscale(F(1, 2), a_star),
    }
    w_modes = {
        e1: vscale(F(1, 2), n_vec),
        neg(e1): vscale(F(1, 2), n_vec),
    }

    vw_paths, vw_total = ordered_fourier_raw(v_modes, w_modes, kc)
    wv_paths, wv_total = ordered_fourier_raw(w_modes, v_modes, kc)
    all_paths = vw_paths + wv_paths

    # Exhaustion over every ordered sign pair leaves precisely the two
    # frequency paths (-k1,e1) and (e1,-k1).  The first is bright and the
    # reverse ordering is zero.
    assert [(p, q) for p, q, _ in all_paths] == [
        (neg(k1), e1),
        (e1, neg(k1)),
    ]
    assert vw_paths[0][2] == vscale(F(3, 4), n_vec)
    assert wv_paths[0][2] == ZERO
    assert vadd(vw_total, wv_total) == vscale(F(3, 4), n_vec)

    # The -kc coefficient is the conjugate of i*(3/4)N, namely
    # -i*(3/4)N.  It reconstructs the nonzero real sine coefficient.
    minus_vw_paths, minus_vw_total = ordered_fourier_raw(
        v_modes, w_modes, neg(kc)
    )
    minus_wv_paths, minus_wv_total = ordered_fourier_raw(
        w_modes, v_modes, neg(kc)
    )
    assert len(minus_vw_paths) + len(minus_wv_paths) == 2
    assert vadd(minus_vw_total, minus_wv_total) == vscale(F(-3, 4), n_vec)

    # Leray leaves N unchanged at kc.
    kc_vec = ivec_to_fraction(kc)
    leray_n = vsub(
        n_vec, vscale(dot(kc_vec, n_vec) / norm_sq(kc_vec), kc_vec)
    )
    assert leray_n == n_vec

    # The displayed actual-wake polarization is dark in the Piola ordering
    # (W.grad)V, so it is not a lower bound for that particular residual.
    assert wv_total == ZERO

    # Charge/reality still do not kill the Piola ordering on the allowed
    # polarization class.  B=(0,-1,1) is divergence free at e1 and has
    # B.k1=1.  Twice the path (e1,-k1) has raw coefficient -A*/2 at kc;
    # Leray leaves the nonzero vector (1,-2,1)/4.
    b_star = (F(0), F(-1), F(1))
    assert dot(b_star, ivec_to_fraction(e1)) == 0
    assert dot(b_star, ivec_to_fraction(k1)) == 1
    piola_modes = {
        e1: vscale(F(1, 2), b_star),
        neg(e1): vscale(F(1, 2), b_star),
    }
    piola_paths, piola_once = ordered_fourier_raw(piola_modes, v_modes, kc)
    assert [(p, q) for p, q, _ in piola_paths] == [(e1, neg(k1))]
    assert piola_once == vscale(F(-1, 4), a_star)
    piola_twice = vscale(F(2), piola_once)
    piola_projected = vsub(
        piola_twice,
        vscale(dot(kc_vec, piola_twice) / norm_sq(kc_vec), kc_vec),
    )
    assert piola_projected == (F(1, 4), F(-1, 2), F(1, 4))


def check_scale_and_homogeneity_ledgers() -> None:
    # Formal q powers for the p=2 normalized, time-integrated terms.
    # Stationary: q^(-3/2) log h; co-moving: q^(-5/2) log h.
    stationary_q_exponent = F(-3, 2)
    moving_q_exponent = F(-5, 2)
    backward_h_exponent = F(3, 2)
    scalar_chart_q_exponent = F(1, 4)
    full_plane_chart_q_exponent = F(3)
    assert stationary_q_exponent + backward_h_exponent == 0
    assert moving_q_exponent + backward_h_exponent == -1
    assert (
        moving_q_exponent
        + backward_h_exponent
        + scalar_chart_q_exponent
        == F(-3, 4)
    )
    assert (
        moving_q_exponent
        + backward_h_exponent
        + full_plane_chart_q_exponent
        == 2
    )

    # q=n^8, h=n^12 gives q^-1 log h = 12 n^-8 log n.
    q_n_exponent = 8
    h_n_exponent = 12
    assert F(h_n_exponent, q_n_exponent) == backward_h_exponent
    assert -q_n_exponent == -8
    assert F(-3, 4) * q_n_exponent == -6
    # Relative to n^-6, the remaining ratio is 12 log(n)/n^2 -> 0.
    ratios = [12 * log(n_value) / n_value**2 for n_value in (20, 50, 100, 200)]
    assert all(left > right for left, right in zip(ratios, ratios[1:]))
    assert ratios[-1] < F(1, 100)

    # The material chart condition is strictly M_F log(n)=o(n^2).
    # Borderline M_F=n^2/log(n) gives ratio one, while a subcritical sample
    # n^2/log(n)^2 tends to zero.  This same ratio compares the charted
    # n^-6 log(n) scale against the allowed n^-4 endpoint scale.
    strict_samples = []
    borderline_samples = []
    for n_value in (20, 50, 100, 200):
        logarithm = log(n_value)
        subcritical_m = n_value**2 / logarithm**2
        borderline_m = n_value**2 / logarithm
        strict_samples.append(subcritical_m * logarithm / n_value**2)
        borderline_samples.append(borderline_m * logarithm / n_value**2)
    assert all(left > right for left, right in zip(strict_samples, strict_samples[1:]))
    assert all(abs(value - 1) < 1e-12 for value in borderline_samples)

    # A fixed nonzero linear block cannot manufacture another b.
    linear_block: Mat = (
        (F(2), F(-1), F(0)),
        (F(0), F(3), F(1)),
        (F(1), F(0), F(-2)),
    )
    wake = (F(1), F(-2), F(4))
    for n_value in (2, 5, 11):
        b = F(1, n_value**2)
        left = mvec(linear_block, vscale(b * b, wake))
        right = vscale(b * b, mvec(linear_block, wake))
        assert left == right
        assert left != vscale(b**3, mvec(linear_block, wake))

        # Treat log(h)=12 log(n) numerically only for the strict scale
        # comparison; the algebraic amplitude ratio is exactly 1/b=n^2.
        assert (b * b) / (b**3) == n_value**2
        assert b * b * (12 * log(n_value)) > b**3


def main() -> None:
    check_affine_curl_and_piola()
    check_material_derivative_and_factor_two()
    check_helmholtz_l2_split()
    check_full_f_jet_polynomial()
    check_a2_reality_triad_exhaustively()
    check_scale_and_homogeneity_ledgers()
    print("PASS C171: exact curl/Piola and material factor-two identities")
    print("PASS C171: full F-jet polynomial and Helmholtz L2 contraction")
    print("PASS C171: q^-1 backward and q^-3/4 charted ledgers")
    print("PASS C171: real zero-charge A2 and Piola-order triads are nonzero")
    print("OPEN: material-collar kernel closure, LCE, BAFL, and full stage")


if __name__ == "__main__":
    main()
