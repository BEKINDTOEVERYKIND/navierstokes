#!/usr/bin/env python3
"""Exact checks for the leading real-child Floquet sideband ladder.

The script verifies the Bloch symbol, the high-sum Wronskian identity,
and the incompatibility between common polarization and an active real
child.  It does not prove a nonlinear WKB error estimate.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def convolution_high_sum(
    a: dict[int, Fraction], b: dict[int, Fraction]
) -> dict[int, Fraction]:
    """Return C_j=sum_{n+m=j-1}(m-n)a_n b_m."""
    result: dict[int, Fraction] = {}
    for n, an in a.items():
        for m, bm in b.items():
            j = n + m + 1
            result[j] = result.get(j, Q(0)) + (m - n) * an * bm
    return {j: value for j, value in result.items() if value}


def polynomial_wronskian(
    a: dict[int, Fraction], b: dict[int, Fraction]
) -> dict[int, Fraction]:
    """Coefficients of A(z) z B'(z)-z A'(z) B(z)."""
    result: dict[int, Fraction] = {}
    for n, an in a.items():
        for m, bm in b.items():
            power = n + m
            result[power] = result.get(power, Q(0)) + (m - n) * an * bm
    return {power: value for power, value in result.items() if value}


def check_high_sum_wronskian() -> None:
    samples = [
        (
            {-2: Q(3, 2), 0: Q(-2), 3: Q(5, 3)},
            {-1: Q(7, 4), 1: Q(2, 5), 4: Q(-3)},
        ),
        (
            {-3: Q(-1), -1: Q(4), 2: Q(9, 7)},
            {-3: Q(5, 2), 0: Q(-8, 3), 2: Q(1, 6)},
        ),
    ]
    for a, b in samples:
        direct = convolution_high_sum(a, b)
        wronskian = polynomial_wronskian(a, b)
        shifted = {power + 1: value for power, value in wronskian.items()}
        assert direct == shifted
    print("positive-positive high sums equal the polynomial Wronskian")


def check_common_polarization_kills_leading_high_sum() -> None:
    a = {-4: Q(2), -1: Q(-3, 2), 0: Q(7), 5: Q(1, 3)}
    rho = Q(-11, 6)
    b = {n: rho * value for n, value in a.items()}
    assert convolution_high_sum(a, b) == {}
    print("common polarization kills every leading e_h high sum")


def check_bloch_symbol() -> None:
    # Store the 2x2 symbol as exact Laurent polynomials in z.  An entry is
    # represented by exponent -> coefficient.  The lower-left entry is
    # gH + QZ(z+z^{-1}).
    g = Q(7, 3)
    h = Q(-5, 4)
    kappa = Q(2, 5)
    qz = Q(11, 6)

    upper_right = g / h
    lower_left = {0: g * h, 1: qz, -1: qz}

    # The determinant of lambda I-M is
    # (lambda+kappa)^2-(g/H)(gH+QZ(z+z^{-1})).
    product = {
        exponent: upper_right * coefficient
        for exponent, coefficient in lower_left.items()
    }
    assert product[0] == g * g
    assert product[1] == g * qz / h
    assert product[-1] == g * qz / h

    # At z=1 and z=-1 this is g(g +/- 2QZ/H).
    at_zero = sum(product.values())
    at_pi = sum((Q(-1) ** exponent) * value
                for exponent, value in product.items())
    assert at_zero == g * (g + 2 * qz / h)
    assert at_pi == g * (g - 2 * qz / h)
    assert h < 0 and at_pi > g * g
    assert kappa < g
    print("Bloch symbol and enhanced theta=pi sector: exact")


def check_ratio_incompatibility() -> None:
    # If B=rho A, its ratio equation at theta=0 and theta=pi differs by
    # 4QZ.  No common rho can satisfy both when the child is active.
    g = Q(13, 7)
    h = Q(-9, 5)
    qz = Q(4, 3)
    rhs_zero = g * h + 2 * qz
    rhs_pi = g * h - 2 * qz
    assert rhs_zero - rhs_pi == 4 * qz
    assert qz != 0
    print("common-polarization ratio equation varies across Bloch angle")


def check_finite_support_leaks() -> None:
    # A finite vector with extreme radial coefficient a_N leaks at once to
    # b_{N+1}; if a_N=0 but b_N!=0, the parent creates a_N first.
    g = Q(5, 2)
    h = Q(-7, 3)
    qz = Q(9, 4)
    a_n = Q(6, 5)
    b_n = Q(-2, 7)
    outward_derivative = qz * a_n
    parent_derivative = (g / h) * b_n
    assert outward_derivative != 0
    assert parent_derivative != 0
    print("every nonzero finite extreme leaks to the next sideband")


def check_palasek_volume_obstruction() -> None:
    # A Bloch packet occupying fraction eta of its slow container needs
    # eta >= Va/Vp=N^{-w}.  Forced 2Lambda subordination needs
    # eta << sqrt(Vc/Va)=N^{-r}.  Compatibility would require w>r,
    # equivalently beta<2/b, contradicting beta>2b.
    b = Q(21, 20)
    beta = Q(12, 5)
    alpha = 1 + Q(3, 2) / b
    assert 1 < b
    assert 2 * b < beta < alpha

    w = (b - 1) * (2 * alpha - beta - 2)
    r = beta * (b - 1) / 2
    assert w == (b - 1) * (Q(3) / b - beta)
    assert w < r
    assert beta > 2 * b > 2 / b

    # At eta=Va/Vp, high-band dissipation / parent work grows like
    # h^2/f=N^{2(r-w)}.
    work_exponent = 2 * (r - w)
    assert work_exponent > 0
    print(
        "Palasek volume obstruction: no eta window; "
        f"h^2/f grows like N^{work_exponent}"
    )


def check_palasek_diffusion_threshold() -> None:
    # epsilon sqrt(G)=N^{b-beta/(2b)} decays only if beta>2b^2.
    b = Q(21, 20)
    beta = Q(12, 5)
    diffusion_decay = beta / (2 * b) - b
    assert beta > 2 * b * b
    assert diffusion_decay > 0
    print(
        "Bloch heat-spread threshold beta>2b^2: exact; "
        f"decay exponent {diffusion_decay}"
    )


def main() -> None:
    check_high_sum_wronskian()
    check_common_polarization_kills_leading_high_sum()
    check_bloch_symbol()
    check_ratio_incompatibility()
    check_finite_support_leaks()
    check_palasek_volume_obstruction()
    check_palasek_diffusion_threshold()
    print("all Floquet sideband-ladder checks passed")


if __name__ == "__main__":
    main()
