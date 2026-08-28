#!/usr/bin/env python3
"""Dependency-free stable-line certificate for the C159/C192 orbit.

The unstable amplitude columns are deliberately not interval-integrated.
Exact half-period reversibility reduces the stable line to one positive
Riccati ratio on half a period.  The C159 tube is then reevaluated on the
same 2048 cells as C192, and scalar comparison plus a directed Taylor
remainder encloses that ratio.

All Decimal endpoint operations in the proof path go through the directed
contexts and interval primitives of zero_drift_cooperative_cone_c159.
In particular, this file never applies Decimal unary minus or abs.
"""

from decimal import Decimal as D
from fractions import Fraction as F

import strong_zero_drift_gain_c192 as gain
import zero_drift_cooperative_cone_c159 as base


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), F(0))


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def exact_symmetry_checks():
    """Check the exact geometric and two-by-two symmetry algebra."""

    delta = F(4, 5)
    r1 = (F(1), F(-1), F(0))
    r2 = (F(0), F(1), F(-1))
    normal = (F(1), F(1), F(1))

    # P reverses the physical coordinates.  It fixes N, has determinant
    # -1, and sends (r1,r2) to (-r2,-r1), inducing the phase reversor
    # (a,b)->(-b,-a), which fixes C152's initial point (A,-A).
    # Consequently P(N x v)=-N x Pv.
    P = (
        (F(0), F(0), F(1)),
        (F(0), F(1), F(0)),
        (F(1), F(0), F(0)),
    )
    assert matvec(P, r1) == tuple(-entry for entry in r2)
    assert matvec(P, r2) == tuple(-entry for entry in r1)
    assert matvec(P, normal) == normal
    columns = tuple(zip(*P))
    assert tuple(tuple(dot(row, column) for column in columns) for row in P) == (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    )
    for vector in ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))):
        left = matvec(P, cross(normal, vector))
        right = tuple(-entry for entry in cross(normal, matvec(P, vector)))
        assert left == right

    # Phase vector field in sine variables.  Checking the three basis
    # inputs proves the linear identities: central inversion is a symmetry,
    # while (a,b)->(-b,-a) is a reversing symmetry.
    def phase_rhs(sa, sb, sab):
        return (3 * (sb + delta * sab), -3 * (sa + delta * sab))

    for values in ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))):
        sa, sb, sab = values
        rhs = phase_rhs(sa, sb, sab)
        central = phase_rhs(-sa, -sb, -sab)
        reversed_phase = phase_rhs(-sb, -sa, -sab)
        assert central == (-rhs[0], -rhs[1])
        assert reversed_phase == (rhs[1], rhs[0])

    # In the C159 formula (2.3), root reflection sends
    # (l.p,p.S.p,p.S.t) to (-l.p,-p.S.p,+p.S.t).  Therefore it flips the
    # two diagonal B entries and fixes both off-diagonal entries.  Test the
    # complete rational formula, including both terms in B12 and B22.
    def abstract_B(lp, psp, pst, m, d_value, q_value, root_c):
        return (
            m * lp / d_value,
            (2 * m * pst + root_c * (d_value - m * m)) / d_value,
            m * m * root_c / (q_value * d_value),
            (2 * psp + m * lp) / d_value,
        )

    for lp, psp, pst in (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
        (F(2), F(-3), F(5)),
    ):
        original = abstract_B(lp, psp, pst, F(7, 3), F(11, 2), F(17, 2), F(13, 5))
        reflected = abstract_B(-lp, -psp, pst, F(7, 3), F(11, 2), F(17, 2), F(13, 5))
        assert reflected == (-original[0], original[1], original[2], -original[3])

    # The reversor gives J M J=M^-1, hence the full M has equal diagonals
    # and determinant one.  Central inversion also gives B(s+1/2)=B(s),
    # so M=A^2 for the half map A.  Since every entry of A is positive,
    # equality of the two diagonals of A^2 forces equality of those of A;
    # det(A)>0 and det(A)^2=1 then give A=[[a,b],[c,a]], a^2-bc=1.
    # The following exact sentinel checks that normal form's eigenlines.
    a, b, c = F(5), F(8), F(3)
    assert a * a - b * c == 1
    # x^2=c/b; avoid irrational arithmetic by checking the eigenline
    # identity after multiplying the second component by x symbolically.
    # A(1,+/- x)=(a+/-sqrt(bc))(1,+/- x) is equivalent to b x^2=c.
    assert b * F(3, 8) == c


def riccati_value(a_value, b_value, c_value, ratio):
    """A+B r-C r^2 with all endpoints outward."""

    return base.I(a_value) + base.I(b_value) * ratio - base.I(c_value) * ratio.square()


def riccati_derivative(b_value, c_value, ratio):
    return base.I(b_value) - base.I(2) * base.I(c_value) * ratio


def point_taylor_step(point, a_value, b_value, c_value, dt_lo, dt_hi):
    """Enclose one constant-Riccati flow from an exact Decimal point.

    A fixed radius 1/100 is first proved self-consistent by the integral
    equation.  Taylor's identity

      r(h)=r(0)+h f(r(0))+int_0^h (h-t) f'(r(t))f(r(t)) dt

    then supplies the directed second-order remainder.
    """

    initial = base.I(point)
    tube = initial.widen(D("0.01"))
    candidate = initial + base.I(0, dt_hi) * riccati_value(
        a_value, b_value, c_value, tube
    )
    assert tube.lo <= candidate.lo and candidate.hi <= tube.hi

    values = riccati_value(a_value, b_value, c_value, tube)
    derivatives = riccati_derivative(b_value, c_value, tube)
    remainder = base.UP.divide(
        base.UP.multiply(
            base.UP.multiply(dt_hi, dt_hi),
            base.UP.multiply(values.abs_upper(), derivatives.abs_upper()),
        ),
        D(2),
    )
    euler = initial + base.I(dt_lo, dt_hi) * riccati_value(
        a_value, b_value, c_value, initial
    )
    return base.I(
        base.DOWN.subtract(euler.lo, remainder),
        base.UP.add(euler.hi, remainder),
    )


def interval_taylor_step(interval, a_value, b_value, c_value, dt_lo, dt_hi):
    """Use order preservation of a scalar autonomous flow at both ends."""

    lower = point_taylor_step(interval.lo, a_value, b_value, c_value, dt_lo, dt_hi)
    upper = point_taylor_step(interval.hi, a_value, b_value, c_value, dt_lo, dt_hi)
    return base.I(lower.lo, upper.hi)


def certify_stable_line():
    base.exact_structural_checks()
    period, beta, sqrt2 = base.parameter_intervals()
    records, *_ = base.generate_reference(period, beta, sqrt2)
    # Reprove the tube; reference generation alone is not a premise.
    base.certify_path(records, period, beta, sqrt2)

    panels = len(records)
    subdivisions = 64
    total_cells = panels * subdivisions
    half_cells = total_cells // 2
    panel_step = D(1) / panels
    dt_lo = base.DOWN.divide(D(1), D(total_cells))
    dt_hi = base.UP.divide(D(1), D(total_cells))
    phase_radius = D("2e-6")
    gamma_radius = D("8e-4")
    sqrt3 = base.sqrt_fraction_bound(base.F(3))

    # r_lower and r_upper enclose the two scalar comparison solutions.
    # For r>=0 and an interval matrix [l_ij,u_ij],
    #
    # l21+(l22-u11)r-u12 r^2 <= r' <=
    # u21+(u22-l11)r-l12 r^2.
    lower_ratio = base.I(0)
    upper_ratio = base.I(0)
    half_column = (D(1), D(0))
    gain_vector = (D(1), base.DOWN.divide(D(3), D(20)))
    maximum_alpha = D(0)
    count = 0

    for coefficients in records:
        polynomials = [base.interval_polynomial(row) for row in coefficients]
        for subcell in range(subdivisions):
            left = panel_step * D(subcell) / subdivisions
            right = panel_step * D(subcell + 1) / subdivisions
            argument = base.I(left, right)
            phase = [
                base.evaluate_polynomial(polynomials[index], argument).widen(phase_radius)
                for index in range(4)
            ]
            gamma = base.evaluate_polynomial(polynomials[5], argument).widen(gamma_radius)
            entries = gain.raw_coefficients(phase, gamma, period, beta, sqrt2, sqrt3)
            b11, b12, b21, b22 = entries

            # C192's full-period positive vector is recomputed here so the
            # stable multiplier conclusion does not call another main().
            g1, g2, alpha = gain.lower_step(
                gain_vector[0], gain_vector[1], entries, dt_lo, dt_hi
            )
            gain_vector = (g1, g2)
            maximum_alpha = max(maximum_alpha, alpha)

            if count < half_cells:
                lower_ratio = interval_taylor_step(
                    lower_ratio,
                    b21.lo,
                    base.DOWN.subtract(b22.lo, b11.hi),
                    b12.hi,
                    dt_lo,
                    dt_hi,
                )
                upper_ratio = interval_taylor_step(
                    upper_ratio,
                    b21.hi,
                    base.UP.subtract(b22.hi, b11.lo),
                    b12.lo,
                    dt_lo,
                    dt_hi,
                )
                assert lower_ratio.lo >= 0
                assert upper_ratio.lo >= 0

                # First column of the half map: its first component is a.
                h1, h2, _ = gain.lower_step(
                    half_column[0], half_column[1], entries, dt_lo, dt_hi
                )
                half_column = (h1, h2)
            count += 1

    assert count == total_cells
    assert maximum_alpha < 5

    # Actual r=A21/A11 lies between the comparison solutions.
    assert lower_ratio.lo > D("0.13")
    assert upper_ratio.hi < D("0.192")  # 24/125
    assert half_column[0] > 29

    # Exact reversibility gives A=[[a,b],[c,a]], a^2-bc=1.  If
    # r=c/a and x=sqrt(c/b), then x=r/sqrt(1-a^-2)>r.  The rational
    # upper calculation avoids evaluating a square root.
    stable_upper_square = F(24, 125) ** 2 / (1 - F(1, 29**2))
    assert stable_upper_square == F(20184, 546875)
    assert stable_upper_square < F(1, 25)

    # The C159 frame is orthogonal rather than orthonormal.  At the return,
    # |E2|/|E1|=|k0|, where |k0|^2=378/25+3 beta^2.  The certified beta
    # interval lies in (21/10,9/4), so the physical orthonormal-frame slope
    # y=x|k0| is safely in (1/2,2).
    k_squared_lower = F(378, 25) + 3 * F(21, 10) ** 2
    k_squared_upper = F(378, 25) + 3 * F(9, 4) ** 2
    assert k_squared_lower == F(567, 20)
    assert k_squared_upper == F(12123, 400)
    assert F(13, 100) ** 2 * k_squared_lower > F(1, 2) ** 2
    assert F(1, 5) ** 2 * k_squared_upper < 2**2

    # For the physical lines (1,+/-y), the oblique projector norm is
    # (1+y^2)/(2y), the angle sine is its reciprocal, and the eigenbasis
    # condition is max(y,1/y).  Both endpoints y=1/2,2 give the safe
    # bounds below.
    projector_endpoint = (1 + F(2) ** 2) / (2 * F(2))
    assert projector_endpoint == F(5, 4)
    assert 1 / projector_endpoint == F(4, 5)
    assert F(2) == max(F(2), 1 / F(1, 2))

    # Recompute the C192 Collatz floor.  Since det M=1, the other real
    # multiplier is the reciprocal of rho(M), hence is <1/3000.
    second_ratio = base.DOWN.divide(gain_vector[1], D(3) / 20)
    assert gain_vector[0] > 3000
    assert second_ratio > 3000
    assert F(1, 3000) < F(1)

    return lower_ratio, upper_ratio, half_column, gain_vector, second_ratio


def main():
    exact_symmetry_checks()
    lower_ratio, upper_ratio, half_column, gain_vector, second_ratio = certify_stable_line()
    print("C159 exact stable-line certificate: PASS")
    print("certified c/a lower:", lower_ratio.lo)
    print("certified c/a upper:", upper_ratio.hi)
    print("half-map first-column lower:", *half_column)
    print("full-period C192 vector lower:", *gain_vector)
    print("second component ratio:", second_ratio)
    print("CONCLUSION: z_s=(1,-x), 13/100<x<1/5, 1/2<x|k0|<2, lambda_s<1/3000")
    print("BOUNDARY: exact central orbit only; no off-ray Floquet-bundle radius")


if __name__ == "__main__":
    main()
