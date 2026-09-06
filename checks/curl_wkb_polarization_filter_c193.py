#!/usr/bin/env python3
"""Exact algebra for the prospective C193 curl-WKB/polarization bridge.

This dependency-free checker verifies three pieces which do not require
numerical ODE integration:

* the leading and first-corrector Kelvin/pressure recurrences, including
  the parallel/perpendicular projection and every sign;
* the rational coefficients in the first-order residual majorant; and
* the fixed-energy hyperbolic polarization filter which converts a
  q=n^8 cocycle gain into a q^(3/8)=n^3 concentration improvement.

It does not certify the off-ray stable bundle, action-angle flow constants,
Fourier tails, viscosity, or nonlinear stage closure.
"""

from __future__ import annotations

from fractions import Fraction as F


Vec = tuple[F, F, F]
Mat = tuple[Vec, Vec, Vec]


def add(left: Vec, right: Vec) -> Vec:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def sub(left: Vec, right: Vec) -> Vec:
    return tuple(left[index] - right[index] for index in range(3))  # type: ignore[return-value]


def scale(scalar: F, vector: Vec) -> Vec:
    return tuple(scalar * entry for entry in vector)  # type: ignore[return-value]


def dot(left: Vec, right: Vec) -> F:
    return sum((left[index] * right[index] for index in range(3)), F(0))


def cross(left: Vec, right: Vec) -> Vec:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def mvec(matrix: Mat, vector: Vec) -> Vec:
    return tuple(dot(row, vector) for row in matrix)  # type: ignore[return-value]


def parallel(k_value: Vec, vector: Vec) -> Vec:
    rho = dot(k_value, k_value)
    return scale(dot(k_value, vector) / rho, k_value)


def perpendicular(k_value: Vec, vector: Vec) -> Vec:
    return sub(vector, parallel(k_value, vector))


def kelvin_generator(k_value: Vec, gradient: Mat, vector: Vec) -> Vec:
    avalue = mvec(gradient, vector)
    return add(scale(F(-1), avalue), scale(2, parallel(k_value, avalue)))


def check_curl_and_pressure_recurrence() -> None:
    cases: tuple[tuple[Vec, Mat, Vec, Vec], ...] = (
        (
            (F(2), F(-1), F(3)),
            ((F(1), F(2), F(-1)), (F(0), F(-2), F(3)), (F(4), F(1), F(1))),
            (F(1), F(2), F(0)),
            (F(3), F(-2), F(5)),
        ),
        (
            (F(1), F(2), F(-2)),
            ((F(-3), F(1), F(2)), (F(2), F(4), F(-1)), (F(0), F(3), F(-1))),
            (F(2), F(-1), F(0)),
            (F(-4), F(7), F(1)),
        ),
    )

    zero: Vec = (F(0), F(0), F(0))
    for k_value, gradient, raw_b, forcing in cases:
        # Project the test amplitude exactly to k^perp.
        bvalue = perpendicular(k_value, raw_b)
        assert dot(k_value, bvalue) == 0
        rho = dot(k_value, k_value)

        # c=-k cross b/|k|^2 gives k cross c=b exactly.  Consequently
        # (iq)^-1 curl(e^{iq phi}c) has leading amplitude b and is an
        # exact curl, hence exactly divergence free.
        cvalue = scale(F(-1, 1) / rho, cross(k_value, bvalue))
        assert cross(k_value, cvalue) == bvalue

        # For either WKB order, b' = Lb-P_perp f.  If pi=i*sigma with
        # sigma=(2 k.Ab+k.f)/|k|^2, then i*k*pi=-k*sigma.  This cancels
        # the complete order equation, not merely its transverse part.
        bprime = sub(
            kelvin_generator(k_value, gradient, bvalue),
            perpendicular(k_value, forcing),
        )
        avalue = mvec(gradient, bvalue)
        sigma = (2 * dot(k_value, avalue) + dot(k_value, forcing)) / rho
        pressure_gradient_leading = scale(-sigma, k_value)
        order_equation = add(add(forcing, add(bprime, avalue)), pressure_gradient_leading)
        assert order_equation == zero

        # The recurrence preserves k.b=0 when k'=-A^T k.  Written without
        # constructing A^T, k'.b=-k.Ab and k.b'=k.Ab exactly.
        assert -dot(k_value, avalue) + dot(k_value, bprime) == 0


def check_residual_majorant_coefficients() -> None:
    # The primitive pointwise/Frobenius estimates use sqrt(2)<3/2.
    # Collecting D_t curl(c), A curl(c), and grad(pi) yields the rounded
    # rational first-order coefficients used in the theorem.
    sqrt2_upper = F(3, 2)
    a1_coefficient = 4 * sqrt2_upper + 2
    d_b1_coefficient = 6 * sqrt2_upper + 2
    d_k1_coefficient = 26 * sqrt2_upper + 6
    assert a1_coefficient == 8
    assert d_b1_coefficient == 11 < F(25, 2)
    assert d_k1_coefficient == 45 < F(99, 2)

    # Thus, with D=||A||, A1=||DA||, K_-=inf|k|,
    # K1=||Dk||_infty, B0=||b||_2, and B1=||Db||_2,
    # ||f0||_2 is bounded by the following deliberately widened row.
    widened_b1 = F(25, 2)
    widened_a1 = F(8)
    widened_k1 = F(99, 2)
    assert widened_b1 > d_b1_coefficient
    assert widened_a1 == a1_coefficient
    assert widened_k1 > d_k1_coefficient

    # The exact A2 global jet bound D=6 lies below C192's first-order
    # exponent threshold by the exact rational margin 8/57.
    assert F(350, 57) - 6 == F(8, 57) > 0


def check_hyperbolic_filter() -> None:
    # Abstract exact principal model.  M has eigenvalues lambda and
    # lambda^-1.  After R returns put G=lambda^R.  C192 plus one extra
    # return gives the strict floor G>3000*n^3 for every n>=2.  The exact
    # product Dirichlet profile with side length 4*n^2 has L2 norm one and
    # Linf=8*n^3; the broad profile psi=1 has both norms equal to one.
    # These coefficient identities prove the Dirichlet normalization and
    # concentration inequalities for every positive integer n; the finite
    # loop below is only an exact-arithmetic diagnostic.
    assert 4**3 == 8**2
    assert F(1) + F(8, 3000) < F(251, 250)
    assert F(7) / F(251, 250) == F(1750, 251) > 1

    # The endpoint norm inequalities are affine in the overlap c.  Their
    # two remainders are 2(c+1)/G and 2(1-c)/G, nonnegative for |c|<=1.
    for overlap_endpoint in (F(-1), F(1)):
        assert 2 * (overlap_endpoint + 1) >= 0
        assert 2 * (1 - overlap_endpoint) >= 0

    for n_value in (2, 3, 5, 11, 101):
        side_length = 4 * n_value**2
        coefficient_count = side_length**3
        normalization_squared = side_length**3
        assert coefficient_count == normalization_squared
        assert (8 * n_value**3) ** 2 == side_length**3
        target = F(n_value**3)  # q=n^8, so q^(3/8)=n^3.
        gain_floor = F(3000 * n_value**3)
        assert gain_floor == 3000 * target

        # The actual gain G is strictly larger than this floor.  Hence
        # the initial local peak is at most 8/3000, while the broad peak
        # is one.  At the endpoint the reverse triangle inequality gives
        # a local peak larger than 8*n^3-1/G>7*n^3.  The common L2
        # normalization cancels in the quotient ratio.
        initial_linf_upper = F(1) + F(8, 3000)
        final_linf_lower = F(7) * target
        assert initial_linf_upper < F(251, 250)
        improvement_lower = final_linf_lower / F(251, 250)
        assert improvement_lower == F(1750, 251) * target
        assert improvement_lower > target

        # Profiles and eigenvectors need not be orthogonal.  If
        # c=Re<phi e_+,psi e_->, |c|<=1, the common endpoint squared norm
        # before normalization is D=1+G^-2+2c/G and obeys
        # (1-G^-1)^2<=D<=(1+G^-1)^2.  Check exact samples including both
        # extremal overlaps; the universal identity is cross-multiplied
        # coefficientwise below.
        gain = gain_floor
        for overlap in (F(-1), F(-2, 7), F(0), F(3, 8), F(1)):
            denominator_squared = 1 + gain**-2 + 2 * overlap / gain
            assert (1 - 1 / gain) ** 2 <= denominator_squared
            assert denominator_squared <= (1 + 1 / gain) ** 2

        # At discrete intermediate returns, x=lambda^(2j) lies in
        # [1,G^2].  The cross coefficient is always 2c/G, so it cancels
        # from endpoint minus intermediate energy.  After multiplication
        # by G^2*x, the universal identity has the following four exact
        # coefficients in the monomials G^2*x, x, x^2, G^2.
        cross_multiplied_left = {
            (2, 1): F(1),
            (0, 1): F(1),
            (0, 2): F(-1),
            (2, 0): F(-1),
        }
        # Independently expand (x-1)(G^2-x).
        cross_multiplied_right = {}
        for (g_power, x_power), coefficient in {
            (0, 1): F(1),
            (0, 0): F(-1),
        }.items():
            for (other_g_power, other_x_power), other_coefficient in {
                (2, 0): F(1),
                (0, 1): F(-1),
            }.items():
                monomial = (g_power + other_g_power, x_power + other_x_power)
                cross_multiplied_right[monomial] = (
                    cross_multiplied_right.get(monomial, F(0))
                    + coefficient * other_coefficient
                )
        assert cross_multiplied_left == cross_multiplied_right

        # Representative exact values additionally check the sign on the
        # closed interval.  The proof uses the displayed factorization,
        # not these samples as a substitute for it.
        for x_value in (F(1), gain, gain**2):
            assert 1 <= x_value <= gain**2
            difference = (1 + gain**-2) - (x_value / gain**2 + 1 / x_value)
            factorized = (x_value - 1) * (gain**2 - x_value) / (gain**2 * x_value)
            assert difference == factorized >= 0

    # The one-extra-return clock changes only the additive constant:
    # T(R_delta+1)<(57/400)log q+2*(76/25).
    assert 2 * F(76, 25) == F(152, 25)


def main() -> None:
    check_curl_and_pressure_recurrence()
    print("C193 curl/Kelvin/pressure recurrence: PASS")
    check_residual_majorant_coefficients()
    print("C193 explicit first-order residual coefficients: PASS")
    check_hyperbolic_filter()
    print("C193 fixed-energy hyperbolic polarization filter: PASS")
    print("OPEN: off-ray stable bundle, finite-frequency composite, band, viscosity, nonlinear stage")


if __name__ == "__main__":
    main()
