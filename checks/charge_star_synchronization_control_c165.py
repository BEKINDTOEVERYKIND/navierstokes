#!/usr/bin/env python3
"""Dependency-free exact checks for C165.

The checker verifies the rigid sin(2 phi) coefficient, finite-sheet source
preweighting, the polynomial identities and fully explicit rational
constants used by the compact-ensemble adiabatic lemma, the rotating-frame
chirp, and the J/b resource ledger.  The tangential full-ladder obstruction
is imported from C164 and is not duplicated here.  This checker does not
assert an unforced Navier--Stokes realization of the prescribed controls.
"""

from fractions import Fraction as F


def kappa(y):
    y = F(y)
    return y * y * (y - 2) / (4 * (y * y + 2 * y + 4))


def return_from_x(x, y):
    """C163 R after x=sin(2 phi), without trigonometric arithmetic."""
    x = F(x)
    value = kappa(y)
    return (1 - x) * (1 + (1 + x) * value) / 2


def check_rigid_rate_coefficient():
    heights = [F(0), F(1, 3), F(2), F(7, 2), F(9)]
    weights = [F(2, 7), F(5, 11), F(3, 5), F(7, 13), F(4, 9)]
    total_weight = sum(weights, F(0))
    height_moment = sum(
        (weight * kappa(height) for weight, height in zip(weights, heights)),
        F(0),
    )

    # Exact coefficients of sum_j w_j R(x,y_j).
    constant = (total_weight + height_moment) / 2
    linear = -total_weight / 2
    quadratic = -height_moment / 2
    assert total_weight > 0
    assert linear < 0

    for x in (F(-1), F(-3, 5), F(-1, 7), F(0), F(2, 9), F(4, 5), F(1)):
        direct = sum(
            (
                weight * return_from_x(x, height)
                for weight, height in zip(weights, heights)
            ),
            F(0),
        )
        polynomial = constant + linear * x + quadratic * x * x
        assert direct == polynomial

        reflected = sum(
            (
                weight * return_from_x(-x, height)
                for weight, height in zip(weights, heights)
            ),
            F(0),
        )
        assert direct - reflected == -total_weight * x

    # The same coefficient identity holds for arbitrary nonnegative rational
    # mixtures, including zero weights and repeated heights.
    test_families = (
        ([F(1)], [F(0)]),
        ([F(0), F(3), F(0)], [F(1), F(2), F(10)]),
        ([F(1, 100), F(99, 100)], [F(1, 5), F(20)]),
    )
    for family_weights, family_heights in test_families:
        weight_sum = sum(family_weights, F(0))
        moment = sum(
            (
                weight * kappa(height)
                for weight, height in zip(family_weights, family_heights)
            ),
            F(0),
        )
        assert -weight_sum / 2 < 0
        for x in (F(-2, 3), F(1, 4), F(3, 4)):
            direct = sum(
                (
                    weight * return_from_x(x, height)
                    for weight, height in zip(family_weights, family_heights)
                ),
                F(0),
            )
            assert direct == (
                (weight_sum + moment) / 2
                - weight_sum * x / 2
                - moment * x * x / 2
            )


def polynomial_add(left, right):
    size = max(len(left), len(right))
    output = [F(0)] * size
    for index in range(size):
        if index < len(left):
            output[index] += left[index]
        if index < len(right):
            output[index] += right[index]
    while len(output) > 1 and output[-1] == 0:
        output.pop()
    return output


def polynomial_scale(poly, scalar):
    return [F(scalar) * coefficient for coefficient in poly]


def polynomial_multiply(left, right):
    output = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i + j] += a * b
    return output


def polynomial_derivative(poly):
    if len(poly) == 1:
        return [F(0)]
    return [F(index) * poly[index] for index in range(1, len(poly))]


def check_adiabatic_algebra_and_bounds():
    # Work in the variable g=2s-1.  Since d/ds=2d/dg,
    # f=1-g^2 and f'=-4g exactly.
    g = [F(0), F(1)]
    one = [F(1)]
    g_squared = polynomial_multiply(g, g)
    f = polynomial_add(one, polynomial_scale(g_squared, -1))
    f_s = polynomial_scale(polynomial_derivative(f), 2)
    g_s = polynomial_scale(polynomial_derivative(g), 2)
    assert f == [F(1), F(0), F(-1)]
    assert f_s == [F(0), F(-4)]
    assert g_s == [F(2)]

    # g f' - f g' = -2(1+g^2), the numerator in theta'.
    angle_numerator = polynomial_add(
        polynomial_multiply(g, f_s),
        polynomial_scale(polynomial_multiply(f, g_s), -1),
    )
    assert angle_numerator == [F(-2), F(0), F(-2)]
    angle_numerator_s = polynomial_scale(
        polynomial_derivative(angle_numerator), 2
    )
    assert angle_numerator_s == [F(0), F(-8)]

    # Instantiate one rational compact interval.  These are explicit rational
    # upper bounds, not sampled numerical premises.
    lambda_minus = F(1, 2)
    lambda_plus = F(2)
    m = min(F(1, 2), 3 * lambda_minus / 4)
    assert m == F(3, 8)
    e_derivative_bound = (2 + 4 * lambda_plus * lambda_plus) / m
    theta_one_bound = 4 * lambda_plus / (m * m)
    theta_two_bound = (
        8 * lambda_plus / (m * m)
        + 8 * lambda_plus * e_derivative_bound / (m**3)
    )
    endpoint_term = theta_one_bound / (2 * m)
    theta_two_term = theta_two_bound / (4 * m)
    gap_derivative_term = (
        theta_one_bound * e_derivative_bound / (4 * m * m)
    )
    feedback_term = theta_one_bound * theta_one_bound / (8 * m)
    constant = (
        endpoint_term
        + theta_two_term
        + gap_derivative_term
        + feedback_term
    )
    assert e_derivative_bound == 48
    assert theta_one_bound == F(512, 9)
    assert theta_two_bound == F(44032, 3)
    assert endpoint_term == F(2048, 27)
    assert theta_two_term == F(88064, 9)
    assert gap_derivative_term == F(131072, 27)
    assert feedback_term == F(262144, 243)
    assert constant == F(3837952, 243)

    # Endpoint eigenstates are exact: f=0 and g=-1,+1.
    for s, expected_f, expected_g in (
        (F(0), F(0), F(-1)),
        (F(1), F(0), F(1)),
    ):
        actual_f = 4 * s * (1 - s)
        actual_g = 2 * s - 1
        assert (actual_f, actual_g) == (expected_f, expected_g)

    # The leading adiabatic dynamical phase cannot be treated as a common
    # phase premise: E_lambda(s)^2 strictly increases with lambda at every
    # interior point where f is nonzero.
    for s in (F(1, 4), F(1, 2), F(3, 4)):
        f_value = 4 * s * (1 - s)
        g_value = 2 * s - 1
        lower_square = g_value * g_value + lambda_minus**2 * f_value**2
        upper_square = g_value * g_value + lambda_plus**2 * f_value**2
        assert f_value > 0
        assert upper_square - lower_square == (
            lambda_plus**2 - lambda_minus**2
        ) * f_value**2
        assert upper_square > lower_square

    # The two-case lower gap proof uses only these exact implications:
    # |g|>=1/2 -> E>=1/2; |g|<=1/2 -> f=1-g^2>=3/4.
    assert 1 - F(1, 2) ** 2 == F(3, 4)
    assert min(F(1, 2), 3 * lambda_minus / 4) == m

    return constant


def cmul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def cscale(value, scalar):
    return (F(scalar) * value[0], F(scalar) * value[1])


def check_finite_sheet_preweighting():
    # Exact two-level rotations with rational sine/cosine pairs.  These stand
    # for finitely many compact-interval rates at one common time.  Preweighting by
    # i/sine makes every daughter exactly the same prescribed complex value,
    # while the source retains i*cot(angle) times that value.
    desired = (F(7, 5), F(-2, 3))
    minus_i = (F(0), F(-1))
    i_unit = (F(0), F(1))
    rotations = (
        (F(3, 5), F(4, 5)),
        (F(5, 13), F(12, 13)),
        (F(8, 17), F(15, 17)),
    )
    inverse_sines = []
    for sine, cosine in rotations:
        assert sine > 0 and cosine > 0
        assert sine * sine + cosine * cosine == 1
        source = cscale(cmul(i_unit, desired), 1 / sine)
        daughter = cscale(cmul(minus_i, source), sine)
        source_remnant = cscale(source, cosine)
        assert daughter == desired
        assert source_remnant == cscale(
            cmul(i_unit, desired), cosine / sine
        )
        inverse_sines.append(1 / sine)
    assert min(inverse_sines) > 0
    assert max(inverse_sines) / min(inverse_sines) < 2


def check_phase_chirp_and_resource_ledger(constant):
    # phi(t)=2J*s(1-s), s=t/J.  Therefore dphi/dt=2-4s=-2g(s),
    # and the rotating-frame diagonal term -phi'/2 is g(s).
    for s in (F(0), F(1, 7), F(1, 2), F(5, 6), F(1)):
        g = 2 * s - 1
        phase_derivative = 2 - 4 * s
        assert phase_derivative == -2 * g
        assert -phase_derivative / 2 == g
    assert 2 * F(0) * (1 - F(0)) == 0
    assert 2 * F(1) * (1 - F(1)) == 0

    # phi=2J*s(1-s) rises to J/2 and returns to zero, so its exact total
    # variation is J.  The rotating and lab frames agree at both endpoints.
    for duration in (F(1), F(7, 3), F(19)):
        phase_left = 2 * duration * F(0) * (1 - F(0))
        phase_mid = 2 * duration * F(1, 2) * (1 - F(1, 2))
        phase_right = 2 * duration * F(1) * (1 - F(1))
        assert phase_left == phase_right == 0
        assert phase_mid == duration / 2
        assert (phase_mid - phase_left) + (phase_mid - phase_right) == duration

    # Integral_0^1 4s(1-s) ds = 2/3.  With duration J this is 2J/3.
    envelope_integral = 4 * (F(1, 2) - F(1, 3))
    assert envelope_integral == F(2, 3)

    for n in (2, 3, 5, 10, 101):
        n = F(n)
        b = n ** -2
        steps = n**2
        assert b * steps == 1
        assert constant / steps == constant * b
        assert envelope_integral * steps == F(2, 3) / b
        # Scaling only the transverse envelope makes the middle half-gap
        # b*lambda.  Scaling it and the chirp-derived detuning together makes
        # the effective adiabatic duration bJ exactly one.
        lambda_value = F(3, 2)
        middle_gap = b * lambda_value
        assert middle_gap == lambda_value / steps
        assert b * steps == 1


def main():
    check_rigid_rate_coefficient()
    check_finite_sheet_preweighting()
    constant = check_adiabatic_algebra_and_bounds()
    check_phase_chirp_and_resource_ledger(constant)
    print("PASS C165: every nonzero nonnegative shared mixture has sin(2phi) coefficient -W/2")
    print(
        "PASS C165: bounded source preweights equalize finite-sheet daughters "
        "but leave source remnants"
    )
    print("PASS C165: the common chirped two-level schedule has a uniform explicit O(1/J) bound")
    print("C165 transfer is only up to a rate-dependent phase; coherent focus remains open")
    print("PASS C165: its normalized envelope area is 2J/3, not the C161 O(1) budget")
    print(
        "Prescribed-control algebra only; unforced NS realization and "
        "physical-energy repair remain open"
    )


if __name__ == "__main__":
    main()
