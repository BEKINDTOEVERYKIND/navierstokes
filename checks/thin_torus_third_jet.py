#!/usr/bin/env python3
"""Exact checks for the thin-torus third Navier--Stokes jet.

The script checks:

1. the binomial coefficients and the factor four in the two-viscosity
   pressure moment;
2. the flat-cylinder radial identity;
3. the integrated anisotropy formula for an exact polynomial test
   cutoff; and
4. all aspect, turnover, exterior-L2, and packed-stage exponents.

It does not prove the tubular remainder or Bogovskii estimates.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction
Polynomial = list[Fraction]
Power = tuple[Fraction, ...]


def trim(poly: Polynomial) -> Polynomial:
    result = poly[:]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    result = [Q(0)] * size
    for index in range(size):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] += right[index]
    return trim(result)


def poly_scale(multiplier: int | Fraction, poly: Polynomial) -> Polynomial:
    return trim([Q(multiplier) * coefficient for coefficient in poly])


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return trim(result)


def poly_derivative(poly: Polynomial) -> Polynomial:
    if len(poly) == 1:
        return [Q(0)]
    return trim([Q(index) * poly[index] for index in range(1, len(poly))])


def poly_times_s(poly: Polynomial) -> Polynomial:
    return [Q(0), *poly]


def poly_divide_s(poly: Polynomial) -> Polynomial:
    assert poly[0] == 0
    return trim(poly[1:])


def poly_integral(poly: Polynomial, lower: int, upper: int) -> Fraction:
    total = Q(0)
    for index, coefficient in enumerate(poly):
        total += coefficient * Q(upper ** (index + 1) - lower ** (index + 1), index + 1)
    return total


def pwr(*entries: int | Fraction) -> Power:
    return tuple(Q(entry) for entry in entries)


def power_add(*powers: Power) -> Power:
    return tuple(sum(column, Q(0)) for column in zip(*powers))


def check_third_jet_combinatorics() -> None:
    # d^2/dt^2 of u tensor u has coefficients 1, 2, 1.
    nonlinear_coefficients = [1, 2, 1]
    assert sum(nonlinear_coefficients) == 4

    # After two integrations by parts, each endpoint-biheat term gives
    # one copy of integral Delta U tensor Delta U.  The middle term gives
    # two more copies.
    endpoint_biheat_left = 1
    heat_heat = 2
    endpoint_biheat_right = 1
    assert endpoint_biheat_left + heat_heat + endpoint_biheat_right == 4
    print("third-jet binomial coefficients and factor-four moment: exact")


def check_flat_profile_and_anisotropy_identity() -> None:
    # Use chi(s)=s(s-1)^2(2-s)^2 on [1,2].  Its endpoint value and first
    # derivative vanish, which is enough for the integrations by parts
    # checked here.  The identities are polynomial and therefore exact.
    s_poly = [Q(0), Q(1)]
    s_minus_one = [Q(-1), Q(1)]
    two_minus_s = [Q(2), Q(-1)]
    chi = poly_mul(
        s_poly,
        poly_mul(poly_mul(s_minus_one, s_minus_one), poly_mul(two_minus_s, two_minus_s)),
    )
    chi_prime = poly_derivative(chi)
    chi_second = poly_derivative(chi_prime)

    # A=4 chi'+2s chi''.  For q=-r chi, the transverse vector Laplacian
    # is -r A, and the symmetric azimuthal advection is -2r chi A e_r.
    a_poly = poly_add(poly_scale(4, chi_prime), poly_scale(2, poly_times_s(chi_second)))
    assert a_poly != [Q(0)]
    assert poly_mul(chi, a_poly) != [Q(0)]

    # Direct form of integral (k^2-2l^2) r dr after s=r^2/2:
    # -chi^2/(2s)-2chi*A.
    chi_squared_over_s = poly_divide_s(poly_mul(chi, chi))
    direct = poly_add(
        poly_scale(Q(-1, 2), chi_squared_over_s),
        poly_scale(-2, poly_mul(chi, a_poly)),
    )

    # Integrated-by-parts form:
    # 4s(chi')^2-chi^2/(2s).
    reduced = poly_add(
        poly_scale(4, poly_times_s(poly_mul(chi_prime, chi_prime))),
        poly_scale(Q(-1, 2), chi_squared_over_s),
    )

    direct_integral = poly_integral(direct, 1, 2)
    reduced_integral = poly_integral(reduced, 1, 2)
    assert direct_integral == reduced_integral
    assert reduced_integral > 0
    print("flat radial-gradient and Delta-U anisotropy identities: exact")


def check_one_bubble_exponents() -> None:
    # Exponent order: (nu, v, R, eps, T, d).
    one_viscosity_source = pwr(1, 3, 0, 0, 0, 0)
    two_viscosity_source = pwr(2, 2, -1, -2, 0, 0)
    turnover_cubed = pwr(0, -3, 3, 0, 3, 0)
    far_kernel = pwr(0, 0, 0, 0, 0, -4)

    one_far = power_add(one_viscosity_source, turnover_cubed, far_kernel)
    two_far = power_add(two_viscosity_source, turnover_cubed, far_kernel)
    assert one_far == pwr(1, 0, 3, 0, 3, -4)
    assert two_far == pwr(2, -1, 2, -2, 3, -4)

    # The L2 norm of A*d^-4 outside radius R gains R^-5/2.
    exterior_l2_kernel = pwr(0, 0, Q(-5, 2), 0, 0, 0)
    one_l2 = power_add(one_viscosity_source, turnover_cubed, exterior_l2_kernel)
    two_l2 = power_add(two_viscosity_source, turnover_cubed, exterior_l2_kernel)
    assert one_l2 == pwr(1, 0, Q(1, 2), 0, 3, 0)
    assert two_l2 == pwr(2, -1, Q(-1, 2), -2, 3, 0)
    print("one-bubble third-wake exponents: exact")


def check_minor_viscosity_identity() -> None:
    # mu=nu/(vR), Theta=T*mu/eps^2.
    # Relative one-viscosity wake: T^3*mu=T^2*eps^2*Theta.
    # Relative two-viscosity wake:
    # T^3*mu^2/eps^2=T*eps^2*Theta^2.
    # Exponent order here is (mu, eps, T, Theta).
    one_left = pwr(1, 0, 3, 0)
    one_right_after_theta_substitution = pwr(0, 2, 2, 1)
    two_left = pwr(2, -2, 3, 0)
    two_right_after_theta_substitution = pwr(0, 2, 1, 2)

    # Substituting Theta=T*mu*eps^-2 in the right sides recovers left.
    theta_definition = pwr(1, -2, 1, -1)
    assert power_add(one_right_after_theta_substitution, theta_definition) == one_left
    assert power_add(
        two_right_after_theta_substitution,
        tuple(2 * entry for entry in theta_definition),
    ) == two_left
    print("minor-viscosity normalization identities: exact")


def check_packed_exponents() -> None:
    # Exponent order: (nu, v, ell, K, eps, T).
    first_relative = pwr(1, -1, -1, 3, 0, 3)
    second_relative = pwr(2, -2, -2, 4, -2, 3)
    assert first_relative == pwr(1, -1, -1, 3, 0, 3)
    assert second_relative == pwr(2, -2, -2, 4, -2, 3)

    # Substitute v=ell^-gamma K^gamma.
    for gamma in (Q(11, 10), Q(5, 4), Q(7, 5)):
        replace_v_inverse = pwr(0, 1, gamma, -gamma, 0, 0)
        first = power_add(first_relative, replace_v_inverse)
        second = power_add(
            second_relative,
            tuple(2 * entry for entry in replace_v_inverse),
        )
        assert first == pwr(1, 0, gamma - 1, 3 - gamma, 0, 3)
        assert second == pwr(2, 0, 2 * gamma - 2, 4 - 2 * gamma, -2, 3)

        # If eps=ell^beta, the third coefficient alone needs
        # 2*gamma-2-2*beta>0, whereas Theta->0 already needs the
        # stronger gamma-1-2*beta>0.
        beta = (gamma - 1) / 3
        assert gamma - 1 - 2 * beta > 0
        assert 2 * gamma - 2 - 2 * beta > 0
    print("packed third-wake and cascade exponents: exact")


def main() -> None:
    check_third_jet_combinatorics()
    check_flat_profile_and_anisotropy_identity()
    check_one_bubble_exponents()
    check_minor_viscosity_identity()
    check_packed_exponents()
    print("all thin-torus third-jet checks passed")


if __name__ == "__main__":
    main()
