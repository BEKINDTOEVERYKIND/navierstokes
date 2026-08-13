#!/usr/bin/env python3
"""Dependency-free exact checks for C169's fixed-plane 2D3C no-go.

The checker verifies with formal polynomial arithmetic that the nonlinear
convective term is triangular and that the pressure source contains no
transverse scalar.  It also checks the real and complex passive-scalar
modulus identities, the sharp finite-mode Cauchy--Schwarz identity, the
signed/half-lattice reality factors, the factor-two full-solution difference
bound, and the C161 norm/energy power ledger.  It does not claim an
in-plane, off-plane, localized, or one-cell result.
"""

from fractions import Fraction as Q


# A tiny exact commutative polynomial algebra.  A monomial is a sorted tuple
# of variable names; a polynomial is a dictionary monomial -> rational.
def var(name):
    return {(name,): Q(1)}


def add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, Q(0)) + coefficient
        if out[monomial] == 0:
            del out[monomial]
    return out


def neg(value):
    return {monomial: -coefficient for monomial, coefficient in value.items()}


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    out = {}
    for monomial_left, coefficient_left in left.items():
        for monomial_right, coefficient_right in right.items():
            monomial = tuple(sorted(monomial_left + monomial_right))
            out[monomial] = out.get(monomial, Q(0)) + coefficient_left * coefficient_right
            if out[monomial] == 0:
                del out[monomial]
    return out


def scale(value, coefficient):
    return {
        monomial: Q(coefficient) * entry
        for monomial, entry in value.items()
        if Q(coefficient) * entry != 0
    }


def substitute_scaled_variable(value, old_name, new_name, factor):
    """Substitute old_name=factor*new_name in an exact polynomial."""
    out = {}
    for monomial, coefficient in value.items():
        count = monomial.count(old_name)
        replaced = tuple(sorted(
            tuple(name for name in monomial if name != old_name)
            + tuple(new_name for _ in range(count))
        ))
        out[replaced] = out.get(replaced, Q(0)) + coefficient * Q(factor) ** count
        if out[replaced] == 0:
            del out[replaced]
    return out


def check_full_nonlinear_triangular_system():
    # Coordinate order is (r,t,z).  Every t derivative is exactly zero.
    velocity = (var("vr"), var("Theta"), var("vz"))
    derivative = {
        (0, 0): var("vr_r"),
        (0, 1): var("Theta_r"),
        (0, 2): var("vz_r"),
        (1, 0): {},
        (1, 1): {},
        (1, 2): {},
        (2, 0): var("vr_z"),
        (2, 1): var("Theta_z"),
        (2, 2): var("vz_z"),
    }

    convection = []
    for component in range(3):
        value = {}
        for direction in range(3):
            value = add(value, mul(velocity[direction], derivative[(direction, component)]))
        convection.append(value)

    assert convection[0] == add(mul(var("vr"), var("vr_r")), mul(var("vz"), var("vr_z")))
    assert convection[2] == add(mul(var("vr"), var("vz_r")), mul(var("vz"), var("vz_z")))
    assert convection[1] == add(
        mul(var("vr"), var("Theta_r")),
        mul(var("vz"), var("Theta_z")),
    )

    # The pressure Poisson source sum_{i,j} d_i u_j d_j u_i has no Theta
    # jet because every term with i=t or j=t contains a t derivative.
    pressure_source = {}
    for i in range(3):
        for j in range(3):
            pressure_source = add(
                pressure_source,
                mul(derivative[(i, j)], derivative[(j, i)]),
            )
    expected = add(
        add(mul(var("vr_r"), var("vr_r")), mul(var("vz_z"), var("vz_z"))),
        scale(mul(var("vz_r"), var("vr_z")), 2),
    )
    assert pressure_source == expected
    assert all("Theta" not in name for monomial in pressure_source for name in monomial)

    # Incompressibility is vr_r+vz_z=0.  After setting vz_z=-vr_r, the
    # pressure source is still purely in-plane and has the familiar 2D form.
    divergence = add(var("vr_r"), var("vz_z"))
    assert divergence == {("vr_r",): Q(1), ("vz_z",): Q(1)}
    expected_divergence_free = scale(
        add(mul(var("vr_r"), var("vr_r")), mul(var("vz_r"), var("vr_z"))),
        2,
    )
    assert substitute_scaled_variable(
        pressure_source,
        "vz_z",
        "vr_r",
        -1,
    ) == expected_divergence_free


def check_complex_passive_modulus_identity():
    # Theta=X+iY.  Insert
    # X_t=nu Delta X-v.grad X, Y_t=nu Delta Y-v.grad Y
    # into L(|Theta|^2), with L=d_t+v.grad-nu Delta.
    nu = var("nu")
    vr = var("vr")
    vz = var("vz")
    x = var("X")
    y = var("Y")
    xr = var("X_r")
    xz = var("X_z")
    yr = var("Y_r")
    yz = var("Y_z")
    lap_x = var("Delta_X")
    lap_y = var("Delta_Y")

    xt = sub(sub(mul(nu, lap_x), mul(vr, xr)), mul(vz, xz))
    yt = sub(sub(mul(nu, lap_y), mul(vr, yr)), mul(vz, yz))
    time_part = scale(add(mul(x, xt), mul(y, yt)), 2)
    advection_part = scale(
        add(
            mul(vr, add(mul(x, xr), mul(y, yr))),
            mul(vz, add(mul(x, xz), mul(y, yz))),
        ),
        2,
    )
    # Delta(X^2+Y^2)=2 X Delta X+2 Y Delta Y
    #                       +2(|grad X|^2+|grad Y|^2).
    diffusion_part = scale(
        mul(
            nu,
            add(
                add(mul(x, lap_x), mul(y, lap_y)),
                add(add(mul(xr, xr), mul(xz, xz)), add(mul(yr, yr), mul(yz, yz))),
            ),
        ),
        2,
    )
    lhs = sub(add(time_part, advection_part), diffusion_part)
    rhs = scale(
        mul(nu, add(add(mul(xr, xr), mul(xz, xz)), add(mul(yr, yr), mul(yz, yz)))),
        -2,
    )
    assert lhs == rhs


def check_real_passive_modulus_and_extremum_signs():
    # For a real Theta, the same exact product-rule computation gives
    # L(Theta^2)=-2 nu |grad Theta|^2.  This is distinct from claiming that
    # a difference of two scalars transported by different v obeys L.
    nu = var("nu")
    vr = var("vr")
    vz = var("vz")
    theta = var("Theta")
    theta_r = var("Theta_r")
    theta_z = var("Theta_z")
    lap_theta = var("Delta_Theta")

    theta_t = sub(
        sub(mul(nu, lap_theta), mul(vr, theta_r)),
        mul(vz, theta_z),
    )
    time_part = scale(mul(theta, theta_t), 2)
    advection_part = scale(
        add(mul(vr, mul(theta, theta_r)), mul(vz, mul(theta, theta_z))),
        2,
    )
    diffusion_part = scale(
        mul(
            nu,
            add(
                mul(theta, lap_theta),
                add(mul(theta_r, theta_r), mul(theta_z, theta_z)),
            ),
        ),
        2,
    )
    lhs = sub(add(time_part, advection_part), diffusion_part)
    rhs = scale(mul(nu, add(mul(theta_r, theta_r), mul(theta_z, theta_z))), -2)
    assert lhs == rhs

    # At a spatial maximum grad Theta=0 and Delta Theta<=0; at a minimum the
    # Laplacian has the opposite sign.  These rational samples check the
    # direction used by the real maximum principle.
    viscosity = Q(5, 7)
    lap_at_max = Q(-3, 2)
    lap_at_min = -lap_at_max
    assert viscosity * lap_at_max < 0
    assert viscosity * lap_at_min > 0


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def abs_sq(a):
    return cmul(a, cconj(a))[0]


def check_sharp_q_mode_bound():
    # Exact Gaussian-rational unit phases and coefficients.  The identity
    # q sum |c_j|^2-|sum c_j z_j|^2
    #   =sum_{j<k}|c_j z_j-c_k z_k|^2
    # proves Cauchy--Schwarz and also shows equality when all c_j z_j agree.
    phases = (
        (Q(1), Q(0)),
        (Q(0), Q(1)),
        (Q(-1), Q(0)),
        (Q(0), Q(-1)),
        (Q(3, 5), Q(4, 5)),
        (Q(3, 5), Q(-4, 5)),
        (Q(5, 13), Q(12, 13)),
        (Q(5, 13), Q(-12, 13)),
    )
    coefficients = (
        (Q(2, 3), Q(1, 7)),
        (Q(-1, 5), Q(3, 11)),
        (Q(4, 9), Q(-2, 13)),
        (Q(1, 4), Q(5, 17)),
        (Q(-3, 8), Q(7, 19)),
        (Q(2, 15), Q(-4, 21)),
        (Q(5, 16), Q(1, 23)),
        (Q(-7, 18), Q(-2, 25)),
    )
    assert all(abs_sq(phase) == 1 for phase in phases)
    values = tuple(cmul(coefficient, phase) for coefficient, phase in zip(coefficients, phases))
    total = (Q(0), Q(0))
    for value in values:
        total = cadd(total, value)
    left = Q(len(values)) * sum((abs_sq(value) for value in values), Q(0)) - abs_sq(total)
    right = Q(0)
    for j, value in enumerate(values):
        for other in values[j + 1 :]:
            right += abs_sq(csub(value, other))
    assert left == right and left > 0

    # Sharpness: choosing c_j=conj(z_j) makes every c_j z_j=1.
    aligned = tuple(cmul(cconj(phase), phase) for phase in phases)
    aligned_sum = (sum((value[0] for value in aligned), Q(0)), Q(0))
    assert aligned_sum == (Q(len(phases)), Q(0))
    assert abs_sq(aligned_sum) == Q(len(phases)) * sum(
        (abs_sq(cconj(phase)) for phase in phases), Q(0)
    )


def check_signed_and_half_lattice_reality_factors():
    # q_s=2h signed nonzero modes.  Equal real coefficients obey the reality
    # relation and saturate the bound at x=0.  G_t is the full signed l2 norm;
    # G_{t,+} is the independent half-lattice norm.  Their squares are
    # energies and differ by exactly two.
    h = 4
    q_signed = 2 * h
    amplitude = Q(3, 11)
    positive = tuple(range(1, h + 1))
    frequencies = tuple(-entry for entry in reversed(positive)) + positive
    coefficients = {frequency: amplitude for frequency in frequencies}
    assert len(coefficients) == q_signed
    assert all(coefficients[-frequency] == coefficients[frequency] for frequency in positive)

    half_norm_sq = sum((coefficients[frequency] ** 2 for frequency in positive), Q(0))
    signed_norm_sq = sum((coefficient**2 for coefficient in coefficients.values()), Q(0))
    point_value = sum(coefficients.values(), Q(0))
    assert signed_norm_sq == 2 * half_norm_sq
    assert point_value**2 == q_signed * signed_norm_sq
    assert point_value**2 == 4 * h * half_norm_sq


def check_full_solution_difference_factor_two():
    # The triangle estimate |Theta^S-Theta^0|<=2||Theta_0||_infty is sharp
    # from the two separate maximum principles.  The values +/-A occur for
    # inviscid translations of cos(z) separated by pi.
    initial_supremum = Q(5, 9)
    first_value = initial_supremum
    second_value = -initial_supremum
    assert abs(first_value) <= initial_supremum
    assert abs(second_value) <= initial_supremum
    assert abs(first_value - second_value) == 2 * initial_supremum


def check_reservoir_lower_bound_conventions():
    # Square the lower bounds to avoid irrational square roots.  The target
    # is L=c*b*q^(3/2), so L^2=c^2*b^2*q^3.  Dividing by q_s is the required
    # full signed transverse coefficient energy G_t^2.  A full-solution
    # difference has an
    # extra factor 2 in amplitude and therefore 4 in energy.
    q = Q(8)
    c_star = Q(3, 5)
    b = Q(2, 7)
    target_sq = c_star**2 * b**2 * q**3

    q_signed_c161 = q
    required_transverse_energy_signed = target_sq / q_signed_c161
    assert required_transverse_energy_signed == c_star**2 * b**2 * q**2

    q_signed_from_half_lattice = 2 * q
    required_transverse_energy_half_convention = target_sq / q_signed_from_half_lattice
    assert required_transverse_energy_half_convention == c_star**2 * b**2 * q**2 / 2
    assert required_transverse_energy_signed / 4 == target_sq / (4 * q_signed_c161)


def check_c161_power_ledger():
    # Store monomials as exponents of n.  q=n^8 and b=n^-2.
    q = 8
    b = -2
    gate_upper = b + q // 2
    required_gate_norm = b + q
    gate_budget_energy = 2 * b
    required_gate_energy = 2 * required_gate_norm
    target = b + 3 * q // 2
    assert gate_upper == 2                 # b sqrt(q)
    assert required_gate_norm == 6         # b q
    assert gate_budget_energy == -4        # b^2
    assert required_gate_energy == 12      # b^2 q^2
    assert target == 10                    # b q^(3/2)
    assert target - gate_upper == q        # exact missing factor q
    assert target - q // 2 == required_gate_norm
    assert required_gate_energy - gate_budget_energy == 2 * q


def main():
    check_full_nonlinear_triangular_system()
    check_complex_passive_modulus_identity()
    check_real_passive_modulus_and_extremum_signs()
    check_sharp_q_mode_bound()
    check_signed_and_half_lattice_reality_factors()
    check_full_solution_difference_factor_two()
    check_reservoir_lower_bound_conventions()
    check_c161_power_ledger()
    print("PASS C169: the full fixed-plane Navier--Stokes system is triangular 2D3C")
    print("PASS C169: real and complex transverse total fields obey exact maximum principles")
    print("PASS C169: signed/half-lattice gate bounds have the exact reality factors")
    print("PASS C169: the two-full-solution difference bound has the sharp factor two")
    print("PASS C169: an O(b) coefficient norm is short of b q^(3/2) by q")
    print("Fixed-plane transverse no-go only; in-plane and off-plane conversion remain open")


if __name__ == "__main__":
    main()
