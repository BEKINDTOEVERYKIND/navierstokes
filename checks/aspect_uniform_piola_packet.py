#!/usr/bin/env python3
"""Exact rational ledgers for aspect-uniform thin-torus Piola transport.

This checks the finite-dimensional first-jet identities in the companion
note. It does not test the analytic pseudomode or Sobolev estimates.
"""

from fractions import Fraction as F


def add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def scale(c, u):
    return tuple(c * a for a in u)


def piola_value(u, h):
    return (u[0] / h, u[1] / h, u[2])


def piola_jet(u, du, epsilon, h):
    """Reference derivatives of (u1/h,u2/h,u3)."""
    out = []
    for component in range(3):
        row = []
        for coordinate in range(3):
            value = du[component][coordinate]
            if component < 2:
                value /= h
                if coordinate == 0:
                    value -= epsilon * u[component] / h**2
            row.append(value)
        out.append(tuple(row))
    return tuple(out)


def directional(direction, derivative):
    return tuple(
        sum(direction[j] * derivative[k][j] for j in range(3))
        for k in range(3)
    )


def covariant(direction, field, derivative, epsilon, h):
    """Physical covariant derivative in the rotating orthonormal frame."""
    out = []
    for component in range(3):
        value = (
            direction[0] * derivative[component][0]
            + direction[1] * derivative[component][1]
            + direction[2] * derivative[component][2] / h
        )
        if component == 0:
            value -= epsilon * direction[2] * field[2] / h
        elif component == 2:
            value += epsilon * direction[2] * field[0] / h
        out.append(value)
    return tuple(out)


def geometric_commutator(c, u, du, epsilon, h):
    dc_u = directional(c, du)
    return (
        -(h - 1) * dc_u[0] / h**2
        - epsilon * c[0] * u[0] / h**3
        - epsilon * c[2] * u[2] / h,
        -(h - 1) * dc_u[1] / h**2
        - epsilon * c[0] * u[1] / h**3,
        -(h - 1) * dc_u[2] / h
        + epsilon * c[2] * u[0] / h**2,
    )


def pressure_defect(gradient, h):
    factor = (h - 1) / h
    return (
        factor * gradient[0],
        factor * gradient[1],
        -factor * gradient[2],
    )


def check_case(case):
    epsilon, y1, a, da, u, du, b, db, pressure_gradient = case
    h = 1 + epsilon * y1
    assert h > 0

    pa = piola_value(a, h)
    pu = piola_value(u, h)
    dpa = piola_jet(a, da, epsilon, h)
    dpu = piola_jet(u, du, epsilon, h)

    # The moving-frame derivative matrix is diag(1,1,h), so both its
    # determinant and its Euclidean volume Jacobian are exactly h.
    derivative_diagonal = (F(1), F(1), h)
    jacobian = (
        derivative_diagonal[0]
        * derivative_diagonal[1]
        * derivative_diagonal[2]
    )
    assert jacobian == h

    # Exact physical L2 density after Piola.
    physical_density = h * sum(value**2 for value in pu)
    reference_weighted_density = (
        (u[0] ** 2 + u[1] ** 2) / h + h * u[2] ** 2
    )
    assert physical_density == reference_weighted_density

    # Exact divergence identity div_x(Pu)=div_y(u)/h.
    curved_divergence = (
        dpu[0][0] + dpu[1][1] + dpu[2][2] / h
        + epsilon * pu[0] / h
    )
    flat_divergence = du[0][0] + du[1][1] + du[2][2]
    assert curved_divergence == flat_divergence / h

    # The hand-written three-component commutator equals its definition.
    commutator_from_definition = sub(
        covariant(pa, pu, dpu, epsilon, h),
        piola_value(directional(a, du), h),
    )
    commutator_from_formula = geometric_commutator(
        a, u, du, epsilon, h
    )
    assert commutator_from_definition == commutator_from_formula

    reverse_from_definition = sub(
        covariant(pu, pa, dpa, epsilon, h),
        piola_value(directional(u, da), h),
    )
    reverse_from_formula = geometric_commutator(
        u, a, da, epsilon, h
    )
    assert reverse_from_definition == reverse_from_formula

    # Scalar pressure transport and Piola transport of the flat gradient
    # differ by the stated signed diagonal multiplier.
    curved_pressure_gradient = (
        pressure_gradient[0],
        pressure_gradient[1],
        pressure_gradient[2] / h,
    )
    piola_pressure_gradient = piola_value(pressure_gradient, h)
    assert sub(
        curved_pressure_gradient, piola_pressure_gradient
    ) == pressure_defect(pressure_gradient, h)

    # Full first-jet identity with an independent base mismatch
    # A=Piola(a)+b.
    physical_base = add(pa, b)
    physical_base_jet = tuple(
        tuple(dpa[k][j] + db[k][j] for j in range(3))
        for k in range(3)
    )
    curved_momentum = add(
        add(
            covariant(
                physical_base, pu, dpu, epsilon, h
            ),
            covariant(
                pu, physical_base, physical_base_jet, epsilon, h
            ),
        ),
        curved_pressure_gradient,
    )

    flat_momentum = add(
        add(directional(a, du), directional(u, da)),
        pressure_gradient,
    )
    decomposed = piola_value(flat_momentum, h)
    decomposed = add(
        decomposed, geometric_commutator(a, u, du, epsilon, h)
    )
    decomposed = add(
        decomposed, geometric_commutator(u, a, da, epsilon, h)
    )
    decomposed = add(decomposed, pressure_defect(pressure_gradient, h))
    decomposed = add(
        decomposed, covariant(b, pu, dpu, epsilon, h)
    )
    decomposed = add(
        decomposed, covariant(pu, b, db, epsilon, h)
    )
    assert curved_momentum == decomposed


def nearest_integer(value):
    """Nearest integer for a positive rational away from half-integers."""
    return (value + F(1, 2)).numerator // (value + F(1, 2)).denominator


def check_winding(epsilon, beta, p):
    target = beta * p / epsilon
    m = nearest_integer(target)
    beta_effective = epsilon * m / p
    assert abs(beta_effective - beta) <= epsilon / (2 * p)
    alpha = epsilon * m
    assert alpha <= abs(beta) * p + epsilon / 2
    return m


def main():
    cases = (
        (
            F(2, 17),
            F(3, 5),
            (F(2, 3), F(-5, 7), F(11, 13)),
            (
                (F(3, 8), F(-2, 9), F(5, 11)),
                (F(-7, 10), F(4, 15), F(9, 14)),
                (F(5, 12), F(-8, 13), F(7, 9)),
            ),
            (F(-4, 9), F(7, 12), F(5, 8)),
            (
                (F(11, 10), F(-3, 7), F(2, 5)),
                (F(6, 13), F(5, 9), F(-7, 11)),
                (F(-2, 15), F(8, 17), F(4, 7)),
            ),
            (F(3, 19), F(-4, 21), F(7, 23)),
            (
                (F(2, 25), F(-3, 26), F(5, 27)),
                (F(-7, 29), F(11, 31), F(13, 33)),
                (F(17, 35), F(-19, 37), F(23, 39)),
            ),
            (F(5, 6), F(-7, 8), F(9, 10)),
        ),
        (
            F(3, 29),
            F(-4, 7),
            (F(-5, 8), F(9, 11), F(7, 10)),
            (
                (F(-2, 7), F(5, 12), F(3, 13)),
                (F(4, 9), F(-7, 16), F(11, 18)),
                (F(13, 20), F(2, 21), F(-5, 22)),
            ),
            (F(8, 15), F(-3, 10), F(11, 14)),
            (
                (F(-9, 17), F(7, 19), F(5, 23)),
                (F(3, 25), F(-11, 27), F(13, 28)),
                (F(2, 31), F(17, 32), F(-19, 34)),
            ),
            (F(-5, 41), F(7, 43), F(11, 45)),
            (
                (F(13, 47), F(-17, 49), F(19, 51)),
                (F(23, 53), F(29, 55), F(-31, 57)),
                (F(-37, 59), F(41, 61), F(43, 63)),
            ),
            (F(-4, 5), F(6, 7), F(-8, 9)),
        ),
    )

    for case in cases:
        check_case(case)

    windings = (
        check_winding(F(2, 101), F(7, 5), 37),
        check_winding(F(3, 211), F(11, 8), 61),
        check_winding(F(5, 307), F(13, 9), 83),
    )
    assert all(isinstance(m, int) and m > 0 for m in windings)

    print("thin-torus Jacobian and weighted L2 density: PASS")
    print("exact Piola divergence identity: PASS")
    print("geometric and pressure commutators: PASS")
    print("full first-jet momentum decomposition: PASS")
    print("aspect-correct axial winding bound: PASS")
    print("all aspect-uniform Piola packet checks passed")


if __name__ == "__main__":
    main()
