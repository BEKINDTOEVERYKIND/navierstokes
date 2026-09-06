#!/usr/bin/env python3
"""Dependency-free exact ledgers for C176.

This checker verifies:

* the action--angle--axial flow and its first/second normal-coordinate
  derivatives for a polynomial test frequency;
* the C154 nilpotent physical and covector return identities;
* exact residence of a correlated model tube and O(q) propagation of its
  reciprocal frequency slab through J returns;
* an explicit axis-aligned representative of the Omega(q^3/J) lattice
  capacity, including q^2 availability on q=n^8;
* the enlarged-collar q/J exponents and their stage asymptotics; and
* the fixed-cocycle signed rewrite endpoint identity and b versus b^3
  budget.

It does not certify the action-angle coordinate-existence theorem, an
actual C159 continuity radius, a compact finite-frequency packet, C125,
MCKC, LCE, BAFL, or an unforced Navier--Stokes stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import ceil, log


Vec3 = tuple[F, F, F]
Mat3 = tuple[Vec3, Vec3, Vec3]
Mat2 = tuple[tuple[F, F], tuple[F, F]]

I3: Mat3 = (
    (F(1), F(0), F(0)),
    (F(0), F(1), F(0)),
    (F(0), F(0), F(1)),
)
I2: Mat2 = ((F(1), F(0)), (F(0), F(1)))


def dot(a: Vec3, b: Vec3) -> F:
    return sum((a[i] * b[i] for i in range(3)), F(0))


def outer(a: Vec3, b: Vec3) -> Mat3:
    return tuple(
        tuple(a[i] * b[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def transpose(a: Mat3) -> Mat3:
    return tuple(
        tuple(a[j][i] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def madd(a: Mat3, b: Mat3) -> Mat3:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mscale(c: F, a: Mat3) -> Mat3:
    return tuple(
        tuple(c * a[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mmul(a: Mat3, b: Mat3) -> Mat3:
    bt = transpose(b)
    return tuple(
        tuple(dot(a[i], bt[j]) for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mvec(a: Mat3, x: Vec3) -> Vec3:
    return tuple(dot(a[i], x) for i in range(3))  # type: ignore[return-value]


def mpow(a: Mat3, exponent: int) -> Mat3:
    result = I3
    factor = a
    power = exponent
    while power:
        if power & 1:
            result = mmul(factor, result)
        factor = mmul(factor, factor)
        power //= 2
    return result


def mmul2(a: Mat2, b: Mat2) -> Mat2:
    return (
        (
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ),
        (
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ),
    )


def mvec2(a: Mat2, x: tuple[F, F]) -> tuple[F, F]:
    return (
        a[0][0] * x[0] + a[0][1] * x[1],
        a[1][0] * x[0] + a[1][1] * x[1],
    )


def check_action_angle_axial_flow() -> None:
    # Use omega(E)=2+3E+5E^2 as an exact representative.  The identities
    # checked here are the algebraic identities used for a general smooth
    # omega; the actual A2 omega is supplied by the regular-level theorem.
    omega = lambda energy: F(2) + 3 * energy + 5 * energy**2
    omega_prime = lambda energy: F(3) + 10 * energy
    omega_second = lambda _energy: F(10)
    sqrt2_symbol = F(7, 5)  # symbolic coefficient; its value is irrelevant

    for energy in (F(-1, 20), F(0), F(1, 17)):
        for tau in (F(0), F(3, 7), F(5), F(19)):
            theta0 = F(4, 9)
            z0 = F(-2, 11)
            endpoint = (
                energy,
                theta0 + omega(energy) * tau,
                z0 - sqrt2_symbol * energy * tau,
            )
            # E is exactly conserved and the two velocities are constant.
            assert endpoint[0] == energy
            assert endpoint[1] - theta0 == omega(energy) * tau
            assert endpoint[2] - z0 == -sqrt2_symbol * energy * tau

            jacobian: Mat3 = (
                (F(1), F(0), F(0)),
                (tau * omega_prime(energy), F(1), F(0)),
                (-tau * sqrt2_symbol, F(0), F(1)),
            )
            assert jacobian[1][0] == tau * omega_prime(energy)
            assert jacobian[2][0] == -tau * sqrt2_symbol
            assert omega_second(energy) * tau == 10 * tau

    # A bounded C^2 coordinate conjugacy gives this exact chain-rule
    # majorant: C2*||Dpsi||^2 + C1*||D2psi||.  It is quadratic in 1+tau.
    c1, c2, cpsi, d2psi = F(7), F(11), F(13), F(10)
    for tau in (F(0), F(1), F(7), F(31)):
        dpsi_bound = cpsi * (1 + tau)
        physical_d2_bound = c2 * dpsi_bound**2 + c1 * d2psi * tau
        quadratic_bound = (c2 * cpsi**2 + c1 * d2psi) * (1 + tau) ** 2
        assert physical_d2_bound <= quadratic_bound


def check_nilpotent_returns_and_tubes() -> None:
    # Axis-aligned model of R=u tensor g with u=e3, g=e1.  This is related
    # to the actual orthogonal (g0, e2, u) frame by a fixed rotation.
    u: Vec3 = (F(0), F(0), F(1))
    g: Vec3 = (F(1), F(0), F(0))
    assert dot(u, g) == 0
    rmat = outer(u, g)
    assert mmul(rmat, rmat) == (
        (F(0), F(0), F(0)),
        (F(0), F(0), F(0)),
        (F(0), F(0), F(0)),
    )

    physical_return = madd(I3, rmat)
    covector_return = madd(I3, mscale(-1, transpose(rmat)))

    for returns in (0, 1, 2, 7, 31):
        expected_physical = madd(I3, mscale(F(returns), rmat))
        expected_covector = madd(
            I3, mscale(F(-returns), transpose(rmat))
        )
        assert mpow(physical_return, returns) == expected_physical
        assert mpow(covector_return, returns) == expected_covector

        # The physical shear is volume preserving.  The J*r^3 tube below
        # is a deliberately correlated envelope/cover, not the volume of
        # the incompressible image of an r-ball.
        assert (
            expected_physical[0][0]
            * expected_physical[1][1]
            * expected_physical[2][2]
        ) == 1

        # Physical tube: |y1|,|y2|<=r and |y3|<=Jr.  For m<=J it stays in
        # |third coordinate|<=2Jr.  Use r=1/q exactly.
        j = max(returns, 1)
        q = F(10_000)
        radius = 1 / q
        test_points = (
            (radius, radius, j * radius),
            (-radius, radius, -j * radius),
            (radius, -radius, -j * radius),
        )
        for point in test_points:
            image = mvec(expected_physical, point)
            assert abs(image[0]) <= radius
            assert abs(image[1]) <= radius
            assert abs(image[2]) <= 2 * j * radius

        # The correlated tube has volume J*r^3.  Relative to an r^3 core,
        # fixed-energy squared amplitude is smaller by exactly 1/J.
        child_volume = radius**3
        tube_volume = F(j) * child_volume
        assert tube_volume / child_volume == j
        fixed_energy_amplitude_squared_ratio = child_volume / tube_volume
        assert fixed_energy_amplitude_squared_ratio == F(1, j)

        # Reciprocal slab: |k1|,|k2|<=q and |k3|<=q/J.  The covector
        # shear changes k1 by at most q for m<=J, hence remains O(q).
        slab_points = (
            (q, q, q / j),
            (-q, q, -q / j),
            (F(0), -q, q / j),
        )
        for wavevector in slab_points:
            image = mvec(expected_covector, wavevector)
            assert abs(image[0]) <= 2 * q
            assert abs(image[1]) <= q
            assert abs(image[2]) <= q / j


def check_lattice_capacity() -> None:
    # In an axis-aligned box with half-widths q,q,floor(q/J), the exact
    # count is (2q+1)^2(2floor(q/J)+1)=Omega(q^3/J).  The note's nearest-
    # integer cube argument makes the exponent rotation independent.
    for n in (2, 3, 5, 10):
        q = n**8
        j = max(1, ceil(12 * log(n)))
        normal_half_width = q // j
        assert normal_half_width >= 1
        count = (2 * q + 1) ** 2 * (2 * normal_half_width + 1)
        assert count >= 4 * q**3 // j
        assert count // 2 >= q**2  # enough unordered reality pairs
        assert q > j

        # C170 comparison: a C/q aperture has O(q) projected capacity,
        # whereas fixed aperture has q^2 projected degrees and q/J normal
        # layers.  Only exponents, not hidden constants, are asserted.
        narrow_projected = q
        fixed_projected = q**2
        normal_layers = q // j
        assert narrow_projected < fixed_projected
        assert fixed_projected * normal_layers >= q**3 // (2 * j)


def check_enlarged_collar_ledgers() -> None:
    # Replacing r by R=Jr contributes J^(5/2) to L2 switch/parent/self
    # terms.  Summing J slices gives J^(7/2).  Square half powers for exact
    # rational checks.
    for n in (2, 3, 5, 9):
        q = F(n**8)
        h = F(n**12)
        j = F(max(1, ceil(12 * log(n))))
        radius_ratio = j / q

        # Switch/parent after h and J slices: J*h*R^(5/2).  Squared.
        actual_switch_squared = (j * h) ** 2 * radius_ratio**5
        claimed_switch_squared = j**7 / q**2
        assert actual_switch_squared == claimed_switch_squared

        # Viscosity: J*h*Re^-1*R^(1/2).  Remove Re^-1 and square.
        actual_viscous_squared = (j * h) ** 2 * radius_ratio
        claimed_viscous_squared = j**3 * q**2
        assert actual_viscous_squared == claimed_viscous_squared

        # Pressure moments: J*h*R^5 and J*h*R^3, respectively.
        parent_pressure = j * h * radius_ratio**5
        viscous_pressure = j * h * radius_ratio**3
        # q=n^8 makes q^(-7/2)=n^-28 and q^(-3/2)=n^-12 exact.
        assert parent_pressure == j**6 / F(n**28)
        assert viscous_pressure == j**4 / F(n**12)

    # q=n^8: q^-1 J^(7/2) is o(n^-6), and one n^2 chart remains o(n^-4).
    raw_ratios: list[float] = []
    chart_ratios: list[float] = []
    for n in (20, 100, 1_000, 100_000, 1_000_000):
        j = 12 * log(n)
        raw = n ** (-8) * j ** 3.5
        raw_ratios.append(raw / n ** (-6))
        chart_ratios.append((n**2 * raw) / n ** (-4))
    assert all(
        abs(raw - chart) <= 1e-14 * max(1.0, abs(raw))
        for raw, chart in zip(raw_ratios, chart_ratios)
    )
    assert all(a > b for a, b in zip(raw_ratios, raw_ratios[1:]))
    assert raw_ratios[-1] < F(1, 1000)


def check_fixed_cocycle_rewrite_identity() -> None:
    # Two exact invertible propagators P1 and P0.
    p0: Mat2 = ((F(2), F(1)), (F(1), F(1)))
    p1: Mat2 = ((F(1), F(-1)), (F(2), F(1)))
    total = mmul2(p1, p0)
    z0 = (F(3, 5), F(-2, 7))
    unrewritten = mvec2(total, z0)

    # Two sources may cancel after terminal weighting.
    xi0 = (F(4, 9), F(-1, 3))
    terminal0 = mvec2(p1, xi0)
    xi1 = (-terminal0[0], -terminal0[1])
    signed_response = (
        terminal0[0] + xi1[0],
        terminal0[1] + xi1[1],
    )
    assert signed_response == (F(0), F(0))
    rewritten = (
        unrewritten[0] + signed_response[0],
        unrewritten[1] + signed_response[1],
    )
    assert rewritten == unrewritten

    # If the desired endpoint has relative mismatch delta in a normalized
    # growing scalar, the signed response is delta*b.  BAFL permits b^3.
    for n in (2, 3, 7, 20):
        b = F(1, n**2)
        delta = F(1)
        response = delta * b
        allowance = b**3
        assert response / allowance == n**4
        assert response > allowance
        # Meeting the allowance forces delta<=b^2.
        critical_delta = b**2
        assert critical_delta * b == allowance


def main() -> None:
    check_action_angle_axial_flow()
    print("PASS C176: exact action-angle-axial flow and polynomial jet ledger")
    check_nilpotent_returns_and_tubes()
    print("PASS C176: nilpotent correlated tube/slab residence")
    check_lattice_capacity()
    print("PASS C176: q^3/J lattice capacity contains q^2 carriers")
    check_enlarged_collar_ledgers()
    print("PASS C176: enlarged-collar factorial schedule")
    check_fixed_cocycle_rewrite_identity()
    print("PASS C176: fixed-cocycle signed rewrites have no free reset")
    print("OPEN: finite-frequency localization, C125, MCKC, LCE, BAFL, stage")


if __name__ == "__main__":
    main()
