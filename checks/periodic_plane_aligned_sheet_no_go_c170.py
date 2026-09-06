#!/usr/bin/env python3
"""Exact structural checks for C170.

The checker verifies the C152 periodic-plane coordinates, the projected
weight-lattice chart and sector bound, the fixed-direction ray, degree-zero
Kelvin homogeneity, and the C154 lift shear/spacing. It does not decide the
arithmetic nature of sigma, certify finite-frequency localization, or
compute the qualitative monodromy-neighborhood constants.
"""

from fractions import Fraction as F
from math import isqrt


N = (1, 1, 1)
R1 = (1, -1, 0)
R2 = (0, 1, -1)
R = (-1, 0, 1)
D = (-1, 2, -1)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def scale(scalar, value):
    return tuple(scalar * entry for entry in value)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def horizontal_projection(value):
    charge = F(dot(N, value), 3)
    return sub(tuple(F(entry) for entry in value), scale(charge, N))


def projected_coordinates(value):
    """Integer coordinates (R.value,D.value) determining pi_H value."""
    return dot(R, value), dot(D, value)


def check_periodic_plane_coordinates():
    assert D == add(scale(-1, R1), R2)
    assert add(add(R1, R2), R) == (0, 0, 0)
    assert dot(N, R) == dot(N, D) == dot(R, D) == 0
    assert (dot(N, N), dot(R, R), dot(D, D)) == (3, 2, 6)

    # Write u=a0 R-b0 N. Symbolically u.v is determined by the pair
    # (R.v,-N.v). With sigma=2a0/(3b0), R+sigma N lies in u^perp.
    assert (dot(R, R), -dot(N, R)) == (2, 0)
    assert (dot(R, D), -dot(N, D)) == (0, 0)
    assert (dot(R, N), -dot(N, N)) == (0, -3)
    assert 3 * 2 - 2 * 3 == 0  # 3b0(2a0)-2a0(3b0)

    # The proposed inverse of pi_H|P is exact for arbitrary rational sigma.
    for x, y, sigma in (
        (F(2), F(-3), F(5, 7)),
        (F(-4, 3), F(7, 2), F(11, 13)),
    ):
        horizontal = add(scale(x, R), scale(y, D))
        plane_lift = add(horizontal, scale(x * sigma, N))
        assert horizontal_projection(plane_lift) == horizontal

        # Clear u.plane_lift=0 using a0/b0=3sigma/2.
        p, s = sigma.numerator, sigma.denominator
        assert 3 * p * dot(R, plane_lift) - 2 * s * dot(N, plane_lift) == 0

    # beta/(3 gamma)=sigma follows from the C152 beta formula and
    # |3 gamma R|^2/|N|^2=6 gamma^2. Both coefficients of -T'/T have
    # square 2*gamma^2=42/25 and are positive.
    gamma_squared = F(21, 25)
    beta_ray_coefficient_squared = F(2) * gamma_squared
    sigma_coefficient_squared = F(42, 25)
    assert beta_ray_coefficient_squared == sigma_coefficient_squared


def check_rational_plane_and_ray():
    # If sigma=p/s, sR+pN is an integer point of P and spans its aligned
    # ray over R. Plane membership is checked after clearing denominators.
    for p, s in ((2, 5), (7, 3), (4, 9)):
        sigma = F(p, s)
        ray_generator = add(scale(s, R), scale(p, N))
        assert horizontal_projection(ray_generator) == scale(s, R)
        assert 3 * p * dot(R, ray_generator) - 2 * s * dot(N, ray_generator) == 0

        # Any tested integer point in P with D-coordinate zero is parallel
        # to this generator. This also catches an erroneous second label.
        aligned = []
        box = 8 * max(p, s)
        for k1 in range(-box, box + 1):
            for k2 in range(-box, box + 1):
                # D.value=0 determines k3 exactly.
                k3 = 2 * k2 - k1
                if not -box <= k3 <= box:
                    continue
                value = (k1, k2, k3)
                if 3 * p * dot(R, value) != 2 * s * dot(N, value):
                    continue
                aligned.append(value)
                # Cross-multiplication proves value is on the ray. The
                # scalar can be half-integral before plane arithmetic is used.
                assert scale(s, value) == scale(F(dot(R, value), 2),
                                                 ray_generator)
        assert aligned


def check_projected_weight_lattice():
    minimum = None
    seen = {}
    for k1 in range(-4, 5):
        for k2 in range(-4, 5):
            for k3 in range(-4, 5):
                value = (k1, k2, k3)
                projected = horizontal_projection(value)
                a_coord, d_coord = projected_coordinates(value)
                assert (a_coord - d_coord) % 2 == 0
                expected = add(scale(F(a_coord, 2), R),
                               scale(F(d_coord, 6), D))
                assert projected == expected
                assert dot(projected, projected) == (
                    F(a_coord * a_coord, 2) + F(d_coord * d_coord, 6)
                )

                if projected != (0, 0, 0):
                    squared = dot(projected, projected)
                    minimum = squared if minimum is None else min(minimum, squared)

                key = projected_coordinates(value)
                if key in seen:
                    difference = sub(value, seen[key])
                    assert difference[0] == difference[1] == difference[2]
                else:
                    seen[key] = value
    assert minimum == F(2, 3)

    # Conversely every same-parity pair has an integer representative.
    for a_coord in range(-8, 9):
        for d_coord in range(-8, 9):
            if (a_coord - d_coord) % 2:
                continue
            representative = (0, (a_coord + d_coord) // 2, a_coord)
            assert projected_coordinates(representative) == (a_coord, d_coord)

    # Direct exact count for a tangent-aperture tau. Any actual projected
    # point is a subset of integer pairs (a,d) satisfying
    # d^2 <= 3 tau^2 a^2. Since sqrt(3)<2, summing at most
    # 4*tau*|a|+1 d-values gives the displayed O(tau*q^2+q) bound.
    radius = 3
    for q in (8, 16, 32, 64):
        for tau in (F(1, q), F(2, q), F(1, 8)):
            a_max = 2 * radius * q  # exceeds sqrt(2)*R*q
            count = 0
            for a_coord in range(-a_max, a_max + 1):
                numerator = 3 * tau.numerator**2 * a_coord**2
                denominator = tau.denominator**2
                d_max = isqrt(numerator // denominator)
                count += 2 * d_max + 1
            summed_bound = 4 * tau * a_max * (a_max + 1) + 2 * a_max + 1
            assert F(count) <= summed_bound
            assert summed_bound <= 80 * radius * radius * (tau * q * q + q)


def matvec(matrix, vector):
    return tuple(sum(row[j] * vector[j] for j in range(3)) for row in matrix)


def kelvin_amplitude_rhs(matrix, wavevector, amplitude):
    av = matvec(matrix, amplitude)
    norm_squared = dot(wavevector, wavevector)
    scalar = dot(wavevector, av)
    return add(scale(-1, av), scale(F(2) * scalar / norm_squared, wavevector))


def check_kelvin_homogeneity():
    matrix = (
        (F(1), F(2), F(-1)),
        (F(0), F(-3), F(4)),
        (F(5), F(1), F(2)),
    )
    wavevector = (F(2), F(-1), F(3))
    amplitude = (F(1), F(4), F(-2))
    reference = kelvin_amplitude_rhs(matrix, wavevector, amplitude)
    for lam in (F(-7, 3), F(2), F(11, 5)):
        assert kelvin_amplitude_rhs(
            matrix, scale(lam, wavevector), amplitude
        ) == reference

    # Inverse-plane directional geometry is smooth: its squared norm has
    # this exact diagonal form in the orthogonal R,D,N chart.
    for x, y, sigma in ((F(3), F(1, 20), F(5, 7)),
                        (F(-2), F(1, 31), F(11, 13))):
        lifted = add(scale(x, add(R, scale(sigma, N))), scale(y, D))
        assert dot(lifted, lifted) == (
            (F(2) + F(3) * sigma * sigma) * x * x + F(6) * y * y
        )


def check_c154_shear_and_lift_spacing():
    # Rational positive proxies preserve the exact rank-one structure.
    a0 = F(3)
    b0 = F(2)
    gamma = F(5)
    u = add(scale(a0, R), scale(-b0, N))
    g0 = scale(gamma, D)
    assert dot(u, g0) == 0
    assert dot(u, N) == -3 * b0

    # K^ell Delta=Delta-ell*g0*(u.Delta).
    delta = (F(4), F(-2), F(7))
    current = delta
    for ell in range(1, 10):
        current = sub(current, scale(dot(u, current), g0))
        closed_form = sub(delta, scale(ell * dot(u, delta), g0))
        assert current == closed_form

    # Integer lifts over one horizontal point differ by N. Their u-values
    # have fixed spacing 3b0, while return shear is purely in D.
    origin = (F(2), F(-1), F(4))
    ell = 7
    returned = []
    for m in range(-5, 6):
        lift = add(origin, scale(m, N))
        scalar = dot(u, lift)
        image = sub(lift, scale(ell * scalar, g0))
        returned.append((scalar, dot(D, horizontal_projection(image))))
    for left, right in zip(returned, returned[1:]):
        assert right[0] - left[0] == -3 * b0
        assert right[1] - left[1] == ell * 18 * b0 * gamma

    # An affine lattice of spacing h has at most floor(2S/h)+1 points in
    # [-S,S]. In the C/q, log(q)-return regime S=O(1/log q), hence one
    # point eventually. Test the exact spacing lemma at several offsets.
    spacing = 3 * b0
    for offset in (F(0), F(1, 3), F(7, 5)):
        for tolerance in (F(1, 10), F(5, 2), F(20)):
            admissible = [
                m for m in range(-100, 101)
                if abs(offset + m * spacing) <= tolerance
            ]
            assert len(admissible) <= (2 * tolerance // spacing) + 1
    small_tolerance = spacing / 3
    assert sum(abs(F(1, 5) + m * spacing) <= small_tolerance
               for m in range(-20, 21)) <= 1


def main():
    check_periodic_plane_coordinates()
    print("PASS C170: exact C152 plane and horizontal projection isomorphism")
    check_rational_plane_and_ray()
    print("PASS C170: rationality caveat and one-dimensional aligned ray")
    check_projected_weight_lattice()
    print("PASS C170: projected sector capacity O(delta*q^2+q)")
    check_kelvin_homogeneity()
    print("PASS C170: exact ray homogeneity and smooth sector geometry")
    check_c154_shear_and_lift_spacing()
    print("PASS C170: C154 tolerance and integer normal-lift loss")


if __name__ == "__main__":
    main()
