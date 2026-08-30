#!/usr/bin/env python3
"""Exact arithmetic diagnostics for C196.

This dependency-free checker verifies representative Fourier curl and
sideband identities, the explicit full-aperture grid inequalities, every
phase-space exponent, the C180/C176 taxes, and the conditional C194 clock
and margins.  The general Fourier/lattice proofs are displayed in C196.

It does not certify a Floquet bundle, a multi-beam FIO estimate, band
retention, viscosity, nonlinear closure, or a singularity.
"""

from fractions import Fraction as F
from math import isqrt


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def curl_and_sideband_identity() -> None:
    # One exact rational instance of the general vector identity in (1.5)
    # and the corrected sideband coefficient in (1.11).
    p = (F(5), F(1), F(0))
    a = (F(0), F(0), F(1))
    s = (F(1), F(2), F(1))
    assert dot(p, a) == 0
    p2 = dot(p, p)
    c = scale(1 / p2, cross(p, a))
    assert scale(-1, cross(p, c)) == a
    corrected = add(a, scale(-1, cross(s, c)))
    assert dot(add(p, s), corrected) == 0
    assert dot(c, c) == 1 / p2
    assert dot(cross(s, c), cross(s, c)) <= dot(s, s) / p2


def explicit_grid_and_concentration() -> None:
    # C196 uses q=m^24, L+=m^8, d=2L++1, K=floor(q/(4d)).
    # Check the smallest allowed m exactly; the displayed proof uses
    # d<=3m^8 and floor(x)>=x/2 to cover every larger m.
    m = 2
    q = m**24
    l_plus = m**8
    l_minus = m**2
    d = 2 * l_plus + 1
    k_count = q // (4 * d)
    assert d <= 3 * m**8
    assert F(q, 4 * d) >= 2
    assert k_count >= m**16 // 24
    modes = k_count**3
    assert 13_824 * modes >= q**2

    # Exponents in (1.18): M~q^2 and envelope mode count L^3.
    carrier_count = F(2)
    stable_l = F(1, 12)
    expanding_l = F(1, 3)
    stable_concentration = (carrier_count + 3 * stable_l) / 2
    expanding_concentration = (carrier_count + 3 * expanding_l) / 2
    assert stable_concentration == F(9, 8)
    assert expanding_concentration == F(3, 2)
    assert expanding_concentration - stable_concentration == F(3, 8)

    # Alpha is a polarization correction, not a curl error.
    assert stable_l - 1 == F(-11, 12)
    assert expanding_l - 1 == F(-2, 3)

    # The coarse rational lower bound e.a_p>=39/40 is weaker than
    # sqrt(1-1/40^2).  Squaring keeps this exact.
    assert F(39, 40) ** 2 < 1 - F(1, 40**2)
    # sqrt(13824)=24*sqrt(24), so the constant in (1.18) is positive.
    assert isqrt(13_824) == 117
    assert 117**2 < 13_824 < 118**2


def support_exponents_and_tax() -> None:
    # Three-coordinate relative q^-1/4 cube.
    three_coordinate_modes = 3 * (1 - F(1, 4))
    assert three_coordinate_modes == F(9, 4)
    assert three_coordinate_modes / 2 == F(9, 8)

    # A projective angular q^-1/4 tube has two narrowed transverse
    # directions but a full q radial direction.
    projective_modes = 1 + 2 * (1 - F(1, 4))
    assert projective_modes == F(5, 2)
    assert projective_modes / 2 == F(5, 4)
    assert F(3, 2) - projective_modes / 2 == F(1, 4)

    # C193's actual one-profile spatial width q^-1/4 has Fourier
    # bandwidth q^1/4 in all three coordinates.
    c193_modes = 3 * F(1, 4)
    assert c193_modes == F(3, 4)
    assert c193_modes / 2 == F(3, 8)

    # C180 Bsharp has N<=C q^3/J^4, so fixed-energy concentration pays J^2.
    assert F(3) / 2 == F(3, 2)
    assert F(-4) / 2 == -2
    # The wider C176 slab pays sqrt(J).
    assert F(-1) / 2 == F(-1, 2)

    multiplicities = [F(3), F(2), F(4), F(1)]
    pairs = sum(multiplicities, F(0))
    outputs = F(len(multiplicities))
    assert sum((m * m for m in multiplicities), F(0)) >= pairs**2 / outputs


def conditional_bridge_window() -> None:
    gain = F(3, 8)
    stable_width = F(1, 12)
    expanding_width = F(1, 3)

    # t<(19g/50)log q+152/25, so e^(6t) costs q^(57g/25).
    clock = F(19, 50) * gain
    energy_exponent = 6 * clock
    assert clock == F(57, 400)
    assert energy_exponent == F(171, 200)
    assert 6 * F(152, 25) == F(912, 25)

    stable_first = 1 - stable_width - energy_exponent
    expanding_first = 1 - expanding_width + gain - energy_exponent
    stable_second = 1 - energy_exponent
    expanding_second = 1 + gain - energy_exponent
    assert stable_first == F(37, 600)
    assert expanding_first == F(14, 75)
    assert stable_second == F(29, 200)
    assert expanding_second == F(13, 25)

    stable_gain_ceiling = F(25, 57) * (1 - stable_width)
    expanding_gain_ceiling = F(25, 32) * (1 - expanding_width)
    assert stable_gain_ceiling == F(275, 684)
    assert expanding_gain_ceiling == F(25, 48)
    assert gain < stable_gain_ceiling < expanding_gain_ceiling
    assert stable_gain_ceiling - gain == F(37, 1368)

    # M~q^2 beams lose q under naive normalized triangle summation.
    assert stable_first - 1 < 0
    assert expanding_first - 1 < 0

    # This is a collision of certified sufficient windows, not a universal
    # no-go: C194's broad-stable majorant gives g<25/57, whereas the
    # principal fixed-band bq entrance cap needs g>=1/2.
    assert F(25, 57) < F(1, 2)
    assert F(3, 2) - gain == F(9, 8)


def main() -> None:
    curl_and_sideband_identity()
    explicit_grid_and_concentration()
    support_exponents_and_tax()
    conditional_bridge_window()
    print("C196 periodic curl/phase-space arithmetic: PASS")
    print("EXACT PROFILE POWERS: q^(9/8) and q^(3/2)")
    print("3-COORDINATE TUBE CEILING: q^(9/8)")
    print("PROJECTIVE TUBE CEILING: q^(5/4)")
    print("C180 Bsharp CEILING: q^(3/2)/J^2")
    print("C194 SINGLE-BEAM MARGINS: 37/600, 14/75")
    print("BOUNDARY: no uniform multi-beam FIO or dynamic retained-band theorem")


if __name__ == "__main__":
    main()
