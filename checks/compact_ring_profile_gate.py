#!/usr/bin/env python3
"""Exact interval ledger for compact-ring profile compatibility gates."""

from fractions import Fraction as F

from ao_batchelor_full_edge_matched import exp_interval
from ao_batchelor_global_bas import I


def main():
    # C65 full-edge-matched Batchelor ring.
    x = I(F(59671214, 10**8), F(59671216, 10**8))
    q_swirl = I(F(8278572, 10**7), F(8278581, 10**7))
    exponential = 1 / exp_interval(x)
    one_minus_exponential = 1 - exponential

    # For the Cao--Zhan generator test, positivity of HH_Psi forces
    # c<=E.  Then A-A' >= xE/(1-E), because the exact difference is
    # (E-c)(x-(1-E))/(1-E)^2 and both factors are positive here.
    x_minus_one_minus_exponential = x - one_minus_exponential
    assert x_minus_one_minus_exponential.lo > F(1473, 10**4)
    generator_gap = (
        x * exponential / one_minus_exponential - q_swirl * q_swirl
    )
    assert generator_gap.subset(F(4579, 10**5), F(4580, 10**5))

    # Exact rational substitutions behind the straight-limit hodograph law
    # s(h^2)' + 2h^2 = 1 and its solution h^2=1/2+C/s^2.
    for s, constant in (
        (F(3, 2), F(7, 5)),
        (F(11, 4), F(-2, 3)),
        (F(5), F(0)),
    ):
        h_squared = F(1, 2) + constant / (s * s)
        h_squared_prime = -2 * constant / (s**3)
        assert s * h_squared_prime + 2 * h_squared == 1

    # Smooth core-through-flow counterexample to the false implication
    # "smooth core => C=0": V=s and W^2=C+s^2/2 obey the law for C>0.
    positive_constant = F(7, 5)
    for s in (F(1, 10), F(3, 4), F(2)):
        velocity_squared = s * s
        axial_squared = positive_constant + s * s / 2
        assert axial_squared / velocity_squared == (
            F(1, 2) + positive_constant / (s * s)
        )

    # Independent modulation of poloidal and toroidal velocities would
    # leave (b^2-a^2)F^2/(2r^2) in the energy.  Two radii on one toroidal
    # level force b^2=a^2 whenever F is nonzero.
    poloidal_scale = F(5, 4)
    toroidal_scale = F(7, 6)
    toroidal_flux = F(9, 5)
    r_one, r_two = F(4, 3), F(7, 4)
    difference = (
        (toroidal_scale**2 - poloidal_scale**2)
        * toroidal_flux**2
        / 2
        * (1 / r_one**2 - 1 / r_two**2)
    )
    assert difference != 0

    print("Gavrilov--Constantin--La profile-law substitutions: PASS")
    print("smooth C>0 core-through-flow escape witness: PASS")
    print("same-foliation modulation nonconstancy witness: PASS")
    print("0.04579 < Cao--Zhan Batchelor generator gap < 0.04580")
    print("all compact-ring profile gate checks passed")


if __name__ == "__main__":
    main()
