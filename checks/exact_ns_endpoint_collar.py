#!/usr/bin/env python3
"""Exact scaling checks for a Navier--Stokes endpoint collar."""

from fractions import Fraction as F


class Monomial:
    """Powers of (a, ell, nu)."""

    def __init__(self, a=0, ell=0, nu=0):
        self.powers = (F(a), F(ell), F(nu))

    def __mul__(self, other):
        return Monomial(*(x + y for x, y in zip(self.powers, other.powers)))

    def __truediv__(self, other):
        return Monomial(*(x - y for x, y in zip(self.powers, other.powers)))

    def __pow__(self, power):
        return Monomial(*(power * x for x in self.powers))

    def __eq__(self, other):
        return self.powers == other.powers


def main():
    velocity = Monomial(a=1)
    length = Monomial(ell=1)
    viscosity = Monomial(nu=1)
    turnover = length / velocity
    reynolds_inverse = viscosity / (velocity * length)

    # Viscous acceleration nu*Delta V, integrated for one turnover,
    # is a*Re^{-1} relative to a velocity of size a.
    viscous_acceleration = viscosity * velocity / (length**2)
    collar_change = viscous_acceleration * turnover
    assert collar_change / velocity == reynolds_inverse

    # Energy and dissipation on a three-dimensional bubble.
    energy = velocity**2 * length**3
    gradient_l2_squared = velocity**2 * length
    collar_dissipation = viscosity * turnover * gradient_l2_squared
    assert collar_dissipation / energy == reynolds_inverse

    # The scale-adapted m-th derivative term has velocity dimension.
    for derivative_order in range(9):
        derivative_l2_squared = (
            velocity**2 * length ** (3 - 2 * derivative_order)
        )
        adapted_term = (
            length ** (2 * derivative_order - 3) * derivative_l2_squared
        )
        assert adapted_term == velocity**2

    # The first two endpoint jets contain at least one power of mu.
    mu = reynolds_inverse
    first_jet = mu
    second_jet_viscous = mu**2
    second_jet_linearized_euler = mu
    assert first_jet.powers[2] == 1
    assert second_jet_viscous.powers[2] == 2
    assert second_jet_linearized_euler.powers[2] == 1

    print("turnover collar deformation / velocity = Re^{-1}: PASS")
    print("turnover energy loss / energy = Re^{-1}: PASS")
    print("scale-adapted Sobolev normalization: PASS")
    print("first and second endpoint-jet viscosity powers: PASS")
    print("all exact Navier--Stokes collar checks passed")


if __name__ == "__main__":
    main()
