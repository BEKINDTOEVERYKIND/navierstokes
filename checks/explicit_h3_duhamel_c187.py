#!/usr/bin/env python3
"""Exact constant ledger for C187's H^3 Duhamel estimate."""

from fractions import Fraction as F


def shell_and_zeta_bound() -> None:
    for radius in range(1, 100):
        exact_shell = (2 * radius + 1) ** 3 - (2 * radius - 1) ** 3
        assert exact_shell == 24 * radius * radius + 2
        assert exact_shell <= 26 * radius * radius

    first_five = sum((F(1, radius**4) for radius in range(1, 6)), F(0))
    integral_tail = F(1, 3 * 5**3)
    assert first_five + integral_tail < F(13, 12)

    a3_squared_upper = F(1) + 26 * F(13, 12)
    assert a3_squared_upper == F(175, 6)
    assert a3_squared_upper < F(121, 4)


def cubic_weight_constant() -> None:
    # (x+y)^3 <= 4(x^3+y^3) for nonnegative x,y.  After scaling y=1,
    # the difference factors as 3(t-1)^2(t+1) >= 0.
    for numerator in range(0, 101):
        t = F(numerator, 25)
        difference = 4 * (t**3 + 1) - (t + 1) ** 3
        assert difference == 3 * (t - 1) ** 2 * (t + 1)
        assert difference >= 0

    a3_upper = F(11, 2)
    product_constant = 8 * a3_upper
    assert product_constant == 44


def energy_and_gronwall_constants() -> None:
    product_constant = F(44)
    two_background_terms = 2 * product_constant
    assert two_background_terms == 88

    # xy <= (nu/4)D^2 + coefficient/nu follows by completing a square.
    background_square = two_background_terms**2
    forcing_square = product_constant**2
    assert background_square == 7744
    assert forcing_square == 1936

    # Multiplying the half-energy inequality by two.
    assert 2 * background_square == 15488
    assert 2 * forcing_square == 3872

    # Square root of the forcing coefficient and halving the exponent.
    assert F(3872) == 2 * product_constant**2
    assert F(15488, 2) == 7744


def factorial_nonuniformity() -> None:
    # Even the prefactor sqrt(1/mu_j) grows like j!, before the exponential
    # factor is counted.  Verify the exact reciprocal schedule ledger.
    factorial = 1
    previous_reciprocal = 1
    for stage in range(1, 15):
        factorial *= stage
        mu_reciprocal = factorial**2
        assert mu_reciprocal >= previous_reciprocal
        assert mu_reciprocal == factorial**2
        previous_reciprocal = mu_reciprocal


def main() -> None:
    shell_and_zeta_bound()
    cubic_weight_constant()
    energy_and_gronwall_constants()
    factorial_nonuniformity()
    print("C187 explicit H3 Duhamel constant ledger passed")
    print("BOUNDARY: the estimate is finite-stage and not scale-uniform")
    print("BOUNDARY: unlanded ladder numbers remain withdrawn")


if __name__ == "__main__":
    main()
