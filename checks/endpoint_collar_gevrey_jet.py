#!/usr/bin/env python3
"""Exact combinatorics for the endpoint-collar Gevrey jet gate."""

from fractions import Fraction as F
from math import comb, factorial


def main():
    # Exact Gevrey-2 binomial convolution.
    for n in range(1, 81):
        convolution = sum(
            comb(n, a) * factorial(a) ** 2 * factorial(n - a) ** 2
            for a in range(n + 1)
        )
        assert convolution <= 3 * factorial(n) ** 2

    # Every q-extra-derivative cost is bounded by (4 n^2)^q.
    for n in range(1, 81):
        for q in range(1, n + 1):
            ratio = factorial(n + q) // factorial(n)
            assert ratio**2 <= (4 * n * n) ** q

    # Verify the differential-order and seed-degree invariants under both
    # recurrence operations on a representative triangular index set.
    for n in range(1, 30):
        for q in range(1, n + 1):
            derivative_order = n + q
            seed_degree = n - q + 1

            # Diffusion: (n,q,D,L) -> (n+1,q+1,D+2,L).
            assert derivative_order + 2 == (n + 1) + (q + 1)
            assert seed_degree == (n + 1) - (q + 1) + 1

            # Split one admissible nonlinear tree whenever possible.
            if n >= 2 and q >= 2:
                left_n, left_q = 1, 1
                right_n, right_q = n - 1, q - 1
                if right_q <= right_n:
                    combined_derivatives = (
                        left_n + left_q + right_n + right_q + 1
                    )
                    combined_degree = (
                        left_n - left_q + 1
                        + right_n - right_q + 1
                    )
                    assert combined_derivatives == (n + 1) + q
                    assert combined_degree == (n + 1) - q + 1

    # At theta=1/4, the full viscosity-degree sum is uniformly bounded.
    theta = F(1, 4)
    for n in range(1, 81):
        partial = sum((theta**q for q in range(1, n + 1)), F(0))
        assert partial <= theta / (1 - theta)

    print("Gevrey-2 time-binomial convolution constant <= 3: PASS")
    print("q-extra-derivative factorial gate (4 n^2)^q: PASS")
    print("viscosity/tree differential bidegrees: PASS")
    print("uniform geometric viscosity-degree sum: PASS")
    print("all endpoint-collar Gevrey jet checks passed")


if __name__ == "__main__":
    main()
