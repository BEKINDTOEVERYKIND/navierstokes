#!/usr/bin/env python3
"""Exact/arithmetic checks for recovered claims C126--C131.

This checker validates only the scalar schedule.  It deliberately does not
claim the conditional localization, pressure-tail, or endpoint-map bounds.
"""

from fractions import Fraction
import math


def factorial_scales(j: int) -> tuple[int, Fraction, int]:
    fact = math.factorial(j)
    ell = Fraction(1, fact**8)
    amp = fact**10
    kappa = fact**8
    return amp, ell, kappa


def exact_ratio_checks() -> None:
    for j in range(1, 12):
        n = j + 1
        amp, ell, kappa = factorial_scales(j)
        amp_next, ell_next, kappa_next = factorial_scales(j + 1)

        q = n**8
        gain = n**10
        seed = Fraction(1, n**2)
        focus = n**12

        assert Fraction(kappa_next, kappa) == q
        assert Fraction(amp_next, amp) == gain
        assert ell_next / ell == Fraction(1, q)
        assert seed * focus == gain

        energy = amp**2 * ell**3
        energy_next = amp_next**2 * ell_next**3
        time = ell / amp
        time_next = ell_next / amp_next
        heat_without_nu = Fraction(kappa, amp)
        heat_next_without_nu = Fraction(kappa_next, amp_next)
        diss_without_nu = energy * heat_without_nu
        diss_next_without_nu = energy_next * heat_next_without_nu

        assert energy == Fraction(1, math.factorial(j) ** 4)
        assert time == Fraction(1, math.factorial(j) ** 18)
        assert heat_without_nu == Fraction(1, math.factorial(j) ** 2)
        assert diss_without_nu == Fraction(1, math.factorial(j) ** 6)
        assert energy_next / energy == Fraction(1, n**4)
        assert time_next / time == Fraction(1, n**18)
        assert heat_next_without_nu / heat_without_nu == Fraction(1, n**2)
        assert diss_next_without_nu / diss_without_nu == Fraction(1, n**6)


def summability_checks() -> None:
    # Numerical partial sums and vanishing ratio bounds complement the exact
    # ratio identities.  Mathematical convergence follows by the ratio and
    # integral tests recorded in the note.
    energy = 0.0
    stage_time = 0.0
    heat = 0.0
    diss = 0.0
    logarithmic_stage_time = 0.0
    logarithmic_heat = 0.0
    logarithmic_diss = 0.0
    pressure = 0.0
    wake = 0.0
    localization = 0.0
    for j in range(1, 200):
        n = j + 1
        log_fact = math.lgamma(j + 1)
        energy += math.exp(-4.0 * log_fact)
        stage_time += math.exp(-18.0 * log_fact)
        heat += math.exp(-2.0 * log_fact)
        diss += math.exp(-6.0 * log_fact)
        log_window = 8.0 * math.log(n)
        logarithmic_stage_time += log_window * math.exp(-18.0 * log_fact)
        logarithmic_heat += log_window * math.exp(-2.0 * log_fact)
        logarithmic_diss += log_window * math.exp(-6.0 * log_fact)
        pressure += 8.0 * math.log(n) / n**1.5
        wake += 1.0 / n**4
        localization += math.exp(-(n**7))

    for value in (
        energy,
        stage_time,
        heat,
        diss,
        logarithmic_stage_time,
        logarithmic_heat,
        logarithmic_diss,
        pressure,
        wake,
        localization,
    ):
        assert math.isfinite(value)

    # Integral-test tail: integral_N^infty log(x)x^(-3/2) dx
    # = (2 log N + 4)/sqrt(N).
    for n in (4, 16, 100, 10_000):
        exact_integral_tail = (2.0 * math.log(n) + 4.0) / math.sqrt(n)
        assert exact_integral_tail > 0.0
    assert (2.0 * math.log(10_000) + 4.0) / 100.0 < 0.23

    # The chart-weighted pressure exponent is exactly
    # n^2 * log(n^8) / n^(7/2) = 8 log(n)/n^(3/2).
    for n in (2, 3, 10, 101):
        lhs = n**2 * math.log(n**8) / n**3.5
        rhs = 8.0 * math.log(n) / n**1.5
        assert math.isclose(lhs, rhs, rel_tol=2e-15, abs_tol=1e-15)


def polynomial_absorption_checks() -> None:
    # Ratio from (5.1) tends to zero for every fixed C and alpha>0.
    for power in (0, 2, 12, 40):
        alpha = 0.5
        ratios = [((1.0 + 1.0 / j) ** power) / ((j + 1) ** alpha)
                  for j in (10_000, 100_000, 1_000_000)]
        assert ratios[2] < ratios[1] < ratios[0]
        assert ratios[2] < 0.01

    # Worst allowed finite-dimensional chart loss times dormant seed.
    for n in range(2, 100):
        chart_loss = n**2
        dormant_seed = Fraction(1, n**2)
        assert chart_loss * dormant_seed == 1


def main() -> None:
    exact_ratio_checks()
    summability_checks()
    polynomial_absorption_checks()
    print("PASS: C126--C131 factorial one-cell scalar ledger")
    print("CONDITIONAL/OPEN: localized pressure bound and shrinking exit map")
    print("OPEN: full Navier--Stokes stage and backward-weighted mode leakage")


if __name__ == "__main__":
    main()
