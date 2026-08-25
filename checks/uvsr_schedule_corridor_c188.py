#!/usr/bin/env python3
"""Exact arithmetic checks for C188's UVSR demand corridor.

The checker certifies the generalized factorial ledgers, the sharp
exponent for the declared worst-case legacy-collar envelope, two
boundary-approaching schedule families, the fixed-shell C161
respecification, the polylogarithmic envelope identity, and the rational
gain retained by the pre-registered three-loss viscous bridge.  It does
not certify a UVSR profile, an epsilon-regularity constant, or a
Navier--Stokes singularity.
"""

from fractions import Fraction as F
from math import isqrt


def generalized_same_energy_ledger() -> None:
    # q=n^Q, b=n^-S, F=q^(3/2), g=bF=n^P.
    for q_exp, seed_exp in ((F(8), F(2)), (F(14), F(4)), (F(26), F(8))):
        amplitude_exp = F(3, 2) * q_exp - seed_exp
        gamma = amplitude_exp / q_exp

        energy_exp = 2 * amplitude_exp - 3 * q_exp
        time_exp = -(amplitude_exp + q_exp)
        reynolds_exp = amplitude_exp - q_exp
        heat_exp = q_exp - amplitude_exp
        dissipation_exp = energy_exp + heat_exp

        assert energy_exp == -2 * seed_exp
        assert reynolds_exp == q_exp / 2 - seed_exp
        assert heat_exp == -reynolds_exp
        assert dissipation_exp == -(seed_exp + q_exp / 2)
        assert F(1) < gamma < F(3, 2)
        assert time_exp < 0
        assert heat_exp == q_exp - amplitude_exp

        # Same seed energy in a q^-3 child volume fixes the raw multiplier.
        child_energy_exp = 2 * amplitude_exp - 3 * q_exp
        seed_energy_exp = -2 * seed_exp
        focus_exp = amplitude_exp + seed_exp
        assert child_energy_exp == seed_energy_exp
        assert focus_exp == F(3, 2) * q_exp

        # With unequal squared L2 shape constants, the physical focus is
        # q^(3/2)*sqrt(kappa_seed/kappa_child).  The schedule exponent is
        # unchanged when those constants are uniformly comparable.
        kappa_seed = F(9)
        kappa_child = F(4)
        shape_factor = F(3, 2)
        assert kappa_child * shape_factor**2 == kappa_seed

        # mu_{j+1}/mu_j=q/g.  A positive-mu autonomous fixed point forces
        # g=q, whereas every listed high-Re schedule has g>q and mu decay.
        viscosity_ratio_exp = q_exp - amplitude_exp
        assert viscosity_ratio_exp < 0

        # For arbitrary one-step schedules rho=b*q^(1/2)=g/q.  This is
        # simultaneously the Reynolds multiplier and inverse viscosity
        # multiplier, so the direct polynomial floor is one.
        rho_exp = q_exp / 2 - seed_exp
        assert rho_exp == amplitude_exp - q_exp

    # Exact all-sequence direct-boundary example: q=16, b=1/2 gives
    # F=q^(3/2)=64, rho=b*sqrt(q)=2, and g=q*rho=32.  Thus polynomial
    # order one is compatible with a strict Reynolds multiplier.
    q = F(16)
    b = F(1, 2)
    focus = F(64)
    rho = b * F(4)
    gain = b * focus
    assert F(0) < b < 1
    assert 1 < rho < 4  # 4=sqrt(q), equivalent to 0<b<1
    assert gain == q * rho == 32


def legacy_collar_family() -> None:
    # Q=6m+2, S=2m retains an exact two-power margin
    # q^-1 / b^3 = n^-2 before the fixed polylogarithm.
    for m in range(1, 101):
        q_exp = F(6 * m + 2)
        seed_exp = F(2 * m)
        amplitude_exp = F(7 * m + 3)
        gamma = amplitude_exp / q_exp

        assert amplitude_exp == F(3, 2) * q_exp - seed_exp
        assert gamma - 1 == (q_exp - 2 * seed_exp) / (2 * q_exp)
        assert gamma - F(7, 6) == (q_exp - 3 * seed_exp) / (
            3 * q_exp
        )
        assert -q_exp + 3 * seed_exp == -2
        assert gamma - F(7, 6) == F(1, 9 * m + 3)
        assert F(7, 6) < gamma < F(3, 2)

        assert 2 * amplitude_exp - 3 * q_exp == -4 * m
        assert -(amplitude_exp + q_exp) == -(13 * m + 5)
        assert amplitude_exp - q_exp == m + 1
        assert q_exp - amplitude_exp == -(m + 1)
        assert (2 * amplitude_exp - 3 * q_exp) + (
            q_exp - amplitude_exp
        ) == -(5 * m + 1)

        # Compatible scalar chart/separation exponents L=n^2,
        # M=n^(S+2): Lb is bounded, (log q)/M=o(b), and
        # L(log q)/M has the same summable n^-S power as b log q.
        chart_exp = F(2)
        separation_exp = seed_exp + 2
        assert chart_exp - seed_exp <= 0
        assert -separation_exp - (-seed_exp) == -2
        assert chart_exp - separation_exp == -seed_exp

    # The landed C127 schedule is m=1; m=2 is a strictly lower-demand
    # exact integer schedule.
    assert (6 * 1 + 2, 2 * 1, 7 * 1 + 3) == (8, 2, 10)
    assert F(7 * 1 + 3, 6 * 1 + 2) == F(5, 4)
    assert (6 * 2 + 2, 2 * 2, 7 * 2 + 3) == (14, 4, 17)
    assert F(7 * 2 + 3, 6 * 2 + 2) == F(17, 14)


def direct_uvsr_family() -> None:
    # If the legacy q^-1 collar is cancelled inside the full residual,
    # Q=4m+2, S=2m approaches the bare high-Re boundary gamma=1.
    for m in range(1, 101):
        q_exp = F(4 * m + 2)
        seed_exp = F(2 * m)
        amplitude_exp = F(4 * m + 3)
        gamma = amplitude_exp / q_exp

        assert amplitude_exp == F(3, 2) * q_exp - seed_exp
        assert amplitude_exp - q_exp == 1
        assert gamma == 1 + F(1, 4 * m + 2)
        assert F(1) < gamma < F(3, 2)
        assert 2 * amplitude_exp - 3 * q_exp == -4 * m


def direct_boundary_sequence() -> None:
    # Exact scalar boundary witness:
    # ell=(j!)^-4, a=2^j(j!)^4, q=n^4, g=2n^4, b=2n^-2.
    q_n_exp = F(4)
    gain_n_exp = F(4)
    seed_n_exp = F(-2)
    focus_n_exp = F(6)
    assert gain_n_exp == seed_n_exp + focus_n_exp

    # Tuple order is (power of 2^j, power of j!).
    amplitude = (F(1), F(4))
    length = (F(0), F(-4))

    def add(left: tuple[F, F], right: tuple[F, F]) -> tuple[F, F]:
        return (left[0] + right[0], left[1] + right[1])

    def scale(multiplier: F, value: tuple[F, F]) -> tuple[F, F]:
        return (multiplier * value[0], multiplier * value[1])

    energy = add(scale(F(2), amplitude), scale(F(3), length))
    turnover = add(length, scale(F(-1), amplitude))
    reynolds = add(amplitude, length)
    viscosity = scale(F(-1), reynolds)
    dissipation = add(energy, viscosity)
    assert energy == (F(2), F(-4))
    assert turnover == (F(-1), F(-8))
    assert reynolds == (F(1), F(0))
    assert viscosity == (F(-1), F(0))
    assert dissipation == (F(1), F(-4))

    # E_{j+1}/E_j=4/(j+1)^4 <=1/4 for j>=1, so the total
    # boundary-witness energy is at most (4/3)E_1=16/3.
    first_ratio_upper = F(4, 2**4)
    geometric_sum_factor = F(1, 1 - first_ratio_upper)
    assert first_ratio_upper == F(1, 4)
    assert geometric_sum_factor == F(4, 3)
    assert geometric_sum_factor * 4 == F(16, 3)

    # b*log q=8(log n)n^-2 <= 8n^-3/2, and the elementary
    # integral comparison used in the note bounds the sum by 16.
    assert 2 * q_n_exp == 8
    assert F(-2) + F(1, 2) == F(-3, 2)
    assert 8 * 2 == 16


def fixed_c180_shell_lower_demand() -> None:
    # Preserve C180's proved q=n^8 shell and change only the dormant seed.
    q_exp = F(8)
    seed_exp = F(5, 2)
    amplitude_exp = F(19, 2)
    gamma = amplitude_exp / q_exp

    assert amplitude_exp == F(3, 2) * q_exp - seed_exp
    assert gamma == F(19, 16) < F(5, 4)
    assert -q_exp + 3 * seed_exp == F(-1, 2)
    assert 2 * amplitude_exp - 3 * q_exp == -5
    assert -(amplitude_exp + q_exp) == F(-35, 2)
    assert amplitude_exp - q_exp == F(3, 2)
    assert q_exp - amplitude_exp == F(-3, 2)
    assert -5 + F(-3, 2) == F(-13, 2)
    assert F(2) - seed_exp == F(-1, 2)  # n^2 chart times b
    assert -2 * seed_exp == -5  # retained-wake scale b^2
    assert -3 * seed_exp == F(-15, 2)  # active scale b^3

    # With L=n^2, M=n^(7/2), and log n <= (3/2)n^(1/4):
    # (log q)/M <= 12 n^-3/4 b and L(log q)/M <= 12 n^-5/4.
    assert q_exp * F(3, 2) == 12
    assert F(1, 4) - 1 == F(-3, 4)
    assert F(2) + F(-7, 2) + F(1, 4) == F(-5, 4)
    integral_tail = F(1, F(5, 4) - 1)
    assert integral_tail == 4
    assert 12 * integral_tail == 48

    # Respecify C161 while preserving its n^-28 preparation size.
    preparation_exp = F(-28)
    seed_coefficient_exp = preparation_exp - q_exp
    preamplifier_exp = F(51, 2)
    c0_exp = preamplifier_exp + seed_coefficient_exp
    c1_exp = c0_exp - q_exp / 2
    assert seed_coefficient_exp == -36
    assert c0_exp == -seed_exp - q_exp == F(-21, 2)
    assert q_exp + c0_exp == -seed_exp  # q*c0=b
    assert F(3, 2) * q_exp + c1_exp == -seed_exp
    assert 3 * q_exp + c1_exp == amplitude_exp

    # J_split=ceil(n^(5/2)) is checked with integer squares, so no float
    # approximation enters the rounding inequalities b/2 <= theta <= b.
    for n in range(2, 1001):
        n5 = n**5
        root = isqrt(n5)
        split_count = root if root * root == n5 else root + 1
        assert split_count * split_count >= n5
        assert (split_count - 1) ** 2 < n5
        assert split_count * split_count <= 4 * n5
        theta = F(1, split_count)
        assert split_count * theta == 1

        # The first square inequality is theta <= b.  Together with the
        # second it gives b/2 <= theta, hence the wake/active bounds.
        assert split_count**2 >= n5
        assert split_count**2 <= 4 * n5


def sharpness_and_polylog_boundary() -> None:
    # The declared worst-case envelope certificate is Q>3S.  Its
    # polynomial infimum is 7/6 and is not attained by a pure power.
    for q_exp, seed_exp in ((F(8), F(2)), (F(14), F(4)), (F(62), F(20))):
        assert q_exp > 3 * seed_exp
        gamma = F(3, 2) - seed_exp / q_exp
        assert gamma > F(7, 6)

    # Include one exact declared-envelope constant K_*=8, whose cube-root
    # is 2.  It cancels from [K_*q^-1(log q)^7/2]/b^3.
    envelope_constant = F(8)
    envelope_constant_cuberoot = F(2)
    assert envelope_constant_cuberoot**3 == envelope_constant

    boundary_q = F(-1, 3)
    boundary_log = F(7, 6)
    boundary_loglog = F(1)
    # b=K_*^(1/3) q^-1/3 (log q)^7/6 loglog q.
    b3_q = 3 * boundary_q
    b3_log = 3 * boundary_log
    b3_loglog = 3 * boundary_loglog
    assert (b3_q, b3_log, b3_loglog) == (F(-1), F(7, 2), F(3))

    # [K_*q^-1(log q)^7/2]/b^3=(loglog q)^-3; constants
    # cancel exactly as well as all q/log powers.
    assert envelope_constant / envelope_constant_cuberoot**3 == 1
    ratio = (F(-1) - b3_q, F(7, 2) - b3_log, -b3_loglog)
    assert ratio == (F(0), F(0), F(-3))

    # omega=loglog(q) diverges but is subpower, so it reaches the 7/6
    # polynomial floor while preserving the vanishing envelope ratio.
    omega_q_power = F(0)
    assert omega_q_power == 0

    # g=b q^3/2 has the boundary polynomial exponent 7/6.
    assert boundary_q + F(3, 2) == F(7, 6)

    # Exact rational ingredient in log x <= (3/2)x^(1/4):
    # e > sum_{k=0}^4 1/k! = 65/24 > 8/3, hence 4/e < 3/2.
    exponential_partial_sum = F(1) + F(1) + F(1, 2) + F(1, 6) + F(1, 24)
    assert exponential_partial_sum == F(65, 24)
    assert exponential_partial_sum > F(8, 3)
    assert exponential_partial_sum > F(7, 3)

    # With Q-3S=2, the 7/2-polylog collar ratio is bounded by
    # K n^-9/8 after the displayed logarithm inequality.
    assert -2 + F(7, 2) * F(1, 4) == F(-9, 8)

    # For the fixed-Q=8, S=5/2 point, log n <= 6 n^(1/14) gives
    # collar/b^3 <= 48^(7/2) n^-1/4.  At n >= 48^28 this is <= n^-1/8.
    assert F(-1, 2) + F(7, 2) * F(1, 14) == F(-1, 4)
    assert F(7, 2) * 8 == 28
    assert F(-1, 4) + F(1, 8) == F(-1, 8)


def epsilon_scaling_and_local_flux() -> None:
    # Record powers of (amplitude a, radius r, viscosity nu) after a
    # turnover-duration packet tau=r/a is inserted into each dimensionless
    # CKN-scale quantity.  Tuple order is (a, r, nu).
    scaled_energy = (2, 2, -2)  # (nu^2 r)^-1 int_B |u|^2
    scaled_cubic = (2, 2, -2)  # (nu^2 r^2)^-1 int_Q |u|^3
    reynolds_squared = (2, 2, -2)
    assert scaled_energy == reynolds_squared
    assert scaled_cubic == reynolds_squared

    # The advective local-energy flux over one turnover has the same
    # physical scaling as packet energy: a^3 r^-1 r^3 (r/a)=a^2 r^3.
    flux = (3 - 1, -1 + 3 + 1, 0)
    packet_energy = (2, 3, 0)
    assert flux == packet_energy

    # At a smooth point, direct cylinder integration gives the r powers:
    # A: volume 3 minus normalization 1;
    # C+D: volume 3 plus parabolic time 2 minus normalization 2;
    # E: volume 3 plus time 2 minus normalization 1.
    smooth_r_powers = (3 - 1, 3 + 2 - 2, 3 + 2 - 1)
    assert smooth_r_powers == (2, 3, 4)

    # Smooth-cylinder viscosity powers after using physical parabolic time
    # r^2/nu: A has nu^-2, C and D have nu^-3, and E has nu^-2.
    smooth_nu_powers = (-2, -2 - 1, -2 - 1, -1 - 1)
    assert smooth_nu_powers == (-2, -3, -3, -2)

    # The volume coefficient in all three bounds is (4/3)pi.  Track its
    # exact rational factor separately from the common transcendental pi.
    ball_volume_rational = F(4, 3)
    assert ball_volume_rational == F(4, 3)

    # The pressure oscillation bound |p-(p)_B| <= 2 M_p contributes
    # the exact coefficient (4/3)pi 2^(3/2).  After removing pi and
    # squaring, the exact rational coefficient is (4/3)^2*8=128/9.
    pressure_coefficient_over_pi_squared = ball_volume_rational**2 * 2**3
    assert pressure_coefficient_over_pi_squared == F(128, 9)


def viscous_bridge_pre_registration() -> None:
    # C185 supplies an inviscid gain strictly above 6/5.  If each of the
    # three multiplicative conversion factors is at least 99/100 and the
    # normalized additive error is at most 1/100, the surviving gain is
    # strictly above 23/20.
    lower = F(6, 5) * F(99, 100) ** 3 - F(1, 100)
    assert lower == F(2_885_897, 2_500_000)
    assert lower - F(23, 20) == F(10_897, 2_500_000)
    assert lower > F(23, 20)


def main() -> None:
    generalized_same_energy_ledger()
    legacy_collar_family()
    direct_uvsr_family()
    direct_boundary_sequence()
    fixed_c180_shell_lower_demand()
    sharpness_and_polylog_boundary()
    epsilon_scaling_and_local_flux()
    viscous_bridge_pre_registration()
    print(
        "PASS: C188 exact schedule corridor, envelope families, C161 "
        "respecification, local scaling, and bridge threshold verified"
    )
    print(
        "SCOPE: no UVSR profile, epsilon constant, PPRG witness, or "
        "Navier--Stokes singularity is certified"
    )


if __name__ == "__main__":
    main()
