#!/usr/bin/env python3
"""Arithmetic checks for the conditional one-cell BAFL reduction.

This script does not prove the PDE estimate (BAFL).
"""

from fractions import Fraction


# Gaussian rationals are represented by (real, imaginary) pairs.  Keeping
# this tiny implementation here makes the A2 Fourier calculation exact and
# dependency-free.
def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gmul_real(c, a):
    return (c * a[0], c * a[1])


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


ZERO = (Fraction(0), Fraction(0))


def vadd(a, b):
    return tuple(gadd(x, y) for x, y in zip(a, b))


def vscale(c, a):
    return tuple(gmul_real(c, x) for x in a)


def dot_real(a, q):
    out = ZERO
    for x, y in zip(a, q):
        out = gadd(out, gmul_real(Fraction(y), x))
    return out


def project(q, v):
    q2 = sum(Fraction(x * x) for x in q)
    coefficient = dot_real(v, q)
    return tuple(gadd(x, gmul_real(-Fraction(y, 1) / q2, coefficient))
                 for x, y in zip(v, q))


def polarization(x, y, n, t):
    return tuple(gadd(gmul(x, (Fraction(ni), Fraction(0))),
                       gmul(y, (Fraction(ti), Fraction(0))))
                 for ni, ti in zip(n, t))


def interaction(k, a, ell, b):
    first = dot_real(a, ell)
    second = dot_real(b, k)
    raw = tuple(gadd(gmul(first, bi), gmul(second, ai))
                for ai, bi in zip(a, b))
    q = tuple(x + y for x, y in zip(k, ell))
    return project(q, raw)


def check_a2_reality_leak() -> None:
    roots = ((1, -1, 0), (0, 1, -1), (-1, 0, 1))
    n = (1, 1, 1)
    tangents = ((1, 1, -2), (-2, 1, 1), (1, -2, 1))
    xs = ((Fraction(2), Fraction(1)),
          (Fraction(-1), Fraction(3)),
          (Fraction(4), Fraction(-2)))
    ys = ((Fraction(1), Fraction(-1)),
          (Fraction(3), Fraction(2)),
          (Fraction(-2), Fraction(1)))
    modes = [polarization(x, y, n, t)
             for x, y, t in zip(xs, ys, tangents)]
    n_gaussian = tuple((Fraction(z), Fraction(0)) for z in n)

    for i, j in ((0, 1), (1, 2), (2, 0)):
        # Intended sum: 3 (y_i x_j - y_j x_i) n.
        sum_out = interaction(roots[i], modes[i], roots[j], modes[j])
        scalar_sum = gadd(gmul(ys[i], xs[j]),
                          gmul_real(-1, gmul(ys[j], xs[i])))
        expected_sum = tuple(gmul(gmul_real(3, scalar_sum), z)
                             for z in n_gaussian)
        assert sum_out == expected_sum

        # Reality companion: coefficient at -k_j is conjugate(a_j).
        minus_root = tuple(-z for z in roots[j])
        conjugate_mode = tuple(gconj(z) for z in modes[j])
        diff_out = interaction(roots[i], modes[i], minus_root,
                               conjugate_mode)
        scalar_diff = gadd(gmul(ys[i], gconj(xs[j])),
                           gmul(gconj(ys[j]), xs[i]))
        expected_diff = tuple(gmul(gmul_real(-3, scalar_diff), z)
                              for z in n_gaussian)
        assert diff_out == expected_diff

    # A clean two-mode paired gate exists: r_1=1, r_2=-1 makes the
    # difference coefficient zero and the intended sum nonzero.
    one = (Fraction(1), Fraction(0))
    minus_one = (Fraction(-1), Fraction(0))
    paired_difference = gadd(gmul(one, gconj(minus_one)),
                             gmul(gconj(one), one))
    paired_sum = gadd(gmul(one, minus_one), gmul_real(-1, gmul(one, one)))
    assert paired_difference == ZERO
    assert paired_sum != ZERO

    # For three nonzero in-plane components, vanishing cyclic difference
    # coefficients forces a common purely imaginary ratio and hence zero
    # intended coefficients.  Check an exact nontrivial representative.
    ratios = ((Fraction(0), Fraction(5)),) * 3
    for i, j in ((0, 1), (1, 2), (2, 0)):
        assert gadd(ratios[i], gconj(ratios[j])) == ZERO
        assert gadd(ratios[j], gmul_real(-1, ratios[i])) == ZERO

    # Quantitative identity behind (3.5): with
    # E_23=r_2+conj(r_3), E_31=r_3+conj(r_1),
    # r_2-r_1=E_23-conj(E_31).
    arbitrary = ((Fraction(2), Fraction(-3)),
                 (Fraction(-5), Fraction(7)),
                 (Fraction(11), Fraction(4)))
    e23 = gadd(arbitrary[1], gconj(arbitrary[2]))
    e31 = gadd(arbitrary[2], gconj(arbitrary[0]))
    lhs = gadd(arbitrary[1], gmul_real(-1, arbitrary[0]))
    rhs = gadd(e23, gmul_real(-1, gconj(e31)))
    assert lhs == rhs


def main() -> None:
    check_a2_reality_leak()

    # BAFL has two distinct channels: an n^-6 pre-chart active response and
    # an n^-4 retained wake.  The worst allowed n^2 chart maps the former to
    # an n^-4 endpoint error.
    active_prechart_exponent = Fraction(-6)
    chart_exponent = Fraction(2)
    endpoint_exponent = active_prechart_exponent + chart_exponent
    assert endpoint_exponent == -4
    retained_wake_exponent = Fraction(-4)
    assert retained_wake_exponent == endpoint_exponent

    # The paired gate's child is n^-2 and its first unavoidable
    # parent-child (second-Picard) output is naturally n^-4.  It therefore
    # saturates the wake allowance but needs one extra n^-2 factor before
    # it may return through the active chart.
    seed_exponent = Fraction(-2)
    second_picard_exponent = 2 * seed_exponent
    assert second_picard_exponent == retained_wake_exponent
    assert active_prechart_exponent - second_picard_exponent == -2

    # Bootstrap: E <= Z + K E^2 preserves E <= 2Z if 4KZ <= 1.
    for K in (Fraction(1, 10), Fraction(1), Fraction(17, 3), Fraction(100)):
        Z = Fraction(1, 8) / K
        assert Z + K * (2 * Z) ** 2 <= 2 * Z

    # The same smallness gives contraction constant 4KZ <= 1/2 if we use
    # the stricter KZ <= 1/8 assumed in the note.
    KZ = Fraction(1, 8)
    assert 4 * KZ <= Fraction(1, 2)

    # Integral-test tail for the retained j^-4 wake/error.  Avoid forming
    # thousands of huge Fraction denominators: for decreasing x^-4,
    # sum_{j=J}^infinity j^-4 <= J^-4 + integral_J^infinity x^-4 dx.
    for J in range(2, 100):
        integral_bound = Fraction(1, J**4) + Fraction(1, 3 * J**3)
        tail_bound = Fraction(1, (J - 1) ** 4) + Fraction(1, 3 * (J - 1) ** 3)
        assert integral_bound < tail_bound

    # The combined norm reserves the worst n^2 chart loss.  An n^-5
    # representative paid residual is negligible relative to the n^-4
    # BAFL target; the note assumes the more general o(n^-4) statement.
    paid_combined_exponent = Fraction(-5)
    assert paid_combined_exponent < endpoint_exponent

    print("PASS: one-cell BAFL reduction and A2 leakage arithmetic")


if __name__ == "__main__":
    main()
