#!/usr/bin/env python3
"""Exact arithmetic checks for C142--C143.

The checker verifies the affine child/wake formulas, the universal
derivative relation behind the selector tradeoff, the curl completion,
and the scale identities in the viscous-envelope/collar audit.  It does
not assert the open finite-energy Lagrangian localization theorem.
"""

from fractions import Fraction as F
from math import factorial, log


N = (F(1), F(1), F(1))
KC = (F(1), F(0), F(-1))
E1 = (F(2), F(-1), F(-1))
E2 = (F(1), F(1), F(-2))


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vscale(c, a):
    return tuple(c * x for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def matvec(a, x):
    return tuple(sum(a[i][j] * x[j] for j in range(3)) for i in range(3))


def transpose(a):
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))


def determinant(a):
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inverse(a):
    det = determinant(a)
    assert det
    cof = (
        (
            a[1][1] * a[2][2] - a[1][2] * a[2][1],
            -(a[1][0] * a[2][2] - a[1][2] * a[2][0]),
            a[1][0] * a[2][1] - a[1][1] * a[2][0],
        ),
        (
            -(a[0][1] * a[2][2] - a[0][2] * a[2][1]),
            a[0][0] * a[2][2] - a[0][2] * a[2][0],
            -(a[0][0] * a[2][1] - a[0][1] * a[2][0]),
        ),
        (
            a[0][1] * a[1][2] - a[0][2] * a[1][1],
            -(a[0][0] * a[1][2] - a[0][2] * a[1][0]),
            a[0][0] * a[1][1] - a[0][1] * a[1][0],
        ),
    )
    return tuple(
        tuple(transpose(cof)[i][j] / det for j in range(3))
        for i in range(3)
    )


def check_geometry_and_kelvin_identity():
    assert vadd(E1, E2) == vscale(3, KC)
    assert dot(N, KC) == 0
    assert dot(N, E1) == 0
    assert dot(N, E2) == 0
    assert dot(KC, KC) == 2
    assert dot(E1, E1) == dot(E2, E2) == 6

    matrices = (
        ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))),
        ((F(2), F(0), F(0)), (F(0), F(1, 2), F(0)), (F(0), F(0), F(1))),
        ((F(1), F(1), F(0)), (F(0), F(1), F(1)), (F(0), F(0), F(1))),
    )
    for g in matrices:
        assert determinant(g) == 1
        f = inverse(transpose(g))
        # F(k x N) = (Gk) x (GN) for det F = 1.
        for k in (KC, E1, E2):
            assert matvec(f, cross(k, N)) == cross(matvec(g, k), matvec(g, N))

        # The vector identity whose triangle inequality is C142 (2.12).
        left = cross(matvec(g, KC), matvec(g, N))
        right = vscale(
            F(1, 3),
            vadd(
                cross(matvec(g, E1), matvec(g, N)),
                cross(matvec(g, E2), matvec(g, N)),
            ),
        )
        assert left == right


def check_selector_formulas():
    for h in (F(1), F(2), F(7, 3), F(11)):
        dcap = (h * h + 1 / (h * h)) / 2
        s = (h + 1 / h) / 2
        d = (h - 1 / h) / 2
        assert s * s - d * d == 1
        assert s * s + d * d == dcap

        child_frequency_sq = F(2)
        child_amplitude_sq = 3 * dcap
        wake_frequency_sq = F(3, 2) * (dcap + 3)
        wake_amplitude_sq = 3 * (3 * dcap + 1) / (dcap + 3)

        # Derive the wake projection rather than merely restating (2.8).
        # In the (u,v,w) frame, Ge_i has v,w part
        # +/-sqrt(6)/2*(s,d), while GN=sqrt(3)*(d,s).
        gn_sq = 3 * dcap
        wake_dot_gn_sq = 18 * s * s * d * d
        assert wake_amplitude_sq == gn_sq - wake_dot_gn_sq / wake_frequency_sq

        relative_child_gain_sq = child_amplitude_sq / 3
        relative_wake_gain_sq = wake_amplitude_sq / 3
        assert relative_child_gain_sq == dcap
        assert 1 <= relative_wake_gain_sq < 3

        x_child_sq = child_frequency_sq * child_amplitude_sq
        x_wake_sq = wake_frequency_sq * wake_amplitude_sq
        assert x_wake_sq == F(9, 4) * x_child_sq + F(9, 2)

        # The transverse child map has eigenvalues h and h^{-1}.
        inverse_norm = h
        condition_number = h * h
        assert inverse_norm * (1 / h) == 1
        assert condition_number == h / (1 / h)


def levi_civita(i, j, k):
    if len({i, j, k}) < 3:
        return F(0)
    return F(1) if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else F(-1)


def check_curl_completion():
    matrices = (
        (
            (F(1), F(2), F(-1)),
            (F(2), F(-3), F(4)),
            (F(-1), F(4), F(2)),
        ),
        (
            (F(0), F(1, 3), F(2)),
            (F(1, 3), F(5), F(-1)),
            (F(2), F(-1), F(-5)),
        ),
    )
    for s_matrix in matrices:
        assert s_matrix == transpose(s_matrix)
        assert sum(s_matrix[i][i] for i in range(3)) == 0

        # A_i = (Sx x x)_i / 3 = sum C[i,m,k] x_m x_k.
        coeff = [
            [[F(0) for _ in range(3)] for _ in range(3)]
            for _ in range(3)
        ]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    eps = levi_civita(i, j, k)
                    for m in range(3):
                        coeff[i][m][k] += eps * s_matrix[j][m] / 3

        # (curl A)_ell = sum_n curl_coeff[ell,n] x_n.
        curl_coeff = [[F(0) for _ in range(3)] for _ in range(3)]
        for ell in range(3):
            for r in range(3):
                for i in range(3):
                    eps = levi_civita(ell, r, i)
                    for n in range(3):
                        curl_coeff[ell][n] += eps * (
                            coeff[i][r][n] + coeff[i][n][r]
                        )
        assert tuple(tuple(row) for row in curl_coeff) == s_matrix


def check_viscous_envelope_obstruction():
    # Treat log(n) as the formal unit 1. Then log(h)=12 for h=n^12.
    for n in (2, 3, 5):
        q = F(n**8)
        h = F(n**12)
        ell = F(1)  # L in the note; avoid shadowing wavelength ell.
        nu_k_sq = F(1)
        log_n = F(1)
        log_h = F(12)
        lam = 3 * nu_k_sq * (h * h - 1 / (h * h)) / (8 * ell * log_n)
        time = log_h / lam

        child_exponent = 2 * nu_k_sq * time
        expected_child = (
            16 * ell * log_n * log_h / (3 * (h * h - 1 / (h * h)))
        )
        assert child_exponent == expected_child

        selected_wake_term = (
            nu_k_sq * 3 * (h * h - 1 / (h * h)) / (8 * lam)
        )
        assert selected_wake_term == ell * log_n
        full_wake_exponent = nu_k_sq * (F(9, 2) * time) + selected_wake_term
        assert full_wake_exponent >= ell * log_n

        envelope_exponent = (
            nu_k_sq * q * q * (h * h - 1) / (2 * lam)
        )
        expected_envelope = (
            4
            * ell
            * q
            * q
            * log_n
            * (h * h - 1)
            / (3 * (h * h - 1 / (h * h)))
        )
        assert envelope_exponent == expected_envelope
        assert envelope_exponent > q * q
        assert h / q == n**4  # launch width is n^4 parent widths.

        # h=q^(3/2) exactly, so a stationary L2 collar becomes log(h)
        # after the backward focus.
        assert h * h == q**3
        assert h * F(1, n**12) == 1  # q^{-3/2}=n^{-12}.

        # Exact p=2 cutoff identities, using an arbitrary dimensionless
        # viscosity mu.  Here q^(1/2)=n^4 and q^(11/2)=n^44.
        mu = F(2, 7)
        lam_normalized = (
            3 * mu * (h * h - 1 / (h * h)) / (8 * ell * log_n)
        )
        epsilon_e_2 = lam_normalized * F(1, n**20) * log_h
        expected_e_2 = (
            F(9, 2)
            * mu
            / ell
            * (F(n**4) - F(1, n**44))
        )
        assert epsilon_e_2 == expected_e_2
        epsilon_nu_2 = mu * F(1, n**4) * log_h
        assert epsilon_nu_2 == 12 * mu / n**4

        # Energy, dissipation, and relative-fidelity ledgers.
        energy_ratio = lam_normalized**2 / n**40
        assert energy_ratio == lam_normalized**2 * q**-5
        dissipation_ratio = mu * lam_normalized * q**-3 * log_h
        expected_dissipation = (
            F(9, 2) * mu * mu / ell * (1 - q**-6)
        )
        assert dissipation_ratio == expected_dissipation
        core_size_p2 = lam_normalized * F(1, n**20)
        assert epsilon_e_2 / core_size_p2 == log_h

        # The stationary L2 collar product bound loses exactly log(h)
        # after the backward factor h=q^(3/2).
        collar_p2 = F(1, n**12) * log_h
        assert h * collar_p2 == log_h


def check_factorial_and_collar_scales():
    # The backward-focused absolute L2 self/viscous errors are factorial
    # small. Check a late exact tail and strict decay of their n^6 weights.
    weighted_self = []
    weighted_viscous = []
    for j in range(16, 23):
        n = j + 1
        mu = 1 / (factorial(j) ** 2)
        # Constants and viscosity are fixed; only the asymptotic powers matter.
        backward_self = mu * n**16
        backward_viscous = 12 * mu * n**8 * log(n)
        weighted_self.append(n**6 * backward_self)
        weighted_viscous.append(n**6 * backward_viscous)
    assert all(x > y for x, y in zip(weighted_self, weighted_self[1:]))
    assert all(x > y for x, y in zip(weighted_viscous, weighted_viscous[1:]))
    assert weighted_self[-1] < 1
    assert weighted_viscous[-1] < 1

    # A curvature gain q^{-1} makes the backward collar scale o(n^{-6}).
    ratios = [12 * log(n) / (n * n) for n in (20, 50, 100, 200)]
    assert all(x > y for x, y in zip(ratios, ratios[1:]))
    assert ratios[-1] < 1

    # Generic wake exposure b^2 log(h) is larger than b^3.
    for n in (2, 5, 20):
        b = 1 / (n * n)
        assert b * b * (12 * log(n)) > b**3


def main():
    check_geometry_and_kelvin_identity()
    check_selector_formulas()
    check_curl_completion()
    check_viscous_envelope_obstruction()
    check_factorial_and_collar_scales()
    print("C142--C143 affine selector/localization checks: PASS")
    print("OPEN: co-moving finite-energy collar slaving and full BAFL")


if __name__ == "__main__":
    main()
