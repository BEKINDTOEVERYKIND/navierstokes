#!/usr/bin/env python3
"""Exact ledger for the one-material-phase pressure-gauge chart."""

from fractions import Fraction as F


def determinant(matrix: list[list[F]]) -> F:
    """Fraction-preserving Gaussian determinant."""
    work = [row[:] for row in matrix]
    n = len(work)
    out = F(1)
    for col in range(n):
        pivot = next((row for row in range(col, n) if work[row][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            out = -out
        pivot_value = work[col][col]
        out *= pivot_value
        for j in range(col, n):
            work[col][j] /= pivot_value
        for row in range(col + 1, n):
            factor = work[row][col]
            if factor:
                for j in range(col, n):
                    work[row][j] -= factor * work[col][j]
    return out


def cosine_cosine_mean(left_frequency: int, right_frequency: int) -> F:
    return F(1, 2) if left_frequency == right_frequency else F(0)


def cosine_sine_mean(left_frequency: int, right_frequency: int) -> F:
    # Every product cos(m theta) sin(n theta) has zero circle mean.
    del left_frequency, right_frequency
    return F(0)


def main() -> None:
    a = F(3, 2)
    b = F(2, 5)

    # Coordinates on Sym(3): (11, 22, 33, 12, 13, 23).  Columns are
    # gauge, q22, q33, q23, kernel rotation toward e2, kernel rotation
    # toward e3.  The last two columns factor as a and b, so checking the
    # unit determinant proves the general determinant is +/- a*b.
    full_columns = [
        [F(1), F(1), F(1), F(0), F(0), F(0)],
        [F(0), F(1), F(0), F(0), F(0), F(0)],
        [F(0), F(0), F(1), F(0), F(0), F(0)],
        [F(0), F(0), F(0), F(0), F(0), F(1)],
        [F(0), F(0), F(0), -a, F(0), F(0)],
        [F(0), F(0), F(0), F(0), -b, F(0)],
    ]
    full = [list(row) for row in zip(*full_columns)]
    full_det = determinant(full)
    assert abs(full_det) == a * b

    unit_columns = [column[:] for column in full_columns]
    unit_columns[-2][3] = F(-1)
    unit_columns[-1][4] = F(-1)
    unit_full = [list(row) for row in zip(*unit_columns)]
    assert abs(determinant(unit_full)) == 1

    # Quotient coordinates are (Q11-Q33, Q22-Q33, Q12, Q13, Q23).
    quotient = [
        [F(0), F(-1), F(0), F(0), F(0)],
        [F(1), F(-1), F(0), F(0), F(0)],
        [F(0), F(0), F(0), -a, F(0)],
        [F(0), F(0), F(0), F(0), -b],
        [F(0), F(0), F(1), F(0), F(0)],
    ]
    quotient_det = determinant(quotient)
    assert abs(quotient_det) == a * b

    # Rational spectral base: R = rho I - Q, Q = diag(0, 1, 1/4).
    rho = F(7, 3)
    q_diag = (F(0), F(1), F(1, 4))
    r_diag = tuple(rho - q for q in q_diag)
    assert r_diag[0] == rho
    assert r_diag[0] > r_diag[2] > r_diag[1]
    recovered_rho = max(r_diag)
    recovered_q = tuple(recovered_rho - r for r in r_diag)
    assert recovered_q == q_diag
    assert recovered_q[1] > 0 and recovered_q[2] > 0

    # Distinct harmonics realize the two covariance columns with no cross
    # covariance and zero averaged helicity. The curl of a cosine harmonic
    # is a transverse sine harmonic, so every helicity pairing has this
    # zero cosine-sine mean.
    assert cosine_cosine_mean(1, 1) == F(1, 2)
    assert cosine_cosine_mean(2, 2) == F(1, 2)
    assert cosine_cosine_mean(1, 2) == 0
    for left_frequency in (1, 2):
        for right_frequency in (1, 2):
            assert cosine_sine_mean(left_frequency, right_frequency) == 0

    # Checked affine design S = diag(-1, -5/4, 9/4).
    strain = (F(-1), F(-5, 4), F(9, 4))
    work = sum(q * s for q, s in zip(q_diag, strain))
    assert work == F(-11, 16)
    assert work < 0

    # k'=-S^T k: k=e1 grows at rate 1.  A'=-SA for the diagonal
    # transverse polarizations: e2 grows at 5/4 and e3 decays at 9/4.
    k_growth = -strain[0]
    e2_amplitude_growth = -strain[1]
    e3_amplitude_growth = -strain[2]
    assert k_growth == F(1)
    assert e2_amplitude_growth == F(5, 4)
    assert e3_amplitude_growth == F(-9, 4)

    # General work gate: -a beta + b(alpha+beta)<0.
    alpha = F(1)
    beta = F(5, 4)
    a0 = F(1)
    b0 = F(1, 4)
    assert F(0) < b0 < a0 * beta / (alpha + beta)
    assert -a0 * beta + b0 * (alpha + beta) == work

    print("one-material-phase pressure-gauge chart: exact checks passed")
    print(f"  full Jacobian determinant at a={a}, b={b}: {full_det}")
    print(f"  quotient determinant at a={a}, b={b}: {quotient_det}")
    print(f"  affine work Q:S = {work}")
    print(f"  carrier/polarization growth = {k_growth}, {e2_amplitude_growth}")
    print("  distinct-harmonic covariance and mean-helicity ledger passed")


if __name__ == "__main__":
    main()
