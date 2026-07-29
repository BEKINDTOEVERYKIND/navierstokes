#!/usr/bin/env python3
"""Exact-regression checks for the two-colour Euler beat calculation.

This script does not simulate Navier--Stokes.  It checks the finite
Fourier algebra used in

    research/2026-07-29-two-colour-endpoint-rank.md.

For two real divergence-free modes at k and l, the Euler coefficient at
n=k+l is, up to the common factor -i,

    P_n ((a.l) b + (b.k) a).

The checks below verify:

* the rank-one collapse for equal carrier radii |k|=|l|;
* rank two on (k+l)^perp for an unequal-radius pair;
* a uniformly nondegenerate explicit integer family;
* the exact O(1/K) difference-sideband leakage ratio; and
* a transverse-polarization control whose derivative at k-l is zero;
* an exact three-beat affine synthesis and rank-five strain chart.

Only the Python standard library is used.
"""

from __future__ import annotations

from fractions import Fraction
from math import isclose, sqrt
from typing import Iterable


Vector = tuple[float, float, float]


def add(a: Vector, b: Vector) -> Vector:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def sub(a: Vector, b: Vector) -> Vector:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def scale(c: float, a: Vector) -> Vector:
    return tuple(c * x for x in a)  # type: ignore[return-value]


def dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm(a: Vector) -> float:
    return sqrt(dot(a, a))


def leray(n: Vector, v: Vector) -> Vector:
    nn = dot(n, n)
    assert nn > 0.0
    return sub(v, scale(dot(v, n) / nn, n))


def interaction(n: Vector, p: Vector, a: Vector, r: Vector, b: Vector) -> Vector:
    """Symmetrized velocity-form Euler interaction at n=p+r."""

    assert all(isclose(x + y, z, abs_tol=1.0e-12) for x, y, z in zip(p, r, n))
    raw = add(scale(dot(a, r), b), scale(dot(b, p), a))
    return leray(n, raw)


def coordinates(v: Vector, e1: Vector, e2: Vector) -> tuple[float, float]:
    return dot(v, e1), dot(v, e2)


def rank_2d(columns: Iterable[tuple[float, float]], tol: float = 1.0e-10) -> int:
    cols = list(columns)
    if not any(abs(x) > tol or abs(y) > tol for x, y in cols):
        return 0
    for i, (x1, y1) in enumerate(cols):
        for x2, y2 in cols[i + 1 :]:
            if abs(x1 * y2 - y1 * x2) > tol:
                return 2
    return 1


def carrier_geometry(k_transverse: float, q_length: float, k_longitudinal: float):
    """Return q,k,l,d and normalized in-plane transverse polarizations."""

    K = k_transverse
    Q = q_length
    kappa = k_longitudinal
    q = (0.0, 0.0, Q)
    k = (K, 0.0, kappa)
    l = (-K, 0.0, Q - kappa)
    d = sub(k, l)
    nk = norm(k)
    nl = norm(l)
    p = (kappa / nk, 0.0, -K / nk)
    r = ((Q - kappa) / nl, 0.0, K / nl)
    e1 = (1.0, 0.0, 0.0)
    e2 = (0.0, 1.0, 0.0)
    assert isclose(dot(k, p), 0.0, abs_tol=1.0e-12)
    assert isclose(dot(l, r), 0.0, abs_tol=1.0e-12)
    return q, k, l, d, p, r, e1, e2, nk, nl


def image_rank(K: float, Q: float, kappa: float) -> int:
    q, k, l, _, p, r, e1, e2, _, _ = carrier_geometry(K, Q, kappa)
    basis_k = (p, e2)
    basis_l = (r, e2)
    columns = [
        coordinates(interaction(q, k, a, l, b), e1, e2)
        for a in basis_k
        for b in basis_l
    ]
    return rank_2d(columns)


def check_equal_radius_collapse() -> None:
    # q=(0,0,2), k=(K,0,1), l=(-K,0,1): |k|=|l|.
    for K in (2.0, 8.0, 64.0):
        assert image_rank(K, 2.0, 1.0) == 1


def check_unequal_radius_family() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for N in (2.0, 4.0, 8.0, 16.0, 64.0, 256.0):
        # Integer family q=(0,0,1), k=(N,0,N), l=(-N,0,1-N).
        q, k, l, d, p, r, e1, e2, nk, nl = carrier_geometry(N, 1.0, N)
        assert image_rank(N, 1.0, N) == 2

        low = interaction(q, k, p, l, r)
        high = interaction(d, k, p, scale(-1.0, l), r)
        low_e1, low_e2 = coordinates(low, e1, e2)
        assert isclose(low_e2, 0.0, abs_tol=1.0e-12)

        predicted_low = N * (2.0 * N - 1.0) / (nk * nl)
        assert isclose(low_e1, predicted_low, rel_tol=1.0e-12, abs_tol=1.0e-12)

        leakage_ratio = norm(high) / abs(low_e1)
        predicted_ratio = 1.0 / norm(d)
        assert isclose(
            leakage_ratio,
            predicted_ratio,
            rel_tol=2.0e-12,
            abs_tol=2.0e-12,
        )

        # Polarization control:
        # da=(nl/(2KQ)) e2 and db=-(nk/(2KQ)) e2.
        # Its low derivative is exactly e2 and its k-l derivative is zero.
        da = scale(nl / (2.0 * N), e2)
        db = scale(-nk / (2.0 * N), e2)
        dlow = add(
            interaction(q, k, da, l, r),
            interaction(q, k, p, l, db),
        )
        dhigh = add(
            interaction(d, k, da, scale(-1.0, l), r),
            interaction(d, k, p, scale(-1.0, l), db),
        )
        assert isclose(dot(dlow, e1), 0.0, abs_tol=2.0e-12)
        assert isclose(dot(dlow, e2), 1.0, rel_tol=2.0e-12, abs_tol=2.0e-12)
        assert norm(dhigh) < 2.0e-12

        # The amplitude-product control and the polarization control give a
        # diagonal 2x2 projected Jacobian diag(low_e1,1).
        determinant = low_e1
        assert determinant > 0.5

        rows.append(
            {
                "N": N,
                "child_parallel": low_e1,
                "jacobian_determinant": determinant,
                "difference_to_child": leakage_ratio,
                "difference_scaled_by_N": N * leakage_ratio,
            }
        )
    return rows


def check_charge_grading() -> None:
    # Labels (r,s) represent r*k+s*l.  With h=r-s and n=s,
    # r*k+s*l = h*k+n*q.  Quadratic addition adds h exactly.
    labels = [
        (1, 0),
        (0, 1),
        (1, 1),
        (1, -1),
        (2, 1),
        (1, 2),
        (-3, 2),
    ]
    q, k, l, _, *_ = carrier_geometry(7.0, 1.0, 7.0)
    for r_label, s_label in labels:
        direct = add(scale(r_label, k), scale(s_label, l))
        h = r_label - s_label
        charged = add(scale(h, k), scale(s_label, q))
        assert norm(sub(direct, charged)) < 1.0e-12

    for a in labels:
        for b in labels:
            charge_a = a[0] - a[1]
            charge_b = b[0] - b[1]
            summed = (a[0] + b[0], a[1] + b[1])
            assert summed[0] - summed[1] == charge_a + charge_b


def check_three_beat_affine_witness() -> None:
    """Three low beat directions exactly synthesize the gamma=5/4 strain.

    The identity is over the rationals.  Each summand c_j tensor q_j is
    trace-free because c_j.q_j=0, as required for a divergence-free sine
    mode c_j sin(q_j.x).
    """

    F = Fraction
    qs = [
        (F(-45), F(-36), F(20)),
        (F(-4), F(-5), F(9)),
        (F(1), F(1), F(1)),
    ]
    cs = [
        (F(1, 13), F(-5, 56), F(9, 728)),
        (F(-4, 13), F(25, 56), F(81, 728)),
        (F(16, 13), F(-125, 56), F(729, 728)),
    ]
    for q, c in zip(qs, cs):
        assert sum(x * y for x, y in zip(q, c)) == 0

    matrix = [
        [sum(c[i] * q[j] for c, q in zip(cs, qs)) for j in range(3)]
        for i in range(3)
    ]
    target = [
        [F(-1), F(0), F(0)],
        [F(0), F(-5, 4), F(0)],
        [F(0), F(0), F(9, 4)],
    ]
    assert matrix == target

    # The same three directions give a full chart of Sym_0(3): each q has
    # two transverse polarizations, hence six columns, and their symmetric
    # gradients have exact rational rank five.
    columns: list[list[Fraction]] = []
    for q in qs:
        transverse = [
            (-q[1], q[0], F(0)),
            (-q[2], F(0), q[0]),
        ]
        for c in transverse:
            symmetric = [
                [F(c[i] * q[j] + q[i] * c[j], 2) for j in range(3)]
                for i in range(3)
            ]
            columns.append(
                [
                    symmetric[0][0],
                    symmetric[1][1],
                    symmetric[0][1],
                    symmetric[0][2],
                    symmetric[1][2],
                ]
            )

    # Fraction-valued Gaussian elimination on the 5-by-6 matrix.
    rows = [[columns[j][i] for j in range(6)] for i in range(5)]

    # The note records this explicit first-five-column minor.  Check its
    # exact value in addition to checking the rank of all six columns.
    minor = [[columns[j][i] for j in range(5)] for i in range(5)]
    determinant = Fraction(1)
    for diagonal in range(5):
        pivot = next(
            row
            for row in range(diagonal, 5)
            if minor[row][diagonal] != 0
        )
        if pivot != diagonal:
            minor[diagonal], minor[pivot] = minor[pivot], minor[diagonal]
            determinant *= -1
        pivot_value = minor[diagonal][diagonal]
        determinant *= pivot_value
        minor[diagonal] = [
            entry / pivot_value for entry in minor[diagonal]
        ]
        for row in range(diagonal + 1, 5):
            multiplier = minor[row][diagonal]
            minor[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(
                    minor[row],
                    minor[diagonal],
                )
            ]
    assert determinant == -1_214_003_700

    rank = 0
    for column in range(6):
        pivot = next(
            (row for row in range(rank, 5) if rows[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [entry / pivot_value for entry in rows[rank]]
        for row in range(5):
            if row == rank:
                continue
            multiplier = rows[row][column]
            rows[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == 5:
            break
    assert rank == 5


def main() -> None:
    check_equal_radius_collapse()
    rows = check_unequal_radius_family()
    check_charge_grading()
    check_three_beat_affine_witness()

    print("equal-radius image rank: 1")
    print("unequal-radius image rank: 2")
    print(
        "       N    child_parallel      det(J)      |difference|/|child|"
        "      N*ratio"
    )
    for row in rows:
        print(
            f"{row['N']:8.0f}"
            f" {row['child_parallel']:17.12f}"
            f" {row['jacobian_determinant']:12.9f}"
            f" {row['difference_to_child']:25.12e}"
            f" {row['difference_scaled_by_N']:12.9f}"
        )
    print("polarization-control difference derivative: exactly zero (checked)")
    print("charge-chain grading: exact (checked)")
    print("three-beat gamma=5/4 affine synthesis: exact over Q (checked)")
    print("three-beat symmetric trace-free chart rank: 5 over Q (checked)")


if __name__ == "__main__":
    main()
