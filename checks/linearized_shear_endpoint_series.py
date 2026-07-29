#!/usr/bin/env python3
"""Rational Taylor audit of the low child map about a real one-phase shear.

The base Fourier modes are +/- K e3 with polarization e1.  A perturbation
starts at q-Ke3, and we propagate the exact linearized Euler Fourier
operator.  The script extracts the q-mode endpoint map through a requested
odd Taylor order.  It uses only the standard library.

This is an algebraic/adversarial checker, not a PDE existence proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import math


Matrix = list[list[F]]


def zeros(rows: int, cols: int) -> Matrix:
    return [[F(0) for _ in range(cols)] for _ in range(rows)]


def eye(n: int) -> Matrix:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = F(1)
    return out


def add(a: Matrix, b: Matrix) -> Matrix:
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a: Matrix, c: F) -> Matrix:
    return [[c * x for x in row] for row in a]


def matmul(a: Matrix, b: Matrix) -> Matrix:
    return [
        [sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
         for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]


def inverse2(a: Matrix) -> Matrix:
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if det == 0:
        raise ZeroDivisionError("singular 2x2 matrix")
    return [
        [a[1][1] / det, -a[0][1] / det],
        [-a[1][0] / det, a[0][0] / det],
    ]


def det2(a: Matrix) -> F:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def determinant(a: Matrix) -> F:
    if len(a) != len(a[0]):
        raise ValueError("determinant requires a square matrix")
    work = [row[:] for row in a]
    out = F(1)
    for col in range(len(work)):
        pivot = next((row for row in range(col, len(work))
                      if work[row][col] != 0), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            out = -out
        value = work[col][col]
        out *= value
        for row in range(col + 1, len(work)):
            factor = work[row][col] / value
            for j in range(col + 1, len(work)):
                work[row][j] -= factor * work[col][j]
    return out


def projector(p: tuple[int, int, int]) -> Matrix:
    den = sum(x * x for x in p)
    out = eye(3)
    for i in range(3):
        for j in range(3):
            out[i][j] -= F(p[i] * p[j], den)
    return out


def source_basis(q: tuple[int, int, int], K: int) -> Matrix:
    alpha, beta, delta = q
    # Uniformly bounded graph basis for (q-Ke3)^perp.
    return [
        [F(1), F(0)],
        [F(0), F(1)],
        [F(alpha, K - delta), F(beta, K - delta)],
    ]


def target_basis(q: tuple[int, int, int]) -> Matrix:
    alpha, beta, delta = q
    return [
        [F(-beta), F(-delta)],
        [F(alpha), F(0)],
        [F(0), F(alpha)],
    ]


def target_coordinates(q: tuple[int, int, int], value: Matrix) -> Matrix:
    c = target_basis(q)
    ct = transpose(c)
    return matmul(matmul(inverse2(matmul(ct, c)), ct), value)


def apply_linearized_shear(
    q: tuple[int, int, int],
    K: int,
    state: dict[int, Matrix],
) -> dict[int, Matrix]:
    """Apply the exact Fourier linearization once.

    If p_m=q+mKe3, interaction with the base sign s sends m to m+s:

        P_{p_{m+s}} [ alpha I + s K e1 tensor e3 ].

    The common harmless Fourier factor -i and real-cosine normalization
    are omitted.  They do not affect rank or K growth.
    """

    alpha, beta, delta = q
    out: dict[int, Matrix] = {}
    for m, value in state.items():
        for sign in (-1, 1):
            mt = m + sign
            p = (alpha, beta, delta + mt * K)
            raw = scale(eye(3), F(alpha))
            raw[0][2] += F(sign * K)
            contribution = matmul(matmul(projector(p), raw), value)
            out[mt] = add(out.get(mt, zeros(3, 2)), contribution)
    return out


def endpoint_coefficients(
    q: tuple[int, int, int],
    K: int,
    max_order: int,
) -> dict[int, Matrix]:
    state = {-1: source_basis(q, K)}
    result: dict[int, Matrix] = {}
    for order in range(1, max_order + 1):
        state = apply_linearized_shear(q, K, state)
        if 0 in state:
            result[order] = scale(
                target_coordinates(q, state[0]),
                F(1, math.factorial(order)),
            )
    return result


def frobenius(a: Matrix) -> float:
    return math.sqrt(sum(float(x) ** 2 for row in a for x in row))


def endpoint_series(coeffs: dict[int, Matrix], time: F) -> Matrix:
    out = zeros(2, 2)
    for order, matrix in coeffs.items():
        out = add(out, scale(matrix, time ** order))
    return out


def strain_columns(
    q: tuple[int, int, int],
    endpoint: Matrix,
) -> list[list[F]]:
    """Return two STF strain columns in 5 independent coordinates."""

    child = matmul(target_basis(q), endpoint)
    columns: list[list[F]] = []
    for control in range(2):
        c = [child[i][control] for i in range(3)]
        alpha, beta, delta = map(F, q)
        # Coordinates (M11,M22,M12,M13,M23) of sym(c tensor q).
        columns.append(
            [
                2 * c[0] * alpha,
                2 * c[1] * beta,
                c[0] * beta + c[1] * alpha,
                c[0] * delta + c[2] * alpha,
                c[1] * delta + c[2] * beta,
            ]
        )
    return columns


def combined_rank_audit(
    max_order: int,
    carriers: list[int],
    time: F,
) -> None:
    qs = ((-45, -36, 20), (-4, -5, 9), (1, 1, 1))
    print(
        f"combined five-coordinate strain minor at t={time} "
        f"(Taylor order {max_order})"
    )
    for K in carriers:
        all_columns: list[list[F]] = []
        for q in qs:
            coeffs = endpoint_coefficients(q, K, max_order)
            all_columns.extend(strain_columns(q, endpoint_series(coeffs, time)))
        # First five of the six columns inherit the exact principal minor.
        minor = [[all_columns[col][row] for col in range(5)]
                 for row in range(5)]
        det = determinant(minor)
        if det == 0:
            raise AssertionError(f"combined rank defect at K={K}")
        normalized = float(det / (time ** 5))
        print(f"  K={K}: det/t^5={normalized:.10g}")


def audit(q: tuple[int, int, int], max_order: int, carriers: list[int]) -> None:
    print(f"q={q}")
    norms: dict[int, list[float]] = {}
    for K in carriers:
        coeffs = endpoint_coefficients(q, K, max_order)
        first_det = det2(coeffs[1])
        if first_det == 0:
            raise AssertionError(f"first child map singular at q={q}, K={K}")
        for order, matrix in coeffs.items():
            norms.setdefault(order, []).append(frobenius(matrix))

    ratio = carriers[-1] / carriers[-2]
    for order, values in sorted(norms.items()):
        exponent = math.log(values[-1] / values[-2]) / math.log(ratio)
        print(
            f"  t^{order}: norm(K={carriers[-1]})={values[-1]:.8g} "
            f"last-scale exponent={exponent:.5f}"
        )
    print("  exact rational first-order determinant: nonzero at all samples")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=7)
    parser.add_argument(
        "--carriers",
        type=int,
        nargs="+",
        default=[128, 256, 512, 1024],
    )
    parser.add_argument(
        "--time-denominator",
        type=int,
        default=100,
        help="audit the truncated endpoint at t=1/N",
    )
    args = parser.parse_args()
    if args.max_order < 1:
        raise SystemExit("--max-order must be positive")
    if len(args.carriers) < 2 or min(args.carriers) <= 20:
        raise SystemExit("supply at least two carriers greater than 20")

    for q in ((-45, -36, 20), (-4, -5, 9), (1, 1, 1)):
        audit(q, args.max_order, args.carriers)
    combined_rank_audit(
        args.max_order,
        args.carriers,
        F(1, args.time_denominator),
    )


if __name__ == "__main__":
    main()
