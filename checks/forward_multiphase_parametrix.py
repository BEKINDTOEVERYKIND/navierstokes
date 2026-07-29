#!/usr/bin/env python3
"""Exact checks for the forward multiphase Kelvin-parametrix audit.

This script checks only finite-dimensional algebra.  It does not simulate
Navier--Stokes and it does not prove blowup.

Checked facts:

* the three-beat charge coordinates p=K h+Q n are injective on the
  order-M band when K>90 M;
* a nondegenerate transparent superposition on the three coordinate
  carrier directions is a common-helicity Beltrami field;
* the corresponding homochiral child chart for the three published low
  directions has rank only three;
* a heterochiral variation creates an order-K unmatched interaction;
* three arbitrary-polarization sidebands sharing one parent recover the
  full rank-five strain chart, with exact rank-two child maps at finite K;
* eight common-shell homochiral sidebands give an exact rank-eight chart
  of sl(3), with determinant -3/10; and
* six of those columns synthesize diag(-1,-5/4,9/4) exactly.

Only the Python standard library is used.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


F = Fraction
RealVector = tuple[Fraction, Fraction, Fraction]
ComplexVector = tuple[complex, complex, complex]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def matrix_rank(columns: list[list[Fraction]], row_count: int) -> int:
    if not columns:
        return 0
    rows = [
        [columns[column][row] for column in range(len(columns))]
        for row in range(row_count)
    ]
    rank = 0
    for column in range(len(columns)):
        pivot = next(
            (
                row
                for row in range(rank, row_count)
                if rows[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [entry / pivot_value for entry in rows[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            multiplier = rows[row][column]
            rows[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def determinant(columns: list[list[Fraction]]) -> Fraction:
    dimension = len(columns)
    assert all(len(column) == dimension for column in columns)
    rows = [
        [columns[column][row] for column in range(dimension)]
        for row in range(dimension)
    ]
    value = F(1)
    for column in range(dimension):
        pivot = next(
            row
            for row in range(column, dimension)
            if rows[row][column] != 0
        )
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            value *= -1
        pivot_value = rows[column][column]
        value *= pivot_value
        rows[column] = [entry / pivot_value for entry in rows[column]]
        for row in range(column + 1, dimension):
            multiplier = rows[row][column]
            rows[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(rows[row], rows[column])
            ]
    return value


def projected_direction(
    w: tuple[int, int, int],
    q: tuple[int, int, int],
):
    qq = sum(x * x for x in q)
    wq = sum(x * y for x, y in zip(w, q))
    return tuple(F(w[i]) - F(wq, qq) * q[i] for i in range(3))


def full_gradient_column(
    w: tuple[int, int, int],
    q: tuple[int, int, int],
) -> list[Fraction]:
    """Coordinates of (P_q w) tensor q in sl(3)."""

    c = projected_direction(w, q)
    matrix = [[c[i] * q[j] for j in range(3)] for i in range(3)]
    assert sum(matrix[i][i] for i in range(3)) == 0
    return [
        matrix[0][0],
        matrix[1][1],
        matrix[0][1],
        matrix[0][2],
        matrix[1][0],
        matrix[1][2],
        matrix[2][0],
        matrix[2][1],
    ]


def symmetric_column(
    w: tuple[int, int, int],
    q: tuple[int, int, int],
) -> list[Fraction]:
    c = projected_direction(w, q)
    matrix = [
        [F(c[i] * q[j] + q[i] * c[j], 2) for j in range(3)]
        for i in range(3)
    ]
    return [
        matrix[0][0],
        matrix[1][1],
        matrix[0][1],
        matrix[0][2],
        matrix[1][2],
    ]


def check_charge_injectivity() -> None:
    # Columns are q_1,q_2,q_3 from the three-beat note.
    q_matrix = (
        (-45, -4, 1),
        (-36, -5, 1),
        (20, 9, 1),
    )
    determinant_q = (
        q_matrix[0][0]
        * (q_matrix[1][1] * q_matrix[2][2] - q_matrix[1][2] * q_matrix[2][1])
        - q_matrix[0][1]
        * (q_matrix[1][0] * q_matrix[2][2] - q_matrix[1][2] * q_matrix[2][0])
        + q_matrix[0][2]
        * (q_matrix[1][0] * q_matrix[2][1] - q_matrix[1][1] * q_matrix[2][0])
    )
    assert determinant_q == 182

    # Exhaust a representative band.  The proof in the note uses
    # |Q(n-n')|_infinity <= 90 M and works for every M.
    M = 2
    K = 181  # K>90M.
    seen: dict[tuple[int, int, int], tuple[tuple[int, ...], tuple[int, ...]]] = {}
    indices = [
        vector
        for vector in product(range(-M, M + 1), repeat=3)
        if sum(abs(x) for x in vector) <= M
    ]
    for h in indices:
        for n in indices:
            if sum(abs(x) for x in h) + sum(abs(x) for x in n) > M:
                continue
            wave = tuple(
                K * h[row] + sum(q_matrix[row][column] * n[column] for column in range(3))
                for row in range(3)
            )
            assert wave not in seen
            seen[wave] = (h, n)


def check_transparent_coordinate_classification() -> None:
    # For
    # A1=(0,a,b), A2=(c,0,d), A3=(e,f,0),
    # cancellation of all e_i +/- e_j interactions gives
    # b/a=r, d/c=-r, f/e=-1/r and r^2=-1.
    # Test both exact roots.
    for r in (1j, -1j):
        a = c = e = 1.0 + 0.0j
        b = r * a
        d = -r * c
        f = -e / r
        vectors: list[ComplexVector] = [
            (0.0j, a, b),
            (c, 0.0j, d),
            (e, f, 0.0j),
        ]

        for i in range(3):
            for j in range(i + 1, 3):
                ei = tuple(1.0 if axis == i else 0.0 for axis in range(3))
                ej = tuple(1.0 if axis == j else 0.0 for axis in range(3))
                ai = vectors[i]
                aj = vectors[j]
                for sign in (1.0, -1.0):
                    bj = tuple(
                        (value if sign > 0 else value.conjugate())
                        for value in aj
                    )
                    output = tuple(
                        ei[axis] + sign * ej[axis] for axis in range(3)
                    )
                    raw = tuple(
                        dot(ai, tuple(sign * x for x in ej)) * bj[axis]
                        + dot(bj, ei) * ai[axis]
                        for axis in range(3)
                    )
                    # Projection vanishes iff raw is parallel to output.
                    cross = (
                        raw[1] * output[2] - raw[2] * output[1],
                        raw[2] * output[0] - raw[0] * output[2],
                        raw[0] * output[1] - raw[1] * output[0],
                    )
                    assert max(abs(value) for value in cross) < 1.0e-12

    # An exact heterochiral variation is not transparent:
    # A1^+=(0,1,i), A2^-=(1,0,i) creates 2 i e3 at e1+e2.
    a1 = (0.0j, 1.0 + 0.0j, 1.0j)
    a2 = (1.0 + 0.0j, 0.0j, 1.0j)
    raw = tuple(
        a1[1] * a2[axis] + a2[0] * a1[axis]
        for axis in range(3)
    )
    # Projection removes the (1,1,0) part and leaves (0,0,2i).
    projected = (
        raw[0] - (raw[0] + raw[1]) / 2,
        raw[1] - (raw[0] + raw[1]) / 2,
        raw[2],
    )
    assert projected == (0.0j, 0.0j, 2.0j)


def check_three_beat_homochiral_rank() -> None:
    qs = [
        (-45, -36, 20),
        (-4, -5, 9),
        (1, 1, 1),
    ]
    coordinate_directions = [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    ]
    columns = [
        symmetric_column(w, q)
        for w, q in zip(coordinate_directions, qs)
    ]
    assert matrix_rank(columns, 5) == 3


def _leray_fraction(q: tuple[int, int, int], v: RealVector) -> RealVector:
    qq = sum(F(x * x) for x in q)
    vq = sum(v[i] * q[i] for i in range(3))
    return tuple(v[i] - vq * q[i] / qq for i in range(3))  # type: ignore[return-value]


def _single_parent_child(
    q: tuple[int, int, int],
    beta: RealVector,
    carrier: int | None,
) -> RealVector:
    """Child from k=K e3, a=e1 and r=q-k.

    If carrier is None, return the K=infinity map on beta in e3^perp.
    Otherwise lift beta to an exactly r-transverse vector before applying
    the exact child symbol.
    """

    a: RealVector = (F(1), F(0), F(0))
    if carrier is None:
        b = beta
    else:
        assert carrier != q[2]
        beta_q = sum(beta[i] * q[i] for i in range(3))
        b = (beta[0], beta[1], beta_q / (carrier - q[2]))
        r = (q[0], q[1], q[2] - carrier)
        assert dot(b, r) == 0
    aq = F(q[0])
    bq = sum(b[i] * q[i] for i in range(3))
    raw = tuple(aq * b[i] + bq * a[i] for i in range(3))
    return _leray_fraction(q, raw)


def check_single_parent_three_sideband_repair() -> None:
    """Three original q's regain rank two around one transverse shear."""

    qs = [
        (-45, -36, 20),
        (-4, -5, 9),
        (1, 1, 1),
    ]
    target_cs = [
        (F(1, 13), F(-5, 56), F(9, 728)),
        (F(-4, 13), F(25, 56), F(81, 728)),
        (F(16, 13), F(-125, 56), F(729, 728)),
    ]
    beta_basis: list[RealVector] = [
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
    ]

    leading_columns: list[list[Fraction]] = []
    for q in qs:
        assert q[0] != 0 and q[2] != 0
        images = [
            _single_parent_child(q, beta, None)
            for beta in beta_basis
        ]
        # Two child vectors span q^perp iff their cross product is nonzero.
        cross = (
            images[0][1] * images[1][2] - images[0][2] * images[1][1],
            images[0][2] * images[1][0] - images[0][0] * images[1][2],
            images[0][0] * images[1][1] - images[0][1] * images[1][0],
        )
        assert any(entry != 0 for entry in cross)
        assert all(dot(image, q) == 0 for image in images)
        for image in images:
            matrix = [
                [F(image[i] * q[j] + q[i] * image[j], 2) for j in range(3)]
                for i in range(3)
            ]
            leading_columns.append(
                [
                    matrix[0][0],
                    matrix[1][1],
                    matrix[0][1],
                    matrix[0][2],
                    matrix[1][2],
                ]
            )

    assert matrix_rank(leading_columns, 5) == 5
    assert determinant(leading_columns[:5]) == F(1_337_394_240_000, 3_721)

    # Exact finite-K transversality and rank persist.  The only possible
    # kernel value is 2 K q_z=|q|^2, which is nonintegral for all three
    # q's.  Check a range and the full strain rank using exact arithmetic.
    for carrier in range(21, 81):
        exact_columns: list[list[Fraction]] = []
        for q in qs:
            assert 2 * carrier * q[2] != sum(x * x for x in q)
            images = [
                _single_parent_child(q, beta, carrier)
                for beta in beta_basis
            ]
            cross = (
                images[0][1] * images[1][2] - images[0][2] * images[1][1],
                images[0][2] * images[1][0] - images[0][0] * images[1][2],
                images[0][0] * images[1][1] - images[0][1] * images[1][0],
            )
            assert any(entry != 0 for entry in cross)
            for image in images:
                matrix = [
                    [F(image[i] * q[j] + q[i] * image[j], 2) for j in range(3)]
                    for i in range(3)
                ]
                exact_columns.append(
                    [
                        matrix[0][0],
                        matrix[1][1],
                        matrix[0][1],
                        matrix[0][2],
                        matrix[1][2],
                    ]
                )
        assert matrix_rank(exact_columns, 5) == 5

    # The published target polarizations are exactly reachable in the
    # limiting map.  These rational beta coordinates are recorded so the
    # synthesis can be reconstructed independently.
    target_beta_coordinates = [
        (F(-183, 104_000), F(61, 40_950)),
        (F(1_301, 11_648), F(-185, 1_456)),
        (F(2_521, 1_456), F(-1_177, 364)),
    ]
    for q, target, coordinates in zip(qs, target_cs, target_beta_coordinates):
        images = [
            _single_parent_child(q, beta, None)
            for beta in beta_basis
        ]
        reconstructed = tuple(
            coordinates[0] * images[0][axis]
            + coordinates[1] * images[1][axis]
            for axis in range(3)
        )
        assert reconstructed == target


def check_common_shell_repair() -> None:
    # All w have |w|^2=2.  Common-helicity waves on this shell form a
    # Beltrami field, so their leading O(K) Euler interaction is a gradient.
    pairs = [
        ((1, 1, 0), (-1, 0, 0)),
        ((1, 1, 0), (0, -1, 0)),
        ((1, 1, 0), (-1, 0, -1)),
        ((1, 1, 0), (-1, 0, 1)),
        ((1, 1, 0), (0, -1, -1)),
        ((1, 1, 0), (0, -1, 1)),
        ((1, 1, 0), (-2, 0, -1)),
        ((1, -1, 0), (-1, 0, -1)),
    ]
    assert all(sum(x * x for x in w) == 2 for w, _ in pairs)
    assert all(dot(w, q) != 0 for w, q in pairs)

    columns = [full_gradient_column(w, q) for w, q in pairs]
    assert matrix_rank(columns, 8) == 8
    assert determinant(columns) == F(-3, 10)

    # The actual homochiral limiting symbol has the scalar
    # -(w.q)/|w|^2 multiplying P_q w.  Its raw determinant differs from
    # the amplitude-normalized chart by the product of those eight
    # nonzero factors.
    physical_columns = []
    for (w, q), column in zip(pairs, columns):
        factor = -F(dot(w, q), dot(w, w))
        physical_columns.append([factor * entry for entry in column])
    assert determinant(physical_columns) == F(-3, 1280)

    # The first six columns synthesize the exact gamma=5/4 affine pump.
    coefficients = [F(-2), F(-5, 2), F(1), F(1), F(5, 4), F(5, 4), F(0), F(0)]
    result = [
        sum(coefficient * column[row] for coefficient, column in zip(coefficients, columns))
        for row in range(8)
    ]
    target = [F(-1), F(-5, 4), F(0), F(0), F(0), F(0), F(0), F(0)]
    assert result == target

    # The first six pairs use one parent w=(1,1,0).  Their symmetric
    # columns already have rank five.  The subset 1,3,4,5,6 (zero-based
    # 0,2,3,4,5) is an exact square chart.
    symmetric_columns = [
        symmetric_column(w, q)
        for w, q in pairs[:6]
    ]
    assert matrix_rank(symmetric_columns, 5) == 5
    symmetric_subset = [
        symmetric_columns[index]
        for index in (0, 2, 3, 4, 5)
    ]
    assert determinant(symmetric_subset) == F(-1, 8)
    symmetric_coefficients = [F(-9, 2), F(1), F(1), F(5, 4), F(5, 4)]
    symmetric_target = [F(-1), F(-5, 4), F(0), F(0), F(0)]
    assert [
        sum(
            coefficient * column[row]
            for coefficient, column in zip(
                symmetric_coefficients,
                symmetric_subset,
            )
        )
        for row in range(5)
    ] == symmetric_target

    # Partner-partner products can make selected low modes.  With
    # l_a=q_a-Kw for the repeated parent direction, l_7-l_1=q_3 and
    # l_7-l_3=q_1 in the table's one-based indexing.
    q1 = pairs[0][1]
    q3 = pairs[2][1]
    q7 = pairs[6][1]
    assert tuple(x - y for x, y in zip(q7, q1)) == q3
    assert tuple(x - y for x, y in zip(q7, q3)) == q1

    # With the single fast direction w=(1,1,0), all possible child
    # gradients lie in a seven-dimensional hyperplane.  The first seven
    # columns span it exactly.  Its missing coordinate is
    # M13-M23-M31+M32, the vorticity component parallel to w.
    fixed_w_columns = columns[:7]
    assert matrix_rank(fixed_w_columns, 8) == 7
    missing_functional = [F(0), F(0), F(0), F(1), F(0), F(-1), F(-1), F(1)]
    assert all(
        sum(coefficient * entry for coefficient, entry in zip(missing_functional, column))
        == 0
        for column in fixed_w_columns
    )

    # The quadratic order-r null estimate costs at most r^2.  Gevrey-2
    # factorials absorb its convolution uniformly:
    #
    # r^2 sum_{p=1}^{r-1} (p!)^2 ((r-p)!)^2 / (r!)^2
    # = r^2 sum binom(r,p)^(-2).
    for order in range(2, 101):
        value = F(order * order) * sum(
            F(1, _binomial(order, part) ** 2)
            for part in range(1, order)
        )
        assert value <= 3


def _binomial(n: int, k: int) -> int:
    k = min(k, n - k)
    value = 1
    for j in range(1, k + 1):
        value = value * (n - k + j) // j
    return value


def main() -> None:
    check_charge_injectivity()
    check_transparent_coordinate_classification()
    check_three_beat_homochiral_rank()
    check_single_parent_three_sideband_repair()
    check_common_shell_repair()
    print("three-beat charge map: injective for K>90M (representative band checked)")
    print("transparent coordinate bath: common helicity r=+/-i (checked)")
    print("heterochiral unmatched symbol: 2 i e3, hence order K (checked)")
    print("three homochiral child-strain columns: rank 3")
    print("one-parent arbitrary-polarization three-sideband strain chart: rank 5")
    print("one-parent exact child maps: rank 2 for every checked integer K>20")
    print("eight common-shell homochiral full-gradient columns: rank 8")
    print("eight-column normalized sl(3) determinant: -3/10")
    print("eight-column raw homochiral-symbol determinant: -3/1280")
    print("first six columns synthesize diag(-1,-5/4,9/4) exactly")
    print("single-parent symmetric chart: rank 5, minor determinant -1/8")
    print("single-fast-direction child-gradient hyperplane: rank 7")
    print("quadratic charge majorant: Gevrey-2 convolution constant <=3")


if __name__ == "__main__":
    main()
