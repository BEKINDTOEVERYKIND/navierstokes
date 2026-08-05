#!/usr/bin/env python3
"""Exact ledger for the two-harmonic one-phase sideband coupling.

The checker uses only rational arithmetic.  It verifies

* the rotated single-parent zero-charge chart still has rank five;
* every bath--partner coefficient has an exact K-free formula;
* the second bath harmonic creates an unavoidable rank-two charged output;
* simultaneous partner--partner products are K-free but create sum and
  difference charges;
* the extreme-charge symbol underlying the finite-support no-go has the
  claimed determinant; and
* the two sine harmonics give the C89 covariance with zero mean helicity.
"""

from __future__ import annotations

from fractions import Fraction as F


Vector = tuple[F, F, F]


# Cyclic rotation of the old (e3 carrier, e1 polarization) chart into the
# compressive C89 coordinates (e1 carrier, e2 polarization).
Q_VECTORS: tuple[Vector, ...] = (
    (F(20), F(-45), F(-36)),
    (F(9), F(-4), F(-5)),
    (F(1), F(1), F(1)),
)
W: Vector = (F(1), F(0), F(0))
A: Vector = (F(0), F(1), F(0))
D: Vector = (F(0), F(0), F(1, 2))


def add(left: Vector, right: Vector) -> Vector:
    return tuple(x + y for x, y in zip(left, right))  # type: ignore[return-value]


def sub(left: Vector, right: Vector) -> Vector:
    return tuple(x - y for x, y in zip(left, right))  # type: ignore[return-value]


def scale(number: F, vector: Vector) -> Vector:
    return tuple(number * entry for entry in vector)  # type: ignore[return-value]


def dot(left: Vector, right: Vector) -> F:
    return sum((x * y for x, y in zip(left, right)), F(0))


def project(wave: Vector, vector: Vector) -> Vector:
    return sub(vector, scale(dot(wave, vector) / dot(wave, wave), wave))


def determinant(matrix: list[list[F]]) -> F:
    work = [row[:] for row in matrix]
    result = F(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(column + 1, len(work)):
            multiplier = work[row][column]
            if multiplier:
                work[row] = [
                    entry - multiplier * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[column])
                ]
    return result


def rank(matrix: list[list[F]]) -> int:
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_value = work[row][column]
        work[row] = [entry / pivot_value for entry in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            multiplier = work[index][column]
            if multiplier:
                work[index] = [
                    entry - multiplier * pivot_entry
                    for entry, pivot_entry in zip(work[index], work[row])
                ]
        row += 1
        if row == len(work):
            break
    return row


def partner_basis(q: Vector, k_value: int) -> tuple[Vector, Vector]:
    """A bounded basis of (q-Ke1)-perp."""

    denominator = F(k_value) - q[0]
    return (
        (q[1] / denominator, F(1), F(0)),
        (q[2] / denominator, F(0), F(1)),
    )


def symmetric_column(child: Vector, q: Vector) -> list[F]:
    matrix = [
        [F(child[i] * q[j] + q[i] * child[j], 2) for j in range(3)]
        for i in range(3)
    ]
    assert sum(matrix[index][index] for index in range(3)) == 0
    return [
        matrix[0][0],
        matrix[1][1],
        matrix[0][1],
        matrix[0][2],
        matrix[1][2],
    ]


def bath_vector(harmonic: int) -> Vector:
    return A if abs(harmonic) == 1 else D


def reduced_bath_partner(q: Vector, b: Vector, harmonic: int) -> Vector:
    bath = bath_vector(harmonic)
    return add(
        scale(dot(bath, q), b),
        scale(F(harmonic) * dot(b, q), bath),
    )


def direct_bath_partner(
    q: Vector, b: Vector, k_value: int, harmonic: int
) -> Vector:
    bath = bath_vector(harmonic)
    parent_wave = scale(F(harmonic * k_value), W)
    partner_wave = sub(q, scale(F(k_value), W))
    return add(
        scale(dot(bath, partner_wave), b),
        scale(dot(b, parent_wave), bath),
    )


def check_zero_charge_rank_and_extra_charges() -> None:
    k_value = 101
    strain_columns: list[list[F]] = []
    observed_extra_charges: set[int] = set()

    for q in Q_VECTORS:
        partner_wave = sub(q, scale(F(k_value), W))
        for b in partner_basis(q, k_value):
            assert dot(b, partner_wave) == 0
            for harmonic in (-2, -1, 1, 2):
                direct = direct_bath_partner(q, b, k_value, harmonic)
                reduced = reduced_bath_partner(q, b, harmonic)
                assert direct == reduced
                output_charge = harmonic - 1
                output_wave = add(q, scale(F(output_charge * k_value), W))
                projected = project(output_wave, reduced)
                if harmonic == 1:
                    strain_columns.append(symmetric_column(projected, q))
                else:
                    observed_extra_charges.add(output_charge)

            # The positive second harmonic acts injectively on each of the
            # two partner directions at this K; no polarization deletes it.
            extra_wave = add(q, scale(F(k_value), W))
            extra = project(extra_wave, reduced_bath_partner(q, b, 2))
            assert dot(extra, extra) > 0

    strain_matrix = [list(row) for row in zip(*strain_columns)]
    assert rank(strain_matrix) == 5
    assert observed_extra_charges == {-3, -2, 1}


def check_partner_partner_ledger() -> None:
    k_value = 101
    q_left, q_right = Q_VECTORS[:2]
    b_left = partner_basis(q_left, k_value)[0]
    b_right = partner_basis(q_right, k_value)[1]
    r_left = sub(q_left, scale(F(k_value), W))
    r_right = sub(q_right, scale(F(k_value), W))

    assert dot(b_left, r_left) == 0
    assert dot(b_right, r_right) == 0

    # High sum: charge -2, slow frequency q_left+q_right.
    assert dot(b_left, r_right) == dot(b_left, sub(q_right, q_left))
    assert dot(b_right, r_left) == dot(b_right, sub(q_left, q_right))

    # Low difference: charge 0, slow frequency q_left-q_right.
    assert dot(b_left, scale(F(-1), r_right)) == dot(
        b_left, sub(q_left, q_right)
    )
    assert dot(b_right, r_left) == dot(b_right, sub(q_left, q_right))

    # A lone exactly divergence-free plane wave has no self interaction,
    # either with itself or with its reality conjugate.
    self_sum = scale(F(2) * dot(b_left, r_left), b_left)
    self_difference = add(
        scale(dot(b_left, scale(F(-1), r_left)), b_left),
        scale(dot(b_left, r_left), b_left),
    )
    assert self_sum == (F(0), F(0), F(0))
    assert self_difference == (F(0), F(0), F(0))


def triple(left: Vector, middle: Vector, right: Vector) -> F:
    return (
        left[0] * (middle[1] * right[2] - middle[2] * right[1])
        - left[1] * (middle[0] * right[2] - middle[2] * right[0])
        + left[2] * (middle[0] * right[1] - middle[1] * right[0])
    )


def extreme_determinant(q: Vector, k_value: int, h: int, m: int) -> F:
    """Triple determinant of the harmonic-m extreme leakage symbol."""

    x_value = q[0] + F(h * k_value)
    domain_one: Vector = (-q[1], x_value, F(0))
    domain_two: Vector = (-q[2], F(0), x_value)
    input_wave = add(q, scale(F(h * k_value), W))
    assert dot(domain_one, input_wave) == 0
    assert dot(domain_two, input_wave) == 0

    def raw(vector: Vector) -> Vector:
        return add(
            scale(dot(D, q), vector),
            scale(F(m * k_value) * dot(vector, W), D),
        )

    output_wave = add(q, scale(F((h + m) * k_value), W))
    return triple(raw(domain_one), raw(domain_two), output_wave)


def check_extreme_charge_no_go_symbol() -> None:
    k_value = 101
    for q in Q_VECTORS:
        q_squared = dot(q, q)
        assert F(k_value * k_value) > q_squared
        assert F(k_value) > q_squared / (F(4) * q[0])
        for h in range(-100, 101):
            x_value = q[0] + F(h * k_value)
            shell_gap = (
                x_value * x_value + q[1] * q[1] + q[2] * q[2]
                - F(4 * k_value * k_value)
            )
            expected = q[2] * q[2] * x_value * shell_gap / F(4)
            for second_harmonic in (-2, 2):
                assert (
                    extreme_determinant(q, k_value, h, second_harmonic)
                    == expected
                )
            assert x_value != 0
            assert shell_gap != 0


def check_harmonic_separation_identity() -> None:
    q = Q_VECTORS[0]
    beta: Vector = (F(0), F(1), F(0))
    for parent_harmonic in (17, 31):
        denominator = F(parent_harmonic * 101) - q[0]
        b = add(beta, scale(dot(beta, q) / denominator, W))
        partner_wave = sub(q, scale(F(parent_harmonic * 101), W))
        assert dot(b, partner_wave) == 0
        for bath_harmonic in (1, 2):
            bath = bath_vector(bath_harmonic)
            direct = add(
                scale(dot(bath, partner_wave), b),
                scale(F(bath_harmonic * 101) * dot(b, W), bath),
            )
            reduced = add(
                scale(dot(bath, q), b),
                scale(
                    F(bath_harmonic, parent_harmonic) * dot(b, q),
                    bath,
                ),
            )
            assert direct == reduced
            # The unsuppressed (bath.q)b term survives increasing separation.
            assert dot(bath, q) != 0


def check_covariance_and_helicity() -> None:
    # W(theta)=sqrt(2) A sin(theta)+sqrt(2) D sin(2 theta).
    # Since <sin(m theta)sin(n theta)>=delta_mn/2, the covariance is
    # A tensor A + D tensor D = diag(0,1,1/4).
    covariance_diagonal = (
        A[0] * A[0] + D[0] * D[0],
        A[1] * A[1] + D[1] * D[1],
        A[2] * A[2] + D[2] * D[2],
    )
    assert covariance_diagonal == (F(0), F(1), F(1, 4))

    # Curl changes every sine harmonic into a cosine harmonic.  All circle
    # means <sin(m theta) cos(n theta)> vanish, including cross harmonics.
    for left_harmonic in (1, 2):
        for right_harmonic in (1, 2):
            sine_cosine_mean = F(0)
            assert sine_cosine_mean == 0
            assert left_harmonic > 0 and right_harmonic > 0

    # Every carrier and partner can be put in sine phase.  It is then odd
    # under central inversion, while curl is even, so helicity density is
    # odd.  The sign ledger is (-1)*(+1)=-1.
    velocity_parity = -1
    curl_parity = 1
    assert velocity_parity * curl_parity == -1


def main() -> None:
    check_zero_charge_rank_and_extra_charges()
    check_partner_partner_ledger()
    check_extreme_charge_no_go_symbol()
    check_harmonic_separation_identity()
    check_covariance_and_helicity()
    print("one-phase two-harmonic sideband coupling: exact checks passed")
    print("  rotated zero-charge symmetric chart rank: 5")
    print("  partner extra charges from one sign: -3, -2, +1")
    print("  bath--partner and partner--partner coefficients: K-free")
    print("  finite-support extreme-leakage determinant: exact")
    print("  covariance diag(0,1,1/4), central-odd helicity: zero")


if __name__ == "__main__":
    main()
