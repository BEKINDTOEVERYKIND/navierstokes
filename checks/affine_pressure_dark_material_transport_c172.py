#!/usr/bin/env python3
"""Dependency-free exact arithmetic for C172.

This checker verifies:

* the oriented compression/cofactor identity in a non-coordinate rational
  orthonormal frame;
* dark-fiber dimensions for invertible, nonsymmetric rank-two, symmetric
  rank-two, and rank-one matrices;
* the compact rank-one curl-kernel algebra and absence of Piola gain;
* brightness, including the full factor 2, of the C140 start-frame
  child/wake polarizations for the normalized C142 selector;
* the linear constraints making Piola (up to a scalar) the unique
  universal zero-order divergence-preserving material push-forward;
* the universal pressure-dark matrix no-go, scalar-Piola rank ledger, and
  inconsistency of any constant symbol with the Kelvin generator.

The measure-zero/L2 localization argument and the open-set polynomial
continuation in the finite-order no-go are proofs in the note, not
numerical claims.  This file does not certify MCKC, LCE, BAFL, or an
unforced stage.
"""

from __future__ import annotations

from fractions import Fraction as F


Vec = tuple[F, F, F]
Mat = tuple[Vec, Vec, Vec]

ZERO: Vec = (F(0), F(0), F(0))
IDENTITY: Mat = (
    (F(1), F(0), F(0)),
    (F(0), F(1), F(0)),
    (F(0), F(0), F(1)),
)


def vec(values) -> Vec:
    return tuple(F(value) for value in values)  # type: ignore[return-value]


def vadd(a: Vec, b: Vec) -> Vec:
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def vsub(a: Vec, b: Vec) -> Vec:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def vscale(c: F, a: Vec) -> Vec:
    return tuple(c * a[i] for i in range(3))  # type: ignore[return-value]


def dot(a: Vec, b: Vec) -> F:
    return sum((a[i] * b[i] for i in range(3)), F(0))


def cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm_sq(a: Vec) -> F:
    return dot(a, a)


def mvec(a: Mat, x: Vec) -> Vec:
    return tuple(dot(a[i], x) for i in range(3))  # type: ignore[return-value]


def transpose(a: Mat) -> Mat:
    return tuple(
        tuple(a[j][i] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mmul(a: Mat, b: Mat) -> Mat:
    bt = transpose(b)
    return tuple(
        tuple(dot(a[i], bt[j]) for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def madd(a: Mat, b: Mat) -> Mat:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mscale(c: F, a: Mat) -> Mat:
    return tuple(
        tuple(c * a[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def outer(a: Vec, b: Vec) -> Mat:
    return tuple(
        tuple(a[i] * b[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def det(a: Mat) -> F:
    return dot(a[0], cross(a[1], a[2]))


def minor2(a: Mat, row: int, column: int) -> F:
    rows = [index for index in range(3) if index != row]
    columns = [index for index in range(3) if index != column]
    return (
        a[rows[0]][columns[0]] * a[rows[1]][columns[1]]
        - a[rows[0]][columns[1]] * a[rows[1]][columns[0]]
    )


def cofactor(a: Mat) -> Mat:
    return tuple(
        tuple(((-1) ** (i + j)) * minor2(a, i, j) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def adjugate(a: Mat) -> Mat:
    return transpose(cofactor(a))


def matrix_rank(rows: list[list[F]]) -> int:
    if not rows:
        return 0
    work = [row[:] for row in rows]
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(row_count):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    work[row][j] - factor * work[pivot_row][j]
                    for j in range(column_count)
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def dark_dimension(a: Mat, k: Vec) -> int:
    """Dimension of {u: k.u=0 and k cross A u=0}."""
    # The cross-product matrix C_k satisfies C_k y = k cross y.
    cross_k: Mat = (
        (F(0), -k[2], k[1]),
        (k[2], F(0), -k[0]),
        (-k[1], k[0], F(0)),
    )
    cross_a = mmul(cross_k, a)
    constraints = [list(k)] + [list(row) for row in cross_a]
    return 3 - matrix_rank(constraints)


def project(k: Vec, value: Vec) -> Vec:
    return vsub(value, vscale(dot(k, value) / norm_sq(k), k))


def compression_det(a: Mat, e: Vec, f: Vec) -> F:
    """Determinant of P_k A on an orthonormal (e,f) plane basis."""
    return (
        dot(e, mvec(a, e)) * dot(f, mvec(a, f))
        - dot(e, mvec(a, f)) * dot(f, mvec(a, e))
    )


def kelvin_generator(a: Mat, k: Vec, u: Vec) -> Vec:
    assert dot(k, u) == 0
    a_u = mvec(a, u)
    return vadd(
        vscale(F(-1), a_u),
        vscale(F(2) * dot(k, a_u) / norm_sq(k), k),
    )


def check_compression_and_dimensions() -> None:
    # The cofactor area identity is the coordinate-free source of (1.3).
    samples = (
        (
            (
                (F(2), F(1), F(0)),
                (F(-1), F(3), F(2)),
                (F(1), F(0), F(-5)),
            ),
            vec((1, 2, -1)),
            vec((0, 1, 3)),
        ),
        (
            (
                (F(0), F(1), F(0)),
                (F(0), F(0), F(1)),
                (F(0), F(0), F(0)),
            ),
            vec((2, -1, 1)),
            vec((1, 0, 4)),
        ),
    )
    for a, e, f in samples:
        assert cross(mvec(a, e), mvec(a, f)) == mvec(
            cofactor(a), cross(e, f)
        )
        # Adjugate and cofactor have the same quadratic form.
        for k in (vec((1, 0, 0)), vec((1, 2, 3)), cross(e, f)):
            assert dot(k, mvec(cofactor(a), k)) == dot(
                k, mvec(adjugate(a), k)
            )

    # Directly test the two-dimensional compression determinant in a
    # non-coordinate rational oriented orthonormal frame.  This catches a
    # cofactor/adjugate transpose error that diagonal samples would miss.
    e = vec((F(1, 9), F(8, 9), F(-4, 9)))
    f = vec((F(8, 9), F(1, 9), F(4, 9)))
    k_hat = vec((F(4, 9), F(-4, 9), F(-7, 9)))
    assert norm_sq(e) == norm_sq(f) == norm_sq(k_hat) == 1
    assert dot(e, f) == dot(e, k_hat) == dot(f, k_hat) == 0
    assert cross(e, f) == k_hat
    nonsymmetric_full: Mat = (
        (F(2), F(-1), F(3)),
        (F(0), F(4), F(1)),
        (F(5), F(2), F(-2)),
    )
    direct_det = compression_det(nonsymmetric_full, e, f)
    assert direct_det == dot(
        k_hat, mvec(cofactor(nonsymmetric_full), k_hat)
    )
    assert direct_det == dot(
        k_hat, mvec(adjugate(nonsymmetric_full), k_hat)
    )

    # Invertible symmetric trace-free example.
    invertible: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(-2)),
    )
    inverse: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(-1, 2)),
    )
    assert det(invertible) == -2
    cone_k = vec((1, 1, 2))
    off_cone_k = vec((1, 0, 1))
    assert dot(cone_k, mvec(inverse, cone_k)) == 0
    assert dark_dimension(invertible, cone_k) == 1
    assert dark_dimension(invertible, off_cone_k) == 0
    cone_u = mvec(inverse, cone_k)
    assert dot(cone_k, cone_u) == 0
    assert cross(mvec(invertible, cone_u), cone_k) == ZERO

    # Nonsymmetric rank-two example with different right and left kernels.
    # Standard adj(A) is c r tensor ell, not its transpose.
    nonsymmetric_rank_two: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(-1), F(0)),
        (F(1), F(0), F(0)),
    )
    right_kernel = vec((0, 0, 1))
    left_kernel = vec((-1, 0, 1))
    assert mvec(nonsymmetric_rank_two, right_kernel) == ZERO
    assert mvec(transpose(nonsymmetric_rank_two), left_kernel) == ZERO
    assert adjugate(nonsymmetric_rank_two) == mscale(
        F(-1), outer(right_kernel, left_kernel)
    )
    assert cofactor(nonsymmetric_rank_two) == mscale(
        F(-1), outer(left_kernel, right_kernel)
    )
    for i in range(-2, 3):
        for j in range(-2, 3):
            for ell_index in range(-2, 3):
                k = vec((i, j, ell_index))
                if k == ZERO:
                    continue
                quadratic = dot(
                    k, mvec(adjugate(nonsymmetric_rank_two), k)
                )
                on_plane_union = (
                    dot(k, right_kernel) == 0
                    or dot(k, left_kernel) == 0
                )
                assert (quadratic == 0) == on_plane_union
                assert (dark_dimension(nonsymmetric_rank_two, k) > 0) == (
                    quadratic == 0
                )

    # Symmetric rank-two selector eigenbasis diag(0,1,-1).
    rank_two: Mat = (
        (F(0), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(-1)),
    )
    assert matrix_rank([list(row) for row in rank_two]) == 2
    assert dark_dimension(rank_two, vec((1, 0, 0))) == 0
    assert dark_dimension(rank_two, vec((0, 1, 0))) == 1
    assert dark_dimension(rank_two, vec((0, 1, 1))) == 2
    assert dark_dimension(rank_two, vec((0, 1, -1))) == 2
    # Its adjugate quadratic is -k_1^2.
    adj = adjugate(rank_two)
    for k in (vec((2, 3, 5)), vec((0, 7, -4))):
        assert dot(k, mvec(adj, k)) == -(k[0] ** 2)

    # Rank-one trace-free shear: dimensions are two on the a and b rays,
    # and one at a generic frequency.
    a_vec = vec((1, 0, 0))
    b_vec = vec((0, 1, 0))
    rank_one = outer(a_vec, b_vec)
    assert dot(a_vec, b_vec) == 0
    assert matrix_rank([list(row) for row in rank_one]) == 1
    assert dark_dimension(rank_one, a_vec) == 2
    assert dark_dimension(rank_one, b_vec) == 2
    assert dark_dimension(rank_one, vec((1, 1, 1))) == 1


def levi(i: int, j: int, k: int) -> int:
    if len({i, j, k}) < 3:
        return 0
    return 1 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1


def check_rank_one_compact_curl_kernel() -> None:
    a = vec((1, 0, 0))
    b = vec((0, 1, 0))
    shear = outer(a, b)
    assert mmul(shear, shear) == mscale(F(0), IDENTITY)

    # For an arbitrary exact gradient g, U=g cross b is in ker A.
    for gradient in (vec((2, -3, 5)), vec((0, 7, -1))):
        u = cross(gradient, b)
        assert dot(b, u) == 0
        assert mvec(shear, u) == ZERO
        # F(t)=I+tA leaves this kernel amplitude fixed.
        for time in (F(-2), F(0), F(7, 3)):
            flow = madd(IDENTITY, mscale(time, shear))
            inverse_flow = madd(IDENTITY, mscale(-time, shear))
            assert mvec(flow, u) == u
            assert det(flow) == 1
            assert mmul(flow, inverse_flow) == IDENTITY

    # The compact kernel exists for a general rank-one matrix as well;
    # nilpotence and the linear flow formula specifically use trace zero.
    general_rank_one = outer(vec((1, 1, 0)), vec((1, 0, 0)))
    assert sum((general_rank_one[i][i] for i in range(3)), F(0)) == 1
    for gradient in (vec((2, -3, 5)), vec((0, 7, -1))):
        u = cross(gradient, vec((1, 0, 0)))
        assert mvec(general_rank_one, u) == ZERO

    # div(grad psi cross b)=epsilon_{i j l} H_{i j} b_l=0 for every
    # symmetric Hessian H.  Check exact generic symmetric samples.
    hessians = (
        (
            (F(1), F(2), F(-1)),
            (F(2), F(3), F(4)),
            (F(-1), F(4), F(5)),
        ),
        (
            (F(0), F(-3), F(2)),
            (F(-3), F(7), F(1)),
            (F(2), F(1), F(-4)),
        ),
    )
    for hessian in hessians:
        divergence = F(0)
        for i in range(3):
            for j in range(3):
                for ell in range(3):
                    divergence += levi(i, j, ell) * hessian[i][j] * b[ell]
        assert divergence == 0


def check_c140_selector_brightness_and_no_gain_if_all_dark() -> None:
    n = vec((1, 1, 1))
    r = vec((1, 0, -1))
    d = vec((-1, 2, -1))
    e1 = vec((2, -1, -1))
    e2 = vec((1, 1, -2))
    selector = madd(outer(d, n), outer(n, d))

    assert transpose(selector) == selector
    assert sum((selector[i][i] for i in range(3)), F(0)) == 0
    assert matrix_rank([list(row) for row in selector]) == 2
    assert dot(n, d) == dot(n, r) == dot(d, r) == 0
    assert e1 == vsub(vscale(F(3, 2), r), vscale(F(1, 2), d))
    assert e2 == vadd(vscale(F(3, 2), r), vscale(F(1, 2), d))
    assert mvec(selector, r) == ZERO
    assert mvec(selector, n) == vscale(F(3), d)
    assert mvec(selector, d) == vscale(F(6), n)
    selector_sq = mmul(selector, selector)
    assert mvec(selector_sq, n) == vscale(F(18), n)
    assert mvec(selector_sq, d) == vscale(F(18), d)
    assert mvec(selector_sq, r) == ZERO

    affine_n = mvec(selector, n)
    assert norm_sq(project(r, affine_n)) == 54
    assert norm_sq(project(e1, affine_n)) == F(81, 2)
    assert norm_sq(project(e2, affine_n)) == F(81, 2)
    # C171's actual affine cross factor is 2 A U, so squared sizes are
    # four times the one-A_* ledger above.
    assert norm_sq(project(r, vscale(F(2), affine_n))) == 216
    assert norm_sq(project(e1, vscale(F(2), affine_n))) == 162
    assert norm_sq(project(e2, vscale(F(2), affine_n))) == 162
    # With A=-rho' A_*/(3 sqrt(2)), these become 12 rho'^2 and
    # 9 rho'^2 without introducing irrational arithmetic.
    rho_prime_sq = F(25, 49)
    sigma_sq = rho_prime_sq / 18
    assert 216 * sigma_sq == 12 * rho_prime_sq
    assert 162 * sigma_sq == 9 * rho_prime_sq

    # A vector parallel to both nonparallel r and e1 must be zero.  Encode
    # cross(w,r)=cross(w,e1)=0 as a rank-three linear system for w.
    def cross_constraints(k: Vec) -> list[list[F]]:
        return [
            [F(0), -k[2], k[1]],
            [k[2], F(0), -k[0]],
            [-k[1], k[0], F(0)],
        ]

    constraints = cross_constraints(r) + cross_constraints(e1)
    assert matrix_rank(constraints) == 3

    # If A N=0 throughout an affine interval, the Kelvin right side for
    # a=N is exactly zero and
    # d/dt(k.N)=-(k.AN)=0.  Verify for generic trace-free examples with
    # N in the kernel.
    plane_strain = madd(
        mscale(F(1, 2), outer(r, r)),
        mscale(F(-1, 6), outer(d, d)),
    )
    kernel_examples: tuple[Mat, ...] = (
        plane_strain,
        (
            (F(0), F(1), F(-1)),
            (F(-1), F(0), F(1)),
            (F(1), F(-1), F(0)),
        ),
    )
    for a in kernel_examples:
        assert sum((a[i][i] for i in range(3)), F(0)) == 0
        assert mvec(a, n) == ZERO
        for k in (r, e1, e2):
            assert dot(k, n) == 0
            a_n = mvec(a, n)
            kelvin_rhs = vadd(
                vscale(F(-1), a_n),
                vscale(
                    F(2) * dot(k, a_n) / norm_sq(k),
                    k,
                ),
            )
            assert kelvin_rhs == ZERO
            assert -dot(k, a_n) == 0


def check_zero_order_uniqueness_and_universal_no_go() -> None:
    # Unknown B is row-major.  The divergence conditions xi.Bv=0 for all
    # v perpendicular to xi have an eight-dimensional row space, leaving
    # precisely the scalar identity line.
    constraints: list[list[F]] = []

    def add_constraint(xi: Vec, value: Vec) -> None:
        constraints.append(
            [xi[i] * value[j] for i in range(3) for j in range(3)]
        )

    axes = (vec((1, 0, 0)), vec((0, 1, 0)), vec((0, 0, 1)))
    for index, xi in enumerate(axes):
        for other, value in enumerate(axes):
            if other != index:
                add_constraint(xi, value)
    add_constraint(vec((1, 1, 0)), vec((1, -1, 0)))
    add_constraint(vec((1, 0, 1)), vec((1, 0, -1)))
    assert matrix_rank(constraints) == 8
    identity_flat = [IDENTITY[i][j] for i in range(3) for j in range(3)]
    assert all(
        sum((row[index] * identity_flat[index] for index in range(9)), F(0))
        == 0
        for row in constraints
    )

    # Reinsert a genuinely nontrivial volume-preserving affine map.  The
    # unique B=cI condition is exactly M=cF, and a perturbed M is detected
    # by an explicit transverse plane wave.
    flow: Mat = (
        (F(1), F(1), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    )
    inverse_flow: Mat = (
        (F(1), F(-1), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    )
    scalar = F(7, 4)
    multiplier = mscale(scalar, flow)
    assert det(flow) == 1
    assert mmul(inverse_flow, flow) == IDENTITY
    assert mmul(inverse_flow, multiplier) == mscale(scalar, IDENTITY)
    plane_wave_pairs = (
        (vec((1, 0, 0)), vec((0, 1, 0))),
        (vec((1, 1, 0)), vec((1, -1, 0))),
        (vec((1, 2, -1)), vec((2, -1, 0))),
    )
    for xi, value in plane_wave_pairs:
        assert dot(xi, value) == 0
        assert dot(xi, mvec(mmul(inverse_flow, multiplier), value)) == 0
    wrong_multiplier = madd(multiplier, outer(axes[0], axes[0]))
    xi = vec((1, 1, 0))
    value = vec((1, -1, 0))
    assert dot(xi, value) == 0
    assert dot(
        xi, mvec(mmul(inverse_flow, wrong_multiplier), value)
    ) != 0

    # If a matrix D maps every k-perpendicular plane into span(k), then
    # all its entries vanish.  Coordinate-axis tests already give rank 9.
    dark_constraints: list[list[F]] = []
    for k_index in range(3):
        perpendicular = [index for index in range(3) if index != k_index]
        for output_index in perpendicular:
            for input_index in perpendicular:
                row = [F(0)] * 9
                row[3 * output_index + input_index] = F(1)
                dark_constraints.append(row)
    assert matrix_rank(dark_constraints) == 9

    # Scalar Piola on eigenvalues (0,lambda,-lambda): the roots that make
    # one eigenvalue of 2A+gamma I vanish are distinct, so rank is never
    # below two.  This exhausts every real gamma.
    lam = F(7, 5)
    zero_locations = {F(0), -2 * lam, 2 * lam}
    assert len(zero_locations) == 3
    for gamma in tuple(zero_locations) + (F(1), F(-3, 2)):
        eigenvalues = (gamma, gamma + 2 * lam, gamma - 2 * lam)
        zero_count = sum(value == 0 for value in eigenvalues)
        if gamma in zero_locations:
            assert zero_count == 1
        else:
            assert zero_count == 0
        assert 3 - zero_count >= 2

    # A constant symbol cannot equal the degree-zero Kelvin generator on
    # transverse data for a nonzero trace-free A.  Axis fibers force
    # C=-A; one tilted fiber makes the augmented exact system inconsistent.
    affine: Mat = (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(-2)),
    )
    fiber_samples = [
        (axes[0], axes[1]),
        (axes[0], axes[2]),
        (axes[1], axes[0]),
        (axes[1], axes[2]),
        (axes[2], axes[0]),
        (axes[2], axes[1]),
    ]
    tilted_k = vec((1, 0, 1))
    tilted_u = vec((1, 0, -1))
    assert kelvin_generator(affine, vscale(F(3), tilted_k), tilted_u) == (
        kelvin_generator(affine, tilted_k, tilted_u)
    )
    fiber_samples.append((tilted_k, tilted_u))

    coefficient_rows: list[list[F]] = []
    augmented_rows: list[list[F]] = []
    for k, u in fiber_samples:
        rhs = kelvin_generator(affine, k, u)
        for output_index in range(3):
            row = [F(0)] * 9
            for input_index in range(3):
                row[3 * output_index + input_index] = u[input_index]
            coefficient_rows.append(row)
            augmented_rows.append(row + [rhs[output_index]])
    assert matrix_rank(coefficient_rows) == 9
    assert matrix_rank(augmented_rows) == 10


def main() -> None:
    check_compression_and_dimensions()
    check_rank_one_compact_curl_kernel()
    check_c140_selector_brightness_and_no_gain_if_all_dark()
    check_zero_order_uniqueness_and_universal_no_go()
    print("PASS C172: exact affine pressure-dark fiber dimensions")
    print("PASS C172: rank-one compact curl exception has no Piola gain")
    print("PASS C172: C140 start-frame fibers are bright, including factor 2")
    print("PASS C172: scalar Piola is the only universal affine multiplier")
    print("PASS C172: no constant symbol matches the open-fiber generator")
    print("OPEN: MCKC, LCE, BAFL, and the full unforced stage")


if __name__ == "__main__":
    main()
