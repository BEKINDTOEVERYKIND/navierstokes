#!/usr/bin/env python3
"""Dependency-free exact checks for C179's planar passive reservoir.

This checker proves finite algebraic identities and ledgers only. It does
not prove terminal-profile preparation, a physical q-star, C125/RIGM,
BAFL, or a one-cell Navier--Stokes stage.
"""

from __future__ import annotations

from fractions import Fraction as F


V = tuple[F, F, F]
M = tuple[V, V, V]

ZERO: V = (F(0), F(0), F(0))
N: V = (F(1), F(1), F(1))
R1: V = (F(1), F(-1), F(0))
R2: V = (F(0), F(1), F(-1))
R3: V = (F(-1), F(0), F(1))
D: V = (F(-1), F(2), F(-1))


def add(a: V, b: V) -> V:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def sub(a: V, b: V) -> V:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def scale(c: F, a: V) -> V:
    return tuple(c * x for x in a)  # type: ignore[return-value]


def dot(a: V, b: V) -> F:
    return sum((x * y for x, y in zip(a, b)), F(0))


def cross(a: V, b: V) -> V:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm_sq(a: V) -> F:
    return dot(a, a)


def leray(k: V, a: V) -> V:
    return sub(a, scale(dot(k, a) / norm_sq(k), k))


def matvec(matrix: M, vector: V) -> V:
    return tuple(dot(row, vector) for row in matrix)  # type: ignore[return-value]


def outer(a: V, b: V) -> M:
    return tuple(tuple(x * y for y in b) for x in a)  # type: ignore[return-value]


def eye() -> M:
    return (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
    )


def matsub(a: M, b: M) -> M:
    return tuple(tuple(x - y for x, y in zip(ar, br))
                 for ar, br in zip(a, b))  # type: ignore[return-value]


def tangent_edge(p: V, g: V, h: V, a: V) -> V:
    """Real symbol P_{p+g}[(h.p)a+(a.g)h], omitting common -i."""
    k = add(p, g)
    return leray(k, add(scale(dot(h, p), a), scale(dot(a, g), h)))


def determinant_numerator(p: V, g: V, h: V, e1: V, e2: V) -> F:
    """|k| times oriented determinant; no square roots are needed."""
    k = add(p, g)
    return dot(k, cross(tangent_edge(p, g, h, e1),
                        tangent_edge(p, g, h, e2)))


def check_triangular_fourier_algebra() -> None:
    # Planar modes have wavevectors perpendicular to N. N-directed scalar
    # modes neither advect one another nor any planar mode.
    planar_modes = (R1, R2, D, scale(F(3), R3))
    for g in planar_modes:
        assert dot(g, N) == 0
        for q in planar_modes:
            assert dot(N, q) == 0
            assert scale(dot(N, q), N) == ZERO

    # A divergence-free planar coefficient v_q is a multiple of N x q.
    # It advects a scalar mode only in the N component.
    for q in planar_modes:
        v = cross(N, q)
        assert dot(v, q) == dot(v, N) == 0
        for g in planar_modes:
            scalar_advection = scale(dot(v, g), N)
            assert dot(scalar_advection, N) == norm_sq(N) * dot(v, g)


def source_basis_with_cross_p(p: V) -> tuple[V, V]:
    """Return rational e1,e2 in p-perp with e1 x e2=p."""
    trial: V = (F(1), F(0), F(0))
    if cross(p, trial) == ZERO:
        trial = (F(0), F(1), F(0))
    e1 = cross(p, trial)
    # p x e1 has norm factor |p|^2 relative to e1 x (p x e1)=|e1|^2 p.
    e2 = scale(F(1) / norm_sq(e1), cross(p, e1))
    assert dot(e1, p) == dot(e2, p) == dot(e1, e2) == 0
    assert cross(e1, e2) == p
    return e1, e2


def check_edge_determinant() -> None:
    # For e1 x e2=p, the note's determinant identity after clearing the
    # unit-area factors becomes exactly the target below.
    samples = (
        ((F(1), F(1), F(2)), R1),
        ((F(2), F(-1), F(3)), R2),
        ((F(3), F(1), F(-1)), D),
        ((F(2), F(3), F(4)), scale(F(2), R3)),
    )
    h = N
    for p, g in samples:
        assert dot(g, h) == 0
        e1, e2 = source_basis_with_cross_p(p)
        k = add(p, g)
        assert k != ZERO
        actual = determinant_numerator(p, g, h, e1, e2)
        expected = dot(h, p) ** 2 * (norm_sq(p) - norm_sq(g))
        assert actual == expected
        # Audit the note's corrected unprojected cross-product identity
        # B1 x B2 = alpha^2 (p-g) in this rational basis e1 x e2=p.
        # (The orthonormal formula has the common 1/|p| factor.)
        alpha = dot(h, p)
        b1 = add(scale(alpha, e1), scale(dot(e1, g), h))
        b2 = add(scale(alpha, e2), scale(dot(e2, g), h))
        assert cross(b1, b2) == scale(alpha * alpha, sub(p, g))

    # Exact normal-charge rank loss: h.p=0 kills the determinant.
    p = R2
    g = D
    e1, e2 = source_basis_with_cross_p(p)
    assert dot(N, p) == 0
    assert determinant_numerator(p, g, N, e1, e2) == 0

    # Exact equal-radius loss with nonzero normal charge and output.
    p = (F(1), F(1), F(0))
    g = R1
    assert norm_sq(p) == norm_sq(g)
    assert dot(N, p) != 0 and add(p, g) != ZERO
    e1, e2 = source_basis_with_cross_p(p)
    assert determinant_numerator(p, g, N, e1, e2) == 0


def check_compact_separated_family_and_normalization() -> None:
    # An explicit finite rational subfamily of a compact cone. This is an
    # exact determinant margin, not a certificate for the full C176 packet.
    gates = (R1, R2, R3, D)
    sources = (
        (F(3), F(2), F(4)),
        (F(4), F(1), F(3)),
        (F(2), F(4), F(3)),
    )
    margins = []
    for p in sources:
        for g in gates:
            k = add(p, g)
            if norm_sq(p) == norm_sq(g) or k == ZERO:
                continue
            numerator_sq = (
                dot(N, p) ** 4
                * (norm_sq(p) - norm_sq(g)) ** 2
            )
            denominator_sq = norm_sq(p) * norm_sq(k)
            margins.append(numerator_sq / denominator_sq)
    assert margins and min(margins) > F(1, 100)

    # q edges of derivative size q, each desired edge b/sqrt(q), require
    # per-mode scalar b/(q sqrt(q)). Squaring avoids irrational sqrt(q).
    for n in (2, 3, 5, 8):
        q = n**8
        b = F(1, n * n)
        theta_sq = b * b / (q**3)
        one_edge_sq = q * q * theta_sq
        assert one_edge_sq == b * b / q
        reservoir_l2_sq = q * theta_sq
        assert reservoir_l2_sq == b * b / (q * q)


def check_coset_translation_arithmetic() -> None:
    # R3 is planar and u=a0 R3-b0 N has u.R3=2a0 != 0.
    assert dot(N, R3) == 0
    assert dot(R3, R3) == 2

    # Pick one representative per R3 coset. Two independent quotient
    # functionals are N and D; both kill R3.
    assert dot(N, R3) == dot(D, R3) == 0
    for q in (2, 4, 6):
        reps = {
            add(scale(F(b), R1), scale(F(2 * q * c), R2))
            for b in range(q) for c in range(q)
        }
        shifts = [scale(F(a), R3)
                  for a in range(-q // 2, q // 2 + 1) if a]
        assert len(reps) == q * q and len(shifts) == q
        # Verify no two representatives differ by an R3 multiple.
        quotient_labels = {(dot(N, s), dot(D, s)) for s in reps}
        assert len(quotient_labels) == q * q
        outputs = {add(s, g) for s in reps for g in shifts}
        assert len(outputs) == q**3

    # Reality completion needs separation from negatives in addition to
    # half-lattice injectivity.  Give one explicit separated translated
    # family; the actual C176-slab selection remains a stated hypothesis.
    for q in (2, 4, 6):
        reps = {
            add(scale(F(20 * q + b), R1), scale(F(50 * q + 2 * q * c), R2))
            for b in range(q) for c in range(q)
        }
        shifts = {scale(F(a), R3)
                  for a in range(-q // 2, q // 2 + 1) if a}
        outputs = {add(s, g) for s in reps for g in shifts}
        negatives = {scale(F(-1), z) for z in reps | shifts | outputs}
        assert outputs.isdisjoint(negatives)


def check_static_vertical_shear_identity() -> None:
    # Unit normal n=e3, planar gradient c. G=I-t c tensor n.
    n: V = (F(0), F(0), F(1))
    gradients = ((F(2), F(-1), F(0)), (F(-3), F(4), F(0)))
    pairs = (
        ((F(2), F(1), F(3)), (F(1), F(1), F(-1))),
        ((F(-1), F(2), F(4)), (F(2), F(-1), F(1, 4))),
    )
    for gradient in gradients:
        assert dot(n, gradient) == 0
        for t in (F(-2), F(0), F(1, 3), F(5)):
            tc_outer: M = tuple(
                tuple(t * x for x in row) for row in outer(gradient, n)
            )  # type: ignore[assignment]
            gmat = matsub(eye(), tc_outer)
            for p, a0 in pairs:
                a = leray(p, a0)
                assert dot(p, a) == 0
                m = dot(n, p)
                assert m != 0
                gp = matvec(gmat, p)
                ga = matvec(gmat, a)
                invariant_left = sub(ga, scale(dot(n, a) / m, gp))
                invariant_right = sub(a, scale(dot(n, a) / m, p))
                assert invariant_left == invariant_right
                assert leray(gp, ga) == leray(gp, invariant_right)


def main() -> None:
    check_triangular_fourier_algebra()
    check_edge_determinant()
    check_compact_separated_family_and_normalization()
    check_coset_translation_arithmetic()
    check_static_vertical_shear_identity()
    print("C179 planar passive reservoir checks passed")


if __name__ == "__main__":
    main()
