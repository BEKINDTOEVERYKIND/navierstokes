#!/usr/bin/env python3
"""Exact bridge arithmetic for C185.

This checker verifies the rational implications from C159's certified cone
inequalities to the explicit e^(1/5) > 6/5 PDE-growth constant.  It does not
rerun C159's interval-Taylor certificate and cannot prove Shvydkoy's
essential-spectrum theorem; those are separately identified premises.
"""

from fractions import Fraction as F


def c159_cone_implication() -> None:
    # C159's rounded interval conclusions are strict.  Here we check that
    # their lower bounds imply B w > (1/5)w for w=(1,3/20).
    w1 = F(1)
    w2 = F(3, 20)
    rate = F(1, 5)

    # First C159 row: B11 + (3/20) B12 > 7/10.
    first_lower = F(7, 10)
    assert first_lower > rate * w1

    # Second C159 row is recorded after division by w2:
    # (20/3) B21 + B22 > 1/5.  Multiplication by w2 gives
    # B21 + (3/20)B22 > (1/5)(3/20).
    normalized_second_lower = F(1, 5)
    second_lower = w2 * normalized_second_lower
    assert second_lower == rate * w2

    # Strictly positive off-diagonal lower bounds make B Metzler.
    assert F(32) > 0
    assert F(9, 10) > 0


def explicit_exponential_constant() -> None:
    # The first three nonnegative Taylor terms already give a strict bound.
    t = F(1, 5)
    exp_partial = F(1) + t + t * t / 2
    assert exp_partial == F(61, 50)
    assert exp_partial > F(6, 5)


def returning_orthogonal_frame_norm() -> None:
    # If coefficients in an orthogonal frame grow componentwise by c>1,
    # the physical squared norm grows by c^2, independently of the two
    # positive frame weights.
    c = F(6, 5)
    for frame_weight_1, frame_weight_2 in (
        (F(1), F(1)),
        (F(7, 5), F(11, 3)),
        (F(1, 101), F(103)),
    ):
        for z1, z2 in ((F(1), F(3, 20)), (F(5), F(2))):
            norm_sq = frame_weight_1 * z1 * z1 + frame_weight_2 * z2 * z2
            grown_sq = (
                frame_weight_1 * (c * z1) ** 2
                + frame_weight_2 * (c * z2) ** 2
            )
            assert grown_sq == c * c * norm_sq


def iteration_ledger() -> None:
    # The theorem gives e^(n/5);  (6/5)^n is the advertised explicit lower
    # constant.  Check exact monotonic iteration on representative n.
    base = F(6, 5)
    previous = F(1)
    for n in range(1, 65):
        current = base**n
        assert current == base * previous
        assert current > previous
        previous = current


def main() -> None:
    c159_cone_implication()
    explicit_exponential_constant()
    returning_orthogonal_frame_norm()
    iteration_ledger()
    print("C185 infinite-dimensional operator-growth bridge arithmetic passed")
    print("BOUNDARY: C159 certificate and Shvydkoy abstract inclusion are premises")
    print("BOUNDARY: essential-radius identification is citation-held by C189")
    print("BOUNDARY: no viscous UVSR, nonlinear stage, or singularity is proved")


if __name__ == "__main__":
    main()
