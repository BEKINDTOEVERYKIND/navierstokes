#!/usr/bin/env python3
"""Dependency-free exact checks for C173's paired Piola multipole.

The checker verifies the curl sign and a nonzero isotropic degree-three
moment, the two-radius weights, A2 strain/reality geometry, the triple-shell
polynomial and its two-jet zeros, material-label versus Eulerian shells,
affine-core preservation, the fixed condition constant, the
first-curvature/exterior-pressure ledger, the fixed-normalization boundary,
and the q,b/order exponents with a sharp univariate Hermite example.
It does not certify a broad-band Duhamel estimate, a material chart, BAFL,
an unforced stage, or a Navier--Stokes singularity.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log


Vec = tuple[F, F, F]
Mat = tuple[Vec, Vec, Vec]


def dot(a: Vec, b: Vec) -> F:
    return sum((a[i] * b[i] for i in range(3)), F(0))


def mvec(a: Mat, x: Vec) -> Vec:
    return tuple(dot(a[i], x) for i in range(3))  # type: ignore[return-value]


def outer(a: Vec, b: Vec) -> Mat:
    return tuple(
        tuple(a[i] * b[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def madd(a: Mat, b: Mat) -> Mat:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def levi(i: int, j: int, k: int) -> F:
    if len({i, j, k}) < 3:
        return F(0)
    return F(1) if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else F(-1)


def delta(i: int, j: int) -> F:
    return F(1) if i == j else F(0)


def poly_mul(a: list[F], b: list[F]) -> list[F]:
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            result[i + j] += left * right
    return result


def poly_diff(a: list[F]) -> list[F]:
    return [F(i) * a[i] for i in range(1, len(a))]


def poly_eval(a: list[F], x: F) -> F:
    result = F(0)
    for coefficient in reversed(a):
        result = result * x + coefficient
    return result


def check_a2_core_strain() -> None:
    n: Vec = (F(1), F(1), F(1))
    d: Vec = (F(-1), F(2), F(-1))
    kc: Vec = (F(1), F(0), F(-1))
    e1: Vec = (F(2), F(-1), F(-1))
    e2: Vec = (F(1), F(1), F(-2))

    # A rational multiple of the symmetric C142 selector B:
    # S=d tensor N + N tensor d. It is symmetric, trace free, and kills
    # the child covector kc.
    s = madd(outer(d, n), outer(n, d))
    assert all(s[i][j] == s[j][i] for i in range(3) for j in range(3))
    assert sum((s[i][i] for i in range(3)), F(0)) == 2 * dot(d, n) == 0
    assert mvec(s, kc) == (F(0), F(0), F(0))
    assert dot(kc, kc) == 2
    assert dot(e1, e1) == dot(e2, e2) == 6
    assert dot(n, kc) == dot(n, e1) == dot(n, e2) == 0

    # The radial integral of the affine potential is proportional to
    # epsilon_npq S_pq, which vanishes componentwise precisely by symmetry.
    axial_skew = tuple(
        sum(
            (levi(axis, p, q) * s[p][q] for p in range(3) for q in range(3)),
            F(0),
        )
        for axis in range(3)
    )
    assert axial_skew == (F(0), F(0), F(0))

    # Re-derive the curl sign coefficientwise from
    # A_n=(1/3) eps_npq S_ps x_s x_q.  This avoids merely assuming the
    # identity curl((Sx)x x/3)=Sx-(tr S)x/3.
    curl_coeff = tuple(
        tuple(
            sum(
                (
                    levi(j, m, axis)
                    * levi(axis, p, q)
                    * s[p][label]
                    * (
                        delta(m, label) * delta(q, coordinate)
                        + delta(m, q) * delta(label, coordinate)
                    )
                    / 3
                    for m in range(3)
                    for axis in range(3)
                    for p in range(3)
                    for q in range(3)
                    for label in range(3)
                ),
                F(0),
            )
            for coordinate in range(3)
        )
        for j in range(3)
    )
    assert curl_coeff == s

    # Use the exact isotropic fourth-moment tensor (with its common positive
    # radial constant normalized to one) to compute a degree-three velocity
    # moment.  It is genuinely allowed/nonzero for this symmetric S, so the
    # paired R^7 cancellation is substantive rather than parity bookkeeping.
    a_second = tuple(
        tuple(
            tuple(
                sum(
                    (
                        levi(component, p, q)
                        * s[p][label]
                        * (
                            delta(left, right) * delta(label, q)
                            + delta(left, label) * delta(right, q)
                            + delta(left, q) * delta(right, label)
                        )
                        / 3
                        for p in range(3)
                        for q in range(3)
                        for label in range(3)
                    ),
                    F(0),
                )
                for right in range(3)
            )
            for left in range(3)
        )
        for component in range(3)
    )
    degree_three = {}
    for component in range(3):
        for first in range(3):
            for second in range(3):
                for third in range(3):
                    value = F(0)
                    for curl_component in range(3):
                        value -= levi(component, first, curl_component) * a_second[curl_component][second][third]
                        value -= levi(component, second, curl_component) * a_second[curl_component][first][third]
                        value -= levi(component, third, curl_component) * a_second[curl_component][first][second]
                    degree_three[(component, first, second, third)] = value
    assert any(value != 0 for value in degree_three.values())


def check_pair_weights_and_moment_orders() -> None:
    alpha = F(128, 127)
    beta = F(-1, 127)

    # Re-derive the unique two-radius weights from core preservation and
    # degree-three cancellation: alpha+beta=1, alpha+2^(3+4) beta=0.
    derived_beta = -F(1, 2 ** 7 - 1)
    derived_alpha = 1 - derived_beta
    assert (alpha, beta) == (derived_alpha, derived_beta)

    # Both collars equal Sx in the common core, so their strains add.
    assert alpha + beta == 1
    assert abs(alpha) + abs(beta) == F(129, 127)

    # A degree-d velocity moment of w_R scales as R^(d+4).
    scale = lambda degree: alpha + beta * F(2) ** (degree + 4)
    assert scale(3) == 0
    assert scale(5) == F(-384, 127)

    # Structural zeros: even degrees by odd parity, degree one by radiality
    # plus symmetry. Pairing adds the degree-three zero, so every moment
    # through degree four vanishes and the first allowed degree is five.
    structural_single_zeros = {0, 1, 2}
    parity_zeros = {0, 2, 4}
    paired_zeros = structural_single_zeros | parity_zeros | {3}
    assert paired_zeros == {0, 1, 2, 3, 4}

    # In dimension three the Riesz/pressure-gradient kernel has exterior
    # order 3; if ghat starts at Fourier order p, the first exterior term is
    # d^-(3+p).  Check paired versus single orders at affine and first
    # curvature level.
    dimension = 3
    constant_g_order = 5
    single_constant_g_order = 3
    curvature_g_order = constant_g_order - 1
    single_curvature_g_order = single_constant_g_order - 1
    assert dimension + constant_g_order == 8
    assert dimension + single_constant_g_order == 6
    assert dimension + curvature_g_order == 7
    assert dimension + single_curvature_g_order == 5

    # At low normalized frequency zeta~epsilon, every first-curvature term
    # retains total exponent five: epsilon*d(w) and
    # epsilon*zeta*d^2(w).
    amplitude_curvature_exponent = 1 + (constant_g_order - 1)
    phase_curvature_exponent = 1 + 1 + (constant_g_order - 2)
    assert amplitude_curvature_exponent == phase_curvature_exponent == 5
    single_amplitude_curvature_exponent = 1 + (single_constant_g_order - 1)
    single_phase_curvature_exponent = 1 + 1 + (single_constant_g_order - 2)
    assert single_amplitude_curvature_exponent == single_phase_curvature_exponent == 3
    assert amplitude_curvature_exponent - single_amplitude_curvature_exponent == 2


def build_a2_notch() -> list[F]:
    q = [F(1)]
    for radius_squared in (F(2), F(6)):
        factor = [F(1), -F(1, radius_squared)]
        for _ in range(3):
            q = poly_mul(q, factor)
    return q


def check_triple_notch_and_condition() -> None:
    q = build_a2_notch()
    expected = [
        F(1),
        F(-2),
        F(19, 12),
        F(-17, 27),
        F(19, 144),
        F(-1, 72),
        F(1, 1728),
    ]
    assert q == expected
    assert len(q) - 1 == 6
    assert 2 * (len(q) - 1) == 12
    assert q[0] == 1
    assert sum((abs(entry) for entry in q), F(0)) == F(343, 64)

    # On the open core w=Sx.  Verify exactly that its componentwise
    # Laplacian is zero, so every positive power of -r^2 Delta vanishes and
    # Q(0)=1 preserves the full affine slope.
    core_s: Mat = (
        (F(1), F(2), F(0)),
        (F(2), F(-1), F(1)),
        (F(0), F(1), F(0)),
    )
    core_point: Vec = (F(2), F(-3), F(5))
    core_value = mvec(core_s, core_point)
    discrete_laplacian = [F(0), F(0), F(0)]
    for axis in range(3):
        plus = tuple(
            core_point[i] + (1 if i == axis else 0) for i in range(3)
        )
        minus = tuple(
            core_point[i] - (1 if i == axis else 0) for i in range(3)
        )
        plus_value = mvec(core_s, plus)  # type: ignore[arg-type]
        minus_value = mvec(core_s, minus)  # type: ignore[arg-type]
        for component in range(3):
            discrete_laplacian[component] += (
                plus_value[component]
                - 2 * core_value[component]
                + minus_value[component]
            )
    assert tuple(discrete_laplacian) == (F(0), F(0), F(0))
    assert tuple(q[0] * entry for entry in core_value) == core_value

    first = poly_diff(q)
    second = poly_diff(first)
    third = poly_diff(second)
    for radius_squared in (F(2), F(6)):
        assert poly_eval(q, radius_squared) == 0
        assert poly_eval(first, radius_squared) == 0
        assert poly_eval(second, radius_squared) == 0
        assert poly_eval(third, radius_squared) != 0

    # For q(xi)=Q(|xi|^2), q, grad q, and Hess q vanish on both shells:
    # grad_i q=2 xi_i Q'; Hess_ij q=4 xi_i xi_j Q''+2 delta_ij Q'.
    shell_vectors: tuple[Vec, ...] = (
        (F(1), F(0), F(-1)),
        (F(2), F(-1), F(-1)),
        (F(1), F(1), F(-2)),
    )
    for xi in shell_vectors:
        radius_squared = dot(xi, xi)
        q0 = poly_eval(q, radius_squared)
        q1 = poly_eval(first, radius_squared)
        q2 = poly_eval(second, radius_squared)
        gradient = tuple(2 * entry * q1 for entry in xi)
        hessian = tuple(
            tuple(
                4 * xi[i] * xi[j] * q2 + (2 * q1 if i == j else 0)
                for j in range(3)
            )
            for i in range(3)
        )
        assert q0 == 0
        assert gradient == (F(0), F(0), F(0))
        assert all(hessian[i][j] == 0 for i in range(3) for j in range(3))

        # Product-rule Fourier jets of Q*w vanish through order two for
        # arbitrary underlying value/first/second jets of w.
        arbitrary_value: Vec = (F(2), F(-3), F(5))
        arbitrary_first = (
            (F(1), F(2), F(3)),
            (F(-2), F(4), F(1)),
            (F(5), F(-1), F(2)),
        )
        filtered_value = tuple(q0 * entry for entry in arbitrary_value)
        assert filtered_value == (F(0), F(0), F(0))
        for axis in range(3):
            filtered_first = tuple(
                gradient[axis] * arbitrary_value[j]
                + q0 * arbitrary_first[axis][j]
                for j in range(3)
            )
            assert filtered_first == (F(0), F(0), F(0))
        # Every second product-rule term contains q0, one gradient entry,
        # or one Hessian entry, all of which are exactly zero.
        for left in range(3):
            for right in range(3):
                assert hessian[left][right] == 0
                assert gradient[left] == gradient[right] == q0 == 0

    # The notch acts on material labels.  Under a nonorthogonal det-one
    # affine chart, physical k=(1,-1,-1) pulls back to kc and is killed,
    # whereas the fixed Eulerian coordinate kc pulls back to a radius-3
    # label and is not killed.
    f_transpose: Mat = (
        (F(1), F(0), F(0)),
        (F(1), F(1), F(0)),
        (F(0), F(0), F(1)),
    )
    physical_k = (F(1), F(-1), F(-1))
    material_kc = shell_vectors[0]
    assert mvec(f_transpose, physical_k) == material_kc
    assert poly_eval(q, dot(material_kc, material_kc)) == 0
    fixed_eulerian_pullback = mvec(f_transpose, material_kc)
    assert dot(fixed_eulerian_pullback, fixed_eulerian_pullback) == 3
    assert poly_eval(q, F(3)) != 0


def check_scale_and_b_ledgers() -> None:
    # C171 L2 parent-cross scale is unchanged by a bounded fixed filter.
    raw_q_power = F(-5, 2)
    backward_h_power = F(3, 2)
    assert raw_q_power + backward_h_power == -1

    # At fixed core normalization the paired transform starts at zeta^5
    # instead of zeta^3, hence q^-2.  If one renormalizes the paired profile
    # by q^2 to hold that suppressed Fourier coefficient fixed, the gain is
    # exactly erased.
    single_fourier_power = -3
    paired_fourier_power = -5
    assert paired_fourier_power - single_fourier_power == -2
    assert paired_fourier_power + 2 == single_fourier_power

    # The paired low-frequency moment adds q^-2 on fixed parent-scale A2
    # fibers. With q=n^8, b=n^-2, h=n^12:
    q_n = 8
    b_n = -2
    h_log_coefficient = 12
    wake_n_power = 2 * b_n - 2 * q_n
    active_budget_n_power = 3 * b_n
    ratio_n_power = wake_n_power - active_budget_n_power
    assert wake_n_power == -20
    assert active_budget_n_power == -6
    assert ratio_n_power == -14

    ratios = [
        h_log_coefficient * log(n) / n**14
        for n in (2, 3, 5, 10, 20)
    ]
    assert all(left > right for left, right in zip(ratios, ratios[1:]))
    assert ratios[-1] < F(1, 10**16)


def check_broad_band_degree_gate() -> None:
    # Restricting a multiplier to the C144 axial line gives q distinct
    # points. Triple zeros force degree at least 3q. The existing q=n^8
    # is already much larger than an n^(7/2) collar window and the
    # j^2/log(n) endpoint jet order; the exact comparison below only uses
    # the weaker n^(7/2) benchmark.
    # Exhaust the distinctness assertion at modest exact q; the formula
    # itself is injective for every q because a/q=b/q implies a=b.
    for q in (4, 8, 16, 32):
        distinct_axial_points = {
            F(4) + F(a, q) for a in range(q)
        }
        assert len(distinct_axial_points) == q
        required_degree = 3 * q
        radial_spatial_order = 6 * q
        assert radial_spatial_order > required_degree

    # Construct the sharp univariate Hermite notch at modest q.  Its value
    # and first two derivatives vanish at every axial point, its third
    # derivative does not, and normalization preserves p(0)=1 because no
    # root is zero.
    for q in (3, 5):
        roots = tuple(F(4) + F(a, q) for a in range(q))
        polynomial = [F(1)]
        for root in roots:
            for _ in range(3):
                polynomial = poly_mul(polynomial, [-root, F(1)])
        assert len(polynomial) - 1 == 3 * q
        assert polynomial[0] != 0
        polynomial = [coefficient / polynomial[0] for coefficient in polynomial]
        assert polynomial[0] == 1
        first = poly_diff(polynomial)
        second = poly_diff(first)
        third = poly_diff(second)
        for root in roots:
            assert poly_eval(polynomial, root) == 0
            assert poly_eval(first, root) == 0
            assert poly_eval(second, root) == 0
            assert poly_eval(third, root) != 0

    for n in (2, 3, 5, 10):
        q = n**8
        required_degree = 3 * q
        assert required_degree**2 > n**7
        gevrey_window = n**2 / log(n)
        assert required_degree > gevrey_window


def main() -> None:
    check_a2_core_strain()
    check_pair_weights_and_moment_orders()
    check_triple_notch_and_condition()
    check_scale_and_b_ledgers()
    check_broad_band_degree_gate()
    print("PASS C173: paired radii retain the affine core and cancel moment 3")
    print("PASS C173: fixed A2 notch has an exact two-jet Fourier kernel")
    print("PASS C173: fixed-fiber b-ledger closes; broad-band BAFL remains OPEN")


if __name__ == "__main__":
    main()
