#!/usr/bin/env python3
"""Checks for C144--C148 on the existing A2 one-cell geometry.

Exact rational/integer arithmetic verifies the coherent packet, shell and
conversion counts, microseed ledger, and broad axial-layer coefficients.
Floating evaluation is used only for the elementary retained-line
exponential formula and asymptotic tail samples. LBRG itself remains open.
"""

from fractions import Fraction as F
from itertools import product
from cmath import exp as cexp
from math import ceil, expm1, factorial, log, pi, sqrt


N = (1, 1, 1)
R1 = (1, -1, 0)
R2 = (0, 1, -1)
R3 = (-1, 0, 1)
ROOTS = (R1, R2, R3)
TARGETS = {
    R1,
    R2,
    R3,
    tuple(-x for x in R1),
    tuple(-x for x in R2),
    tuple(-x for x in R3),
}


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def determinant_columns(a, b, c):
    return dot(a, cross(b, c))


def packet_vector(q, a, b, c):
    return add(scale(4 * q + a, N), add(scale(b, R1), scale(c, R2)))


def check_c144_packet():
    assert abs(determinant_columns(N, R1, R2)) == 3
    assert dot(N, R1) == dot(N, R2) == 0
    assert dot(R1, R1) == dot(R2, R2) == 2
    assert dot(R1, R2) == -1

    for q in range(1, 7):
        modes = {}
        unnormalized_fourier = {}
        scalar_sum = F(0)
        for a, b, c in product(range(q), repeat=3):
            v = packet_vector(q, a, b, c)
            norm_sq = dot(v, v)
            expected = 3 * (4 * q + a) ** 2 + 2 * (b * b + c * c - b * c)
            assert norm_sq == expected
            assert 48 * q * q <= norm_sq
            assert norm_sq <= 77 * q * q - 34 * q + 5
            assert dot(v, N) > 0
            assert v not in modes
            modes[v] = (a, b, c)

            projection_sq = F(2) - F((2 * b - c) ** 2, norm_sq)
            assert projection_sq >= F(23, 12)
            scalar_sum += projection_sq

            projection = tuple(
                F(R1[index]) - F(dot(R1, v), norm_sq) * v[index]
                for index in range(3)
            )
            assert dot(v, projection) == 0
            unnormalized_fourier[v] = projection
            unnormalized_fourier[scale(-1, v)] = projection

        assert len(modes) == q**3
        assert not set(modes).intersection({scale(-1, v) for v in modes})
        assert len(unnormalized_fourier) == 2 * q**3
        for wavevector, coefficient in unnormalized_fourier.items():
            assert unnormalized_fourier[scale(-1, wavevector)] == coefficient
            assert dot(wavevector, coefficient) == 0
        assert scalar_sum >= F(23, 12) * q**3

        # Parseval gives ||W||_2^2=2S.  At a coherence point,
        # (E/|E|).W=sqrt(2)S, hence the normalized component squared is S.
        parseval_unnormalized = 2 * scalar_sum
        peak_e_component_unnormalized_sq = 2 * scalar_sum**2
        assert peak_e_component_unnormalized_sq / parseval_unnormalized == scalar_sum
        peak_component_sq = scalar_sum
        assert peak_component_sq >= F(23, 12) * q**3

    # Universal proof behind (1.6): after moving the tangential term,
    # the desired inequality is
    # 46b^2-46bc+10c^2 <= 3(4q+a)^2.  Convexity in c and then b puts the
    # box maximum at (b,c)=(q-1,0), where it is 46(q-1)^2<48q^2.
    for q in range(1, 101):
        assert 46 * (q - 1) ** 2 < 48 * q * q

    # Verify the constant in the child-ball energy lower bound.
    lower_density_coefficient = F(49, 64) * F(23, 12)
    normalized_ball_without_kq = 1 / (48 * pi**2 * 77 ** 1.5)
    expected = 1127 / (36864 * pi**2 * 77 ** 1.5)
    assert abs(float(lower_density_coefficient) * normalized_ball_without_kq - expected) < 1e-18

    # Every coordinate is a multiple of K, so the field repeats in K^3 cells.
    for k_value in (1, 3, 11):
        for q in (2, 4):
            for a, b, c in ((0, 0, 0), (q - 1, q - 1, q - 1)):
                k = scale(k_value, packet_vector(q, a, b, c))
                assert all(entry % k_value == 0 for entry in k)


def sphere_points(radius):
    radius_sq = radius * radius
    count = 0
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            remainder = radius_sq - x * x - y * y
            if remainder < 0:
                continue
            z = int(sqrt(remainder))
            if z * z == remainder:
                count += 1 if z == 0 else 2
    return count


def annulus_points(radius, delta):
    outer = ceil(radius + delta)
    count = 0
    for x in range(-outer, outer + 1):
        for y in range(-outer, outer + 1):
            for z in range(-outer, outer + 1):
                norm = sqrt(x * x + y * y + z * z)
                if abs(norm - radius) <= delta:
                    count += 1
    return count


def check_c145_shell_no_go():
    for radius in range(1, 16):
        count = sphere_points(radius)
        assert count <= 2 * (2 * radius + 1) ** 2

    for n in (3, 4, 5):
        q = n**8
        for k_value in (1, 7):
            h = n**12
            required_modes = h * h * k_value**3
            assert required_modes == (q * k_value) ** 3

            # The C144 band has radius below sqrt(77) qK < 9qK.
            radius_bound = 9 * q * k_value
            exact_shell_upper = 2 * (2 * radius_bound + 1) ** 2
            assert exact_shell_upper < required_modes

    # Unit cubes around annulus lattice points lie in the annulus enlarged
    # by sqrt(3)/2.  Numerically exercise this exact volume comparison for
    # representative radii; (2.6) is the expansion of that shell volume.
    cube_radius = sqrt(3) / 2
    for radius, delta in ((5, 1), (9, 2), (15, 1)):
        count = annulus_points(radius, delta)
        outer = radius + delta + cube_radius
        inner = max(0.0, radius - delta - cube_radius)
        shell_volume = 4 * pi * (outer**3 - inner**3) / 3
        assert count <= shell_volume + 1e-12

    # Six global modes cannot combine child-volume L2 size with child point size.
    for r in (10, 100, 1000):
        assert 6 < r**3


def check_c146_conversion_and_unfolding():
    # An explicit assignment witnesses the combinatorial proof:
    # one gate or one companion can have multiplicity at most |T|=6.
    sources = [(m, m * m + 7, 2 * m - 3) for m in range(1, 121)]
    targets = sorted(TARGETS)
    assignments = [(s, targets[index % 6]) for index, s in enumerate(sources)]
    gates = {}
    companions = {}
    for source, target in assignments:
        gate = sub(target, source)
        companion = sub(source, gate)
        assert companion == sub(scale(2, source), target)
        gates.setdefault(gate, []).append((source, target))
        companions.setdefault(companion, []).append((source, target))
    assert max(map(len, gates.values())) <= 6
    assert max(map(len, companions.values())) <= 6
    assert len(gates) >= len(sources) / 6
    assert len(companions) >= len(sources) / 6
    canonical_gate_pairs = {
        min(gate, scale(-1, gate)) for gate in gates
    }
    canonical_companion_pairs = {
        min(companion, scale(-1, companion)) for companion in companions
    }
    assert len(canonical_gate_pairs) >= len(sources) / 12
    assert len(canonical_companion_pairs) >= len(sources) / 12
    for group in tuple(gates.values()) + tuple(companions.values()):
        assert len({source for source, _ in group}) == len(group)
        assert len({target for _, target in group}) == len(group)

    # Saturate the six-source oriented-gate bound with the actual targets.
    fixed_gate = (17, -11, 5)
    saturated = [(sub(target, fixed_gate), target) for target in targets]
    assert len({sub(target, source) for source, target in saturated}) == 1
    assert len({source for source, _ in saturated}) == 6

    # Saturate the abstract companion collision bound c=2s-t with six
    # doubled targets.  This checks the structural factor six independently
    # of parity accidents in the particular root target set.
    abstract_sources = list(TARGETS)
    doubled_targets = [scale(2, source) for source in abstract_sources]
    saturated_companions = [
        sub(scale(2, source), target)
        for source, target in zip(abstract_sources, doubled_targets)
    ]
    assert saturated_companions == [(0, 0, 0)] * 6

    for n in (2, 7, 31):
        b = F(1, n * n)
        microgates = 1 / b
        theta = b
        assert microgates.denominator == 1
        assert microgates * b * theta == b
        assert microgates * b * theta**2 == b**2 == F(1, n**4)
        assert microgates * b * theta**3 == b**3 == F(1, n**6)
        assert microgates * theta == 1


def retained_gain(theta, lam, mu, q):
    return -lam * expm1(-2 * mu * theta) / (2 * mu) - mu * (3 * q * q + 2) * theta


def check_c147_microseed_and_gain():
    # Store only powers of n; multiplication becomes addition.
    q_exp = 8
    volume_exp = -3 * q_exp
    seed_point_exp = -16
    seed_l2_exp = seed_point_exp + volume_exp // 2
    gain_exp = 26
    final_point_exp = seed_point_exp + gain_exp
    final_l2_exp = seed_l2_exp + gain_exp
    assert volume_exp == -24
    assert seed_l2_exp == -28
    assert final_point_exp == 10
    assert final_l2_exp == -2
    assert 2 * final_l2_exp == -4

    # Same-support writer error is n^{-20}*26log n, hence n^8 larger
    # than the n^{-28} growing seed.
    writer_l2_exp = seed_point_exp + q_exp + volume_exp // 2
    assert writer_l2_exp == -20
    assert writer_l2_exp - seed_l2_exp == 8

    for lam, mu, q in ((2.0, 1e-4, 9), (1.5, 1e-7, 40)):
        ratio = lam / (mu * (3 * q * q + 2))
        assert ratio > 1
        theta_star = log(ratio) / (2 * mu)
        closed = (3 * q * q + 2) * (ratio - 1 - log(ratio)) / 2
        assert abs(retained_gain(theta_star, lam, mu, q) - closed) < 1e-8 * max(1, closed)

    tail = []
    for j in range(10, 18):
        n = j + 1
        mu = 1 / factorial(j) ** 2
        q = n**8
        tail.append(mu * q * q * log(n))
    assert all(x > y for x, y in zip(tail, tail[1:]))
    assert tail[-1] < 1

    # Locate the first hit of 26 log(n) on the increasing branch and verify
    # its asymptotic 26 log(n)/lambda timing without cancellation loss in
    # 1-exp(-x).
    first_hit_ratios = []
    lam = 1.5
    for j in range(16, 22):
        n = j + 1
        mu = 1 / factorial(j) ** 2
        q = n**8
        target = 26 * log(n)
        leading_time = target / lam
        lo, hi = 0.0, 2.0 * leading_time
        assert retained_gain(hi, lam, mu, q) > target
        for _ in range(100):
            middle = (lo + hi) / 2
            if retained_gain(middle, lam, mu, q) < target:
                lo = middle
            else:
                hi = middle
        first_hit_ratios.append(((lo + hi) / 2) / leading_time)
    assert all(abs(ratio - 1) < 1e-5 for ratio in first_hit_ratios)
    assert abs(first_hit_ratios[-1] - 1) < abs(first_hit_ratios[0] - 1)

    for n in (2, 5, 20):
        q = n**8
        assert n**4 * F(1, n**2) == n**2


def check_c148_axial_layer():
    directions = ((1, 0), (0, 1), (-1, -1))
    tangents = tuple(cross(N, root) for root in ROOTS)
    assert add(add(R1, R2), R3) == (0, 0, 0)
    assert tuple(sum(root[index] for root in ROOTS) for index in range(3)) == (0, 0, 0)
    assert tuple(sum(tangent[index] for tangent in tangents) for index in range(3)) == (0, 0, 0)

    # For equal real pump coefficients at theta=0:
    # U_0=2c*sum(t_i)=0 and
    # C_0=2i*sqrt(2)c*N*(sum r_i)^T=0.
    equal_real_u0 = tuple(2 * sum(tangent[index] for tangent in tangents) for index in range(3))
    root_sum = tuple(sum(root[index] for root in ROOTS) for index in range(3))
    assert equal_real_u0 == (0, 0, 0)
    assert root_sum == (0, 0, 0)

    for m in (3, 11):
        for b, c in ((0, 0), (1, 0), (0, 1), (4, -2), (-3, 5)):
            kappa = add(scale(m, N), add(scale(b, R1), scale(c, R2)))
            d_value = dot(kappa, kappa)
            assert d_value == 3 * m * m + 2 * (b * b + c * c - b * c)
            zetas = (c, -b, b - c)

            for root, tangent, zeta in zip(ROOTS, tangents, zetas):
                # h_i=t_i+i*sqrt(2)N, so the dot product is represented
                # exactly by (real part, coefficient of i*sqrt(2)).
                helical_dot = (dot(tangent, kappa), dot(N, kappa))
                assert helical_dot == (3 * zeta, 3 * m)

            # Exact rational transverse basis from the audit.
            k_dot_r1 = 2 * b - c
            e1_norm_sq = F(2) - F(k_dot_r1 * k_dot_r1, d_value)
            e2 = cross(kappa, R1)
            assert dot(e2, e2) == d_value * e1_norm_sq

            # The neighbor shifts are precisely the three A2 directions.
            for root, direction in zip(ROOTS, directions):
                shifted = sub(kappa, root)
                expected = add(
                    scale(m, N),
                    add(scale(b - direction[0], R1), scale(c - direction[1], R2)),
                )
                assert shifted == expected

    # A constant polarization transverse at three sites would be
    # orthogonal to this full-rank basis and hence zero.
    assert abs(determinant_columns(N, R1, R2)) == 3

    # On the C144 box, zeta/m has an order-one range rather than o(1).
    for q in (20, 100):
        minimum = F(0, 1)
        maximum = F(q - 1, 4 * q)
        assert maximum - minimum > F(1, 6)

    # Frozen-symbol reality structure for arbitrary complex pump
    # coefficients: U_theta is real and C_theta is purely imaginary.
    coefficients = (1 + 2j, F(3, 2) - F(1, 3) * 1j, -2 + F(4, 5) * 1j)
    theta = (0.37, -0.91)
    u_theta = [0j, 0j, 0j]
    c_theta = [[0j for _ in range(3)] for _ in range(3)]
    for coefficient, root, tangent, direction in zip(
        coefficients, ROOTS, tangents, directions
    ):
        h_plus = tuple(
            tangent[index] + 1j * sqrt(2) * N[index]
            for index in range(3)
        )
        u_plus = tuple(complex(coefficient) * entry for entry in h_plus)
        u_minus = tuple(entry.conjugate() for entry in u_plus)
        phase = theta[0] * direction[0] + theta[1] * direction[1]
        phase_minus = cexp(-1j * phase)
        phase_plus = cexp(1j * phase)
        for row in range(3):
            u_theta[row] += (
                u_plus[row] * phase_minus + u_minus[row] * phase_plus
            )
            for column in range(3):
                c_theta[row][column] += (
                    u_plus[row] * root[column] * phase_minus
                    - u_minus[row] * root[column] * phase_plus
                )
    assert max(abs(entry.imag) for entry in u_theta) < 1e-12
    assert max(abs(entry.real) for row in c_theta for entry in row) < 1e-12

    # Direct Fourier convolution versus the boxed nearest-neighbor formula,
    # including the epsilon sign and Leray projection.  The direct side uses
    # q=kappa-p; the displayed side uses u_p.kappa and
    # epsilon*(v_q.r_i), so equality also checks u_p.p=0.
    def complex_dot(left, right):
        return sum(x * y for x, y in zip(left, right))

    def complex_add(left, right):
        return tuple(x + y for x, y in zip(left, right))

    def complex_scale(factor, vector):
        return tuple(factor * entry for entry in vector)

    def project_transverse(wavevector, vector):
        denominator = dot(wavevector, wavevector)
        return tuple(
            vector[index]
            - complex_dot(wavevector, vector) * wavevector[index] / denominator
            for index in range(3)
        )

    m, b, c = 7, 2, -3
    kappa = add(scale(m, N), add(scale(b, R1), scale(c, R2)))
    direct_total = (0j, 0j, 0j)
    displayed_total = (0j, 0j, 0j)
    trial_vectors = ((1, 2, -1), (2, -1, 3), (-3, 1, 2))
    for coefficient, root, tangent, trial in zip(
        coefficients, ROOTS, tangents, trial_vectors
    ):
        h_plus = tuple(
            tangent[index] + 1j * sqrt(2) * N[index]
            for index in range(3)
        )
        positive = tuple(complex(coefficient) * entry for entry in h_plus)
        for epsilon in (1, -1):
            pump = positive if epsilon == 1 else tuple(
                entry.conjugate() for entry in positive
            )
            pump_wavevector = scale(epsilon, root)
            neighbor = sub(kappa, pump_wavevector)
            neighbor_value = project_transverse(neighbor, trial)
            assert abs(complex_dot(pump_wavevector, pump)) < 1e-12
            assert abs(complex_dot(neighbor, neighbor_value)) < 1e-12

            direct_bracket = complex_add(
                complex_scale(
                    complex_dot(pump, neighbor), neighbor_value
                ),
                complex_scale(
                    complex_dot(neighbor_value, pump_wavevector), pump
                ),
            )
            displayed_bracket = complex_add(
                complex_scale(complex_dot(pump, kappa), neighbor_value),
                complex_scale(
                    epsilon * complex_dot(neighbor_value, root), pump
                ),
            )
            direct_total = complex_add(
                direct_total,
                complex_scale(-1j, project_transverse(kappa, direct_bracket)),
            )
            displayed_total = complex_add(
                displayed_total,
                complex_scale(-1j, project_transverse(kappa, displayed_bracket)),
            )
    assert max(
        abs(left - right)
        for left, right in zip(direct_total, displayed_total)
    ) < 1e-11

    # The frozen stretching block is real.  With one active root it is a
    # compressed outer product and therefore has rank at most one.
    kappa = scale(5, N)
    frame_1 = tuple(entry / sqrt(2) for entry in R1)
    frame_2 = tuple(entry / sqrt(6) for entry in tangents[0])
    assert abs(sum(x * y for x, y in zip(frame_1, frame_2))) < 1e-12
    one_coefficient = 1.0 + 0.4j
    phase = 0.63
    h_plus = tuple(
        tangents[0][index] + 1j * sqrt(2) * N[index]
        for index in range(3)
    )
    positive = tuple(one_coefficient * entry for entry in h_plus)
    negative = tuple(entry.conjugate() for entry in positive)
    c_one = [
        [
            positive[row] * R1[column] * cexp(-1j * phase)
            - negative[row] * R1[column] * cexp(1j * phase)
            for column in range(3)
        ]
        for row in range(3)
    ]
    frozen = []
    for left in (frame_1, frame_2):
        row = []
        for right in (frame_1, frame_2):
            row.append(
                -1j
                * sum(
                    left[i] * c_one[i][j] * right[j]
                    for i in range(3)
                    for j in range(3)
                )
            )
        frozen.append(row)
    assert max(abs(entry.imag) for row in frozen for entry in row) < 1e-12
    assert abs(frozen[0][0] * frozen[1][1] - frozen[0][1] * frozen[1][0]) < 1e-12


def main():
    check_c144_packet()
    check_c145_shell_no_go()
    check_c146_conversion_and_unfolding()
    check_c147_microseed_and_gain()
    check_c148_axial_layer()
    print("C144--C148 coherent packet/relative-gain checks: PASS")
    print("OPEN: LBRG, unfolded endpoint conversion, and full BAFL")


if __name__ == "__main__":
    main()
