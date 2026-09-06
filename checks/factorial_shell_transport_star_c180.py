#!/usr/bin/env python3
"""Dependency-free checks for C180's factorial-shell transport star.

The script checks the finite arithmetic, tight-frame, Fourier-symbol,
collision, shrinking-slab, and principal-composition ledgers.  It does not
prove FFCC, endpoint coherence, C125/RIGM/BAFL, or a one-cell stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import permutations, product
from math import ceil, cos, factorial, log, pi, sin, sqrt


IV = tuple[int, int, int]
RV = tuple[float, float, float]
CV = tuple[complex, complex, complex]


def iadd(a: IV, b: IV) -> IV:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def isub(a: IV, b: IV) -> IV:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def ineg(a: IV) -> IV:
    return tuple(-x for x in a)  # type: ignore[return-value]


def idot(a: IV, b: IV) -> int:
    return sum(x * y for x, y in zip(a, b))


def rdot(a: RV, b: RV) -> float:
    return sum(x * y for x, y in zip(a, b))


def rnorm(a: RV) -> float:
    return sqrt(rdot(a, a))


def rscale(c: float, a: RV) -> RV:
    return tuple(c * x for x in a)  # type: ignore[return-value]


def rsub(a: RV, b: RV) -> RV:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def rcross(a: RV, b: RV) -> RV:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def cdot_bilinear(a: CV, b: RV) -> complex:
    return sum((x * y for x, y in zip(a, b)), 0j)


def hdot(a: CV, b: CV) -> complex:
    return sum((x.conjugate() * y for x, y in zip(a, b)), 0j)


def cnorm(a: CV) -> float:
    return sqrt(hdot(a, a).real)


def cadd(a: CV, b: CV) -> CV:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def csub(a: CV, b: CV) -> CV:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def cscale(c: complex, a: CV) -> CV:
    return tuple(c * x for x in a)  # type: ignore[return-value]


def real_as_complex(a: RV) -> CV:
    return tuple(complex(x) for x in a)  # type: ignore[return-value]


def cross_real_complex(a: RV, b: CV) -> CV:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def leray(k: RV, a: CV) -> CV:
    kk = rdot(k, k)
    assert kk > 0
    alpha = cdot_bilinear(a, k) / kk
    return csub(a, tuple(alpha * x for x in k))  # type: ignore[arg-type]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def factor_integer(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def factorial_valuation(j: int, p: int) -> int:
    total = 0
    power = p
    while power <= j:
        total += j // power
        power *= p
    return total


def r2_of_square_formula(m: int) -> int:
    """r_2(m^2)=4 product_{p=1 mod 4}(2 v_p(m)+1)."""
    answer = 4
    for p, exponent in factor_integer(m).items():
        if p % 4 == 1:
            answer *= 2 * exponent + 1
    return answer


def r2_of_square_bruteforce(m: int) -> int:
    target = m * m
    return sum(
        1
        for a in range(-m, m + 1)
        for b in range(-m, m + 1)
        if a * a + b * b == target
    )


def check_factorial_shell_count() -> None:
    # Independently compare the divisor formula with literal ordered signed
    # representation counts on several composite radii.
    for m in (1, 2, 5, 10, 13, 25, 65, 85):
        assert r2_of_square_formula(m) == r2_of_square_bruteforce(m)

    fixed_primes = (5, 13, 17, 29, 37, 41, 53, 61, 73)
    assert len(fixed_primes) == 9
    assert all(is_prime(p) and p % 4 == 1 for p in fixed_primes)

    # Formula (1.2), restricted to nine factors, already beats q at a
    # representative large stage.  No enumeration of the enormous shell is
    # performed and no prime-number theorem is used.
    j = 50_000
    restricted_r2 = 4
    elementary_lower = 4
    for p in fixed_primes:
        valuation = factorial_valuation(j, p)
        assert valuation >= F(j, 2 * p)
        restricted_r2 *= 16 * valuation + 1
        elementary_lower *= F(8 * j, p)
    q = (j + 1) ** 8
    target = q
    assert restricted_r2 >= elementary_lower > target

    # Check the exact exponent conversion for small factorials without ever
    # materializing K^2 in a shell loop.
    for small_j in (2, 3, 4, 5, 7, 10):
        k = factorial(small_j) ** 8
        from_factorization = r2_of_square_formula(k)
        from_valuations = 4
        for p in range(2, small_j + 1):
            if is_prime(p) and p % 4 == 1:
                from_valuations *= 16 * factorial_valuation(small_j, p) + 1
        assert from_factorization == from_valuations


def octahedral_orbit(seed: IV) -> list[IV]:
    return sorted({
        tuple(signs[i] * perm[i] for i in range(3))
        for perm in set(permutations(seed))
        for signs in product((-1, 1), repeat=3)
    })  # type: ignore[return-value]


def positive_lex(k: IV) -> bool:
    return k > ineg(k)


def helical_palette(gates: list[IV]) -> dict[IV, CV]:
    palette: dict[IV, CV] = {}
    for k in gates:
        if not positive_lex(k):
            continue
        kr: RV = tuple(float(x) for x in k)  # type: ignore[assignment]
        radius = rnorm(kr)
        khat = rscale(1.0 / radius, kr)
        refs: tuple[RV, ...] = (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
        ref = min(refs, key=lambda e: abs(rdot(khat, e)))
        tangent = rcross(khat, ref)
        tangent = rscale(1.0 / rnorm(tangent), tangent)
        second = rcross(khat, tangent)
        h: CV = tuple(
            (tangent[i] + 1j * second[i]) / sqrt(2.0)
            for i in range(3)
        )  # type: ignore[assignment]
        palette[k] = h
        palette[ineg(k)] = tuple(z.conjugate() for z in h)  # type: ignore[assignment]
    assert set(palette) == set(gates)
    return palette


def check_octahedral_tight_frame() -> tuple[list[IV], dict[IV, CV]]:
    gates = octahedral_orbit((1, 2, 2))
    qg = len(gates)
    radius_sq = idot(gates[0], gates[0])
    assert qg > 0 and radius_sq == 9 and all(idot(g, g) == 9 for g in gates)
    assert set(gates) == {ineg(g) for g in gates}

    second_moment = [[0 for _ in range(3)] for _ in range(3)]
    for g in gates:
        for i in range(3):
            for j in range(3):
                second_moment[i][j] += g[i] * g[j]
    for i in range(3):
        for j in range(3):
            expected = qg * radius_sq // 3 if i == j else 0
            assert second_moment[i][j] == expected

    palette = helical_palette(gates)
    radius = 3.0
    for g, h in palette.items():
        gr: RV = tuple(float(x) for x in g)  # type: ignore[assignment]
        assert abs(cdot_bilinear(h, gr)) < 2e-14
        assert abs(cnorm(h) - 1.0) < 2e-14
        curl_h = cscale(1j, cross_real_complex(gr, h))
        assert cnorm(csub(curl_h, cscale(radius, h))) < 5e-14
        assert all(
            abs(x - y.conjugate()) < 2e-14
            for x, y in zip(h, palette[ineg(g)])
        )

    # This probes the helical modulus identity itself, not merely the real
    # second-moment matrix above.
    for p in ((2.0, -3.0, 5.0), (7.0, 1.0, -4.0), (1.0, 1.0, 1.0)):
        actual = sum(abs(cdot_bilinear(h, p)) ** 2 for h in palette.values())
        expected = qg * rdot(p, p) / 3.0
        assert abs(actual - expected) < 2e-12 * expected
    return gates, palette


def tangent_edge(p: RV, g: IV, h: CV, a: RV) -> CV:
    gr: RV = tuple(float(x) for x in g)  # type: ignore[assignment]
    k: RV = tuple(x + y for x, y in zip(p, gr))  # type: ignore[assignment]
    ac = real_as_complex(a)
    raw = cadd(
        cscale(cdot_bilinear(h, p), ac),
        cscale(rdot(a, gr), h),
    )
    return leray(k, raw)


def source_basis(p: RV) -> tuple[RV, RV]:
    phat = rscale(1.0 / rnorm(p), p)
    refs: tuple[RV, ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    ref = min(refs, key=lambda e: abs(rdot(phat, e)))
    e1 = rcross(phat, ref)
    e1 = rscale(1.0 / rnorm(e1), e1)
    e2 = rcross(phat, e1)
    return e1, e2


def hermitian_two_by_two_error(gram: list[list[complex]]) -> float:
    a = gram[0][0].real - 1.0
    d = gram[1][1].real - 1.0
    b = gram[0][1]
    trace_half = (a + d) / 2.0
    radius = sqrt(((a - d) / 2.0) ** 2 + abs(b) ** 2)
    return max(abs(trace_half + radius), abs(trace_half - radius))


def check_dominant_symbol_and_star(gates: list[IV], palette: dict[IV, CV]) -> None:
    p0: RV = (4.0, 5.0, 6.0)
    previous = 1.0
    for q in (20, 50, 100, 250, 500):
        p = rscale(float(q), p0)
        e1, e2 = source_basis(p)
        basis = (e1, e2)
        common = sqrt(3.0) / (sqrt(len(gates)) * rnorm(p))
        images: list[list[CV]] = [[], []]
        max_edge_remainder = 0.0
        for g in gates:
            h = palette[g]
            for index, a in enumerate(basis):
                exact = tangent_edge(p, g, h, a)
                leading = cscale(cdot_bilinear(h, p), real_as_complex(a))
                remainder = cnorm(csub(exact, leading))
                max_edge_remainder = max(max_edge_remainder, remainder)
                images[index].append(cscale(common, exact))
        # The absolute subprincipal remainder is uniformly shell-sized.
        assert max_edge_remainder < 9.0

        gram = [[0j, 0j], [0j, 0j]]
        for i in range(2):
            for j in range(2):
                gram[i][j] = sum(
                    hdot(images[i][r], images[j][r])
                    for r in range(len(gates))
                )
        error = hermitian_two_by_two_error(gram)
        assert error < previous
        assert q * error < 6.0
        previous = error
    assert previous < 0.01


def sphere_shell(kscale: int) -> list[IV]:
    target = 2 * kscale * kscale
    bound = ceil(sqrt(target))
    return [
        (x, y, z)
        for x in range(-bound, bound + 1)
        for y in range(-bound, bound + 1)
        for z in range(-bound, bound + 1)
        if x * x + y * y + z * z == target
    ]


def check_collision_pruning() -> None:
    kscale = 5
    gates = sphere_shell(kscale)
    assert gates and set(gates) == {ineg(g) for g in gates}

    differences: set[IV] = set()
    for g in gates:
        for h in gates:
            delta = isub(h, g)
            if all(x % kscale == 0 for x in delta):
                d = tuple(x // kscale for x in delta)  # type: ignore[assignment]
                differences.add(d)
                assert all(-2 <= x <= 2 for x in d)
    nonzero = differences - {(0, 0, 0)}
    assert len(nonzero) <= 124

    candidates = [
        (100 + x, 200 + y, 300 + z)
        for x in range(9) for y in range(9) for z in range(9)
    ]
    chosen: list[IV] = []
    chosen_set: set[IV] = set()
    for k in candidates:
        if all(iadd(k, d) not in chosen_set for d in nonzero):
            chosen.append(k)
            chosen_set.add(k)
    assert len(chosen) * 125 >= len(candidates)

    positive = chosen[:24]
    sources = positive + [ineg(k) for k in positive]
    source_set = set(sources)
    assert len(source_set) == 2 * len(positive)
    for k in sources:
        for d in nonzero:
            assert iadd(k, d) not in source_set

    outputs = {
        tuple(kscale * k[i] + g[i] for i in range(3))
        for k in sources for g in gates
    }
    assert len(outputs) == len(sources) * len(gates)


def check_shrinking_slab_ledger() -> None:
    # The nearest-cube proof may erode every half-width by 1 because
    # sqrt(3)/2<1.  At a representative stage the resulting rigorous volume
    # lower bound already supplies 125 q^2 candidates.
    n = 100
    q = n**8
    window = ceil(12 * log(n))
    delta = F(1, 100)
    w1 = delta * q / window
    w2 = delta * q / window
    w3 = delta * q / (window * window)
    assert min(w1, w2, w3) > 1
    eroded_volume = 8 * (w1 - 1) * (w2 - 1) * (w3 - 1)
    assert eroded_volume > 125 * q * q

    # In normalized C176 frame coordinates, an initial normal deviation
    # q/J^2 accumulates only q/J over J nilpotent returns.  This is the
    # exponent arithmetic behind (3.6); geometry constants are external.
    xi1 = delta * q / window
    xi2 = delta * q / window
    xi3 = delta * q / (window * window)
    shear_shift = window * xi3
    assert shear_shift == delta * q / window
    assert xi1 + shear_shift == 2 * delta * q / window
    assert xi2 == delta * q / window

    # Reciprocal widths have relative volume J^4 and filled-envelope
    # point-amplitude tax J^2.
    reciprocal_volume_ratio = window**4
    point_tax = window**2
    assert point_tax * point_tax == reciprocal_volume_ratio


def check_transport_composition_boundary() -> None:
    # W=(cos y,0) has the exact area-preserving flow
    # X_t(a,b)=(a+t cos b,b).  Evaluate the transported scalar at X_t and
    # verify that all point values, not only an energy average, are retained.
    def datum(x: float, y: float) -> float:
        return sin(2.0 * x - y) + 0.37 * cos(3.0 * y)

    def transported(t: float, x: float, y: float) -> float:
        return datum(x - t * cos(y), y)

    for t in (-2.3, -0.1, 0.0, 0.7, 4.2):
        initial_values = []
        transported_values = []
        for ix in range(37):
            for iy in range(31):
                a = 2.0 * pi * ix / 37
                b = 2.0 * pi * iy / 31
                x = a + t * cos(b)
                y = b
                initial_values.append(datum(a, b))
                transported_values.append(transported(t, x, y))
        assert max(
            abs(x - y) for x, y in zip(initial_values, transported_values)
        ) < 3e-15
        assert abs(
            max(map(abs, initial_values))
            - max(map(abs, transported_values))
        ) < 3e-15

    # The first Fourier sideband matrix of the same scalar transport is
    # skew-Hermitian on every symmetric finite block.  It has bright
    # off-diagonal edges but zero instantaneous L2-energy derivative.
    m = 7
    modes = list(range(-8, 9))
    size = len(modes)
    matrix = [[0j for _ in modes] for _ in modes]
    for col, source_n in enumerate(modes):
        for shift in (-1, 1):
            output_n = source_n + shift
            if output_n in modes:
                row = modes.index(output_n)
                matrix[row][col] += -1j * m / 2.0
    for i in range(size):
        for j in range(size):
            assert abs(matrix[i][j] + matrix[j][i].conjugate()) < 1e-15
    vector = [complex((3 * i + 2) % 11 - 5, (5 * i + 1) % 13 - 6)
              for i in range(size)]
    av = [sum(matrix[i][j] * vector[j] for j in range(size))
          for i in range(size)]
    energy_derivative = 2.0 * sum(
        vector[i].conjugate() * av[i] for i in range(size)
    ).real
    assert abs(energy_derivative) < 1e-12

    # Common tight-star normalization has vanishing strain even though the
    # high source sees order-one collective Fourier action.
    previous = 1.0
    for n in (10, 30, 100, 300):
        q = n**8
        qg = q + 47
        rho = sqrt(3.0) / (q * sqrt(qg))
        strain = sqrt(2.0) * rho * qg
        assert strain <= 4.0 / sqrt(q)
        assert strain < previous
        previous = strain


def main() -> None:
    check_factorial_shell_count()
    gates, palette = check_octahedral_tight_frame()
    check_dominant_symbol_and_star(gates, palette)
    check_collision_pruning()
    check_shrinking_slab_ledger()
    check_transport_composition_boundary()
    print("C180 factorial-shell transport-star checks passed")
    print("BOUNDARY: FFCC and physical endpoint focus are not proved")


if __name__ == "__main__":
    main()
