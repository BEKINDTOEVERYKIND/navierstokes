#!/usr/bin/env python3
"""Dependency-free exact checks for C177's same-curl reservoir.

The checker uses Gaussian rationals for exact Fourier/Beltrami algebra and
plain integer/Fraction arithmetic for the stage ledgers.  It does not prove
the full-polarization symbol, C125, RIGM, BAFL, or a one-cell stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import exp, factorial, log


C = tuple[F, F]
CV = tuple[C, C, C]
K = tuple[int, int, int]
Field = dict[K, CV]

CZ: C = (F(0), F(0))
CO: C = (F(1), F(0))
CI: C = (F(0), F(1))
CVZ: CV = (CZ, CZ, CZ)


def ca(a: C, b: C) -> C:
    return (a[0] + b[0], a[1] + b[1])


def cn(a: C) -> C:
    return (-a[0], -a[1])


def cs(x: F, a: C) -> C:
    return (x * a[0], x * a[1])


def cm(a: C, b: C) -> C:
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def cc(a: C) -> C:
    return (a[0], -a[1])


def cabs2(a: C) -> F:
    return a[0] * a[0] + a[1] * a[1]


def va(a: CV, b: CV) -> CV:
    return tuple(ca(x, y) for x, y in zip(a, b))  # type: ignore[return-value]


def vn(a: CV) -> CV:
    return tuple(cn(x) for x in a)  # type: ignore[return-value]


def vs(x: F, a: CV) -> CV:
    return tuple(cs(x, y) for y in a)  # type: ignore[return-value]


def vc(a: CV) -> CV:
    return tuple(cc(x) for x in a)  # type: ignore[return-value]


def dot_real(a: CV, k: K) -> C:
    out = CZ
    for x, y in zip(a, k):
        out = ca(out, cs(F(y), x))
    return out


def cross(a: CV, b: CV) -> CV:
    return (
        ca(cm(a[1], b[2]), cn(cm(a[2], b[1]))),
        ca(cm(a[2], b[0]), cn(cm(a[0], b[2]))),
        ca(cm(a[0], b[1]), cn(cm(a[1], b[0]))),
    )


def cross_real(k: K, a: CV) -> CV:
    kr: CV = tuple((F(x), F(0)) for x in k)  # type: ignore[assignment]
    return cross(kr, a)


def times_i(a: CV) -> CV:
    return tuple(cm(CI, x) for x in a)  # type: ignore[return-value]


def leray(k: K, a: CV) -> CV:
    norm = sum(x * x for x in k)
    assert norm
    kdota = dot_real(a, k)
    correction: CV = tuple(cs(F(x, norm), kdota) for x in k)  # type: ignore[assignment]
    return va(a, vn(correction))


def negk(k: K) -> K:
    return tuple(-x for x in k)  # type: ignore[return-value]


def kadd(a: K, b: K) -> K:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def reality_mode(k: K, h: CV) -> Field:
    return {k: h, negk(k): vc(h)}


def field_add(*fields: Field) -> Field:
    out: Field = {}
    for field in fields:
        for k, value in field.items():
            out[k] = va(out.get(k, CVZ), value)
    return {k: v for k, v in out.items() if v != CVZ}


def adv_cross(left: Field, right: Field) -> Field:
    """(left.grad)right + (right.grad)left in Fourier variables."""
    out: Field = {}
    for p, u in left.items():
        for q, r in right.items():
            first = vs(F(1), tuple(cm(dot_real(u, q), x) for x in r))
            second = tuple(cm(dot_real(r, p), x) for x in u)
            value = times_i(va(first, second))
            k = kadd(p, q)
            out[k] = va(out.get(k, CVZ), value)
    return {k: v for k, v in out.items() if v != CVZ}


def cross_convolution(left: Field, right: Field) -> Field:
    out: Field = {}
    for p, u in left.items():
        for q, r in right.items():
            k = kadd(p, q)
            out[k] = va(out.get(k, CVZ), cross(u, r))
    return {k: v for k, v in out.items() if v != CVZ}


def projected(field: Field) -> Field:
    out: Field = {}
    for k, value in field.items():
        if k == (0, 0, 0):
            # The advective zero mode is zero in the examples below.  A
            # constant cross product is killed after multiplication by the
            # curl-eigenvalue difference only when the identity requires it.
            assert value == CVZ
            continue
        pv = leray(k, value)
        if pv != CVZ:
            out[k] = pv
    return out


def field_scale(x: F, field: Field) -> Field:
    return {k: vs(x, v) for k, v in field.items() if vs(x, v) != CVZ}


def field_norm_sq(field: Field) -> F:
    return sum((sum(cabs2(x) for x in v) for v in field.values()), F(0))


def check_beltrami_and_cross_curl() -> None:
    # Unit shell, positive helicity: h=t+i k cross t.
    k1 = (1, 0, 0)
    hp1: CV = (CZ, CO, CI)
    k2 = (0, 1, 0)
    hp2: CV = (CI, CZ, CO)
    hm2: CV = (cn(CI), CZ, CO)

    p = reality_mode(k1, hp1)
    g = reality_mode(k2, hp2)
    r = reality_mode(k2, hm2)

    for field, kappa in ((p, 1), (g, 1), (r, -1)):
        for k, value in field.items():
            assert dot_real(value, k) == CZ
            assert times_i(cross_real(k, value)) == vs(F(kappa), value)

    # Each coefficient in the epsilon-polynomial for P+epsilon*G is dark:
    # pump self, cross tangent, and reservoir self.  Hence no cancellation
    # between different epsilon orders is being used.
    assert projected(adv_cross(p, p)) == {}
    assert projected(adv_cross(p, g)) == {}
    assert projected(adv_cross(g, g)) == {}

    # Any same-curl sum is projected-Euler dark, including every cross path.
    same = field_add(p, g)
    assert projected(adv_cross(same, same)) == {}
    for eps_test in (F(-3, 5), F(0), F(2, 7), F(11, 4)):
        background = field_add(p, field_scale(eps_test, g))
        assert projected(adv_cross(background, background)) == {}

    # Cross-curl identity for the opposite-helicity field.
    lhs = projected(adv_cross(p, r))
    rhs_raw = cross_convolution(p, r)
    # Constant Fourier output of U x R has zero coefficient here.
    rhs = projected(field_scale(F(2), rhs_raw))
    assert lhs == rhs
    assert lhs  # The generic boundary is genuinely bright in this sample.
    # Swapping fields reverses both the cross product and the curl gap, so
    # the symmetric advective cross term is unchanged.
    lhs_swapped = projected(adv_cross(r, p))
    rhs_swapped = projected(
        field_scale(F(-2), cross_convolution(r, p))
    )
    assert lhs_swapped == lhs == rhs_swapped

    # Orthogonal supports give the exact Pythagorean reservoir energy.
    assert set(p).isdisjoint(g)
    eps = F(1, 7)
    combined = field_add(p, field_scale(eps, g))
    assert field_norm_sq(combined) == field_norm_sq(p) + eps * eps * field_norm_sq(g)

    # A common positive heat scalar preserves all coefficient ratios,
    # phases, and polarizations on the same shell.  Both helicities at one
    # wavevector receive the same multiplier because the rate uses |k|^2.
    heat_scalar = F(5, 13)
    assert vs(heat_scalar, hp2)[0] == cm((heat_scalar, F(0)), hp2[0])
    assert vs(heat_scalar, hp2)[2] == cm((heat_scalar, F(0)), hp2[2])
    assert field_norm_sq(field_scale(heat_scalar, g)) == heat_scalar**2 * field_norm_sq(g)
    assert field_norm_sq(field_scale(heat_scalar, r)) == heat_scalar**2 * field_norm_sq(r)


def check_action_and_stage_ledgers() -> None:
    # Heat action bounds for x in [0,1].
    eps = 0.01
    duration = 2.0
    for x in (0.0, 0.01, 0.2, 1.0):
        ratio = 1.0 if x == 0.0 else -__import__("math").expm1(-x) / x
        action = eps * duration * ratio
        assert eps * duration / 2.0 <= action <= eps * duration * (1.0 + 1e-14)

    # Exact exponent ledger: b gate, b source norm, J=b^-1 windows.
    for n in (4, 7, 12):
        b = F(1, n * n)
        j_windows = n * n
        per_window_work = b ** 3
        assert j_windows * per_window_work == b ** 2
        # Energy of J mutually orthogonal, individually normalized
        # reservoirs.  Without orthogonality only the norm of their sum is
        # intrinsic, as stated in the note.
        assert j_windows * b ** 2 == b
        assert b * j_windows == 1       # One static gate has order-one action.

    # Conditional C125 relative charge b log n tends to zero and never
    # contains the n^28 inverse-seed factor.
    previous = 1.0
    for n in (10, 100, 1000, 10000):
        charge = log(n) / (n * n)
        assert charge < previous
        previous = charge
    assert previous < 1e-6


def check_factorial_heat_flatness() -> None:
    # Exact C127 mu=nu/((n-1)!)^2.  For representative fixed polynomial
    # bands/windows and every tested M, n^M x_n rapidly decreases to zero.
    c_freq = 4
    d_time = 3
    for m_power in (0, 2, 8, 16):
        values = []
        for n in range(20, 31):
            numerator = n ** (2 * c_freq + d_time + m_power)
            denominator = factorial(n - 1) ** 2
            values.append(F(numerator, denominator))
        assert all(values[i + 1] < values[i] for i in range(len(values) - 1))
        assert values[-1] < F(1, 10**12)

    # Pure-normal vectors a*N have a fixed nonzero radius only for +/-a.
    nvec_norm_sq = 3
    radius_sq = 12
    allowed = [a for a in range(-20, 21) if a and a * a * nvec_norm_sq == radius_sq]
    assert allowed == [-2, 2]


def main() -> None:
    check_beltrami_and_cross_curl()
    check_action_and_stage_ledgers()
    check_factorial_heat_flatness()
    print("C177 same-curl unforced reservoir checks passed")


if __name__ == "__main__":
    main()
