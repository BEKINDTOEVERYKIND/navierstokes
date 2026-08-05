#!/usr/bin/env python3
"""Dependency-free checks for multicolour_core_buffer_audit.md."""

import math
import random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def mul(c, a):
    return tuple(c * x for x in a)


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def outer(a, b):
    return tuple(tuple(x * y for y in b) for x in a)


def madd(A, B):
    return tuple(tuple(x + y for x, y in zip(ar, br)) for ar, br in zip(A, B))


def mscale(c, A):
    return tuple(tuple(c * x for x in row) for row in A)


def sym_outer(a, b):
    return mscale(0.5, madd(outer(a, b), outer(b, a)))


def proj(k, v):
    return add(v, mul(-dot(k, v) / dot(k, k), k))


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def interaction(k, u, ell, v):
    out = add(mul(dot(u, ell), v), mul(dot(v, k), u))
    return proj(add(k, ell), out)


def mat_close(A, B, tol=2e-11):
    return max(abs(A[i][j] - B[i][j]) for i in range(3) for j in range(3)) < tol


def main():
    e1, e2, n = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    root = (-5 - math.sqrt(22)) / 3
    f = 1.5 * root + 3 + 0.5 / root
    assert abs(-0.5 * f - 1) < 2e-14 and root < -2

    frames = []
    for m, t in ((e1, e2), (e2, mul(-1, e1))):
        assert norm(add(cross(t, n), mul(-1, m))) < 1e-14
        er = mul(1 / math.sqrt(2), add(n, mul(-1, m)))
        eh = mul(-1 / math.sqrt(2), add(n, m))
        assert abs(dot(er, eh)) < 1e-14
        assert norm(add(add(er, mul(-1, eh)), mul(-math.sqrt(2), n))) < 1e-14
        atom = sym_outer(eh, er)
        target = mscale(0.5, madd(outer(m, m), mscale(-1, outer(n, n))))
        assert mat_close(atom, target)
        frames.append((m, t, er, eh))

    # Two atoms synthesize an arbitrary diagonal STF matrix in its eigenframe.
    for _ in range(100):
        l1, l2 = random.uniform(-3, 3), random.uniform(-3, 3)
        S = ((l1, 0.0, 0.0), (0.0, l2, 0.0), (0.0, 0.0, -l1 - l2))
        made = ((0.0, 0.0, 0.0),) * 3
        for lam, (_, _, er, eh) in zip((l1, l2), frames):
            made = madd(made, mscale(2 * lam, sym_outer(eh, er)))
        assert mat_close(made, S)

    # Exact within-colour incompressibility, pressure high sum, and low atom.
    Lam, Q = 11.0, 0.7
    delta = Q / (2 * Lam)
    modes = []
    for m, t, er, eh in frames:
        kp = add(mul(Lam, t), mul(Q / 2, er))
        km = add(mul(Lam, t), mul(-Q / 2, er))
        ap = add(mul(math.sqrt(2), n), mul(-delta, t))
        am = add(mul(math.sqrt(2), n), mul(delta, t))
        assert abs(dot(kp, ap)) < 1e-12 and abs(dot(km, am)) < 1e-12
        assert norm(interaction(kp, ap, km, am)) < 1e-11
        low = interaction(kp, ap, mul(-1, km), am)
        assert norm(add(low, mul(2 * Q, eh))) < 1e-11
        modes.append(((kp, ap), (km, am)))

    # Simultaneous colours have unavoidable O(Q) mixed outputs.
    leading_sum = []
    leading_diff = []
    for si, (k1, a1) in zip((1, -1), modes[0]):
        for tj, (k2, a2) in zip((1, -1), modes[1]):
            ps = interaction(k1, a1, k2, a2)
            pd = interaction(k1, a1, mul(-1, k2), a2)
            leading_sum.append((si, tj, norm(ps)))
            leading_diff.append((si, tj, norm(pd)))
    for si, tj, value in leading_sum:
        if si == tj:
            assert value > 0.5 * Q
    for si, tj, value in leading_diff:
        if si != tj:
            assert value > 0.5 * Q

    # Nonlinear relative-phase cell identities at arbitrary psi', psi''.
    # Use one of the orthonormal frames and a generic H (the identities do
    # not require H=-1).
    _, t, er, eh = frames[0]
    H = -1.0
    w = add(er, mul(H, eh))
    for _ in range(100):
        Lam = random.uniform(5.0, 30.0)
        psip = random.uniform(-4.0, 4.0)
        psipp = random.uniform(-3.0, 3.0)
        delta = psip / (2 * Lam)
        deltap = psipp / (2 * Lam)
        kp = add(mul(Lam, t), mul(psip / 2, er))
        km = add(mul(Lam, t), mul(-psip / 2, er))
        ap = add(w, mul(-delta, t))
        am = add(w, mul(delta, t))

        # Exact WKB divergence, including div(a_pm).
        div_ap = dot(mul(-deltap, t), er)
        div_am = dot(mul(deltap, t), er)
        assert abs(div_ap + dot(ap, kp)) < 1e-11
        assert abs(div_am + dot(am, km)) < 1e-11

        # Cross polarization derivatives cancel.
        deriv_cross = add(mul(dot(ap, er), mul(deltap, t)),
                          mul(dot(am, er), mul(-deltap, t)))
        assert norm(deriv_cross) < 1e-11

        # Cross high phase term is exactly parallel to 2 Lambda t.
        phase_high = add(mul(dot(ap, km), am), mul(dot(am, kp), ap))
        assert norm(proj(mul(2 * Lam, t), phase_high)) < 1e-11

        # Low phase term is 2 psi' (e_r+H e_h); its e_r part is gradient.
        phase_low = add(mul(dot(ap, mul(-1, km)), am), mul(dot(am, kp), ap))
        expected = mul(2 * psip, w)
        assert norm(add(phase_low, mul(-1, expected))) < 1e-11
        low_projected = proj(mul(psip, er), phase_low) if abs(psip) > 1e-8 else mul(2 * psip * H, eh)
        assert norm(add(low_projected, mul(-2 * psip * H, eh))) < 1e-10

        # The projected self-polarization derivative has the claimed size.
        selfp = proj(mul(2, kp), mul(-deltap, t))
        exact_norm = abs(deltap * delta) / math.sqrt(1 + delta * delta)
        assert abs(norm(selfp) - exact_norm) < 2e-11

    # The envelope-dark scalar equation does not make the full common high
    # sum pressure.  Take A=B=exp(i mu m.x), where m is orthogonal to the
    # common polarization n and carrier direction t.  Both envelope
    # directional derivatives vanish, but the product envelope tilts the
    # output wavevector away from t and leaves a nonzero Leray projection.
    m, t, er, eh = frames[0]
    Lam, Q, mu = 17.0, 0.8, 0.3
    delta = Q / (2 * Lam)
    w = add(er, mul(-1.0, eh))
    ap = add(w, mul(-delta, t))
    am = add(w, mul(delta, t))
    kp = add(mul(Lam, t), mul(Q / 2, er))
    km = add(mul(Lam, t), mul(-Q / 2, er))
    assert abs(dot(ap, m)) < 1e-14 and abs(dot(am, m)) < 1e-14
    leak = interaction(
        add(kp, mul(mu, m)), ap, add(km, mul(mu, m)), am
    )
    expected = Q * Q / Lam * mu / math.sqrt(Lam * Lam + mu * mu)
    assert abs(norm(leak) - expected) < 2e-11
    assert norm(leak) > 0

    print("PASS: dictionary, mixed-band scope, chirp identities, and envelope-pressure counterexample")


if __name__ == "__main__":
    main()
