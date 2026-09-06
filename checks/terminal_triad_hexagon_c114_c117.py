#!/usr/bin/env python3
"""Exact, dependency-free checks for reconstructed C114--C117.

All helical algebra is performed in Q(sqrt(2)) + i Q(sqrt(2)); no floating
point arithmetic or external symbolic package is used.
"""

from fractions import Fraction as F


class Q2:
    """a + b*sqrt(2), with exact rational a,b."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = F(a)
        self.b = F(b)

    def __add__(self, other):
        other = q2(other)
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-q2(other))

    def __rsub__(self, other):
        return q2(other) - self

    def __mul__(self, other):
        other = q2(other)
        return Q2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = q2(other)
        den = other.a * other.a - 2 * other.b * other.b
        assert den != 0
        return Q2(
            (self.a * other.a - 2 * self.b * other.b) / den,
            (self.b * other.a - self.a * other.b) / den,
        )

    def __eq__(self, other):
        other = q2(other)
        return self.a == other.a and self.b == other.b

    def __repr__(self):
        return f"Q2({self.a}, {self.b})"


def q2(x):
    return x if isinstance(x, Q2) else Q2(x)


SQRT2 = Q2(0, 1)
INV_SQRT2 = Q2(0, F(1, 2))


class C2:
    """Exact complexification of Q(sqrt(2))."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = q2(re)
        self.im = q2(im)

    def __add__(self, other):
        other = c2(other)
        return C2(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return C2(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-c2(other))

    def __rsub__(self, other):
        return c2(other) - self

    def __mul__(self, other):
        other = c2(other)
        return C2(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = c2(other)
        den = other.re * other.re + other.im * other.im
        assert den != Q2(0)
        return C2(
            (self.re * other.re + self.im * other.im) / den,
            (self.im * other.re - self.re * other.im) / den,
        )

    def conj(self):
        return C2(self.re, -self.im)

    def __eq__(self, other):
        other = c2(other)
        return self.re == other.re and self.im == other.im

    def __repr__(self):
        return f"C2({self.re}, {self.im})"


def c2(x):
    return x if isinstance(x, C2) else C2(x)


ZERO = C2(0)
I = C2(0, 1)


def vadd(x, y):
    return tuple(a + b for a, b in zip(x, y))


def vscale(a, x):
    a = c2(a)
    return tuple(a * z for z in x)


def vconj(x):
    return tuple(z.conj() for z in x)


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), ZERO)


def hdot(x, y):
    return sum((a.conj() * b for a, b in zip(x, y)), ZERO)


def cross(x, y):
    return (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )


def as_cvec(k):
    return tuple(C2(x) for x in k)


def idot(k, l):
    return sum(a * b for a, b in zip(k, l))


def iadd(k, l):
    return tuple(a + b for a, b in zip(k, l))


def ineg(k):
    return tuple(-a for a in k)


def helical(k, sign):
    """Unnormalised H_s(k)=n+i*s*(k cross n)/sqrt(2)."""
    n = (C2(1), C2(1), C2(1))
    kxn = cross(as_cvec(k), n)
    return vadd(n, vscale(C2(0, sign * INV_SQRT2), kxn))


def project(k, v):
    kc = as_cvec(k)
    return vadd(v, vscale(-dot(kc, v) / idot(k, k), kc))


def interaction(p, a, q, b):
    """The symmetrised coefficient -i P_(p+q)[(a.q)b+(b.p)a]."""
    k = iadd(p, q)
    raw = vadd(vscale(dot(a, as_cvec(q)), b), vscale(dot(b, as_cvec(p)), a))
    return vscale(-I, project(k, raw))


def ordered_interaction(p, a, q, b):
    """One ordered Fourier-convolution term -i P_(p+q)[(a.q)b]."""
    k = iadd(p, q)
    raw = vscale(dot(a, as_cvec(q)), b)
    if k == (0, 0, 0):
        # Here q=-p and transversality makes the raw coefficient zero, so
        # the undefined zero-mode Leray formula is never needed.
        assert raw == (ZERO, ZERO, ZERO)
        return raw
    return vscale(-I, project(k, raw))


def convolution(left, right):
    """Exact ordered projected convolution between two Fourier dictionaries."""
    outputs = {}
    for p, ap in left.items():
        for q, aq in right.items():
            k = iadd(p, q)
            term = ordered_interaction(p, ap, q, aq)
            outputs[k] = vadd(outputs.get(k, (ZERO, ZERO, ZERO)), term)
    return outputs


def norm_sq(v):
    z = hdot(v, v)
    assert z.im == Q2(0)
    return z.re


K1 = (1, -1, 0)
K2 = (0, 1, -1)
K3 = (-1, 0, 1)
ROOTS = (K1, K2, K3, ineg(K1), ineg(K2), ineg(K3))
N = (C2(1), C2(1), C2(1))


def check_c114_geometry():
    assert iadd(iadd(K1, K2), K3) == (0, 0, 0)
    for k in (K1, K2, K3):
        assert idot(k, k) == 2
    assert idot(K1, K2) == idot(K2, K3) == idot(K3, K1) == -1

    diffs = (iadd(K1, ineg(K2)), iadd(K2, ineg(K3)), iadd(K3, ineg(K1)))
    assert diffs == ((1, -2, 1), (1, 1, -2), (-2, 1, 1))
    assert all(idot(d, d) == 6 for d in diffs)


def check_helical_definitions():
    for k in ROOTS:
        for sign in (-1, 1):
            H = helical(k, sign)
            assert dot(as_cvec(k), H) == ZERO
            lhs = vscale(I, cross(as_cvec(k), H))
            rhs = vscale(sign * SQRT2, H)
            assert lhs == rhs
            assert norm_sq(H) == Q2(6)
            assert helical(ineg(k), sign) == tuple(z.conj() for z in H)


def check_c115_terminal_coefficients():
    expected_plus = vscale(3 * SQRT2, N)
    expected_minus = vscale(-3 * SQRT2, N)

    for p, q, third in ((K1, K2, K3), (K2, K3, K1), (K3, K1, K2)):
        target = ineg(third)
        for sign in (-1, 1):
            got = interaction(p, helical(p, sign), q, helical(q, sign))
            assert got == (ZERO, ZERO, ZERO)

        mixed_1 = interaction(p, helical(p, -1), q, helical(q, 1))
        mixed_2 = interaction(p, helical(p, 1), q, helical(q, -1))
        assert mixed_1 == expected_plus
        assert mixed_2 == expected_minus
        assert dot(as_cvec(target), mixed_1) == ZERO
        assert norm_sq(mixed_1) == Q2(54)

        # Unit input polarisations divide the bilinear output by 6.
        normalized_output = vscale(F(1, 6), mixed_1)
        assert norm_sq(normalized_output) == Q2(F(3, 2))

        # Squared projection onto either unit output helicity:
        # |H_out^* . (output/6)|^2 / |H_out|^2 = 3/4.
        for out_sign in (-1, 1):
            Hout = helical(target, out_sign)
            numerator = norm_sq((hdot(Hout, normalized_output),))
            assert numerator / Q2(6) == Q2(F(3, 4))


def unordered_root_pairs_with_sum(target):
    pairs = []
    for i, p in enumerate(ROOTS):
        for q in ROOTS[i:]:
            if iadd(p, q) == target:
                pairs.append((p, q))
    return pairs


def mode(k, x, y):
    return vadd(vscale(x, helical(k, 1)), vscale(y, helical(k, -1)))


def check_c116_paired_helicity_gate():
    d12 = iadd(K1, ineg(K2))
    expected_minus = vscale(-3 * SQRT2, N)
    expected_plus = vscale(3 * SQRT2, N)

    leak_1 = interaction(K1, helical(K1, -1), ineg(K2), helical(ineg(K2), 1))
    leak_2 = interaction(K1, helical(K1, 1), ineg(K2), helical(ineg(K2), -1))
    assert leak_1 == expected_minus
    assert leak_2 == expected_plus
    assert norm_sq(leak_1) == Q2(54)

    terminal = interaction(K1, helical(K1, -1), K2, helical(K2, 1))
    assert norm_sq(leak_1) == norm_sq(terminal)

    # Only {k1,-k2} inside the six roots sums to d12.  Hence no independent
    # internal quadratic channel can cancel it.
    pairs = unordered_root_pairs_with_sum(d12)
    assert len(pairs) == 1
    assert set(pairs[0]) == {K1, ineg(K2)}

    # Two helicity channels on the same wavevector pair cancel d12 while
    # retaining the terminal root.  (x1,y1)=(i,1), (x2,y2)=(-i,1).
    a1 = mode(K1, I, C2(1))
    a2 = mode(K2, -I, C2(1))
    a2_reality = tuple(z.conj() for z in a2)
    terminal_paired = interaction(K1, a1, K2, a2)
    difference_paired = interaction(K1, a1, ineg(K2), a2_reality)
    assert difference_paired == (ZERO, ZERO, ZERO)
    assert terminal_paired != (ZERO, ZERO, ZERO)

    # The terminal vector is parallel to n, hence has equal +/- helicities.
    assert cross(terminal_paired, N) == (ZERO, ZERO, ZERO)

    # Exhaust the full ordered convolution of the real four-mode input.
    # This verifies that no output omitted from the prose survives: only
    # the terminal pair +/- (k1+k2) is nonzero at first Picard order.
    support = {
        K1: a1,
        K2: a2,
        ineg(K1): vconj(a1),
        ineg(K2): vconj(a2),
    }
    outputs = convolution(support, support)
    nonzero = {k: value for k, value in outputs.items() if value != (ZERO, ZERO, ZERO)}
    kc = iadd(K1, K2)
    assert set(nonzero) == {kc, ineg(kc)}
    assert nonzero[ineg(kc)] == vconj(nonzero[kc])


def check_c117_second_picard_leakage():
    a1 = mode(K1, I, C2(1))
    a2 = mode(K2, -I, C2(1))
    kc = iadd(K1, K2)  # -k3
    child = N  # H_+(kc)+H_-(kc)=2n: equal amplitudes of both helicities.
    assert vadd(helical(kc, 1), helical(kc, -1)) == vscale(2, N)

    e1 = iadd(K1, kc)
    e2 = iadd(K2, kc)
    assert e1 == (2, -1, -1)
    assert e2 == (1, 1, -2)
    assert idot(e1, e1) == idot(e2, e2) == 6
    assert len(unordered_root_pairs_with_sum(e1)) == 1
    assert len(unordered_root_pairs_with_sum(e2)) == 1

    leak1 = interaction(K1, a1, kc, child)
    leak2 = interaction(K2, a2, kc, child)
    assert leak1 != (ZERO, ZERO, ZERO)
    assert leak2 != (ZERO, ZERO, ZERO)

    # Establish the general coefficient formula, rather than only checking
    # the explicit phase choice.  Same-helicity pairs are dark, and for
    # each parent the two mixed-helicity basis vectors are nonzero opposites.
    # Hence a parent (x_i,y_i) against an equal-amplitude child (c,c) has
    # coefficient c*(y_i-x_i) times a fixed nonzero vector.
    for parent in (K1, K2):
        same_plus = interaction(parent, helical(parent, 1), kc, helical(kc, 1))
        same_minus = interaction(parent, helical(parent, -1), kc, helical(kc, -1))
        mixed_yx = interaction(parent, helical(parent, -1), kc, helical(kc, 1))
        mixed_xy = interaction(parent, helical(parent, 1), kc, helical(kc, -1))
        assert same_plus == same_minus == (ZERO, ZERO, ZERO)
        assert mixed_yx != (ZERO, ZERO, ZERO)
        assert mixed_xy == vscale(-1, mixed_yx)
        assert cross(mixed_yx, N) == (ZERO, ZERO, ZERO)

    # Full second-Picard check for the explicit real gate.  The first
    # convolution has only +/-kc; differentiating the quadratic map gives
    # Q(u0,q1)+Q(q1,u0), and both named positive leaks survive exactly.
    support = {
        K1: a1,
        K2: a2,
        ineg(K1): vconj(a1),
        ineg(K2): vconj(a2),
    }
    first = {
        k: value
        for k, value in convolution(support, support).items()
        if value != (ZERO, ZERO, ZERO)
    }
    assert set(first) == {kc, ineg(kc)}
    left = convolution(support, first)
    right = convolution(first, support)
    second = {}
    for k in set(left) | set(right):
        second[k] = vadd(
            left.get(k, (ZERO, ZERO, ZERO)),
            right.get(k, (ZERO, ZERO, ZERO)),
        )
    assert second[e1] != (ZERO, ZERO, ZERO)
    assert second[e2] != (ZERO, ZERO, ZERO)

    # General coefficient identity.  With L_i=y_i-x_i,
    # T=y1*x2-x1*y2 = y2*L1-y1*L2.  Thus L1=L2=0 forces T=0.
    samples = (
        (F(2), F(5), F(-3), F(4)),
        (F(5, 2), F(-1, 3), F(7), F(11, 4)),
        (F(-9), F(-9), F(13, 2), F(13, 2)),
    )
    for x1, y1, x2, y2 in samples:
        L1, L2 = y1 - x1, y2 - x2
        T = y1 * x2 - x1 * y2
        assert T == y2 * L1 - y1 * L2
        if L1 == 0 and L2 == 0:
            assert T == 0


def main():
    check_c114_geometry()
    check_helical_definitions()
    check_c115_terminal_coefficients()
    check_c116_paired_helicity_gate()
    check_c117_second_picard_leakage()
    print("PASS C114: exact integer A2 shell and difference geometry")
    print("PASS C115: exact helical selection and terminal coefficients")
    print("PASS C116: paired-helicity cancellation retains terminal output")
    print("PASS C117: two-helicity child forces named second-Picard leaks")
    print("No bare-gate closure or Navier-Stokes stage map is claimed")


if __name__ == "__main__":
    main()
