#!/usr/bin/env python3
"""Dependency-free exact ledgers for C174's time-sliced rigid collar.

The checker verifies:

* exact divergence/core preservation under a collar ramp and its signed
  integrated switching identity;
* pressure-exact interpolation of symmetric affine cores;
* all q, h, b powers in the switching, residual, pressure, and rewrite
  ledgers;
* telescoping invariance of a sliced cocycle under exact recharts;
* the C154 nilpotent return powers and a finite-radius residence inequality;
* the distinction between a sum of rewrite norms and a possibly cancelling
  signed total.

It does not prove a finite-tube second-flow-jet bound, packet coherence,
the dynamic pressure/wake estimate, MCKC, LCE, BAFL, or an unforced stage.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log


Vec3 = tuple[F, F, F]
Mat3 = tuple[Vec3, Vec3, Vec3]
Mat2 = tuple[tuple[F, F], tuple[F, F]]

I3: Mat3 = (
    (F(1), F(0), F(0)),
    (F(0), F(1), F(0)),
    (F(0), F(0), F(1)),
)
I2: Mat2 = ((F(1), F(0)), (F(0), F(1)))


def dot(a: Vec3, b: Vec3) -> F:
    return sum((a[i] * b[i] for i in range(3)), F(0))


def vadd(a: Vec3, b: Vec3) -> Vec3:
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def vsub(a: Vec3, b: Vec3) -> Vec3:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def vscale(c: F, a: Vec3) -> Vec3:
    return tuple(c * a[i] for i in range(3))  # type: ignore[return-value]


def outer(a: Vec3, b: Vec3) -> Mat3:
    return tuple(
        tuple(a[i] * b[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def transpose3(a: Mat3) -> Mat3:
    return tuple(
        tuple(a[j][i] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mscale3(c: F, a: Mat3) -> Mat3:
    return tuple(
        tuple(c * a[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def madd3(a: Mat3, b: Mat3) -> Mat3:
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mmul3(a: Mat3, b: Mat3) -> Mat3:
    bt = transpose3(b)
    return tuple(
        tuple(dot(a[i], bt[j]) for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def mvec3(a: Mat3, x: Vec3) -> Vec3:
    return tuple(dot(a[i], x) for i in range(3))  # type: ignore[return-value]


def mpow3(a: Mat3, exponent: int) -> Mat3:
    result = I3
    for _ in range(exponent):
        result = mmul3(a, result)
    return result


def mmul2(a: Mat2, b: Mat2) -> Mat2:
    return (
        (
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ),
        (
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ),
    )


def inverse2(a: Mat2) -> Mat2:
    determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    assert determinant
    return (
        (a[1][1] / determinant, -a[0][1] / determinant),
        (-a[1][0] / determinant, a[0][0] / determinant),
    )


def trace2(a: Mat2) -> F:
    return a[0][0] + a[1][1]


def det2(a: Mat2) -> F:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def check_ramp_divergence_and_variation() -> None:
    # Two exact Fourier collar coefficients at one nonzero wavevector.
    # Both are transverse, so every convex ramp is divergence free and has
    # no zero Fourier mode.
    k = (F(2), F(-1), F(1))
    u0 = (F(1), F(1), F(-1))
    u1 = (F(3), F(1), F(-5))
    assert dot(k, u0) == 0
    assert dot(k, u1) == 0
    difference = vsub(u1, u0)
    assert dot(k, difference) == 0
    for theta in (F(0), F(1, 5), F(2, 3), F(1)):
        ramped = vadd(vscale(1 - theta, u0), vscale(theta, u1))
        assert dot(k, ramped) == 0

    # On theta(t)=t for t in [0,1], integral theta'(u1-u0) dt is exactly
    # u1-u0.  A slower/faster monotone affine ramp has the same identity.
    for duration in (F(1, 7), F(1), F(11, 3)):
        theta_prime = 1 / duration
        signed_integral = vscale(duration * theta_prime, difference)
        assert signed_integral == difference

    # Monotonicity is load-bearing for the unit total-variation bound.
    # A backtracking ramp has the same signed endpoint variation but a
    # strictly larger sum-of-norms budget.
    monotone_values = (F(0), F(1, 5), F(2, 3), F(1))
    backtracking_values = (F(0), F(1), F(1, 2), F(1))
    total_variation = lambda values: sum(
        (abs(right - left) for left, right in zip(values, values[1:])),
        F(0),
    )
    assert total_variation(monotone_values) == 1
    assert total_variation(backtracking_values) == 2
    assert backtracking_values[-1] - backtracking_values[0] == 1

    # Discrete profile bookkeeping: adjacent completions agree on the
    # first two "core" entries and differ only in three collar entries.
    profile0 = (F(4), F(-2), F(1), F(3), F(-5))
    profile1 = (F(4), F(-2), F(7), F(-1), F(2))
    profile_difference = tuple(
        profile1[index] - profile0[index] for index in range(5)
    )
    assert profile_difference[:2] == (F(0), F(0))
    assert profile_difference[2:] != (F(0), F(0), F(0))

    # Leray leaves a divergence-free switch coefficient unchanged.
    projection = vsub(
        difference,
        vscale(dot(k, difference) / dot(k, k), k),
    )
    assert projection == difference


def check_symmetric_affine_ramp_is_pressure_exact() -> None:
    s0: Mat3 = (
        (F(1), F(2), F(0)),
        (F(2), F(-1), F(1)),
        (F(0), F(1), F(0)),
    )
    s1: Mat3 = (
        (F(-2), F(0), F(3)),
        (F(0), F(1), F(-1)),
        (F(3), F(-1), F(1)),
    )
    for s in (s0, s1):
        assert transpose3(s) == s
        assert sum((s[i][i] for i in range(3)), F(0)) == 0
    s_prime = madd3(s1, mscale3(-1, s0))
    for theta in (F(0), F(1, 4), F(3, 5), F(1)):
        s_theta = madd3(mscale3(1 - theta, s0), mscale3(theta, s1))
        affine_force = madd3(s_prime, mmul3(s_theta, s_theta))
        assert transpose3(affine_force) == affine_force
        # A linear field Mx is a gradient exactly when M is symmetric.


def check_scale_powers_and_rewrite_cost() -> None:
    # q powers after one full backward h=q^(3/2).
    h_q = F(3, 2)
    assert F(-5, 2) + h_q == -1       # switch/parent/self
    assert F(-1, 2) + h_q == 1        # viscosity
    assert -5 + h_q == F(-7, 2)       # separated pressure
    assert -3 + h_q == F(-3, 2)       # separated viscous pressure

    # q=n^8, h=n^12, b=n^-2.
    assert F(12, 8) == h_q
    assert 8 * (-1) == -8
    assert 8 * F(-7, 2) == -28
    assert 8 * F(-3, 2) == -12

    # q^-1 log h is o(n^-6); its ratio to n^-6 is 12 log(n)/n^2.
    ratios = [12 * log(n) / n**2 for n in (20, 50, 100, 200)]
    assert all(left > right for left, right in zip(ratios, ratios[1:]))
    assert ratios[-1] < F(1, 100)

    for n in (2, 3, 7):
        q = F(n**8)
        h = F(n**12)
        b = F(1, n**2)
        # Avoid fractional Fraction powers: h^2=q^3 is the exact identity.
        assert h * h == q * q * q

        # Current child size is delta*b*g/h because q^(3/2)=h.
        for gain in (F(1), F(n**3), h):
            delta = F(2, 5)
            current = delta * b * gain / h
            remaining_gain = h / gain
            assert current * remaining_gain == delta * b

        # One order-one rewrite is b, versus the b^3 allowance.
        assert b / (b**3) == n**4

        deltas = (F(1, 3), F(2, 5), F(1, 7))
        sum_of_norms = b * sum(deltas, F(0))
        assert sum_of_norms == sum((b * delta for delta in deltas), F(0))
        # Signed responses may cancel even though the sum of norms is
        # positive, so the latter is not a lower bound.
        signed = b * F(1, 3) - b * F(1, 3)
        assert signed == 0
        assert 2 * b * F(1, 3) > 0


def check_dimensional_slice_ledgers() -> None:
    # Square every half-integer L2 scale so the complete dimensional
    # derivation can be checked with exact rational arithmetic. Set
    # a=ell=1; then r=q^-1, Lambda=lambda, and Re^-1=nu.
    for n in (2, 3, 5):
        q = F(n**8)
        h = F(n**12)
        radius = 1 / q
        lam = F(3, 7)
        reynolds_inverse = F(2, 11)
        slice_time = 1 / lam

        # Switch: h*lambda*r^(5/2), squared.
        switch_back_squared = h**2 * lam**2 * radius**5
        assert switch_back_squared == lam**2 / q**2

        # Parent: (a*lambda*r^(5/2))*lambda^-1, then h.
        parent_back_squared = (
            h**2
            * (lam**2 * radius**5)
            * slice_time**2
        )
        assert parent_back_squared == 1 / q**2

        # Self: (lambda^2*r^(5/2))*lambda^-1, then h.
        self_back_squared = (
            h**2
            * (lam**4 * radius**5)
            * slice_time**2
        )
        assert self_back_squared == lam**2 / q**2

        # Viscosity: (nu*lambda*r^(1/2))*lambda^-1, then h.
        viscous_back_squared = (
            h**2
            * (reynolds_inverse**2 * lam**2 * radius)
            * slice_time**2
        )
        assert viscous_back_squared == reynolds_inverse**2 * q**2

        # Instantaneous-zero-mean pressure moments integrated over a slice
        # and evaluated at separation d=ell=1.
        parent_point = radius**5
        self_point = lam * radius**5
        viscous_point = reynolds_inverse * radius**3
        assert parent_point == q**-5
        assert self_point == lam * q**-5
        assert viscous_point == reynolds_inverse * q**-3
        assert (h * parent_point) ** 2 == q**-7
        assert (h * self_point) ** 2 == lam**2 * q**-7
        assert (h * viscous_point) ** 2 == reynolds_inverse**2 * q**-3


def check_zero_mean_kernel_gain() -> None:
    # A one-dimensional scalar analogue of subtracting a degree -3 kernel
    # at the center. Equal and opposite point sources at +/-r have zero
    # mean. Their far field gains one inverse power: d^4*tail stays bounded
    # while d^3*tail tends to zero. The Leray tensor uses the same mean-value
    # subtraction componentwise.
    radius = F(1)
    scaled_fourth = []
    scaled_third = []
    for distance in (10, 20, 40, 80):
        d = F(distance)
        tail = 1 / (d - radius) ** 3 - 1 / (d + radius) ** 3
        assert tail > 0
        scaled_fourth.append(d**4 * tail)
        scaled_third.append(d**3 * tail)
    assert all(F(6) < value < F(7) for value in scaled_fourth)
    assert all(
        left > right for left, right in zip(scaled_third, scaled_third[1:])
    )
    assert scaled_third[-1] < F(1, 10)


def check_rechart_telescoping() -> None:
    p0: Mat2 = ((F(2), F(1)), (F(1), F(1)))
    p1: Mat2 = ((F(1), F(2)), (F(0), F(1)))
    p2: Mat2 = ((F(3), F(-1)), (F(1), F(0)))
    q0: Mat2 = ((F(1), F(1)), (F(0), F(1)))
    q1: Mat2 = ((F(2), F(0)), (F(1), F(1)))
    q2: Mat2 = ((F(1), F(-1)), (F(1), F(2)))
    q3: Mat2 = ((F(3), F(1)), (F(2), F(1)))
    ps = (p0, p1, p2)
    qs = (q0, q1, q2, q3)

    sliced = I2
    for index in range(3):
        tilde = mmul2(inverse2(qs[index + 1]), mmul2(ps[index], qs[index]))
        sliced = mmul2(tilde, sliced)

    physical = mmul2(p2, mmul2(p1, p0))
    expected = mmul2(inverse2(q3), mmul2(physical, q0))
    assert sliced == expected

    # Only a periodic endpoint chart gives similarity and hence preserves
    # trace/determinant (the two-dimensional characteristic polynomial).
    periodic = mmul2(inverse2(q0), mmul2(physical, q0))
    assert trace2(periodic) == trace2(physical)
    assert det2(periodic) == det2(physical)

    # Arbitrary endpoint charts can change coordinate norms and even the
    # characteristic polynomial; telescoping only preserves the represented
    # physical map. The identity physical map is enough to exhibit this.
    stretch: Mat2 = ((F(7), F(0)), (F(0), F(1)))
    two_sided_identity = mmul2(inverse2(stretch), I2)
    assert two_sided_identity == ((F(1, 7), F(0)), (F(0), F(1)))
    assert trace2(two_sided_identity) != trace2(I2)


def check_nilpotent_return_and_residence() -> None:
    u: Vec3 = (F(1), F(1), F(0))
    g: Vec3 = (F(1), F(-1), F(2))
    assert dot(g, u) == 0
    r_matrix = outer(u, g)
    assert mmul3(r_matrix, r_matrix) == mscale3(F(0), I3)
    f = madd3(I3, r_matrix)
    inverse_transpose = madd3(I3, mscale3(-1, transpose3(r_matrix)))

    for count in (0, 1, 2, 5, 11):
        return_matrix = madd3(I3, mscale3(F(count), r_matrix))
        assert mpow3(f, count) == return_matrix
        assert mpow3(inverse_transpose, count) == madd3(
            I3, mscale3(F(-count), transpose3(r_matrix))
        )
        infinity_norm = max(
            sum((abs(entry) for entry in row), F(0))
            for row in return_matrix
        )
        assert infinity_norm <= 1 + 4 * count

    # A J^-1 radius shrink pays J^-3 in volume and J in reciprocal width.
    for slice_count in (2, 5, 13):
        radius_ratio = F(1, slice_count)
        assert radius_ratio**3 == F(1, slice_count**3)
        assert 1 / radius_ratio == slice_count

    # Exact scalar version of the finite-radius sufficient condition.
    for n in (2, 5, 11):
        q = F(n**8)
        ell = F(1)
        radius = ell / q
        for distortion in (F(1), F(3), F(20)):
            initial_radius = radius / (4 * distortion)
            hessian_bound = 4 * distortion * distortion * q
            linear_part = distortion * initial_radius
            nonlinear_part = hessian_bound * initial_radius**2 / ell
            assert linear_part == radius / 4
            assert nonlinear_part == radius / 4
            assert linear_part + nonlinear_part <= radius


def main() -> None:
    check_ramp_divergence_and_variation()
    check_symmetric_affine_ramp_is_pressure_exact()
    check_scale_powers_and_rewrite_cost()
    check_dimensional_slice_ledgers()
    check_zero_mean_kernel_gain()
    check_rechart_telescoping()
    check_nilpotent_return_and_residence()
    print("PASS C174: exact rigid-collar ramp and backward scale ledgers")
    print("PASS C174: recharts telescope to the same physical child cocycle")
    print("PASS C174: C154 shear survives recharting; residence gate recorded")
    print("PASS C174: rewrite sum-of-norms costs b per relative unit")
    print("OPEN: finite-tube residence, MCKC(ii), LCE, BAFL, full stage")


if __name__ == "__main__":
    main()
