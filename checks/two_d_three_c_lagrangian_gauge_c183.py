#!/usr/bin/env python3
"""Dependency-free checks for C183's exact 2D3C gauge/return identities.

The checker certifies finite-dimensional algebra and power ledgers only.
It does not turn the Kelvin principal system into a finite-frequency PDE
estimate or construct the common Floquet frame required by CFFC.
"""

from __future__ import annotations

from fractions import Fraction as F
import cmath
import math


V2 = tuple[F, F]
V3 = tuple[F, F, F]
M2 = tuple[tuple[F, F], tuple[F, F]]
M3 = tuple[tuple[F, F, F], tuple[F, F, F], tuple[F, F, F]]


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def scale(scalar, vector):
    return tuple(scalar * value for value in vector)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def cross(left: V3, right: V3) -> V3:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def transpose(matrix):
    return tuple(tuple(matrix[j][i] for j in range(len(matrix)))
                 for i in range(len(matrix[0])))


def mat_vec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def mat_mul(left, right):
    right_t = transpose(right)
    return tuple(tuple(dot(row, column) for column in right_t) for row in left)


def mat_add(left, right):
    return tuple(tuple(a + b for a, b in zip(lrow, rrow))
                 for lrow, rrow in zip(left, right))


def mat_scale(scalar, matrix):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def identity(size: int):
    return tuple(tuple(F(int(i == j)) for j in range(size)) for i in range(size))


def inv2(matrix: M2) -> M2:
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    assert determinant != 0
    return (
        (matrix[1][1] / determinant, -matrix[0][1] / determinant),
        (-matrix[1][0] / determinant, matrix[0][0] / determinant),
    )


def matrix_power(matrix, exponent: int):
    assert exponent >= 0
    result = identity(len(matrix))
    factor = matrix
    power = exponent
    while power:
        if power & 1:
            result = mat_mul(result, factor)
        factor = mat_mul(factor, factor)
        power //= 2
    return result


def passive_covector_reduction() -> None:
    """Check d(D^T p)/dt=-m g0 at arbitrary exact states."""
    cases: tuple[tuple[M2, M2, V2, V2, F, F], ...] = (
        (
            ((F(2), F(-1)), (F(3), F(-2))),
            ((F(2), F(1)), (F(1), F(1))),
            (F(3), F(-2)), (F(1), F(4)), F(5, 3), F(7, 4),
        ),
        (
            ((F(-1), F(5, 2)), (F(-3), F(1))),
            ((F(3), F(-2)), (F(-1), F(1))),
            (F(-2), F(5)), (F(7), F(-1)), F(-4, 5), F(3, 2),
        ),
    )
    for m_planar, deformation, g0, p0, charge, time in cases:
        assert m_planar[0][0] + m_planar[1][1] == 0
        inverse_transpose = transpose(inv2(deformation))
        g = mat_vec(inverse_transpose, g0)
        label_covector = sub(p0, scale(charge * time, g0))
        p = mat_vec(inverse_transpose, label_covector)

        p_dot = sub(scale(F(-1), mat_vec(transpose(m_planar), p)),
                    scale(charge, g))
        deformation_dot = mat_mul(m_planar, deformation)
        lhs = add(mat_vec(transpose(deformation_dot), p),
                  mat_vec(transpose(deformation), p_dot))
        assert lhs == scale(-charge, g0)
        assert mat_vec(transpose(deformation), p) == label_covector


def gauge_plane_wave_algebra() -> None:
    """Check the covariant frequency and Piola divergence identities."""
    deformation: M2 = ((F(2), F(1)), (F(1), F(1)))  # determinant one
    inverse = inv2(deformation)
    inverse_transpose = transpose(inverse)
    p0 = (F(7), F(-3))
    g0 = (F(2), F(5))
    charge = F(3)
    time = F(4, 7)
    covariant_label = sub(p0, scale(charge * time, g0))
    physical_p = mat_vec(inverse_transpose, covariant_label)
    horizontal_amplitude = (F(5), F(-2))
    pulled_horizontal = mat_vec(inverse, horizontal_amplitude)
    vertical_amplitude = -dot(physical_p, horizontal_amplitude) / charge

    # The covariant Piola divergence equals the physical Fourier divergence.
    assert dot(covariant_label, pulled_horizontal) == dot(
        physical_p, horizontal_amplitude
    )
    assert dot(covariant_label, pulled_horizontal) + charge * vertical_amplitude == 0

    # The removed scalar phase has exactly unit modulus.
    for theta0 in (-7.5, -1.0, 0.0, 2.25, 11.0):
        phase = cmath.exp(-1j * float(charge) * float(time) * theta0)
        assert abs(abs(phase) - 1.0) < 2.0e-15


def kelvin_constraint_and_vorticity() -> None:
    matrices: tuple[M3, ...] = (
        (
            (F(1), F(2), F(-3)),
            (F(4), F(-2), F(5)),
            (F(-1), F(7), F(1)),
        ),
        (
            (F(-3), F(1), F(2)),
            (F(5), F(4), F(-1)),
            (F(6), F(-2), F(-1)),
        ),
    )
    for gradient in matrices:
        assert sum(gradient[i][i] for i in range(3)) == 0
        for k, seed in (
            ((F(2), F(-1), F(3)), (F(1), F(4), F(-2))),
            ((F(-3), F(5), F(2)), (F(2), F(-1), F(6))),
        ):
            amplitude = cross(k, seed)
            assert dot(k, amplitude) == 0
            r2 = dot(k, k)
            aa = mat_vec(gradient, amplitude)
            sigma = dot(k, aa)
            k_dot = scale(F(-1), mat_vec(transpose(gradient), k))
            a_dot = add(scale(F(-1), aa), scale(2 * sigma / r2, k))
            assert dot(k_dot, amplitude) + dot(k, a_dot) == 0

            vorticity = cross(k, amplitude)
            vorticity_dot = add(cross(k_dot, amplitude), cross(k, a_dot))
            curl_base = (
                gradient[2][1] - gradient[1][2],
                gradient[0][2] - gradient[2][0],
                gradient[1][0] - gradient[0][1],
            )
            corrected_cauchy = add(
                mat_vec(gradient, vorticity),
                scale(dot(curl_base, k), amplitude),
            )
            assert vorticity_dot == corrected_cauchy


def connection_trace_identity() -> None:
    """Verify tr(B)=k.A.k/|k|^2=-d log|k| without choosing a frame."""
    gradients: tuple[M3, ...] = (
        (
            (F(2), F(-1), F(4)),
            (F(3), F(-5), F(2)),
            (F(1), F(6), F(3)),
        ),
        (
            (F(-1), F(2), F(0)),
            (F(4), F(3), F(-2)),
            (F(5), F(1), F(-2)),
        ),
    )
    for gradient in gradients:
        assert sum(gradient[i][i] for i in range(3)) == 0
        for k in ((F(2), F(-1), F(3)), (F(-4), F(5), F(1))):
            r2 = dot(k, k)
            k_a = tuple(dot(k, tuple(gradient[i][j] for i in range(3)))
                        for j in range(3))
            kelvin = tuple(
                tuple(-gradient[i][j] + 2 * k[i] * k_a[j] / r2
                      for j in range(3))
                for i in range(3)
            )
            scalar = dot(k, mat_vec(gradient, k))
            trace_full = sum(kelvin[i][i] for i in range(3))
            normal_quadratic = dot(k, mat_vec(kelvin, k)) / r2
            trace_plane = trace_full - normal_quadratic
            assert trace_full == 2 * scalar / r2
            assert normal_quadratic == scalar / r2
            assert trace_plane == scalar / r2
            r2_dot = -2 * scalar
            assert trace_plane == -r2_dot / (2 * r2)


def return_resonance() -> None:
    """Check the general return formula and exact unipotent drift ledger."""
    for shear in (F(1, 3), F(2), F(-5)):
        nilpotent: M2 = ((F(0), shear), (F(0), F(0)))
        assert mat_mul(nilpotent, nilpotent) == ((F(0), F(0)), (F(0), F(0)))
        return_map = mat_add(identity(2), nilpotent)
        inverse_return = mat_add(identity(2), mat_scale(F(-1), nilpotent))
        assert mat_mul(return_map, inverse_return) == identity(2)

        p0 = (F(7), F(-4))
        g0 = (F(3), F(5))
        charge = F(2)
        period = F(7, 3)
        for count in range(0, 12):
            inverse_power = matrix_power(inverse_return, count)
            direct_g = mat_vec(inverse_power, g0)
            direct_p = mat_vec(
                inverse_power,
                sub(p0, scale(charge * count * period, g0)),
            )
            formula_g = sub(g0, scale(count, mat_vec(nilpotent, g0)))
            formula_p = add(
                sub(p0, scale(count, add(
                    mat_vec(nilpotent, p0), scale(charge * period, g0)
                ))),
                scale(charge * count * count * period,
                      mat_vec(nilpotent, g0)),
            )
            assert direct_g == formula_g
            assert direct_p == formula_p

        # Construct exact resonance: N g=0 and N p=-mT g.
        invariant_g = (shear, F(0))
        resonant_p = (F(11), -charge * period)
        assert mat_vec(nilpotent, invariant_g) == (F(0), F(0))
        assert mat_vec(nilpotent, resonant_p) == scale(
            -charge * period, invariant_g
        )
        for count in range(0, 20):
            inverse_power = matrix_power(inverse_return, count)
            assert mat_vec(inverse_power, invariant_g) == invariant_g
            assert mat_vec(
                inverse_power,
                sub(resonant_p, scale(charge * count * period, invariant_g)),
            ) == resonant_p


def c152_stationary_scalar_specialization() -> None:
    """Check that Theta=H(f) reduces return resonance to one scalar plane."""
    gamma_f: V2 = (F(3), F(-2))
    # u_h is tangent to the level loop.
    u_h: V2 = (F(4), F(6))
    assert dot(gamma_f, u_h) == 0
    period = F(7, 3)
    scalar_derivative = F(-5, 4)
    charge = F(2)
    scalar_gradient = scale(scalar_derivative, gamma_f)

    planar_return = mat_add(identity(2), (
        (u_h[0] * gamma_f[0], u_h[0] * gamma_f[1]),
        (u_h[1] * gamma_f[0], u_h[1] * gamma_f[1]),
    ))
    return_transpose = transpose(planar_return)
    assert mat_vec(return_transpose, scalar_gradient) == scalar_gradient

    # Choose p with u_h.p=-m*T*c.  A gamma_f component is freely allowed.
    target = -charge * period * scalar_derivative
    tangent_coefficient = target / dot(u_h, u_h)
    p0 = add(scale(tangent_coefficient, u_h), scale(F(11, 7), gamma_f))
    assert dot(u_h, p0) == target
    assert mat_vec(return_transpose, p0) == sub(
        p0, scale(charge * period, scalar_gradient)
    )

    # The full return F=I+(u_h+Tc*n) outer (gamma_f,0) fixes k=(p,m)
    # by exactly the same scalar condition.
    u_full: V3 = (u_h[0], u_h[1], period * scalar_derivative)
    gamma_full: V3 = (gamma_f[0], gamma_f[1], F(0))
    rank_one: M3 = tuple(
        tuple(u_full[i] * gamma_full[j] for j in range(3))
        for i in range(3)
    )
    assert mat_mul(rank_one, rank_one) == mat_scale(F(0), identity(3))
    full_return = mat_add(identity(3), rank_one)
    wavevector: V3 = (p0[0], p0[1], charge)
    assert mat_vec(transpose(full_return), wavevector) == wavevector


def floquet_schedule_ledger() -> None:
    # epsilon=q^-1/2, N=O(log q): epsilon*N -> 0 and cannot be q^1/2.
    previous = None
    for exponent in range(4, 41):
        q = F(2) ** exponent
        # Use N=exponent, proportional to log_2(q).
        epsilon_times_returns = F(exponent, 2 ** (exponent // 2))
        if exponent >= 8 and previous is not None and exponent % 2 == 0:
            assert epsilon_times_returns < previous
        if exponent % 2 == 0:
            previous = epsilon_times_returns
        assert float(math.exp(float(epsilon_times_returns))) < math.sqrt(float(q))
    assert previous is not None and previous < F(1, 1000)


def main() -> None:
    passive_covector_reduction()
    gauge_plane_wave_algebra()
    kelvin_constraint_and_vorticity()
    connection_trace_identity()
    return_resonance()
    c152_stationary_scalar_specialization()
    floquet_schedule_ledger()
    print("C183 exact 2D3C Lagrangian gauge checks passed")


if __name__ == "__main__":
    main()
