#!/usr/bin/env python3
"""Directed cellwise Metzler certificate for the strengthened C159 gain.

The landed C159 checker certifies the phase/covector path and a global
cooperative cone.  This checker reruns that certificate, reevaluates all
four raw transverse-generator entries on 2048 cells, and propagates a
rigorous componentwise lower solution without integrating an unstable
amplitude column.

On a cell, let L be the matrix of entrywise lower endpoints.  Both the
actual B(t) and L are Metzler.  With alpha=max(0,-L_11,-L_22)+10^-60 and
C=L+alpha I >= 0,

  exp(L dt) = exp(-alpha dt) exp(C dt)
            >= (1-alpha dt) sum_{r=0}^6 (C dt)^r/r!.

All entries propagated below are nonnegative, so directed-down Decimal
rounding preserves a valid lower vector.  The positive Duhamel comparison
then puts the true time-dependent cell flow above this vector.
"""

from decimal import Decimal as D
from fractions import Fraction as F

import zero_drift_cooperative_cone_c159 as c


def raw_coefficients(phase, gamma, period, beta, sqrt2, sqrt3):
    ca, sa, cb, sb = phase
    sab = sa * cb + ca * sb
    cab = ca * cb - sa * sb
    fa = -sa - c.I(c.DELTA) * sab
    fb = -sb - c.I(c.DELTA) * sab
    hessian = [[
        -ca * c.R1[i] * c.R1[j]
        - cb * c.R2[i] * c.R2[j]
        - c.I(c.DELTA) * cab * c.RS[i] * c.RS[j]
        for j in range(3)
    ] for i in range(3)]
    gradient = [fa * c.R1[i] + fb * c.R2[i] for i in range(3)]
    velocity = c.vector_cross((c.I(1), c.I(1), c.I(1)), gradient)
    h_value = c.square_norm(gradient)
    p_vector = [c.I(c.CARRIER / 3) * velocity[i] / h_value + gamma * gradient[i] for i in range(3)]
    normal = tuple(c.I(1) / sqrt3 for _ in range(3))
    m_value = sqrt3 * beta
    d_value = c.square_norm(p_vector)
    q_value = d_value + m_value * m_value
    tangent = c.vector_cross(p_vector, normal)
    velocity_gradient = [[
        sum((c.I(c.CN[i][inner]) * hessian[inner][j] for inner in range(3)), c.I(0))
        - sqrt2 * gradient[j]
        for j in range(3)
    ] for i in range(3)]
    p_s_t = c.vector_dot(p_vector, c.matrix_vector(velocity_gradient, tangent))
    l_p = c.vector_dot(normal, c.matrix_vector(velocity_gradient, p_vector))
    p_s_p = c.vector_dot(p_vector, c.matrix_vector(velocity_gradient, p_vector))
    b11 = m_value * l_p / d_value
    b22 = (2 * p_s_p + m_value * l_p) / d_value
    b21 = m_value * m_value * sqrt2 * c.I(c.CARRIER) / (q_value * d_value)
    b12 = (2 * m_value * p_s_t + sqrt2 * c.I(c.CARRIER) * (d_value - m_value * m_value)) / d_value
    return tuple(period * value for value in (b11, b12, b21, b22))


def lower_step(z1, z2, entries, dt_lo, dt_hi):
    l11, l12, l21, l22 = (entry.lo for entry in entries)
    assert l12 > 0 and l21 > 0
    alpha = c.UP.add(
        max(D(0), c.UP.minus(l11), c.UP.minus(l22)),
        D("1e-60"),
    )
    x = c.UP.multiply(alpha, dt_hi)
    assert x < 1
    scalar = c.DOWN.subtract(D(1), x)
    c11 = c.DOWN.add(l11, alpha)
    c22 = c.DOWN.add(l22, alpha)
    assert c11 >= 0 and c22 >= 0, (l11, l22, alpha, c11, c22)
    cdt = (
        (c.DOWN.multiply(dt_lo, c11), c.DOWN.multiply(dt_lo, l12)),
        (c.DOWN.multiply(dt_lo, l21), c.DOWN.multiply(dt_lo, c22)),
    )

    def mmul(a, b):
        return tuple(tuple(
            c.DOWN.add(c.DOWN.multiply(a[i][0], b[0][j]), c.DOWN.multiply(a[i][1], b[1][j]))
            for j in range(2)
        ) for i in range(2))

    series = ((D(1), D(0)), (D(0), D(1)))
    term = series
    for order in range(1, 7):
        raw = mmul(term, cdt)
        term = tuple(tuple(c.DOWN.divide(raw[i][j], D(order)) for j in range(2)) for i in range(2))
        series = tuple(tuple(c.DOWN.add(series[i][j], term[i][j]) for j in range(2)) for i in range(2))
    q11 = c.DOWN.multiply(scalar, series[0][0])
    q12 = c.DOWN.multiply(scalar, series[0][1])
    q21 = c.DOWN.multiply(scalar, series[1][0])
    q22 = c.DOWN.multiply(scalar, series[1][1])
    next1 = c.DOWN.add(c.DOWN.multiply(q11, z1), c.DOWN.multiply(q12, z2))
    next2 = c.DOWN.add(c.DOWN.multiply(q21, z1), c.DOWN.multiply(q22, z2))
    return next1, next2, alpha


def main():
    c.exact_structural_checks()
    period, beta, sqrt2 = c.parameter_intervals()
    records, *_ = c.generate_reference(period, beta, sqrt2)
    c.certify_path(records, period, beta, sqrt2)

    panels = len(records)
    subdivisions = 64
    step = D(1) / panels
    total_cells = panels * subdivisions
    dt_lo = c.DOWN.divide(D(1), D(total_cells))
    dt_hi = c.UP.divide(D(1), D(total_cells))
    phase_radius = D("2e-6")
    gamma_radius = D("8e-4")
    sqrt3 = c.sqrt_fraction_bound(c.F(3))
    z1, z2 = D(1), D(3) / 20
    maximum_alpha = D(0)
    minima = [D("1e100")] * 4
    for coefficients in records:
        polys = [c.interval_polynomial(row) for row in coefficients]
        for subcell in range(subdivisions):
            left = step * D(subcell) / subdivisions
            right = step * D(subcell + 1) / subdivisions
            argument = c.I(left, right)
            phase = [c.evaluate_polynomial(polys[index], argument).widen(phase_radius) for index in range(4)]
            gamma = c.evaluate_polynomial(polys[5], argument).widen(gamma_radius)
            entries = raw_coefficients(phase, gamma, period, beta, sqrt2, sqrt3)
            for index, entry in enumerate(entries):
                minima[index] = min(minima[index], entry.lo)
            z1, z2, alpha = lower_step(z1, z2, entries, dt_lo, dt_hi)
            maximum_alpha = max(maximum_alpha, alpha)

    ratio2 = c.DOWN.divide(z2, D(3) / 20)
    assert maximum_alpha < 5
    assert z1 > 3000
    assert ratio2 > 3000

    # A rational exponential tail proves exp(8)<3000 without floats.
    # After term 8^32/32!, every successive ratio is at most 8/33<1/4.
    term = F(1)
    exponential_upper = F(1)
    for order in range(1, 33):
        term *= F(8, order)
        exponential_upper += term
    exponential_upper += term * F(8, 33) * F(4, 3)
    assert exponential_upper < 3000

    # Exact clock coefficients on q=n^8 from rho(M)>exp(8).
    assert F(3, 8) * 8 == 3                 # q^(3/8)=n^3
    assert F(57, 16) * 8 == F(57, 2)       # H q^(3/8)=n^(57/2)
    assert F(76, 25) * F(3, 8) == F(57, 50)
    assert F(57, 50) / 8 == F(57, 400)
    assert F(76, 25) * F(57, 16) == F(1083, 100)
    assert F(1083, 100) / 8 == F(1083, 800)
    assert 40 * F(57, 400) == F(57, 10)
    assert 40 * F(1083, 800) == F(1083, 20)

    # C188's log(n)<=6 n^(1/14), the shortened collars, and heat charge.
    assert 6 * F(57, 16) + 1 < 23
    assert F(7, 2) * F(1, 14) - F(1, 2) == F(-1, 4)
    assert F(3, 2) * F(1, 14) + 8 + F(15, 2) == F(437, 28)
    assert 2 * F(76, 25) == F(152, 25)
    assert 1 / (1 - F(1, 100)) == F(100, 99)

    # Relative to q^(3/8), the q^(-1/2) and q^(-1) remainder thresholds.
    first_order_threshold = F(7, 8) / F(57, 400)
    order_minus_one_threshold = F(11, 8) / F(57, 400)
    assert first_order_threshold == F(350, 57)
    assert order_minus_one_threshold == F(550, 57)
    assert -F(7, 8) + 6 * F(57, 400) == F(-1, 50)
    assert -F(11, 8) + 6 * F(57, 400) == F(-13, 25)

    print("C192 strong zero-drift gain certificate: PASS")
    print("raw entry lower bounds:", *(str(value) for value in minima))
    print("maximum shift alpha:", maximum_alpha)
    print("componentwise lower endpoint:", z1, z2)
    print("ratios to w:", z1, ratio2)
    print("CONCLUSION: M w > 3000 w and log rho(M) > 8 per period")
    print("BOUNDARY: principal cocycle/operator norm only; no finite-band or viscous bridge")


if __name__ == "__main__":
    main()
