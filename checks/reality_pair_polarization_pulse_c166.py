#!/usr/bin/env python3
"""Dependency-free exact checker for C166's reality-pair pulse.

The checker verifies:

* the signed-pair Lie-bracket coefficient as a rational identity;
* the complete y=1 six-coordinate radial/tangential matrices, rederived
  from the three-dimensional projected Euler symbol;
* the positive radial energy metric and square-zero active operator;
* the exact radial-quarter/tangential-unit endpoint;
* physical coefficient-energy increase and the coherent reality-completed
  point-gain identity, both before and after equal-energy rescaling.

This is a limiting-symbol, first-neighbour finite-dimensional check.  It
does not certify invariance of the full charge ladder or a Navier--Stokes
stage.
"""

from fractions import Fraction as Q


RADICANDS = (2, 3, 7)


class R:
    """Q(sqrt(2),sqrt(3),sqrt(7)) in the square-free mask basis."""

    def __init__(self, terms=()):
        if isinstance(terms, (int, Q)):
            terms = {0: Q(terms)}
        elif isinstance(terms, dict):
            terms = {mask: Q(value) for mask, value in terms.items()}
        else:
            terms = dict(terms)
        self.terms = {mask: value for mask, value in terms.items() if value}

    def __add__(self, other):
        other = as_r(other)
        out = dict(self.terms)
        for mask, value in other.terms.items():
            out[mask] = out.get(mask, Q(0)) + value
            if not out[mask]:
                del out[mask]
        return R(out)

    __radd__ = __add__

    def __neg__(self):
        return R({mask: -value for mask, value in self.terms.items()})

    def __sub__(self, other):
        return self + (-as_r(other))

    def __rsub__(self, other):
        return as_r(other) - self

    def __mul__(self, other):
        other = as_r(other)
        out = {}
        for left_mask, left_value in self.terms.items():
            for right_mask, right_value in other.terms.items():
                common = left_mask & right_mask
                factor = 1
                for bit, radicand in enumerate(RADICANDS):
                    if common & (1 << bit):
                        factor *= radicand
                mask = left_mask ^ right_mask
                out[mask] = out.get(mask, Q(0)) + (
                    left_value * right_value * factor
                )
        return R(out)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = Q(scalar)
        return R({mask: value / scalar for mask, value in self.terms.items()})

    def __eq__(self, other):
        return self.terms == as_r(other).terms

    def __repr__(self):
        return f"R({self.terms})"


def as_r(value):
    return value if isinstance(value, R) else R(value)


ONE = R(1)
ZERO = R(0)
SQRT2 = R({1: Q(1)})
SQRT3 = R({2: Q(1)})
SQRT7 = R({4: Q(1)})
SQRT14 = SQRT2 * SQRT7
SQRT42 = SQRT2 * SQRT3 * SQRT7


def zero_matrix(rows, columns=None):
    if columns is None:
        columns = rows
    return [[ZERO for _ in range(columns)] for _ in range(rows)]


def matmul(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    assert len(left[0]) == middle
    out = zero_matrix(rows, columns)
    for i in range(rows):
        for j in range(columns):
            value = ZERO
            for k in range(middle):
                value += left[i][k] * right[k][j]
            out[i][j] = value
    return out


def matsub(left, right):
    return [
        [a - b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def transpose(value):
    return [list(row) for row in zip(*value)]


def matvec(matrix, vector):
    return [
        sum((entry * value for entry, value in zip(row, vector)), ZERO)
        for row in matrix
    ]


def diagonal(entries):
    out = zero_matrix(len(entries))
    for index, value in enumerate(entries):
        out[index][index] = as_r(value)
    return out


def vadd(left, right):
    return [a + b for a, b in zip(left, right)]


def vneg(value):
    return [-entry for entry in value]


def vscale(scalar, value):
    return [scalar * entry for entry in value]


def vdot(left, right):
    return sum((a * b for a, b in zip(left, right)), ZERO)


def projected_symbol(k, a, q, b):
    """The exact symmetric Euler symbol P_{k+q}((a.q)b+(b.k)a)."""
    total = vadd(k, q)
    raw = vadd(vscale(vdot(a, q), b), vscale(vdot(b, k), a))
    norm_square = vdot(total, total)
    assert set(norm_square.terms) <= {0}
    projection_coefficient = vdot(total, raw) / norm_square.terms[0]
    projected = vadd(raw, vneg(vscale(projection_coefficient, total)))
    assert vdot(total, projected) == ZERO
    return projected


def matrix_block(row_basis, column_basis, operator):
    return [
        [vdot(row, operator(column)) for column in column_basis]
        for row in row_basis
    ]


def check_general_bracket_identity():
    # After a common denominator, the signed-pair coefficient is
    # y[(4-y)D_+-(y+4)D_-].  Expand it and compare with 2y^2(4-y^2).
    # Polynomials are low-to-high tuples of rational coefficients.
    def pmul(left, right):
        out = [Q(0)] * (len(left) + len(right) - 1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                out[i + j] += a * b
        return tuple(out)

    def psub(left, right):
        out = [Q(0)] * max(len(left), len(right))
        for i, value in enumerate(left):
            out[i] += value
        for i, value in enumerate(right):
            out[i] -= value
        return tuple(out)

    y = (Q(0), Q(1))
    d_plus = (Q(4), Q(2), Q(1))
    d_minus = (Q(4), Q(-2), Q(1))
    numerator = pmul(
        y,
        psub(
            pmul((Q(4), Q(-1)), d_plus),
            pmul((Q(4), Q(1)), d_minus),
        ),
    )
    claimed = (Q(0), Q(0), Q(8), Q(0), Q(-2))
    assert numerator == claimed
    assert pmul(d_plus, d_minus) == (Q(16), Q(0), Q(4), Q(0), Q(1))


def claimed_y_one_matrices():
    # Coordinate order:
    # source sigma, source t, daughter + radial, daughter + t,
    # daughter - radial, daughter - t.
    hr = zero_matrix(6)
    ht = zero_matrix(6)

    a_plus = -3 * SQRT7 / 14
    c_plus = -3 * SQRT7 / 7
    b_plus = R(Q(-1, 2))
    d_plus = -SQRT7 / 7

    a_minus = -SQRT3 / 2
    c_minus = -SQRT3 / 3
    b_minus = R(Q(1, 2))
    d_minus = SQRT3 / 3

    hr[0][2], hr[2][0] = c_plus, a_plus
    hr[0][4], hr[4][0] = c_minus, a_minus
    hr[1][3] = hr[3][1] = ONE
    hr[1][5] = hr[5][1] = ONE

    ht[1][2], ht[3][0] = d_plus, b_plus
    ht[1][4], ht[5][0] = d_minus, b_minus
    return hr, ht


def derive_y_one_matrices_from_symbol():
    """Independently build every entry from the 3D projected symbol."""
    e_r = [ONE, ZERO, ZERO]
    e_t = [ZERO, ONE, ZERO]
    e_z = [ZERO, ZERO, ONE]
    e_sigma = [ONE / 2, ZERO, -SQRT3 / 2]
    p = [ONE, ZERO, SQRT3 / 3]

    assert vdot(e_sigma, e_sigma) == ONE
    assert vdot(e_t, e_t) == ONE
    assert vdot(p, e_sigma) == ZERO
    assert vdot(p, e_t) == ZERO

    source_basis = [e_sigma, e_t]
    geometry = {}
    hr = zero_matrix(6)
    ht = zero_matrix(6)

    for eta, offset in ((1, 2), (-1, 4)):
        g = vscale(R(Q(eta)) * SQRT3 / 3, e_z)
        daughter_wavevector = vadd(p, g)
        if eta == 1:
            e_perp = [-2 * SQRT7 / 7, ZERO, SQRT3 * SQRT7 / 7]
        else:
            e_perp = list(e_z)
        daughter_basis = [e_perp, e_t]

        assert vdot(e_perp, e_perp) == ONE
        assert vdot(daughter_wavevector, e_perp) == ZERO
        assert vdot(daughter_wavevector, e_t) == ZERO

        def forward(gate):
            return matrix_block(
                daughter_basis,
                source_basis,
                lambda source: projected_symbol(p, source, g, gate),
            )

        def reverse(gate):
            return matrix_block(
                source_basis,
                daughter_basis,
                lambda daughter: projected_symbol(
                    daughter_wavevector, daughter, vneg(g), gate
                ),
            )

        for assembled, gate in ((hr, e_r), (ht, e_t)):
            f_block = forward(gate)
            r_block = reverse(gate)
            for row in range(2):
                for column in range(2):
                    assembled[offset + row][column] = f_block[row][column]
                    assembled[row][offset + column] = r_block[row][column]

        geometry[eta] = (daughter_wavevector, e_perp)

    claimed_hr, claimed_ht = claimed_y_one_matrices()
    assert hr == claimed_hr
    assert ht == claimed_ht
    return hr, ht, geometry


def check_matrices_metric_and_bracket():
    hr, ht, geometry = derive_y_one_matrices_from_symbol()

    # The radial metric has daughter-radial weights 2 and 2/3.
    metric = diagonal((1, 1, 2, 1, Q(2, 3), 1))
    assert matmul(metric, hr) == matmul(transpose(hr), metric)

    ht_squared = matmul(ht, ht)
    assert ht_squared == zero_matrix(6)
    assert ht != zero_matrix(6)

    bracket = matsub(matmul(hr, ht), matmul(ht, hr))
    assert bracket[1][0] == R(Q(2, 7))

    source = [ONE, ZERO, ZERO, ZERO, ZERO, ZERO]
    bracket_source = matvec(bracket, source)
    assert bracket_source[1] == R(Q(2, 7))
    assert all(
        value == ZERO for index, value in enumerate(bracket_source)
        if index != 1
    )
    return hr, ht, metric, geometry


def qform(metric, vector):
    weighted = matvec(metric, vector)
    return sum((a * b for a, b in zip(vector, weighted)), ZERO)


def euclidean_square(vector):
    return sum((value * value for value in vector), ZERO)


def reality_completed_value(real_part, imaginary_part, cosine, sine):
    """Return 2 Re((real_part+i imaginary_part)(cosine+i sine))."""
    return vscale(
        2,
        vadd(vscale(cosine, real_part), vscale(-sine, imaginary_part)),
    )


def check_exact_two_pulse_endpoint(hr, ht, metric):
    # lambda=sqrt(8/7)=2sqrt(14)/7 and 1/lambda=sqrt(14)/4.
    lam = 2 * SQRT14 / 7
    inv_lam = SQRT14 / 4
    assert lam * lam == R(Q(8, 7))
    assert lam * inv_lam == ONE

    source = [ONE, ZERO, ZERO, ZERO, ZERO, ZERO]
    radial_velocity = matvec(hr, source)
    assert matvec(hr, radial_velocity) == [
        R(Q(8, 7)), ZERO, ZERO, ZERO, ZERO, ZERO
    ]

    # At the radial quarter turn, source=0 and daughter coordinates are
    # -i*a_eta/lambda.  Store the real magnitudes after factoring common i,
    # but derive them from H_r rather than inserting the endpoint formula.
    x_plus = 3 * SQRT2 / 8
    x_minus = SQRT42 / 8
    assert x_plus * x_plus == R(Q(9, 32))
    assert x_minus * x_minus == R(Q(21, 32))

    radial_real = vscale(-inv_lam, radial_velocity)
    assert radial_real == [ZERO, ZERO, x_plus, ZERO, x_minus, ZERO]
    assert euclidean_square(radial_real) == R(Q(15, 16))
    assert qform(metric, radial_real) == ONE

    # The tangential unit pulse adds -i H_t to the common-i radial state,
    # hence a real source-t component H_t*radial_real.
    mixed = matvec(ht, radial_real)
    assert mixed[1] == SQRT14 / 14
    assert all(
        value == ZERO for index, value in enumerate(mixed)
        if index != 1
    )

    endpoint_energy = euclidean_square(radial_real) + euclidean_square(mixed)
    assert endpoint_energy == R(Q(113, 112))
    endpoint_metric = qform(metric, radial_real) + qform(metric, mixed)
    assert endpoint_metric == R(Q(15, 14))

    # Exact finite-ordering defect: H_r H_t e_sigma=0, while
    # H_t H_r e_sigma=-(2/7)e_t.
    hr_ht_source = matvec(hr, matvec(ht, source))
    ht_hr_source = matvec(ht, matvec(hr, source))
    assert hr_ht_source == [ZERO] * 6
    assert ht_hr_source[1] == R(Q(-2, 7))

    # At the radial quarter turn and s=1, U_t U_r has real part `mixed`
    # and imaginary part `radial_real`.  In the reverse order H_t e_sigma
    # is H_r-dark, so U_r U_t has zero real part and imaginary part
    # radial_real-H_t e_sigma.  This directly verifies the finite ordering
    # defect, not merely its mixed derivative at the origin.
    tangential_daughters = matvec(ht, source)
    assert matvec(hr, tangential_daughters) == [ZERO] * 6
    reverse_imaginary = [
        radial - tangent
        for radial, tangent in zip(radial_real, tangential_daughters)
    ]
    ordering_real = mixed
    ordering_imaginary = [
        forward - reverse
        for forward, reverse in zip(radial_real, reverse_imaginary)
    ]
    assert ordering_real[1] == SQRT14 / 14
    assert ordering_imaginary == tangential_daughters

    return radial_real, mixed, endpoint_energy


def check_physical_point_gain(geometry, radial_real, mixed, endpoint_energy):
    # At z=0, endpoint radial daughters sum to
    # i[-3sqrt(14)/28 e_r + 5sqrt(42)/28 e_z].
    e_plus = geometry[1][1]
    e_minus = geometry[-1][1]
    coherent = vadd(
        vscale(radial_real[2], e_plus),
        vscale(radial_real[4], e_minus),
    )
    radial = -3 * SQRT14 / 28
    vertical = 5 * SQRT42 / 28
    assert coherent == [radial, ZERO, vertical]
    gain_squared = vdot(coherent, coherent)
    assert gain_squared == R(Q(3, 2))

    # At z=0 take the common horizontal phase exp(i theta)=-i.  Derive the
    # reality-completed value from its real/imaginary coefficient parts.
    e_t = [ZERO, ONE, ZERO]
    zero_vector = [ZERO, ZERO, ZERO]
    source_coefficient = vscale(mixed[1], e_t)
    source_point = reality_completed_value(
        source_coefficient, zero_vector, ZERO, -ONE
    )
    daughter_point = reality_completed_value(
        zero_vector, coherent, ZERO, -ONE
    )
    endpoint_point = vadd(source_point, daughter_point)
    assert source_point == zero_vector
    assert endpoint_point == vscale(2, coherent)
    assert mixed[1] == SQRT14 / 14
    endpoint_point_square = euclidean_square(endpoint_point)

    # A real unit source coefficient reaches its physical supremum at
    # phase one.  This also checks the source basis normalization used in
    # the gain ratio.
    e_sigma = [ONE / 2, ZERO, -SQRT3 / 2]
    initial_point = reality_completed_value(
        e_sigma, zero_vector, ONE, ZERO
    )
    initial_supremum_square = euclidean_square(initial_point)
    assert initial_supremum_square == R(4)
    assert (
        endpoint_point_square / initial_supremum_square.terms[0]
        == gain_squared
    )
    assert gain_squared == R(Q(3, 2))

    # The raw endpoint has coefficient energy 113/112, so the displayed
    # gain is not free energy.  After rescaling the entire endpoint back to
    # the initial coefficient energy, the exhibited point-gain squared is
    # still (3/2)/(113/112)=168/113>1.
    assert endpoint_energy == R(Q(113, 112))
    normalized_gain_squared = Q(3, 2) / Q(113, 112)
    assert normalized_gain_squared == Q(168, 113)
    assert normalized_gain_squared > 1

    # Certify the action bound without floating point.  Machin's identity
    # pi=16 atan(1/5)-4 atan(1/239), together with alternating-series
    # upper/lower truncations, gives the rational upper enclosure below.
    upper_atan_fifth = Q(1, 5) - Q(1, 5) ** 3 / 3 + Q(1, 5) ** 5 / 5
    lower_atan_239 = Q(1, 239) - Q(1, 239) ** 3 / 3
    pi_upper = 16 * upper_atan_fifth - 4 * lower_atan_239
    assert pi_upper < Q(22, 7)
    assert Q(14) < Q(15, 4) ** 2
    assert pi_upper * Q(15, 4) / 8 < Q(3, 2)


def main():
    check_general_bracket_identity()
    hr, ht, metric, geometry = check_matrices_metric_and_bracket()
    radial_real, mixed, endpoint_energy = check_exact_two_pulse_endpoint(
        hr, ht, metric
    )
    check_physical_point_gain(geometry, radial_real, mixed, endpoint_energy)
    print("PASS C166: the signed-pair radial/tangential bracket is exact")
    print("PASS C166: the radial metric is positive and H_t is square-zero")
    print("PASS C166: the y=1 two-pulse endpoint has active branch mixing")
    print("PASS C166: raw point-gain squared is 3/2; equal-energy is 168/113")


if __name__ == "__main__":
    main()
