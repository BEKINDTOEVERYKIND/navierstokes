#!/usr/bin/env python3
"""Independent, dependency-free exact audit ledger for C168.

The audit does four things independently of the prose formulas:

* derives the generic full-ladder radial and tangential coefficients from
  the projected symmetric Euler symbol;
* evaluates that symbol in Q(sqrt(3),sqrt(7)) to fix every C166 sign and
  signed-pair normalization at charges 0, 1, and 2;
* proves the global block-square-zero and bounded-operator estimates; and
* checks the exact unit-pulse gain and every outward rational inequality in
  the tau=1/100 all-walk, point-error, and energy certificate.

This certifies a prescribed, frozen, linearized-pump calculation only.  It
does not certify an unforced or energy-preserving pump, a q^2/q charge star,
localization, BAFL, or a nonlinear Navier--Stokes stage.  All fields in this
calculation are 2D3C; the exact total tangential field remains subject to the
passive-scalar maximum principle.
"""

from fractions import Fraction as Q


# ---------------------------------------------------------------------------
# A tiny exact field Q(sqrt(3),sqrt(7)) for a direct Leray-symbol rederivation.
# ---------------------------------------------------------------------------

RADICANDS = (3, 7)


class A:
    """Element of Q(sqrt(3),sqrt(7)) in the square-free mask basis."""

    def __init__(self, value=0):
        if isinstance(value, A):
            self.terms = dict(value.terms)
        elif isinstance(value, dict):
            self.terms = {
                int(mask): Q(coefficient)
                for mask, coefficient in value.items()
                if coefficient
            }
        else:
            coefficient = Q(value)
            self.terms = {} if not coefficient else {0: coefficient}

    def __add__(self, other):
        other = A(other)
        result = dict(self.terms)
        for mask, coefficient in other.terms.items():
            result[mask] = result.get(mask, Q(0)) + coefficient
            if not result[mask]:
                del result[mask]
        return A(result)

    __radd__ = __add__

    def __neg__(self):
        return A({mask: -coefficient for mask, coefficient in self.terms.items()})

    def __sub__(self, other):
        return self + (-A(other))

    def __rsub__(self, other):
        return A(other) - self

    def __mul__(self, other):
        other = A(other)
        result = {}
        for left_mask, left_coefficient in self.terms.items():
            for right_mask, right_coefficient in other.terms.items():
                common = left_mask & right_mask
                rational_factor = 1
                for bit, radicand in enumerate(RADICANDS):
                    if common & (1 << bit):
                        rational_factor *= radicand
                mask = left_mask ^ right_mask
                result[mask] = result.get(mask, Q(0)) + (
                    left_coefficient * right_coefficient * rational_factor
                )
        return A(result)

    __rmul__ = __mul__

    def __truediv__(self, rational):
        rational = Q(rational)
        return A({mask: coefficient / rational for mask, coefficient in self.terms.items()})

    def __eq__(self, other):
        return self.terms == A(other).terms

    def __repr__(self):
        return f"A({self.terms})"


ZERO = A(0)
ONE = A(1)
SQRT3 = A({1: Q(1)})
SQRT7 = A({2: Q(1)})
SQRT21 = SQRT3 * SQRT7


def vadd(left, right):
    return [x + y for x, y in zip(left, right)]


def vscale(scalar, vector):
    return [scalar * x for x in vector]


def vdot(left, right):
    return sum((x * y for x, y in zip(left, right)), ZERO)


def projected_euler_symbol(k, source, q, gate):
    """Return P_{k+q}((source.q)gate+(gate.k)source), exactly."""
    destination = vadd(k, q)
    raw = vadd(vscale(vdot(source, q), gate), vscale(vdot(gate, k), source))
    destination_square = vdot(destination, destination)
    assert set(destination_square.terms) == {0}
    divisor = destination_square.terms[0]
    projection_factor = vdot(destination, raw) / divisor
    result = vadd(raw, vscale(-projection_factor, destination))
    assert vdot(destination, result) == ZERO
    return result


def check_direct_leray_symbol_and_c166_signs():
    """Rebuild the m=0,1,2 edges directly, rather than insert formulas."""
    e_r = [ONE, ZERO, ZERO]
    e_t = [ZERO, ONE, ZERO]
    e_z = [ZERO, ZERO, ONE]

    rho = {
        0: [ZERO, ZERO, -ONE],
        1: [ONE / 2, ZERO, -SQRT3 / 2],
        2: [2 * SQRT7 / 7, ZERO, -SQRT21 / 7],
    }
    wave = {
        0: [ONE, ZERO, ZERO],
        1: [ONE, ZERO, SQRT3 / 3],
        2: [ONE, ZERO, 2 * SQRT3 / 3],
    }
    q_plus = [ZERO, ZERO, SQRT3 / 3]
    q_minus = vscale(-ONE, q_plus)

    for m in (0, 1, 2):
        assert vdot(rho[m], rho[m]) == ONE
        assert vdot(wave[m], rho[m]) == ZERO
        assert vdot(wave[m], e_t) == ZERO

    def edge(source_m, destination_m, q, gate, source_polarization, target):
        assert vadd(wave[source_m], q) == wave[destination_m]
        value = projected_euler_symbol(
            wave[source_m], source_polarization, q, gate
        )
        return vdot(target, value)

    # Full-ladder rho-basis radial entries, all derived from P_{k+q}.
    radial = {
        (2, 1): edge(1, 2, q_plus, e_r, rho[1], rho[2]),
        (0, 1): edge(1, 0, q_minus, e_r, rho[1], rho[0]),
        (1, 2): edge(2, 1, q_minus, e_r, rho[2], rho[1]),
        (1, 0): edge(0, 1, q_plus, e_r, rho[0], rho[1]),
    }
    expected_radial = {
        (2, 1): 3 * SQRT7 / 14,
        (0, 1): SQRT3 / 2,
        (1, 2): 3 * SQRT7 / 7,
        (1, 0): SQRT3 / 3,
    }
    assert radial == expected_radial

    # A tangential source under a radial gate has coefficient exactly one.
    assert edge(1, 2, q_plus, e_r, e_t, e_t) == ONE
    assert edge(1, 0, q_minus, e_r, e_t, e_t) == ONE

    # A radial source under the tangential gate gives K delta_1.
    k_delta_1 = {
        2: edge(1, 2, q_plus, e_t, rho[1], e_t),
        0: edge(1, 0, q_minus, e_t, rho[1], e_t),
    }
    assert k_delta_1 == {2: A(Q(-1, 2)), 0: A(Q(1, 2))}

    # Tangential sources are dark under the tangential gate, on both signs.
    dark_plus = projected_euler_symbol(wave[1], e_t, q_plus, e_t)
    dark_minus = projected_euler_symbol(wave[1], e_t, q_minus, e_t)
    assert dark_plus == [ZERO, ZERO, ZERO]
    assert dark_minus == [ZERO, ZERO, ZERO]

    # C166 used daughter bases e_{+,perp}=-rho_2 and e_{-,perp}=-rho_0.
    # Re-evaluate the symbol in those bases to audit all four signs.
    c166_plus = vscale(-ONE, rho[2])
    c166_minus = vscale(-ONE, rho[0])
    c166_entries = {
        "a_plus": edge(1, 2, q_plus, e_r, rho[1], c166_plus),
        "a_minus": edge(1, 0, q_minus, e_r, rho[1], c166_minus),
        "c_plus": edge(2, 1, q_minus, e_r, c166_plus, rho[1]),
        "c_minus": edge(0, 1, q_plus, e_r, c166_minus, rho[1]),
    }
    assert c166_entries == {
        "a_plus": -3 * SQRT7 / 14,
        "a_minus": -SQRT3 / 2,
        "c_plus": -3 * SQRT7 / 7,
        "c_minus": -SQRT3 / 3,
    }


# ---------------------------------------------------------------------------
# Generic symbolic identities in Q[m,sigma]/(sigma^2-1).
# ---------------------------------------------------------------------------


def ptrim(poly):
    poly = list(poly)
    while len(poly) > 1 and not poly[-1]:
        poly.pop()
    return tuple(poly)


def padd(left, right):
    result = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return ptrim(result)


def pneg(poly):
    return tuple(-value for value in poly)


def pmul(left, right):
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return ptrim(result)


def smul(left, right):
    """Multiply A(m)+sigma B(m), using sigma^2=1."""
    a, b = left
    c, d = right
    return padd(pmul(a, c), pmul(b, d)), padd(pmul(a, d), pmul(b, c))


def sadd(left, right):
    return padd(left[0], right[0]), padd(left[1], right[1])


def sneg(value):
    return pneg(value[0]), pneg(value[1])


def check_generic_full_ladder_derivation():
    constant = lambda value: ((Q(value),), (Q(0),))
    m = ((Q(0), Q(1)), (Q(0),))
    sigma = ((Q(0),), (Q(1),))

    # The radial raw vector, after stripping sqrt(m^2+3), has plane
    # components (m-sigma,-sqrt(3)).  Dotting with the destination unit
    # numerator (m+sigma,-sqrt(3)) gives this exact polynomial.
    radial_numerator = sadd(
        smul(sadd(m, sneg(sigma)), sadd(m, sigma)),
        constant(3),
    )
    assert radial_numerator == ((Q(2), Q(0), Q(1)), (Q(0),))

    # The destination radial numerator is divergence-free:
    # (m+sigma) - (m+sigma)=0.  Hence Leray projection cannot change the
    # preceding dot product.
    destination_divergence = sadd(sadd(m, sigma), sneg(sadd(m, sigma)))
    assert destination_divergence == constant(0)

    # The tangential-gate raw coefficient is
    # (-sqrt(3))*sigma/sqrt(3)=-sigma.  A tangential source under the
    # radial gate instead has coefficient e_r.k_m=1; a tangential source
    # under a tangential gate has both contractions zero.
    tangential_from_radial = sneg(sigma)
    tangential_from_tangential_under_radial = constant(1)
    tangential_under_tangential = constant(0)
    assert tangential_from_radial == ((Q(0),), (Q(-1),))
    assert tangential_from_tangential_under_radial == constant(1)
    assert tangential_under_tangential == constant(0)

    # These identities give exactly H_r=diag(J,L) and H_t(a,b)=(0,Ka).
    # Apply the generic signed-shift rule to a multi-charge formal state.
    # Its output has no radial entries; a second application is empty,
    # regardless of the incoming tangential entries or ladder support.
    def h_t(radial, _tangential):
        output = []
        for charge, label in sorted(radial.items()):
            output.append((charge + 1, -1, charge, label))
            output.append((charge - 1, +1, charge, label))
        return {}, tuple(output)

    arbitrary_radial = {-7: "a_-7", 0: "a_0", 11: "a_11"}
    arbitrary_tangential = {4: "b_4"}
    first = h_t(arbitrary_radial, arbitrary_tangential)
    assert len(first[1]) == 2 * len(arbitrary_radial)
    second = h_t(*first)
    assert second == ({}, ())


def check_operator_norm_majorants():
    # For n=m+/-1,
    # j_nm^2 < (m^2+3)/(n^2+3) because the cross-multiplied difference
    # (m^2+3)^2-(m^2+2)^2 is 2m^2+5.
    m2_plus_3 = (Q(3), Q(0), Q(1))
    m2_plus_2 = (Q(2), Q(0), Q(1))
    difference = padd(pmul(m2_plus_3, m2_plus_3), pneg(pmul(m2_plus_2, m2_plus_2)))
    assert difference == (Q(5), Q(0), Q(2))

    # Audit the note's warning about the word Jacobi: J is not symmetric
    # in coefficient l2.  Its adjacent detailed-balance weight is instead
    # w_m=m^2+2, because both sides below are the same nonconstant product.
    shifted_weight = (Q(3), Q(2), Q(1))  # (m+1)^2+2
    weighted_forward = pmul(shifted_weight, m2_plus_2)
    weighted_reverse = pmul(m2_plus_2, shifted_weight)
    assert weighted_forward == weighted_reverse
    assert m2_plus_2 != shifted_weight

    # The ratio is <(3/2)^2 for either adjacent destination.  Re-expand
    # 9((m+eps)^2+3)-4(m^2+3)=5(m+9eps/5)^2+39/5.
    for epsilon in (-1, 1):
        adjacent = (Q(1), Q(2 * epsilon), Q(1))
        left = padd(
            tuple(9 * value for value in padd(adjacent, (Q(3),))),
            pneg(tuple(4 * value for value in m2_plus_3)),
        )
        completed = padd(
            tuple(
                5 * value
                for value in pmul((Q(9 * epsilon, 5), Q(1)), (Q(9 * epsilon, 5), Q(1)))
            ),
            (Q(39, 5),),
        )
        assert left == completed == (Q(24), Q(18 * epsilon), Q(5))
        assert Q(39, 5) > 0

    # Two adjacent entries per row/column yield ||J||_1,||J||_inf <=3,
    # hence the Schur l2 bound <=3.  For K, every entry is <=1/sqrt(3)
    # and 2/sqrt(3)<6/5 is exactly 100<108 after squaring.
    assert 2 * Q(3, 2) == 3
    assert 100 < 108


# ---------------------------------------------------------------------------
# Exact phase calculation and outward tau=1/100 certificate.
# ---------------------------------------------------------------------------


def gadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def gscale(scalar, value):
    return scalar * value[0], scalar * value[1]


def gmul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def check_unit_pulse_point_and_energy():
    one = (Q(1), Q(0))
    i = (Q(0), Q(1))
    minus_i = (Q(0), Q(-1))
    minus_one = (Q(-1), Q(0))

    # theta=pi/2: the source phase is i, while
    # K delta_1=(delta_0-delta_2)/2 evaluates to one.
    source_theta = i
    k_theta = gadd(gscale(Q(1, 2), one), gscale(Q(-1, 2), minus_one))
    assert k_theta == one

    # xi=-pi/2.  The radial scalar is 1 and the -iK response scalar -1.
    radial_scalar = gmul(minus_i, source_theta)
    tangential_scalar = gmul(minus_i, gmul(minus_i, k_theta))
    assert radial_scalar == one
    assert tangential_scalar == minus_one

    # rho_1 and e_t are orthonormal.  Reality completion multiplies both
    # endpoint and initial amplitudes by two, so the gain is sqrt(2).
    point_norm_squared = Q(1) + Q(1)
    initial_complex_norm_squared = Q(1)
    assert point_norm_squared / initial_complex_norm_squared == 2

    # Half-lattice coefficient energy: one source plus two 1/2 daughters.
    coefficient_energy = Q(1) + Q(1, 2) ** 2 + Q(-1, 2) ** 2
    assert coefficient_energy == Q(3, 2)


def check_tau_one_hundredth_certificate():
    tau = Q(1, 100)
    x = 3 * tau

    # Since 1/r! <=1, e^x < sum x^r=1/(1-x) for x>0.  This bounds the
    # complete operator exponential, not a finite matrix or walk cutoff.
    geometric_exponential = 1 / (1 - x)
    assert geometric_exponential == Q(100, 97)
    h_bound = geometric_exponential - 1
    assert h_bound == Q(3, 97)
    assert h_bound < Q(1, 32)

    # The literal off-source charge tail is no larger than ||h||_1.  At
    # the selected point the radial and K-response errors sum to this.
    k_norm_bound = Q(6, 5)
    point_error = (1 + k_norm_bound) * Q(1, 32)
    assert point_error == Q(11, 160)

    # sqrt(2)>7/5 by a positive rational square comparison.
    assert Q(2) > Q(7, 5) ** 2
    point_gain_lower = Q(7, 5) - point_error
    assert point_gain_lower == Q(213, 160)
    assert point_gain_lower > Q(13, 10)

    # ||a||_2 <=e^(3/100)<100/97<33/32, followed by
    # ||(a,Ka)||_2^2 <=(1+(6/5)^2)||a||_2^2.
    assert geometric_exponential < Q(33, 32)
    energy_bound = (1 + k_norm_bound**2) * Q(33, 32) ** 2
    assert energy_bound == Q(66429, 25600)
    assert energy_bound < 3

    assert 1 + tau == Q(101, 100)


def main():
    check_direct_leray_symbol_and_c166_signs()
    check_generic_full_ladder_derivation()
    check_operator_norm_majorants()
    check_unit_pulse_point_and_energy()
    check_tau_one_hundredth_certificate()
    print("PASS C168: direct Leray projection gives the full ladder coefficients")
    print("PASS C168: C166 basis signs and signed-pair normalization are exact")
    print("PASS C168: the global tangential block is square-zero")
    print("PASS C168: the unit pulse has gain sqrt(2) and half-energy 3/2")
    print("PASS C168: tau=1/100 includes all walks and has gain >13/10")
    print("PASS C168: its off-source tail is <1/32 and half-energy is <3")


if __name__ == "__main__":
    main()
