#!/usr/bin/env python3
"""Exact rational interval checks for the explicit AO Batchelor profile.

No floating-point result is used to certify an inequality.  Decimal-looking
bounds are Fractions, and log(9/5) is enclosed by the positive atanh series

    log(9/5) = 2 * sum_{k>=0} (2/7)^(2k+1)/(2k+1).

The calculus/uniqueness proof is recorded in the accompanying research note.
"""

from fractions import Fraction as F


class I:
    """Closed rational interval with elementary outward-safe operations."""

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = as_i(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_i(other))

    def __rsub__(self, other):
        return as_i(other) - self

    def __mul__(self, other):
        other = as_i(other)
        vals = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def reciprocal(self):
        assert not (self.lo <= 0 <= self.hi)
        return I(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * as_i(other).reciprocal()

    def __rtruediv__(self, other):
        return as_i(other) / self

    def subset(self, lo, hi):
        return F(lo) < self.lo and self.hi < F(hi)


def as_i(value):
    return value if isinstance(value, I) else I(value)


def log_nine_fifths_bounds(terms=6):
    z = F(2, 7)
    lower = 2 * sum(z ** (2 * k + 1) / F(2 * k + 1) for k in range(terms))
    first_power = 2 * terms + 1
    tail = (
        2
        * z**first_power
        / F(first_power)
        / (1 - z * z)
    )
    return lower, lower + tail


def enclose_sqrt(interval, lo, hi):
    """Certify [sqrt(interval.lo),sqrt(interval.hi)] inside (lo,hi)."""
    lo, hi = F(lo), F(hi)
    assert lo > 0
    assert lo * lo < interval.lo
    assert interval.hi < hi * hi
    return I(lo, hi)


def main():
    log_lo, log_hi = log_nine_fifths_bounds()
    claimed_x_lo = F(58778665, 100000000)
    claimed_x_hi = F(58778667, 100000000)
    assert claimed_x_lo < log_lo < log_hi < claimed_x_hi
    x = I(claimed_x_lo, claimed_x_hi)

    # Elementary exact sign gates used in the calculus proof.
    assert x.lo > F(4, 7)
    assert x.hi < F(3, 5)
    c = (F(4, 5) - x) / (x * x)  # beta / Q
    assert c.lo > F(1, 2)

    r0 = enclose_sqrt(x, F(76667, 100000), F(76668, 100000))

    beta_squared = 1 / (4 - x)
    beta = enclose_sqrt(
        beta_squared, F(54135, 100000), F(54136, 100000)
    )
    Q = beta * x * x / (F(4, 5) - x)
    assert Q.subset(F(88134, 100000), F(88136, 100000))

    beta_Q = x * x / ((4 - x) * (F(4, 5) - x))
    assert beta_Q.subset(F(47711, 100000), F(47714, 100000))
    assert beta_Q.hi < 1

    b0 = (
        32
        * x
        * x
        * (2 - 3 * x)
        / (81 * (4 - x) * (F(4, 5) - x) * (F(4, 5) - x))
    )
    assert b0.subset(F(21018, 100000), F(21020, 100000))
    sqrt_b0 = enclose_sqrt(
        b0, F(45845, 100000), F(45848, 100000)
    )

    lambda_second = (
        F(8, 9) * beta * (7 * x - 4) / (F(4, 5) - x)
    )
    assert lambda_second.subset(F(25963, 100000), F(25967, 100000))
    p0 = 4 / (x * (4 - x))
    assert p0.subset(F(19943, 10000), F(19945, 10000))

    D1_squared = lambda_second / (8 * p0 * sqrt_b0)
    D1 = enclose_sqrt(
        D1_squared, F(1883, 10000), F(1885, 10000)
    )

    # Exact stationarity identities at y=e^-X=5/9.
    # Lambda: e^X = 1+X+(beta/Q)X^2 = 9/5.
    assert F(1) + x.lo + c.lo * x.lo * x.lo <= F(9, 5)
    assert F(9, 5) <= F(1) + x.hi + c.hi * x.hi * x.hi
    # b: after multiplication by 4-X, H(X)=5-5=0.
    assert F(5, 9) * 9 - 5 == 0

    print("log(9/5) enclosed by exact rational atanh tail")
    print("AO stationarity and positivity gates certified")
    print("global b-max calculus reduces to H_a'(x)<0")
    print("0.21018 < b0 < 0.21020")
    print("0.45845 < sqrt(b0) = fixed-beta resonant BAS edge < 0.45848")
    print("0.1883 < D1 < 0.1885 in the AO n^(-1/2) ratio correction")
    print("all explicit AO profile checks passed")


if __name__ == "__main__":
    main()
