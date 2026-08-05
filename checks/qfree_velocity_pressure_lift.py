#!/usr/bin/env python3
"""Exact Gaussian-rational checks for the q-free velocity-pressure lift."""

from fractions import Fraction as F


class C:
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        if isinstance(re, C):
            self.re, self.im = re.re, re.im
        else:
            self.re, self.im = F(re), F(im)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, C) else C(value)

    def __add__(self, other):
        other = self.coerce(other)
        return C(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return C(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    __rmul__ = __mul__

    def conjugate(self):
        return C(self.re, -self.im)

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.re**2 + other.im**2
        if denominator == 0:
            raise ZeroDivisionError
        return self * other.conjugate() * C(1 / denominator)

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, exponent):
        if exponent == 2:
            return self * self
        raise ValueError("only square is used")

    def __eq__(self, other):
        other = self.coerce(other)
        return self.re == other.re and self.im == other.im

    def __repr__(self):
        return f"C({self.re}, {self.im})"


I = C(0, 1)
ZERO = C(0)


def check_case(case):
    (
        radius,
        beta,
        n_value,
        omega_spectral,
        velocity,
        velocity_prime,
        velocity_second,
        axial,
        axial_prime,
        axial_second,
        h_value,
        h_prime,
        h_second,
    ) = case

    alpha = beta * n_value
    omega_flow = velocity / radius
    omega_flow_prime = velocity_prime / radius - velocity / radius**2
    gamma_prime_value = velocity + radius * velocity_prime
    gamma_second_value = 2 * velocity_prime + radius * velocity_second
    lambda_value = beta * axial - omega_flow
    lambda_prime = beta * axial_prime - omega_flow_prime
    gamma = n_value * lambda_value - omega_spectral
    gamma_derivative = n_value * lambda_prime

    d_value = 1 + beta**2 * radius**2
    d_derivative = 2 * beta**2 * radius
    numerator = beta * radius**2 * axial_prime - gamma_prime_value
    numerator_derivative = (
        beta * (2 * radius * axial_prime + radius**2 * axial_second)
        - gamma_second_value
    )
    denominator = radius * d_value
    denominator_derivative = d_value + radius * d_derivative
    coefficient = numerator / denominator
    coefficient_derivative = (
        numerator_derivative * denominator
        - numerator * denominator_derivative
    ) / denominator**2
    a_value = radius * coefficient_derivative
    b_value = (
        -2
        * beta
        * velocity
        * (axial_prime + beta * gamma_prime_value)
        / d_value
    )

    u_r = h_value / radius
    pressure = (
        -I * gamma * radius * h_prime / (n_value**2 * d_value)
        + I * numerator * h_value / (n_value * radius * d_value)
    )
    u_theta = (
        I * gamma_prime_value * u_r / (radius * gamma)
        + n_value * pressure / (radius * gamma)
    )
    u_z = I * axial_prime * u_r / gamma - alpha * pressure / gamma

    # Differentiate P exactly from the three input two-jets.
    first_coefficient = -I * gamma * radius / (n_value**2 * d_value)
    first_coefficient_derivative = -I / n_value**2 * (
        (gamma_derivative * radius + gamma) / d_value
        - gamma * radius * d_derivative / d_value**2
    )
    second_coefficient = I * numerator / (n_value * radius * d_value)
    second_coefficient_derivative = (
        I * coefficient_derivative / n_value
    )
    pressure_derivative = (
        first_coefficient_derivative * h_prime
        + first_coefficient * h_second
        + second_coefficient_derivative * h_value
        + second_coefficient * h_prime
    )

    radial_residual = (
        I * gamma * u_r
        - 2 * omega_flow * u_theta
        + pressure_derivative
    )
    theta_residual = (
        I * gamma * u_theta
        + gamma_prime_value * u_r / radius
        - I * n_value * pressure / radius
    )
    axial_residual = (
        I * gamma * u_z + axial_prime * u_r + I * alpha * pressure
    )
    divergence = (
        h_prime / radius
        - I * n_value * u_theta / radius
        + I * alpha * u_z
    )

    radial_weight = radius / d_value
    radial_weight_derivative = (
        d_value - radius * d_derivative
    ) / d_value**2
    scalar_residual = (
        radial_weight_derivative * h_prime
        + radial_weight * h_second
        - n_value**2
        / radius
        * (1 + a_value / (n_value * gamma) + b_value / gamma**2)
        * h_value
    )

    assert theta_residual == ZERO
    assert axial_residual == ZERO
    assert divergence == ZERO
    assert I * n_value**2 * radial_residual / gamma == scalar_residual
    assert radial_residual == gamma * scalar_residual / (I * n_value**2)


def main():
    cases = (
        (
            F(7, 5), F(4, 3), 11, C(F(3, 7), F(5, 4)),
            F(5, 6), F(-7, 4), F(13, 9),
            F(2, 5), F(9, 8), F(-11, 6),
            F(3, 2), F(-5, 7), F(17, 10),
        ),
        (
            F(9, 4), F(2, 5), 17, C(F(-2, 3), F(11, 7)),
            F(-4, 9), F(13, 6), F(-5, 3),
            F(7, 8), F(-8, 9), F(19, 12),
            F(-5, 4), F(7, 11), F(23, 13),
        ),
        (
            F(5, 3), F(7, 8), 23, C(F(5, 8), F(-7, 6)),
            F(8, 7), F(2, 3), F(5, 11),
            F(-3, 10), F(13, 10), F(-17, 9),
            F(11, 6), F(-9, 14), F(29, 15),
        ),
    )
    for case in cases:
        check_case(case)
    print("azimuthal/axial momentum reconstruction: PASS")
    print("Fourier divergence reconstruction: PASS")
    print("scalar residual to radial force factor gamma/(i*n^2): PASS")
    print("all q-free velocity-pressure lift checks passed")


if __name__ == "__main__":
    main()
