#!/usr/bin/env python3
"""Checks for C181's static vertical-shear Kelvin propagator.

The exact algebra is Fraction-valued.  Floating calculations are only
diagnostics for the elementary arctangent closed form, direct ODE evolution,
and asymptotic orientation; they are not a PDE or finite-frequency proof.
"""

from fractions import Fraction as F
import math


def exact_kelvin_reduction() -> None:
    # Symbols are represented by many exact rational substitutions.
    values = (F(-3), F(-1), F(1, 2), F(2), F(5))
    for x in values:
        for y in values:
            for m in (F(1, 2), F(2), F(5)):
                for g in (F(1, 3), F(2), F(7)):
                    for u in (F(-2), F(1, 2), F(3)):
                        for v in (F(-1), F(2)):
                            r2 = x * x + y * y + m * m
                            w = -(x * u + y * v) / m
                            h = g * u
                            x_dot = -m * g
                            u_dot = 2 * m * x * h / r2
                            v_dot = 2 * m * y * h / r2
                            w_dot_constraint = -(
                                x_dot * u + x * u_dot + y * v_dot
                            ) / m
                            # Direct third Kelvin component:
                            # -g*u + 2*m*(m*g*u)/R^2.
                            w_dot_kelvin = -g * u + 2 * m * m * h / r2
                            assert w_dot_constraint == w_dot_kelvin
                            # (g*u*R^2)'=0.
                            r2_dot = 2 * x * x_dot
                            h_dot = g * u_dot
                            assert h_dot * r2 + h * r2_dot == 0
                            assert x * u + y * v + m * w == 0


def antiderivative(x: float, d: float) -> float:
    return (
        x / (2.0 * d * d * (x * x + d * d))
        + math.atan(x / d) / (2.0 * d**3)
    )


def closed_form(
    x0: float, y: float, m: float, g: float,
    u0: float, v0: float, time: float,
) -> tuple[float, float, float, float]:
    x = x0 - m * g * time
    d = math.hypot(y, m)
    r0 = x0 * x0 + d * d
    r = x * x + d * d
    u = u0 * r0 / r
    v = v0 - 2.0 * y * u0 * r0 * (
        antiderivative(x, d) - antiderivative(x0, d)
    )
    w = -(x * u + y * v) / m
    return x, u, v, w


def elementary_closed_form_checks() -> None:
    # Central differences verify I_d'=R^-4 with a wide safety margin.
    for d in (0.3, 1.0, math.sqrt(2.0), 5.0):
        for x in (-8.0, -1.0, 0.0, 0.7, 11.0):
            step = 1.0e-6 * max(1.0, abs(x), d)
            derivative = (
                antiderivative(x + step, d)
                - antiderivative(x - step, d)
            ) / (2.0 * step)
            target = 1.0 / (x * x + d * d) ** 2
            assert abs(derivative - target) <= 2.0e-7 * max(1.0, target)

    # Direct RK4 comparison for the three independent Kelvin components.
    cases = (
        (-3.0, 1.0, 1.0, 0.7, 0.4, -0.2, 4.0),
        (2.0, -1.3, 0.8, 1.7, -0.6, 0.9, 1.2),
        (-7.0, 2.0, 3.0, 0.4, 0.2, 0.1, 8.0),
    )
    for x0, y, m, g, u0, v0, terminal in cases:
        state = [x0, u0, v0, -(x0 * u0 + y * v0) / m]
        steps = 20000
        dt = terminal / steps

        def rhs(z: list[float]) -> list[float]:
            x, u, v, _w = z
            r2 = x * x + y * y + m * m
            scalar = m * g * u / r2
            return [
                -m * g,
                2.0 * x * scalar,
                2.0 * y * scalar,
                -g * u + 2.0 * m * scalar,
            ]

        for _ in range(steps):
            k1 = rhs(state)
            k2 = rhs([state[i] + dt * k1[i] / 2 for i in range(4)])
            k3 = rhs([state[i] + dt * k2[i] / 2 for i in range(4)])
            k4 = rhs([state[i] + dt * k3[i] for i in range(4)])
            state = [
                state[i] + dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6
                for i in range(4)
            ]
        target = closed_form(x0, y, m, g, u0, v0, terminal)
        assert max(abs(left - right) for left, right in zip(state, target)) < 2e-9


def fixed_cone_grid() -> None:
    # Check the deliberately crude K_C bound over an exhaustive rational
    # input-cone grid and a much wider terminal-x grid. Scaling sets m=1;
    # g/time enter only via final x.
    for cone in (1, 2, 4):
        kc = 4 * (1 + cone) * (1 + 2 * cone * cone) * (1 + math.pi * cone)
        input_grid = [F(j, 4) for j in range(-4 * cone, 4 * cone + 1)]
        terminal_grid = [F(j * cone, 2) for j in range(-64, 65)]
        for x0 in input_grid:
            for x1 in terminal_grid:
                for y in input_grid:
                    d = math.hypot(float(y), 1.0)
                    r0 = float(x0 * x0 + y * y + 1)
                    r1 = float(x1 * x1 + y * y + 1)
                    au = r0 / r1
                    c = -2 * float(y) * r0 * (
                        antiderivative(float(x1), d)
                        - antiderivative(float(x0), d)
                    )
                    # Columns for (u0,v0), including reconstructed w1.
                    col_u = (au, c, -(float(x1) * au + float(y) * c))
                    col_v = (0.0, 1.0, -float(y))
                    # Frobenius norm bounds the operator from horizontal
                    # input; |a0| >= |(u0,v0)|.
                    bound = math.sqrt(sum(z * z for z in col_u + col_v))
                    assert bound <= kc


def long_sweep_gain() -> None:
    ratios = []
    for x in (32.0, 64.0, 128.0, 256.0, 512.0):
        u0 = 1.0 / math.sqrt(1.0 + x * x)
        # Choose g=-1 and T=2X so x0=-X is swept to x1=X.
        _xf, u, v, w = closed_form(-x, 1.0, 1.0, -1.0, u0, 0.0, 2 * x)
        gain = math.sqrt(u * u + v * v + w * w)
        ratios.append(gain / x)
    assert all(left < right for left, right in zip(ratios, ratios[1:]))
    assert abs(ratios[-1] - math.pi / 2) < 0.003
    assert ratios[-1] > 1.56


def schedule_ledger() -> None:
    # q=n^8, mu=nu/((n-1)!)^2. For fixed nu and polynomial q,T,
    # mu*q^2*T is superpolynomially small. Check representative ratios.
    previous = None
    factorial = 1
    for n in range(2, 21):
        if n > 2:
            factorial *= n - 1
        q = n**8
        heat = q * q * n**3 / (factorial * factorial)
        if n >= 12 and previous is not None:
            assert heat < previous
        previous = heat
    assert previous is not None and previous < 1e-8


def main() -> None:
    exact_kelvin_reduction()
    elementary_closed_form_checks()
    fixed_cone_grid()
    long_sweep_gain()
    schedule_ledger()
    print("C181 static vertical-shear Kelvin checks passed")


if __name__ == "__main__":
    main()
