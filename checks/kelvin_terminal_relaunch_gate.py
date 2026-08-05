#!/usr/bin/env python3
"""Exact arithmetic ledger for the C96 long-window conditioning gate."""

from fractions import Fraction as F


def determinant_2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def trace_2(matrix):
    return matrix[0][0] + matrix[1][1]


def diagonal_congruence(left, matrix):
    return [
        [left[row] * matrix[row][column] * left[column]
         for column in range(2)]
        for row in range(2)
    ]


def main():
    # C89 rates S=diag(-alpha,-beta,alpha+beta).
    alpha = F(1)
    beta = F(5, 4)
    assert -alpha - beta + (alpha + beta) == 0

    wave_number_log_rate = alpha
    kelvin_log_rates = (beta, -(alpha + beta))
    assert sum(kelvin_log_rates) == -wave_number_log_rate

    phase_2_forward_rate = beta - alpha
    phase_3_forward_rate = -(2 * alpha + beta)
    covariance_22_forward_rate = 2 * beta
    covariance_33_forward_rate = -2 * (alpha + beta)
    covariance_23_forward_rate = -alpha
    kernel_12_forward_rate = 3 * beta - alpha
    kernel_13_forward_rate = -(4 * alpha + 3 * beta)

    # Product of the five physical chart columns; pressure has rate zero.
    jacobian_log_rate = sum((
        covariance_22_forward_rate,
        covariance_33_forward_rate,
        covariance_23_forward_rate,
        kernel_12_forward_rate,
        kernel_13_forward_rate,
    ))
    assert jacobian_log_rate == -8 * alpha

    assert -phase_3_forward_rate == F(13, 4)
    assert -covariance_33_forward_rate == F(9, 2)
    assert -kernel_13_forward_rate == F(31, 4)

    # Make every exponential an exact integer power.  Put u=exp(G/4)=2.
    # Then exp(alpha G)=u^4, exp(beta G)=u^5, and
    # exp((alpha+beta)G)=u^9.
    u = F(2)
    wave_number_gain = u**4
    kelvin = (u**5, u**-9)
    assert kelvin[0] * kelvin[1] == 1 / wave_number_gain

    # A unit terminal transverse covariance pulls back exactly.  Its launch
    # determinant grows by s^2, and its trace obeys the invariant lower bound
    # trace(Q0) >= 2*s*lambda_min(QT).
    terminal_covariance = [[F(1), F(0)], [F(0), F(1)]]
    inverse_kelvin = (1 / kelvin[0], 1 / kelvin[1])
    launch_covariance = diagonal_congruence(
        inverse_kelvin, terminal_covariance
    )
    assert launch_covariance == [
        [F(1, 1024), F(0)],
        [F(0), F(262144)],
    ]
    assert determinant_2(launch_covariance) == wave_number_gain**2
    assert determinant_2(terminal_covariance) == (
        determinant_2(launch_covariance) / wave_number_gain**2
    )
    assert trace_2(launch_covariance) >= 2 * wave_number_gain

    # Labelled viscous harmonics multiply the inviscid chart determinant by
    # (delta_1*delta_2)^2*c_D, delta_j=d_j^2.  A common floor d_j>=d_*
    # gives the crude but uniform factor d_*^10.
    damping_1 = F(1, 2)
    damping_2 = F(3, 4)
    damping_floor = min(damping_1, damping_2)
    delta_1 = damping_1**2
    delta_2 = damping_2**2
    sqrt_a = F(2)
    sqrt_b = F(3)
    c_d = (delta_1 * sqrt_a + delta_2 * sqrt_b) / (
        sqrt_a + sqrt_b
    )
    viscous_factor = (delta_1 * delta_2) ** 2 * c_d
    assert c_d >= damping_floor**2
    assert viscous_factor >= damping_floor**10

    # Exact current-to-launch coordinate pullback.
    # theta_3(T)=u^-13*y(0), B_33(T)=u^-18*B_33(0), and
    # C_23(T)=u^-4*C_23(0).
    theta_3_terminal = F(1, 10)
    theta_3_launch = u**13 * theta_3_terminal
    assert u**-13 * theta_3_launch == theta_3_terminal

    b_33_terminal = F(3, 7)
    b_33_launch = u**18 * b_33_terminal
    assert u**-18 * b_33_launch == b_33_terminal

    c_23_terminal = F(2, 9)
    c_23_launch = u**4 * c_23_terminal
    assert u**-4 * c_23_launch == c_23_terminal

    # The fixed-launch-base Jacobian factor is exp(-8G)=u^-32.
    assert u**(-32) == F(1, 2**32)

    # A constant terminal flow-feedback pulse of duration delta with
    # exp(delta/4)=u has exact phase kernels.  Unit delta-S_13 produces
    # -(1-exp(-(2alpha+beta)delta))/(2alpha+beta), while unit delta-S_12
    # produces -(exp((beta-alpha)delta)-1)/(beta-alpha).
    feedback_3 = -F(4, 13) * (1 - u**-13)
    feedback_2 = -F(4) * (u - 1)
    assert feedback_3 != 0
    assert feedback_2 == -4

    # If the same pulse is separated from the endpoint by delta, its hard
    # phase-3 effect acquires the exact suppression u^-13.
    separated_feedback_3 = u**-13 * feedback_3
    assert abs(separated_feedback_3) < abs(feedback_3)

    # Endpoint order doubling: amplitudes s^(M+1) produce covariance order
    # 2M+2.
    for order in range(1, 20):
        amplitude_power = order + 1
        covariance_power = 2 * amplitude_power
        assert covariance_power == 2 * order + 2

    print("Kelvin transverse area identity: PASS")
    print(f"passive wave-number gain: {wave_number_gain}")
    print(f"launch covariance trace: {trace_2(launch_covariance)}")
    print(f"viscous determinant floor factor: {viscous_factor}")
    print(f"diagonal chart log-determinant rate: {jacobian_log_rate}")
    print(f"hard phase inverse rate: {-phase_3_forward_rate}")
    print(f"hard covariance inverse rate: {-covariance_33_forward_rate}")
    print(f"bounded terminal phase-feedback kernel: {feedback_3}")
    print("all Kelvin terminal-relaunch checks passed")


if __name__ == "__main__":
    main()
