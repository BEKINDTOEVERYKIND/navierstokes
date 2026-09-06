#!/usr/bin/env python3
"""C198 exact rational scaling, telescoping, and residual arithmetic.

The analytic energy and circulation PDE arguments are in the claim note;
this checker does not represent numerical samples as PDE certification.
"""
from fractions import Fraction as F


def main():
    # Rational directed bounds for all parameters in the registered interval.
    q = F(6, 5)
    gamma_lo, gamma_hi = F(7, 6), F(37, 25)
    assert 2 * gamma_hi - 3 == -F(1, 25)
    assert gamma_lo - 1 == F(1, 6)
    energy_ratio_cap, circulation_cap = F(993, 1000), F(971, 1000)
    assert energy_ratio_cap ** 25 > 1 / q
    assert circulation_cap ** 6 > 1 / q
    assert circulation_cap ** 3 == F(915498611, 10**9)
    residual = F(1, 100000)
    energy_loss_floor = 1 - energy_ratio_cap * (1 + residual)**2
    assert energy_loss_floor == F(69801399007, 10**13)
    assert energy_loss_floor > F(69, 10000)
    assert 1 - energy_ratio_cap == F(7, 1000)
    assert 1 - circulation_cap == F(29, 1000)

    # Local dilation checked on two nonconstant polynomial fields. These
    # polynomials test chart identities only, not finite-energy existence.
    # Swirl u_theta=r(2+r^2+3z^2); divergence-free poloidal streamfunction
    # psi=r^4 z gives u_r=-r^3, u_z=4r^2 z and eta=-8z.
    for radius, height, gain in ((F(2), F(3), F(5, 4)),
                                 (F(1, 3), F(-2, 5), F(7, 6))):
        old_radius, old_height = radius/q, height/q
        old_swirl = old_radius * (2 + old_radius**2 + 3*old_height**2)
        circulation_after_velocity_chart = radius * old_swirl / gain
        circulation_before = old_radius * old_swirl
        assert circulation_after_velocity_chart == q/gain * circulation_before
        # Curl computed after chart: partial_z u_r - partial_r u_z.
        omega_after_chart = -8 * radius * height / (gain*q**3)
        eta_after_chart = omega_after_chart / radius
        assert eta_after_chart == (-8*old_height) / (gain*q**2)

    # Telescoping checked for several nonuniform-energy exact orbits.
    for k in (F(101, 100), F(11, 10), F(3, 2)):
        for initial in (F(1), F(7, 3)):
            energy, weighted_loss = initial, F(0)
            for j in range(7):
                loss = energy / (j + 4)
                weighted_loss += k**(-j) * loss
                energy = k * (energy - loss)
                assert weighted_loss == initial - k**(-(j + 1)) * energy
        # Zero dissipation explicitly permits finite energy growth.
        assert k**(-3) * (k**3 * initial) == initial

    # Three-stage imperfect recurrence, including varying measured energies.
    k = F(9, 8)
    energies = [F(2)]
    losses = [F(1, 10), F(1, 11), F(1, 12)]
    defects = [F(1, 1000), -F(1, 2000), F(1, 3000)]
    for loss, defect in zip(losses, defects):
        energies.append(k * (energies[-1] - loss) + defect)
    discrepancy = (sum(k**(-j) * losses[j] for j in range(3))
                   - energies[0] + k**(-3) * energies[3])
    assert discrepancy == sum(k**(-(j + 1)) * defects[j] for j in range(3))
    assert abs(discrepancy) <= sum(k**(-(j + 1)) * abs(defects[j])
                                  for j in range(3))

    # Reverse-triangle interval in (2.1), including e > output amplitude.
    # Choose rational exact amplitudes to avoid uncertified square roots.
    for output_amplitude in (F(0), F(1, 10), F(1), F(7, 3)):
        for e in (F(0), F(1, 100), F(1, 2), F(3)):
            lower = max(F(0), output_amplitude - e)
            upper = output_amplitude + e
            for predicted in (lower, (lower + upper) / 2, upper):
                entrance_energy = F(10) + predicted**2 / k
                loss = entrance_energy - predicted**2 / k
                assert entrance_energy - upper**2 / k <= loss
                assert loss <= entrance_energy - lower**2 / k
                defect = output_amplitude**2 - predicted**2
                assert abs(defect) <= 2 * e * output_amplitude + e**2

    # Exact and perturbed circulation iterations, with nonconstant residuals.
    rho, amplitude = F(19, 20), F(7, 3)
    start = amplitude
    errors = (F(1, 1000), F(1, 700), F(1, 800))
    for error in errors:
        amplitude = rho * amplitude + error
    assert amplitude == rho**3 * start + sum(rho**(2-j) * errors[j]
                                             for j in range(3))
    assert amplitude <= circulation_cap**3 * start + sum(
        circulation_cap**(2-j) * errors[j] for j in range(3))

    # Full-trajectory weighted-Sobolev energy bound: the scalar Fourier
    # inequality is a square, including zero and the equality frequency.
    for sqrt_mu in (F(1, 100), F(1, 3), F(2)):
        for frequency in (F(0), F(1, 3), F(1), 1/sqrt_mu, F(100)):
            mu = sqrt_mu**2
            difference = frequency + mu*frequency**3 - 2*sqrt_mu*frequency**2
            assert difference == frequency*(sqrt_mu*frequency-1)**2 >= 0
            assert 2*mu*frequency**2 <= sqrt_mu*(frequency+mu*frequency**3)

    # Bounded-energy contradiction factor and positive limit dissipation.
    for k in (F(101, 100), F(11, 10), F(3, 2)):
        assert k * (1 - (1 - 1/k) / 2) == (k + 1) / 2 > 1
        target_energy = F(7, 3)
        required_loss = (1 - 1/k) * target_energy
        assert k * (target_energy - required_loss) == target_energy
        assert required_loss > 0

    print('C198 full-wake energy/circulation arithmetic: PASS')
    print('EQUAL-L2 COROLLARY: relative residual <= 1/100000 requires loss > 69/10000')
    print('THREE EXACT RETURNS: circulation factor <= 915498611/1000000000')
    print('BOUNDARY: no nonzero regular finite-energy Euler fixed face; singular graphs remain open')


if __name__ == '__main__':
    main()
