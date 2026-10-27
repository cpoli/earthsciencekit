"""Physical identities, limits and dynamic checks for earthsciencekit.oceanography."""

import numpy as np
import pytest
from scipy.integrate import trapezoid

from earthsciencekit.oceanography.systems.circulation import (
    ekman_spiral_velocity,
    ekman_transport,
    simulate_thermohaline_boxes,
    stommel_streamfunction,
    sverdrup_transport,
    thermohaline_forcing,
)
from earthsciencekit.oceanography.systems.hydrography import abyssal_temperature, linear_density, simplified_seawater_density
from earthsciencekit.oceanography.systems.inertial import simulate_inertial_current
from earthsciencekit.oceanography.systems.waves import gerstner_orbit_radius


def test_inertial_current_direction_and_damping():
    f, damping = 0.1, 0.01
    time = np.linspace(0, 60, 101)
    result = simulate_inertial_current(f, time, damping=damping)
    np.testing.assert_allclose(result.values[:, 0], np.exp(-damping * time) * np.cos(f * time), atol=3e-8)
    np.testing.assert_allclose(result.values[:, 1], -np.exp(-damping * time) * np.sin(f * time), atol=3e-8)
    undamped = simulate_inertial_current(f, time)
    np.testing.assert_allclose(np.sum(undamped.values**2, axis=1), 1, atol=3e-8)
    southern = simulate_inertial_current(-f, time)
    assert southern.values[1, 1] > 0 > undamped.values[1, 1]
    for coriolis in [0, [f, f]]:
        with pytest.raises(ValueError):
            simulate_inertial_current(coriolis, time)


def test_linear_density_rejects_invalid_reference_linearization():
    with pytest.raises(ValueError):
        linear_density(1000, 1e5, 0, 0.001, 0.001)


def test_ekman_spiral_integrates_to_ekman_transport():
    depth = np.linspace(0, 1000, 100001)
    for coriolis in [1e-4, -1e-4]:
        velocity = ekman_spiral_velocity(depth, 0.1, 1025, 0.05, coriolis)
        assert np.isclose(trapezoid(velocity[:, 0], depth), 0, atol=1e-6)
        assert np.isclose(trapezoid(velocity[:, 1], depth), ekman_transport(0.1, 1025, coriolis), rtol=1e-6)
        assert np.isclose(abs(np.arctan2(velocity[0, 1], velocity[0, 0])), np.pi / 4)


def test_stommel_interior_is_sverdrup_and_beta_free_gyre_is_symmetric():
    width, length, tau, rho, beta, friction = 6e6, 5e6, 0.1, 1025, 2e-11, 1e-7
    x = np.array([0.3, 0.5, 0.8]) * width
    psi = stommel_streamfunction(x, length / 2, tau, width, length, rho, beta, friction)
    curl = -tau * np.pi / length
    np.testing.assert_allclose(psi, -sverdrup_transport(curl, rho, beta) * (width - x), rtol=0.02)
    symmetric = stommel_streamfunction([0.2 * width, 0.8 * width], length / 3, tau, width, length, rho, 0, friction)
    assert np.isclose(symmetric[0], symmetric[1])
    edges = stommel_streamfunction([0, width, width / 2], [length / 2, length / 2, 0], tau, width, length, rho, beta, friction)
    np.testing.assert_allclose(edges, 0, atol=1e-6 * psi[0])


def test_thermohaline_boxes_reach_both_steady_states():
    thermal, haline, forcing, volume = 3e7, 6e6, 3e7, 1e16
    time = np.linspace(0, 3e10, 61)
    for initial, sign in [(0.0, 1), (8.0, -1)]:
        result = simulate_thermohaline_boxes(thermal, haline, forcing, volume, time, initial_salinity_difference=initial)
        flow = result.values[-1, 1]
        assert np.sign(flow) == sign
        assert np.isclose(thermohaline_forcing(flow, thermal, haline), forcing, rtol=1e-6)
        np.testing.assert_allclose(result.values[:, 1], thermal - haline * result.values[:, 0])
    flows = np.linspace(0, thermal, 1001)
    assert np.isclose(thermohaline_forcing(flows, thermal, haline).max(), thermal**2 / (4 * haline))
    for bad in [[thermal, thermal], 0]:
        with pytest.raises(ValueError):
            simulate_thermohaline_boxes(bad, haline, forcing, volume, time)


def test_abyssal_profile_solves_advection_diffusion():
    z = np.linspace(0, 4000, 4001)
    temperature = abyssal_temperature(z, 274, 284, 4000, 1e-4, 1e-7)
    gradient, curvature = np.gradient(temperature, z), np.gradient(np.gradient(temperature, z), z)
    np.testing.assert_allclose(1e-7 * gradient[5:-5], 1e-4 * curvature[5:-5], rtol=1e-4)
    linear = abyssal_temperature(z, 274, 284, 4000, 1e-4, 1e-12)
    np.testing.assert_allclose(linear, 274 + 10 * z / 4000, atol=1e-3)


def test_simplified_density_cabbeling_and_linear_limit():
    from scipy.optimize import brentq

    cold = (275.15, 34.6)
    rho_cold = simplified_seawater_density(*cold, 0)
    warm = (293.15, brentq(lambda s: simplified_seawater_density(293.15, s, 0) - rho_cold, 30, 45, xtol=1e-14))
    mixed = simplified_seawater_density((cold[0] + warm[0]) / 2, (cold[1] + warm[1]) / 2, 0)
    assert mixed > rho_cold
    small = np.array([1e-4, -1e-4])
    nonlinear = simplified_seawater_density(283.15 + small, 35 + small[::-1], 0)
    linear = linear_density(1026, small, small[::-1], 0.16550 / 1026, 0.76554 / 1026)
    np.testing.assert_allclose(nonlinear, linear, atol=1e-8)


def test_gerstner_orbits_decay_one_wavelength_down():
    k = 2 * np.pi / 100
    ratio = gerstner_orbit_radius(100, 1.0, k) / gerstner_orbit_radius(0, 1.0, k)
    assert np.isclose(ratio, np.exp(-2 * np.pi))


@pytest.mark.parametrize(
    "absolute_salinity,celsius,decibar,density",
    [
        # GSW-C gsw_check_data.h (TEOS-10 v3.06.16), the computed check values of gsw_rho for the test casts.
        (34.468236430490606, 27.996436412058213, 0.0, 1021.8863044505447),
        (34.67238066266718, 27.30397518400217, 40.0, 1022.4287988266693),
        (34.86377002371946, 1.0291291795169686, 5098.0, 1050.5884824364505),
        (10.389468455026284, 4.611397992798779, 101.0, 1008.6875042532482),
    ],
)
def test_teos10_density_matches_teos10_check_values(absolute_salinity, celsius, decibar, density):
    from earthsciencekit.oceanography import teos10_density

    np.testing.assert_allclose(teos10_density(celsius + 273.15, absolute_salinity, decibar * 1e4), density, rtol=1e-13)


def test_teos10_density_brackets_the_fits_and_is_thermobaric():
    from earthsciencekit.oceanography import teos10_density

    kelvin, salinity = np.meshgrid(np.linspace(272, 300, 15), np.linspace(33, 37, 9))
    exact = teos10_density(kelvin, salinity, 0)
    # Roquet et al. (2015) tune the S-EOS to TEOS-10 gradients; its absolute density stays within a few kg/m^3.
    assert np.max(np.abs(simplified_seawater_density(kelvin, salinity, 0) - exact)) < 3
    cold, warm = (teos10_density(t, 35, [0, 4e7]) for t in (275.15, 285.15))
    # Thermal expansion grows with pressure (thermobaricity), and seawater is compressed by about 4.5 % per 100 MPa.
    assert cold[0] - warm[0] < cold[1] - warm[1]
    assert 0.015 < cold[1] / cold[0] - 1 < 0.02
