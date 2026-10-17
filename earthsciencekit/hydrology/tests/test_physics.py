"""Hydrology identities, limits and the Green-Ampt infiltration simulation."""

import numpy as np
import pytest
from scipy.integrate import quad

from earthsciencekit.atmosphere import saturation_vapor_pressure
from earthsciencekit.constants import CP_DRY_AIR, LATENT_HEAT_VAPORIZATION
from earthsciencekit.hydrology import (
    chezy_from_manning,
    chezy_velocity,
    cooper_jacob_drawdown,
    cooper_jacob_transmissivity,
    darcy_flux,
    dupuit_discharge,
    dupuit_head,
    gardner_conductivity,
    green_ampt_rate,
    green_ampt_time,
    horton_cumulative_infiltration,
    horton_infiltration_rate,
    hydraulic_radius,
    manning_velocity,
    muskingum_coefficients,
    muskingum_routing,
    nash_unit_hydrograph,
    penman_evaporation,
    penman_monteith_evaporation,
    psychrometric_constant,
    richards_steady_head,
    saturation_slope,
    simulate_infiltration,
    stream_number,
    theis_drawdown,
    topmodel_baseflow,
    topmodel_deficit,
    unit_hydrograph_convolution,
    unsaturated_darcy_flux,
    water_balance,
    well_function,
)


def test_water_balance_closes_over_a_cycle():
    rain = np.array([3.0, 2.0, 1.0, 0.5])
    evaporation = np.array([0.5, 1.0, 1.5, 1.0])
    runoff = np.array([1.0, 1.0, 0.5, 0.0])
    storage = np.cumsum(water_balance(rain, evaporation, runoff))
    assert np.isclose(storage[-1], 0) and storage.max() > 0


def test_manning_is_chezy_with_a_sixth_root_coefficient():
    radius = hydraulic_radius(np.array([2.0, 10.0, 1000.0]), 1.0)
    assert np.isclose(radius[-1], 1, rtol=3e-3)
    np.testing.assert_allclose(manning_velocity(0.03, radius, 1e-3), chezy_velocity(chezy_from_manning(0.03, radius), radius, 1e-3))


def test_dupuit_profile_carries_the_dupuit_discharge():
    x = np.linspace(0, 100, 2001)
    head = dupuit_head(x, 10, 4, 100)
    discharge = head * darcy_flux(1e-4, np.gradient(head, x))
    np.testing.assert_allclose(discharge[1:-1], dupuit_discharge(1e-4, 10, 4, 100), rtol=1e-5)


def test_richards_profile_carries_a_constant_flux():
    height = np.linspace(0, 0.5, 4001)
    for flux in (-4e-6, 0.0, 1e-6):
        head = richards_steady_head(height, flux, 1e-5, 4)
        conductivity = gardner_conductivity(head, 1e-5, 4)
        result = unsaturated_darcy_flux(conductivity, np.gradient(head, height))
        np.testing.assert_allclose(result[1:-1], flux, atol=1e-11)
    np.testing.assert_allclose(richards_steady_head(height, 0.0, 1e-5, 4), -height, atol=1e-15)
    deep = richards_steady_head(50.0, -4e-6, 1e-5, 4)
    assert np.isclose(gardner_conductivity(deep, 1e-5, 4), 4e-6)


def test_green_ampt_time_inverts_its_rate():
    depth, step = np.linspace(1e-3, 0.2, 40), 1e-7
    slope = (green_ampt_time(depth + step, 1e-6, 0.1, 0.3) - green_ampt_time(depth - step, 1e-6, 0.1, 0.3)) / (2 * step)
    np.testing.assert_allclose(1 / slope, green_ampt_rate(depth, 1e-6, 0.1, 0.3), rtol=1e-6)
    early = green_ampt_time(1e-6, 1e-6, 0.1, 0.3)
    assert np.isclose(early, 1e-12 / (2 * 1e-6 * 0.03), rtol=1e-4)


def test_unit_hydrograph_conserves_volume_and_superposes():
    area, step = 50e6, 3600.0
    ordinates = area * nash_unit_hydrograph(np.arange(0, 200) * step, 3, 5 * step)
    assert np.isclose(quad(lambda t: nash_unit_hydrograph(t, 3, 1), 0, np.inf)[0], 1)
    storm = np.array([0.01, 0.02, 0.005])
    flood = unit_hydrograph_convolution(storm, ordinates)
    assert np.isclose(flood.sum() * step, storm.sum() * area, rtol=1e-3)
    np.testing.assert_allclose(flood, sum(unit_hydrograph_convolution(np.eye(3)[i] * storm[i], ordinates) for i in range(3)))
    assert nash_unit_hydrograph(0, 1, 10) == 0.1 and nash_unit_hydrograph(0, 2, 10) == 0


def test_horton_cumulative_is_the_integral_of_the_rate():
    total = quad(lambda t: horton_infiltration_rate(t, 3e-5, 3e-6, 2e-3), 0, 7200)[0]
    assert np.isclose(horton_cumulative_infiltration(7200, 3e-5, 3e-6, 2e-3), total, rtol=1e-12)
    assert np.isclose(horton_infiltration_rate(1e6, 3e-5, 3e-6, 2e-3), 3e-6)


def test_theis_tends_to_cooper_jacob_and_its_slope():
    time = np.logspace(4, 7, 50)
    theis = theis_drawdown(0.01, 1e-3, 1e-4, 5, time)
    jacob = cooper_jacob_drawdown(0.01, 1e-3, 1e-4, 5, time)
    u = 25 * 1e-4 / (4 * 1e-3 * time)
    assert np.all(np.abs(theis - jacob)[u < 0.01] < 0.01 * theis[u < 0.01])
    slope = theis[-1] - theis_drawdown(0.01, 1e-3, 1e-4, 5, time[-1] / 10)
    assert np.isclose(cooper_jacob_transmissivity(0.01, slope), 1e-3, rtol=1e-3)
    assert np.isclose(well_function(1e-8), -np.euler_gamma - np.log(1e-8), rtol=1e-8)


def test_muskingum_conserves_volume_and_attenuates():
    coefficients = muskingum_coefficients([3600.0, 7200.0], [0.1, 0.3], 1800.0)
    np.testing.assert_allclose(coefficients.sum(axis=-1), 1)
    time = np.arange(120)
    inflow = 10 + 90 * np.exp(-(((time - 20) / 5.0) ** 2))
    outflow = muskingum_routing(inflow, 7200, 0.2, 3600)
    assert np.isclose(outflow.sum(), inflow.sum(), rtol=1e-6)
    assert outflow.max() < inflow.max() and outflow.argmax() > inflow.argmax()
    np.testing.assert_allclose(muskingum_routing(np.full(5, 7.0), 7200, 0.2, 3600), 7)


def test_stream_numbers_fall_geometrically():
    counts = stream_number(np.arange(1, 6), 5, 4)
    np.testing.assert_allclose(counts[:-1] / counts[1:], 4)
    assert counts[-1] == 1


def test_penman_monteith_reduces_to_penman_and_open_water():
    temperature, pressure, density, resistance, deficit, energy = 293.15, 101325.0, 1.2, 50.0, 1200.0, 300.0
    gamma = psychrometric_constant(pressure)
    aerodynamic = density * CP_DRY_AIR * deficit / (gamma * LATENT_HEAT_VAPORIZATION * resistance)
    open_water = penman_evaporation(energy, temperature, aerodynamic, pressure)
    assert np.isclose(penman_monteith_evaporation(energy, temperature, deficit, density, resistance, 0, pressure), open_water)
    canopy = penman_monteith_evaporation(energy, temperature, deficit, density, resistance, [70, 500], pressure)
    assert np.all(np.diff(np.r_[open_water, canopy]) < 0)
    step = 1e-3
    numeric = (saturation_vapor_pressure(temperature + step) - saturation_vapor_pressure(temperature - step)) / (2 * step)
    assert np.isclose(saturation_slope(temperature), numeric, rtol=1e-6)
    assert np.isclose(penman_evaporation(energy, temperature, 0, pressure) * LATENT_HEAT_VAPORIZATION, energy * 0.7, rtol=0.05)


def test_topmodel_saturates_the_wettest_ground_first():
    index = np.array([4.0, 6.0, 8.0, 10.0])
    deficit = topmodel_deficit(0.04, index, index.mean(), 0.02)
    assert np.isclose(deficit.mean(), 0.04) and deficit[-1] < 0 < deficit[0]
    ratio = topmodel_baseflow(1e-3, 7, 0.06, 0.02) / topmodel_baseflow(1e-3, 7, 0.04, 0.02)
    assert np.isclose(ratio, np.exp(-1))


def test_infiltration_matches_mein_larson_after_ponding():
    conductivity, suction, deficit, rain = 1e-6, 0.1, 0.3, 5e-6
    time = np.linspace(600, 600 + 6 * 3600, 721)
    result = simulate_infiltration(conductivity, suction, deficit, rain, time)
    depth, rate, runoff = result.values.T
    storage = suction * deficit
    ponding = conductivity * storage / (rain * (rain - conductivity))
    elapsed = time - time[0]
    before = elapsed <= ponding
    np.testing.assert_allclose(depth[before], rain * elapsed[before], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(runoff[before], 0)
    after = elapsed > ponding + 60
    ponded_depth = rain * ponding
    exact = ponding + (depth[after] - ponded_depth - storage * np.log((storage + depth[after]) / (storage + ponded_depth))) / conductivity
    np.testing.assert_allclose(exact, elapsed[after], rtol=1e-5)
    np.testing.assert_allclose(rate[after], green_ampt_rate(depth[after], conductivity, suction, deficit))
    np.testing.assert_allclose(rate + runoff, rain)
    assert result.names == ("cumulative_infiltration", "infiltration_rate", "runoff_rate") and result.time[0] == 600


def test_light_rain_never_ponds():
    result = simulate_infiltration(1e-6, 0.1, 0.3, 5e-7, np.linspace(0, 86400, 50))
    np.testing.assert_allclose(result.values[:, 0], 5e-7 * result.time, rtol=1e-8)
    np.testing.assert_allclose(result.values[:, 2], 0)
    with pytest.raises(ValueError):
        simulate_infiltration(1e-6, [0.1, 0.2], 0.3, 5e-6, [0, 1])
