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


DAY = 86400.0


@pytest.mark.parametrize("flux", [-4e-6, 0.0, 1e-7])
def test_richards_column_relaxes_to_the_steady_gardner_profile(flux):
    from earthsciencekit.hydrology import richards_steady_head, simulate_richards_column

    result = simulate_richards_column(2.0, flux, 1e-5, 2.0, 0.4, 0.05, np.linspace(0, 400 * DAY, 5), nodes=80)
    height = (np.arange(80) + 0.5) * 2.0 / 80
    np.testing.assert_allclose(result.values[-1, 2:], richards_steady_head(height, flux, 1e-5, 2.0), atol=1e-6)
    np.testing.assert_allclose(result.values[-1, 1], flux, atol=1e-6 * 1e-5)
    assert result.names[:3] == ("storage", "water_table_flux", "pressure_head_0") and result.units[2] == "m"


def test_richards_column_conserves_water_and_wets_from_the_top():
    from scipy.integrate import trapezoid

    from earthsciencekit.hydrology import exponential_water_content, simulate_richards_column

    time = np.linspace(0, 20 * DAY, 2001)
    result = simulate_richards_column(3.0, -5e-6, 1e-5, 3.0, 0.42, 0.04, time, nodes=60)
    storage, bottom = result.values[:, 0], result.values[:, 1]
    # Storage gains what enters at the surface minus what drains to the water table.
    np.testing.assert_allclose(storage - storage[0], trapezoid_cumulative(bottom, time) + 5e-6 * time, atol=1e-5)
    head = result.values[:, 2:]
    assert head[-1, -1] > head[0, -1] and head[200, -1] > head[200, 30]
    assert trapezoid(bottom, time) < 0
    initial = exponential_water_content(-(np.arange(60) + 0.5) * 0.05, 0.42, 0.04, 3.0)
    np.testing.assert_allclose(storage[0], np.sum(initial - 0.04) * 0.05, rtol=1e-12)


def trapezoid_cumulative(values, time):
    return np.concatenate([[0.0], np.cumsum((values[1:] + values[:-1]) / 2 * np.diff(time))])


@pytest.mark.parametrize(
    "args,kwargs",
    [
        ((2.0, -2e-5, 1e-5, 2.0, 0.4, 0.05), {}),
        ((2.0, 0.0, 1e-5, 2.0, 0.04, 0.05), {}),
        (([2.0, 1.0], 0.0, 1e-5, 2.0, 0.4, 0.05), {}),
        ((2.0, 0.0, 1e-5, 2.0, 0.4, 0.05), {"nodes": 2}),
        ((2.0, 0.0, 1e-5, 2.0, 0.4, 0.05), {"nodes": 5.0}),
        ((2.0, 0.0, 1e-5, 2.0, 0.4, 0.05), {"initial_head": [-1.0, -1.0]}),
        ((2.0, 0.0, 1e-5, 2.0, 0.4, 0.05), {"nodes": 3, "initial_head": [-1.0, 0.5, -1.0]}),
    ],
)
def test_richards_column_rejects_bad_setup(args, kwargs):
    from earthsciencekit.hydrology import simulate_richards_column

    with pytest.raises(ValueError):
        simulate_richards_column(*args, [0, 1], **kwargs)


def test_kinematic_wave_is_steady_conserves_water_and_travels_at_its_celerity():
    from earthsciencekit.hydrology import kinematic_wave_celerity, simulate_kinematic_wave

    hour = 3600.0
    time = np.linspace(0, 48 * hour, 481)
    steady = simulate_kinematic_wave(np.full(time.size, 20.0), 1e4, 30, 0.03, 1e-3, time, nodes=20)
    np.testing.assert_allclose(steady.values[:, 2:], 20.0, rtol=1e-8)
    rating = np.sqrt(1e-3) / (0.03 * 30 ** (2 / 3))
    np.testing.assert_allclose(steady.values[:, 1], 1e4 * (20 / rating) ** 0.6, rtol=1e-8)
    inflow = 10 + 90 * np.exp(-(((time - 8 * hour) / (3 * hour)) ** 2))
    flood = simulate_kinematic_wave(inflow, 5e4, 50, 0.035, 5e-4, time, nodes=200)
    outflow, storage = flood.values[:, 0], flood.values[:, 1]
    gained = np.concatenate([[0.0], np.cumsum(((inflow - outflow)[1:] + (inflow - outflow)[:-1]) / 2 * np.diff(time))])
    np.testing.assert_allclose(storage - storage[0], gained, atol=1e-3 * storage[0])
    velocity = 100 / (100 / (np.sqrt(5e-4) / (0.035 * 50 ** (2 / 3)))) ** 0.6
    delay = time[np.argmax(outflow)] - 8 * hour
    assert 0.9 < delay / (5e4 / kinematic_wave_celerity(velocity, 5 / 3)) < 1.25
    assert 85 < outflow.max() < 100


@pytest.mark.parametrize(
    "inflow,kwargs",
    [(np.zeros(3), {}), (np.ones(2), {}), (np.ones(3), {"nodes": 1}), (np.ones(3), {"nodes": True})],
)
def test_kinematic_wave_rejects_bad_setup(inflow, kwargs):
    from earthsciencekit.hydrology import simulate_kinematic_wave

    with pytest.raises(ValueError):
        simulate_kinematic_wave(inflow, 1e3, 10, 0.03, 1e-3, [0, 1, 2], **kwargs)
    with pytest.raises(ValueError):
        simulate_kinematic_wave(np.ones(3), [1e3, 2e3], 10, 0.03, 1e-3, [0, 1, 2])


def test_snowpack_conserves_water_and_melts_by_degree_days():
    from earthsciencekit.hydrology import degree_day_melt, snowpack_water_equivalent

    rng = np.random.default_rng(5)
    days = np.arange(240)
    temperature = 273.15 - 8 * np.cos(2 * np.pi * (days - 20) / 365) + rng.normal(0, 3, days.size) + 2
    precipitation = rng.exponential(0.003, days.size)
    pack = snowpack_water_equivalent(temperature, precipitation, 3e-3 / DAY, 273.15, DAY)
    swe, released = pack.T
    np.testing.assert_allclose(swe[-1] + released.sum(), precipitation.sum(), rtol=1e-12)
    assert np.all(swe >= 0) and swe.max() > 0.05
    melting = (temperature > 273.15) & (np.r_[0, swe[:-1]] > 0.02)
    np.testing.assert_allclose(released[melting], degree_day_melt(temperature[melting], 3e-3 / DAY, 273.15) * DAY + precipitation[melting])
