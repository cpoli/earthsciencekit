"""Geomorphology identities, limits and the river-profile simulation."""

import numpy as np
import pytest
from scipy.integrate import quad

from earthsciencekit.constants import LAMBDA_BE10, SECONDS_PER_YEAR, VON_KARMAN
from earthsciencekit.geomorphology import (
    bagnold_sand_flux,
    bed_shear_stress,
    bedload_flux,
    chi_coordinate,
    cosmogenic_concentration,
    cosmogenic_production_at_depth,
    critical_shields_number,
    diffusive_hillslope_flux,
    dimensionless_grain_size,
    downstream_hydraulic_geometry,
    dune_migration_rate,
    effective_discharge,
    ferguson_church_settling_velocity,
    flint_slope,
    hack_length,
    hilltop_curvature,
    kawamura_sand_flux,
    keulegan_mean_velocity,
    knickpoint_celerity,
    meyer_peter_muller_transport,
    morphologic_age,
    nonlinear_hillslope_flux,
    rouse_concentration,
    rouse_number,
    scarp_elevation,
    scarp_maximum_gradient,
    shear_velocity,
    shields_number,
    simulate_river_profile,
    soil_production_rate,
    steady_state_erosion_rate,
    steady_state_slope,
    strahler_order,
    stream_power_erosion,
    threshold_shear_velocity,
)

YEAR = SECONDS_PER_YEAR
RHO_S, RHO_W, NU, G = 2650.0, 1000.0, 1e-6, 9.81


def test_ferguson_church_tends_to_stokes_for_fine_grains_and_to_drag_for_coarse():
    from earthsciencekit.geomorphology import stokes_settling_velocity

    fine = 1e-6
    stokes = stokes_settling_velocity(fine, RHO_S, RHO_W, RHO_W * NU, G)
    assert np.isclose(ferguson_church_settling_velocity(fine, RHO_S, RHO_W, NU, G), stokes, rtol=1e-3)
    coarse = np.array([0.1, 0.4])
    velocity = ferguson_church_settling_velocity(coarse, RHO_S, RHO_W, NU, G)
    assert np.isclose(velocity[1] / velocity[0], 2.0, rtol=1e-2)
    assert np.isclose(velocity[1], np.sqrt(1.65 * G * 0.4 / 0.75), rtol=1e-2)


def test_shields_number_of_a_uniform_flow_and_the_threshold_curve():
    stress = bed_shear_stress(1.5, 2e-3, RHO_W, G)
    assert np.isclose(shear_velocity(stress, RHO_W), np.sqrt(G * 1.5 * 2e-3))
    assert np.isclose(shields_number(stress, 0.02, RHO_S, RHO_W, G), 1.5 * 2e-3 / (1.65 * 0.02))
    size = dimensionless_grain_size(np.geomspace(1e-5, 0.1, 200), RHO_S, RHO_W, NU, G)
    threshold = critical_shields_number(size)
    assert size[0] < 1 < 1e3 < size[-1]
    lowest = np.argmin(threshold)
    assert 5 < size[lowest] < 50 and threshold[0] > threshold[lowest] < threshold[-1]
    assert np.isclose(critical_shields_number(1e6), 0.055, rtol=1e-5)


def test_keulegan_velocity_is_the_depth_average_of_the_log_profile():
    shear, depth, roughness = 0.08, 1.2, 0.01
    bottom = roughness / 30
    mean = quad(lambda z: shear / VON_KARMAN * np.log(z / bottom), bottom, depth)[0] / depth
    assert np.isclose(keulegan_mean_velocity(shear, depth, roughness), mean, rtol=2e-3)


def test_rouse_profile_balances_settling_and_turbulent_mixing():
    depth, reference, settling, shear = 2.0, 0.1, 0.01, 0.05
    number = rouse_number(settling, shear)
    height = np.linspace(0.2, 1.8, 801)
    concentration = rouse_concentration(height, depth, reference, 3.0, number)
    diffusivity = VON_KARMAN * shear * height * (1 - height / depth)
    flux = settling * concentration + diffusivity * np.gradient(concentration, height)
    assert np.all(np.abs(flux[1:-1]) < 1e-3 * settling * concentration[1:-1])
    assert rouse_concentration(reference, depth, reference, 3.0, number) == 3.0
    assert rouse_concentration(depth, depth, reference, 3.0, number) == 0.0
    np.testing.assert_allclose(rouse_concentration(height, depth, reference, 3.0, 0), 3.0)


def test_meyer_peter_muller_threshold_and_bedload_scaling():
    assert meyer_peter_muller_transport(0.03, 0.047) == 0.0
    einstein = meyer_peter_muller_transport(0.1, 0.047)
    assert np.isclose(einstein, 8 * 0.053**1.5)
    flux = bedload_flux(einstein, np.array([0.01, 0.04]), RHO_S, RHO_W, G)
    assert np.isclose(flux[1] / flux[0], 8.0)


def test_nonlinear_creep_is_linear_on_gentle_slopes_and_diverges_near_critical():
    assert np.isclose(nonlinear_hillslope_flux(0.01, 1e-4, 1.25), diffusive_hillslope_flux(0.01, 1e-4), rtol=1e-8)
    assert nonlinear_hillslope_flux(0.01, -1.24, 1.25) > 50 * diffusive_hillslope_flux(0.01, -1.24)


def test_scarp_solves_the_diffusion_equation_and_its_age_inverts():
    kappa, age, half, slope = 1e-2 / YEAR, 2e4 * YEAR, 3.0, 0.05
    x = np.linspace(-40, 40, 4001)
    dt = 1e-3 * age
    rate = (scarp_elevation(x, half, kappa, age + dt, slope) - scarp_elevation(x, half, kappa, age - dt, slope)) / (2 * dt)
    curvature = np.gradient(np.gradient(scarp_elevation(x, half, kappa, age, slope), x), x)
    np.testing.assert_allclose(rate[2:-2], kappa * curvature[2:-2], atol=1e-6 * half / age * 1e3)
    gradient = np.gradient(scarp_elevation(x, half, kappa, age, slope), x)
    assert np.isclose(gradient.max(), scarp_maximum_gradient(half, kappa, age, slope), rtol=1e-5)
    measured = scarp_maximum_gradient(half, kappa, age, slope)
    assert np.isclose(morphologic_age(half, measured, slope), kappa * age)
    np.testing.assert_allclose(scarp_elevation([-1, 0, 1], half, 0, age, 0), [-half, 0, half])


def test_soil_production_and_hilltop_curvature_balance_steady_lowering():
    assert np.isclose(soil_production_rate(0.5, 2e-12, 0.5), 2e-12 / np.e)
    erosion, kappa, rock, soil = 1e-4 / YEAR, 5e-3 / YEAR, 2600, 1300
    curvature = hilltop_curvature(erosion, kappa, rock, soil)
    x = np.linspace(-20, 20, 81)
    flux = diffusive_hillslope_flux(kappa, curvature * x)
    np.testing.assert_allclose(np.gradient(flux, x), rock * erosion / soil)


def test_stream_power_steady_state_flint_law_and_knickpoint_speed():
    area = np.geomspace(1e5, 1e10, 6)
    uplift, k, m, n = 1e-3 / YEAR, 2e-5 / YEAR, 0.45, 1.5
    slope = steady_state_slope(uplift, area, k, m, n)
    np.testing.assert_allclose(stream_power_erosion(area, slope, k, m, n), uplift)
    np.testing.assert_allclose(flint_slope(area, (uplift / k) ** (1 / n), m / n), slope)
    np.testing.assert_allclose(knickpoint_celerity(area, slope, k, m, n) * slope, uplift)


def test_chi_reduces_to_distance_at_the_reference_area_and_hack_inverts():
    distance = np.linspace(0, 1e4, 11)
    np.testing.assert_allclose(chi_coordinate(distance, np.full(11, 3e6), 3e6, 0.5), distance)
    linear = chi_coordinate([0, 1, 2, 3], [1, 1, 1, 1], 4, 0.5)
    np.testing.assert_allclose(linear, [0, 2, 4, 6])
    area = 2.5e7
    assert np.isclose(hack_length(area, 0.32, 0.6) ** (1 / 0.6) / 0.32 ** (1 / 0.6), area)


def test_strahler_order_of_balanced_trees_and_side_tributaries():
    # A complete binary tree of depth 3: eight sources, the outlet is order 4.
    downstream = [-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    order = strahler_order(downstream)
    assert order[0] == 4 and order[1] == order[2] == 3 and np.all(order[7:] == 1)
    # A first-order tributary joining a second-order stream leaves it second order; two outlets are allowed.
    np.testing.assert_allclose(strahler_order([3, 3, 4, 4, -1, -1]), [1, 1, 1, 2, 2, 1])
    np.testing.assert_allclose(strahler_order([4, 3, 3, 4, -1]), [1, 1, 1, 2, 2])
    with pytest.raises(ValueError, match="cycles"):
        strahler_order([1, 2, 0, -1])


def test_hydraulic_geometry_conserves_discharge_and_effective_discharge_peaks():
    discharge = np.geomspace(1, 1e4, 5)
    geometry = downstream_hydraulic_geometry(discharge, 4.0, 0.5, 0.3, 0.4)
    np.testing.assert_allclose(geometry.prod(axis=-1), discharge)
    np.testing.assert_allclose(geometry[:, 2], discharge**0.1 / 1.2)
    median, sigma, b = 40.0, 0.9, 2.2
    q = np.geomspace(1, 1e4, 200001)
    work = q**b * np.exp(-((np.log(q / median)) ** 2) / (2 * sigma**2)) / q
    assert np.isclose(q[np.argmax(work)], effective_discharge(median, sigma, b), rtol=1e-4)
    assert effective_discharge(median, sigma, 1) == median


def test_blown_sand_threshold_cube_law_and_dune_speed():
    threshold = threshold_shear_velocity(2.5e-4, RHO_S, 1.2, G, 0.1)
    assert 0.2 < threshold < 0.3
    assert kawamura_sand_flux(threshold, threshold, 1.2, G, 2.78) == 0.0
    assert np.isclose(kawamura_sand_flux(0.6, 0, 1.2, G, 1.8), bagnold_sand_flux(0.6, 2.5e-4, 1.2, G, 1.8))
    speeds = dune_migration_rate(0.01, np.array([2.0, 8.0]), 1600)
    assert np.isclose(speeds[0] / speeds[1], 4.0)


def test_cosmogenic_accumulation_saturation_and_erosion_inversion():
    production, density, attenuation = 4.0 / 1e-3, 2650.0, 1600.0  # atoms/(kg year) with years as the time unit
    assert np.isclose(LAMBDA_BE10 * 1.387e6 * YEAR, np.log(2))
    assert np.isclose(cosmogenic_production_at_depth(production, attenuation / density, density, attenuation), production / np.e)
    np.testing.assert_allclose(cosmogenic_concentration([0, 1e3], production, 0, 0, density, attenuation), [0, production * 1e3])
    erosion = 2e-5  # m/year
    saturated = cosmogenic_concentration(1e8, production, 5e-7, erosion, density, attenuation)
    assert np.isclose(saturated, production / (5e-7 + density * erosion / attenuation))
    assert np.isclose(steady_state_erosion_rate(saturated, production, 5e-7, density, attenuation), erosion)
    young = cosmogenic_concentration(10.0, production, 5e-7, erosion, density, attenuation)
    assert np.isclose(young, production * 10.0, rtol=1e-3)


def _channel(nodes=201):
    distance = np.linspace(0, 50e3, nodes)
    area = ((50e3 - distance + 300) / 0.32) ** (1 / 0.6)
    return distance, area


def test_river_profile_reaches_the_discrete_steady_state_close_to_chi():
    distance, area = _channel()
    uplift, k = 1e-3 / YEAR, 5.4e-6 / YEAR
    steady = np.concatenate([[0.0], np.cumsum(steady_state_slope(uplift, area[1:], k, 0.5, 1) * np.diff(distance))])
    run = simulate_river_profile(distance, area, uplift, k, np.linspace(0, 5e6, 26) * YEAR)
    assert run.names[:3] == ("relief", "mean_elevation", "elevation_0") and run.units[0] == "m"
    np.testing.assert_allclose(run.values[-1, 2:], steady, atol=1e-3)
    assert np.all(run.values[:, 2] == 0) and np.all(np.diff(run.values[:, 0]) > -1e-3)
    chi = chi_coordinate(distance, area, 1e6, 0.5)
    np.testing.assert_allclose(run.values[-1, 2:], uplift / (k * 1e3) * chi, rtol=0.06, atol=1.0)


def test_knickpoint_travels_at_the_stream_power_celerity():
    distance, area = _channel()
    uplift, k = 1e-3 / YEAR, 5.4e-6 / YEAR
    steady = np.concatenate([[0.0], np.cumsum(steady_state_slope(uplift, area[1:], k, 0.5, 1) * np.diff(distance))])
    years = np.array([0.0, 2.5e5, 5e5, 1e6])
    run = simulate_river_profile(distance, area, 2 * uplift, k, years * YEAR, initial_elevation=steady)
    celerity = knickpoint_celerity(area[1:], 1.0, k, 0.5, 1)
    travel = np.concatenate([[0.0], np.cumsum(np.diff(distance) / celerity)])
    slope = np.diff(steady) / np.diff(distance)
    for row, year in zip(run.values[1:, 2:], years[1:]):
        exact = np.minimum(2 * steady, steady + uplift * year * YEAR)
        assert np.max(np.abs(row - exact)) < 0.05 * uplift * year * YEAR
        knick = distance[1:][np.flatnonzero(np.diff(row) / np.diff(distance) / slope < 1.5)[0]]
        assert abs(knick - np.interp(year * YEAR, travel, distance)) <= 2 * (distance[1] - distance[0])


def test_relief_decays_without_uplift_and_ever_more_slowly_for_a_steeper_slope_exponent():
    distance, area = _channel(101)
    k = 5.4e-6 / YEAR
    steady = np.concatenate([[0.0], np.cumsum(steady_state_slope(1e-3 / YEAR, area[1:], k, 0.5, 1) * np.diff(distance))])
    time = np.linspace(0, 3e6, 31) * YEAR
    linear = simulate_river_profile(distance, area, 0, k, time, initial_elevation=steady).values[:, 0]
    # The same erosion rate at the outlet's initial slope, so both start losing relief alike.
    k_steep = k / (np.diff(steady)[0] / np.diff(distance)[0])
    steep = simulate_river_profile(distance, area, 0, k_steep, time, slope_exponent=2, initial_elevation=steady).values[:, 0]
    assert np.all(np.diff(linear) < 0) and linear[-1] < 0.05 * linear[0]
    assert np.all(np.diff(steep) < 0) and np.all(np.diff(steep, 2) > 0) and steep[-1] > 3 * linear[-1]


@pytest.mark.parametrize(
    "args,kwargs",
    [
        (([0, 0, 1], [1, 1, 1], 0, 1, [0, 1]), {}),
        (([0, 1], [1, 1, 1], 0, 1, [0, 1]), {}),
        (([0, 1], [1, 1], [0, 0], 1, [0, 1]), {}),
        (([0, 1], [1, 1], 0, 1, [0, 1]), {"initial_elevation": [0, 0, 0]}),
        (([0, 1], [1, 1], 0, 0, [0, 1]), {}),
        (([0, 1], [1, 1], 0, 1, [1, 0]), {}),
    ],
)
def test_river_profile_rejects_bad_setup(args, kwargs):
    with pytest.raises(ValueError):
        simulate_river_profile(*args, **kwargs)
