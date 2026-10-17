"""Glaciology identities, limits and the flowline glacier simulation."""

import numpy as np
import pytest
from scipy.integrate import quad

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.glaciology import (
    dansgaard_delta,
    dansgaard_johnsen_age,
    dansgaard_johnsen_downward_velocity,
    dansgaard_temperature,
    driving_stress,
    effective_viscosity,
    flotation_thickness,
    glen_rate_factor,
    glen_strain_rate,
    glen_stress,
    glen_surface_velocity,
    glen_velocity_profile,
    haskell_viscosity,
    ice_shelf_strain_rate,
    linear_mass_balance,
    nye_age,
    nye_crevasse_depth,
    nye_plastic_profile,
    rebound_uplift_rate,
    remaining_depression,
    schoof_grounding_line_flux,
    sea_ice_freeboard,
    shallow_ice_diffusivity,
    shallow_ice_flux,
    simulate_flowline_glacier,
    stefan_growth_rate,
    stefan_ice_thickness,
    vialov_dome_height,
    vialov_profile,
    viscous_relaxation_time,
    viscous_slab_velocity,
    water_filled_crevasse_depth,
)

A, N, RHO, G = 2.4e-24, 3.0, 917.0, 9.81


def test_glen_with_n_one_is_forbes_viscous_slab():
    height = np.linspace(0, 300, 7)
    viscosity = 1e14
    newtonian = glen_velocity_profile(height, 300, 0.05, 1 / (2 * viscosity), 1, RHO, G)
    np.testing.assert_allclose(newtonian, viscous_slab_velocity(height, 300, 0.05, viscosity, RHO, G), rtol=1e-12)
    assert glen_velocity_profile(300, 300, 0.05, A, N, RHO, G) == glen_surface_velocity(300, 0.05, A, N, RHO, G)


def test_glen_stress_viscosity_and_rate_factor_are_consistent():
    stress = np.array([2e4, 1e5, 3e5])
    rate = glen_strain_rate(stress, A, N)
    np.testing.assert_allclose(glen_stress(rate, A, N), stress, rtol=1e-12)
    np.testing.assert_allclose(2 * effective_viscosity(stress, A, N) * rate, stress, rtol=1e-12)
    assert glen_rate_factor(263.15, A, 6e4, 263.15) == A and glen_rate_factor(240, A, 0, 263.15) == A
    assert glen_rate_factor(243.15, A, 6e4, 263.15) < A / 5


def test_shallow_ice_flux_is_the_integral_of_the_glen_profile():
    thickness, slope = 800.0, 0.01
    integral = quad(lambda z: glen_velocity_profile(z, thickness, slope, A, N, RHO, G), 0, thickness)[0]
    assert np.isclose(shallow_ice_flux(thickness, slope, A, N, RHO, G), integral, rtol=1e-10)
    flux = shallow_ice_diffusivity(thickness, slope, A, N, RHO, G) * slope
    assert np.isclose(flux, shallow_ice_flux(thickness, slope, A, N, RHO, G))
    stress = driving_stress(thickness, slope, RHO, G)
    assert np.isclose(shallow_ice_flux(thickness, slope, A, N, RHO, G), 2 * A * stress**N * thickness**2 / (N + 2))


def test_nye_profile_carries_the_yield_stress():
    distance = np.geomspace(1e3, 5e5, 4001)
    height = nye_plastic_profile(distance, 1e5, RHO, G)
    np.testing.assert_allclose(driving_stress(height, np.gradient(height, distance), RHO, G)[1:-1], 1e5, rtol=1e-4)


def test_vialov_profile_carries_the_accumulated_snowfall():
    half_width, accumulation = 5e5, 0.2 / SECONDS_PER_YEAR
    dome = vialov_dome_height(half_width, accumulation, A, N, RHO, G)
    distance = np.linspace(0, half_width, 20001)
    height = vialov_profile(distance, half_width, dome, N)
    slope = -np.gradient(height, distance)
    flux = shallow_ice_flux(height, slope, A, N, RHO, G)
    interior = (distance > 0.05 * half_width) & (distance < 0.9 * half_width)
    np.testing.assert_allclose(flux[interior], accumulation * distance[interior], rtol=1e-3)
    assert height[0] == dome and height[-1] == 0


def test_crevasses_deepen_with_stretching_and_water():
    stress = glen_stress(1e-10, A, N)
    assert np.isclose(nye_crevasse_depth(1e-10, A, N, RHO, G), 2 * stress / (RHO * G))
    assert water_filled_crevasse_depth(1e-10, A, N, 0, RHO, 1000, G) == nye_crevasse_depth(1e-10, A, N, RHO, G)
    assert water_filled_crevasse_depth(1e-10, A, N, 20, RHO, 1000, G) > nye_crevasse_depth(1e-10, A, N, RHO, G) + 20


def test_stefan_thickness_grows_at_the_stefan_rate():
    difference, conductivity, latent = 25.0, 2.2, 3.34e5
    time = np.linspace(86400, 100 * 86400, 200)
    thickness = stefan_ice_thickness(difference * time, conductivity, RHO, latent)
    step = 1.0
    later, earlier = (stefan_ice_thickness(difference * (time + d), conductivity, RHO, latent) for d in (step, -step))
    rate = stefan_growth_rate(thickness, difference, conductivity, RHO, latent)
    np.testing.assert_allclose((later - earlier) / (2 * step), rate, rtol=1e-6)
    freeboard = sea_ice_freeboard(thickness, RHO, 1025)
    np.testing.assert_allclose(1025 * (thickness - freeboard), RHO * thickness)


def test_rebound_relaxes_haskell_style():
    tau = viscous_relaxation_time(3e6, 1e21, 3300, G)
    assert np.isclose(haskell_viscosity(tau, 3e6, 3300, G), 1e21)
    assert np.isclose(viscous_relaxation_time(1.5e6, 1e21, 3300, G), 2 * tau)
    uplift = quad(lambda t: rebound_uplift_rate(t, 250, tau), 0, 20 * tau)[0]
    assert np.isclose(uplift, 250 - remaining_depression(20 * tau, 250, tau))
    assert np.isclose(remaining_depression(tau, 250, tau), 250 / np.e)


def test_isotope_thermometer_inverts():
    temperature = np.array([220.0, 250.0, 270.0])
    np.testing.assert_allclose(dansgaard_temperature(dansgaard_delta(temperature, 0.69, -13.6), 0.69, -13.6), temperature)


def test_dansgaard_johnsen_age_integrates_its_velocity_and_tends_to_nye():
    thickness, kink, accumulation = 1370.0, 400.0, 0.35 / SECONDS_PER_YEAR
    depth = np.linspace(0, 1360, 2001)
    age = dansgaard_johnsen_age(depth, thickness, kink, accumulation)
    velocity = dansgaard_johnsen_downward_velocity(depth, thickness, kink, accumulation)
    step = 1e-3
    deeper = dansgaard_johnsen_age(depth[1:] + step, thickness, kink, accumulation)
    shallower = dansgaard_johnsen_age(depth[1:] - step, thickness, kink, accumulation)
    np.testing.assert_allclose((deeper - shallower) / (2 * step), 1 / velocity[1:], rtol=1e-6)
    assert dansgaard_johnsen_downward_velocity(0, thickness, kink, accumulation) == accumulation
    assert dansgaard_johnsen_downward_velocity(thickness, thickness, kink, accumulation) == 0
    nye = nye_age(depth[:-1], thickness, accumulation)
    np.testing.assert_allclose(dansgaard_johnsen_age(depth[:-1], thickness, 1e-6, accumulation), nye, rtol=1e-6)
    assert np.all(age[1:] > nye_age(depth[1:], thickness, accumulation))


def test_grounding_line_flux_rises_steeply_with_depth():
    shelf = ice_shelf_strain_rate(500, A, N, RHO, 1028, G)
    assert np.isclose(shelf, glen_strain_rate(RHO * G * 500 * (1 - RHO / 1028) / 4, A, N))
    assert np.isclose(flotation_thickness(800, RHO, 1028) * RHO, 800 * 1028)
    ratio = schoof_grounding_line_flux(2000, A, N, 7.6e6, 1 / 3, RHO, 1028, G) / schoof_grounding_line_flux(
        1000, A, N, 7.6e6, 1 / 3, RHO, 1028, G
    )
    assert np.isclose(ratio, 2**4.75)


BED = 3000 - 0.1 * np.arange(0, 20e3 + 1, 200.0)
GRADIENT = 0.007 / SECONDS_PER_YEAR


def test_flowline_glacier_reaches_a_steady_flux_balance():
    time = np.linspace(0, 1500, 31) * SECONDS_PER_YEAR
    result = simulate_flowline_glacier(BED, 200.0, 2600.0, GRADIENT, time)
    thickness = result.values[-1, 2:]
    assert result.values[0, 1] == 0 and np.isclose(result.values[-1, 1], result.values[-2, 1], rtol=1e-6)
    balance = 200.0 * np.cumsum(linear_mass_balance(BED + thickness, 2600.0, GRADIENT))
    slope = -np.diff(BED + thickness) / 200.0
    flux = shallow_ice_flux(0.5 * (thickness[1:] + thickness[:-1]), slope, 2.4e-24, 3, 917, 9.80665)
    interior = thickness[1:] > 20
    np.testing.assert_allclose(flux[interior], balance[:-1][interior], rtol=1e-3, atol=1e-9)
    assert result.names[:3] == ("length", "area", "thickness_0") and result.units[1] == "m^2"


def test_raising_the_equilibrium_line_shrinks_the_glacier():
    time = np.linspace(0, 1000, 21) * SECONDS_PER_YEAR
    steady = simulate_flowline_glacier(BED, 200.0, 2600.0, GRADIENT, time).values[-1, 2:]
    warmer = simulate_flowline_glacier(BED, 200.0, 2700.0, GRADIENT, time, initial_thickness=steady)
    assert np.all(np.diff(warmer.values[:, 1]) <= 1e-8 * warmer.values[0, 1]) and warmer.values[-1, 0] < warmer.values[0, 0]
    bare = simulate_flowline_glacier(BED, 200.0, 3500.0, GRADIENT, time[:3])
    np.testing.assert_array_equal(bare.values, 0)


@pytest.mark.parametrize(
    "args,kwargs",
    [
        ((BED[:2], 200.0, 2600.0, GRADIENT, [0, 1]), {}),
        ((BED, [200.0, 100.0], 2600.0, GRADIENT, [0, 1]), {}),
        ((BED, 200.0, 2600.0, GRADIENT, [0, 1]), {"initial_thickness": np.zeros(3)}),
        ((BED, 200.0, 2600.0, 0.0, [0, 1]), {}),
    ],
)
def test_flowline_glacier_rejects_bad_setups(args, kwargs):
    with pytest.raises(ValueError):
        simulate_flowline_glacier(*args, **kwargs)
