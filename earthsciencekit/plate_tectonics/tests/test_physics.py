"""Physical identities, limits and dynamic checks for earthsciencekit.plate_tectonics."""

import numpy as np
import pytest

from earthsciencekit.plate_tectonics.systems.cooling import diffuse_temperature, halfspace_heat_flux, mantle_rayleigh_number, ocean_depth
from earthsciencekit.plate_tectonics.systems.isostasy import flexural_deflection, flexural_parameter
from earthsciencekit.plate_tectonics.systems.kinematics import (
    angular_speed,
    compose_euler_vectors,
    euler_velocity_x,
    euler_velocity_y,
    euler_velocity_z,
    hotspot_age,
    plate_speed,
)


def test_periodic_diffusion_conservation_and_eigenmode():
    n, spacing, diffusivity = 16, 2.0, 0.1
    initial = 300 + 10 * np.cos(2 * np.pi * np.arange(n) / n)
    time = np.linspace(10, 110, 31)
    result = diffuse_temperature(initial, time, diffusivity, spacing)
    np.testing.assert_allclose(result.values.mean(axis=1), 300, atol=1e-12)
    rate = 4 * diffusivity * np.sin(np.pi / n) ** 2 / spacing**2
    expected = 300 + (initial - 300)[None, :] * np.exp(-rate * (time - time[0]))[:, None]
    np.testing.assert_allclose(result.values, expected, atol=1e-12)
    assert result.values.min() >= initial.min() and result.values.max() <= initial.max()
    for initial, diffusion, spacing in [([300, 310], 1, 1), ([[300] * 3], 1, 1), ([300] * 3, [1, 2], 1), ([300] * 3, 0, 1)]:
        with pytest.raises(ValueError):
            diffuse_temperature(initial, [0, 1], diffusion, spacing)


def test_plate_rotation_tangent_and_speed():
    omega = np.array([1.0, 2.0, 3.0]) * 1e-9
    r = np.array([3.0, 4.0, 5.0]) * 1e6
    v = np.array(
        [
            euler_velocity_x(omega[1], omega[2], r[1], r[2]),
            euler_velocity_y(omega[2], omega[0], r[2], r[0]),
            euler_velocity_z(omega[0], omega[1], r[0], r[1]),
        ]
    )
    np.testing.assert_allclose(v, np.cross(omega, r))
    assert abs(np.dot(v, r)) < 1e-10
    distance = np.arccos(np.dot(omega, r) / (np.linalg.norm(omega) * np.linalg.norm(r)))
    assert np.isclose(plate_speed(np.linalg.norm(omega), np.linalg.norm(r), distance), np.linalg.norm(v))


def test_flexure_balances_load_and_places_bulge():
    rigidity, contrast, gravity, load = 1e23, 2300.0, 9.8, 2e12
    alpha = flexural_parameter(rigidity, contrast, gravity)
    x = np.linspace(-25 * alpha, 25 * alpha, 200001)
    w = flexural_deflection(x, load, rigidity, contrast, gravity)
    assert np.isclose(np.trapezoid(contrast * gravity * w, x), load, rtol=1e-6)
    w0 = flexural_deflection(0, load, rigidity, contrast, gravity)
    assert np.isclose(flexural_deflection(3 * np.pi * alpha / 4, load, rigidity, contrast, gravity), 0, atol=1e-12 * w0)
    assert np.isclose(flexural_deflection(np.pi * alpha, load, rigidity, contrast, gravity), -np.exp(-np.pi) * w0)


def test_ocean_depth_matches_heat_flux_and_parsons_sclater():
    year = 365.25 * 86400
    age = np.array([1, 4, 25]) * 1e6 * year
    depth = ocean_depth(age, 2500, 3300, 1000, 3e-5, 1e-6, 273, 1553)
    np.testing.assert_allclose((depth - 2500) / np.sqrt(age / (1e6 * year)), 350, rtol=0.01)
    np.testing.assert_allclose(np.diff(depth - 2500) / (depth - 2500)[:-1], [1, 1.5])
    # Subsidence is the thermal contraction of the heat lost since the ridge: alpha * int(q dt) / (c_p (rho_m - rho_w)).
    flux = halfspace_heat_flux(age, 1e-6, 3.3, 273, 1553)
    heat_capacity = 3.3 / (3300 * 1e-6)
    np.testing.assert_allclose(depth - 2500, 3e-5 * 2 * flux * age / (heat_capacity * (3300 - 1000)))


def test_euler_circuit_closes_and_velocities_add():
    omega_ab = np.array([1.0, -2.0, 0.5]) * 1e-9
    omega_bc = np.array([-0.3, 0.8, 2.0]) * 1e-9
    omega_ac = compose_euler_vectors(omega_ab, omega_bc)
    np.testing.assert_allclose(compose_euler_vectors(omega_ac, -omega_ac), 0)
    r = np.array([3.0, 4.0, 5.0]) * 1e6
    np.testing.assert_allclose(np.cross(omega_ac, r), np.cross(omega_ab, r) + np.cross(omega_bc, r))


def test_hotspot_age_is_rotation_angle_over_rate():
    omega, radius, colatitude = 2e-15, 6.371e6, np.array([0.3, 1.0, 2.5])
    angle = 0.2
    ages = hotspot_age(angle * radius * np.sin(colatitude), omega, radius, colatitude)
    np.testing.assert_allclose(ages, angle / omega)
    assert np.isclose(angular_speed(angle, ages[0]), omega)


def test_mantle_rayleigh_number_cubic_in_thickness():
    args = (3e-5, 4000, 9.8, 2500, 1e6, 1e-6, 1e21)
    ratio = mantle_rayleigh_number(*args[:4], 2 * args[4], *args[5:]) / mantle_rayleigh_number(*args)
    assert np.isclose(ratio, 8)
    assert mantle_rayleigh_number(3e-5, 4000, 9.8, 0, 1e6, 1e-6, 1e21) == 0
