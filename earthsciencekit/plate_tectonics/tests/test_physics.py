"""Physical identities, limits and dynamic checks for earthsciencekit.plate_tectonics."""

import numpy as np
import pytest
from scipy.integrate import trapezoid

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
    for bad_initial, diffusion, spacing in [([300, 310], 1, 1), ([[300] * 3], 1, 1), ([300] * 3, [1, 2], 1), ([300] * 3, 0, 1)]:
        with pytest.raises(ValueError):
            diffuse_temperature(bad_initial, [0, 1], diffusion, spacing)


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
    assert np.isclose(trapezoid(contrast * gravity * w, x), load, rtol=1e-6)
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


def _rotations(seed, size):
    rng = np.random.default_rng(seed)
    return np.arcsin(rng.uniform(-1, 1, size)), rng.uniform(-np.pi, np.pi, size), rng.uniform(-3, 3, size)


def test_finite_rotation_matrix_is_a_proper_rotation_about_its_pole():
    from earthsciencekit.plate_tectonics import finite_rotation_matrix, rotate_position

    latitude, longitude, angle = _rotations(1, 20)
    matrix = finite_rotation_matrix(latitude, longitude, angle)
    np.testing.assert_allclose(matrix @ np.swapaxes(matrix, -1, -2), np.broadcast_to(np.eye(3), matrix.shape), atol=1e-14)
    np.testing.assert_allclose(np.linalg.det(matrix), 1)
    pole = np.stack([np.cos(latitude) * np.cos(longitude), np.cos(latitude) * np.sin(longitude), np.sin(latitude)], axis=-1)
    np.testing.assert_allclose(np.einsum("nij,nj->ni", matrix, pole), pole, atol=1e-14)
    np.testing.assert_allclose(np.arccos(np.clip((np.trace(matrix, axis1=1, axis2=2) - 1) / 2, -1, 1)), np.abs(angle), atol=1e-7)
    point_lat, point_lon = _rotations(2, 20)[:2]
    moved = rotate_position(point_lat, point_lon, latitude, longitude, angle)
    before = np.stack([np.cos(point_lat) * np.cos(point_lon), np.cos(point_lat) * np.sin(point_lon), np.sin(point_lat)], axis=-1)
    after = np.stack([np.cos(moved[:, 0]) * np.cos(moved[:, 1]), np.cos(moved[:, 0]) * np.sin(moved[:, 1]), np.sin(moved[:, 0])], axis=-1)
    np.testing.assert_allclose(after, np.einsum("nij,nj->ni", matrix, before), atol=1e-14)
    # A rotation keeps every point at its angular distance from the pole: points move along small circles.
    np.testing.assert_allclose(np.sum(after * pole, axis=-1), np.sum(before * pole, axis=-1), atol=1e-14)


def test_composed_rotations_are_matrix_products_and_do_not_commute():
    from earthsciencekit.plate_tectonics import compose_finite_rotations, finite_rotation_matrix

    first, second = _rotations(3, 30), _rotations(4, 30)
    composite = compose_finite_rotations(*first, *second)
    product = finite_rotation_matrix(*second) @ finite_rotation_matrix(*first)
    np.testing.assert_allclose(finite_rotation_matrix(*np.moveaxis(composite, -1, 0)), product, atol=1e-13)
    assert np.all((composite[:, 2] >= 0) & (composite[:, 2] <= np.pi))
    reversed_order = compose_finite_rotations(*second, *first)
    assert not np.allclose(composite, reversed_order)
    identity = compose_finite_rotations(*first, first[0], first[1], -first[2])
    np.testing.assert_allclose(identity, np.broadcast_to([np.pi / 2, 0, 0], identity.shape), atol=1e-7)


def test_interpolated_rotation_turns_steadily_about_the_stage_pole():
    from earthsciencekit.plate_tectonics import compose_finite_rotations, finite_rotation_matrix, interpolate_finite_rotation

    first, second = (0.8, -0.6, 0.3), (0.2, 0.4, 1.1)
    for fraction, expected in ((0, first), (1, second)):
        rotation = interpolate_finite_rotation(*first, *second, fraction)
        np.testing.assert_allclose(finite_rotation_matrix(*rotation), finite_rotation_matrix(*expected), atol=1e-13)
    stage = compose_finite_rotations(first[0], first[1], -first[2], *second)
    fractions = np.linspace(0, 1, 11)
    path = interpolate_finite_rotation(*first, *second, fractions)
    # The stage pole is fixed in the reference frame, so each step is the same rotation about it.
    turned = compose_finite_rotations(first[0], first[1], -first[2], *np.moveaxis(path, -1, 0))
    np.testing.assert_allclose(turned[1:, 2], fractions[1:] * stage[2], atol=1e-12)
    np.testing.assert_allclose(turned[1:, :2], np.broadcast_to(stage[:2], (10, 2)), atol=1e-9)
    same_pole = interpolate_finite_rotation(0.5, 1.0, 0.2, 0.5, 1.0, 0.6, 0.25)
    np.testing.assert_allclose(same_pole, [0.5, 1.0, 0.3], atol=1e-12)
