"""Physical identities, limits and dynamic checks for earthsciencekit.geophysics."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from earthsciencekit.constants import MU0
from earthsciencekit.geophysics.systems.deformation import equivalent_water_height, love_number_h2, love_number_k2
from earthsciencekit.geophysics.systems.dynamo import rikitake_fixed_point, simulate_rikitake_dynamo
from earthsciencekit.geophysics.systems.figure import (
    centrifugal_ratio,
    clairaut_gravity,
    degree_geoid_height,
    flattening_from_j2,
    nodal_precession_rate,
    normal_gravity,
    stokes_function,
)
from earthsciencekit.geophysics.systems.gravity import bouguer_anomaly, bouguer_correction, free_air_correction, sphere_gravity
from earthsciencekit.geophysics.systems.magnetism import (
    block_magnetic_anomaly,
    dipole_field,
    dipole_inclination,
    direction_angle,
    field_elements,
    paleolatitude,
    virtual_geomagnetic_pole,
    zonal_field,
)


@pytest.mark.parametrize(
    "function,args",
    [
        (sphere_gravity, [0.5, 0, 1, 1000]),
        (nodal_precession_rate, [1.5, 0.5, 0, 1e-3, 1, 1]),
        (field_elements, [0, 0, 0]),
        (block_magnetic_anomaly, [0, 1, 1, 1, 2, 1]),
        (block_magnetic_anomaly, [0, -1, 1, 2, 2, 1]),
    ],
)
def test_cross_parameter_rejections(function, args):
    with pytest.raises(ValueError):
        function(*args)


def test_sphere_gravity_peak_and_half_width():
    depth = 1000.0
    field = sphere_gravity(np.array([-depth * 0.766, 0, depth * 0.766]), depth, 200, 500)
    assert np.isclose(field[0, 1], field[1, 1] / 2, rtol=1e-3)
    assert field[0, 0] > 0 > field[2, 0] and field[1, 0] == 0


def test_bouguer_anomaly_recovers_buried_body_on_topography():
    height = np.array([0.0, 250.0, 800.0])
    radius, density, buried = 6.371e6, 2670.0, 3e-5
    normal = normal_gravity(0.7)
    observed = normal - free_air_correction(height, normal, radius) + bouguer_correction(density, height) + buried
    np.testing.assert_allclose(bouguer_anomaly(observed, normal, height, density, radius), buried, rtol=1e-9)


def test_clairaut_matches_somigliana_to_first_order():
    latitude = np.linspace(0, np.pi / 2, 19)
    flattening, ratio = 1 / 298.257222101, 0.00344978600308
    np.testing.assert_allclose(clairaut_gravity(latitude, normal_gravity(0), flattening, ratio), normal_gravity(latitude), rtol=5e-5)


def test_satellite_precession_returns_flattening():
    gm, radius, j2 = 3.986004418e14, 6378137.0, 1.08263e-3
    rate = nodal_precession_rate(7.3e6, 0.1, np.radians(65.3), j2, radius, gm)
    recovered = rate / nodal_precession_rate(7.3e6, 0.1, np.radians(65.3), 1, radius, gm)
    assert np.isclose(recovered, j2)
    assert nodal_precession_rate(7.3e6, 0, np.pi / 2, j2, radius, gm) == pytest.approx(0, abs=1e-20)
    ratio = centrifugal_ratio(7.292115e-5, radius, 9.7803267715)
    assert np.isclose(1 / flattening_from_j2(j2, ratio), 298.257, rtol=2e-3)


@pytest.mark.parametrize("degree", [2, 3, 5])
def test_stokes_integral_of_a_zonal_anomaly(degree):
    def integrand(psi):
        return eval_legendre(degree, np.cos(psi)) * stokes_function(psi) * np.sin(psi)

    integral = quad(integrand, 0, np.pi, limit=200)[0] / 2
    assert np.isclose(integral * 6.371e6 / 9.8, degree_geoid_height(1, degree, 6.371e6, 9.8), rtol=1e-6)


def test_dipole_directions_poles_and_paleolatitude():
    latitude = np.linspace(-1.4, 1.4, 15)
    field = dipole_field(latitude, 6.371e6, 8e22)
    elements = field_elements(field[..., 0], 0, field[..., 1])
    np.testing.assert_allclose(elements[..., 1], dipole_inclination(latitude))
    np.testing.assert_allclose(paleolatitude(dipole_inclination(latitude)), latitude)
    poles = dipole_field([np.pi / 2, 0], 6.371e6, 8e22)
    assert np.isclose(poles[0, 1], 2 * poles[1, 0])
    pole = virtual_geomagnetic_pole(latitude, 2.0, 0, dipole_inclination(latitude))
    np.testing.assert_allclose(pole[..., 0], np.pi / 2)
    reversed_pole = virtual_geomagnetic_pole(0.6, 2.0, np.pi, -dipole_inclination(0.6))
    assert np.isclose(reversed_pole[0], -np.pi / 2)


def test_zonal_dipole_term_is_the_axial_dipole():
    latitude, distance, g10, a = np.linspace(-1.5, 1.5, 7), 8e6, -29404.8e-9, 6371.2e3
    moment = -4 * np.pi * a**3 * g10 / MU0
    np.testing.assert_allclose(zonal_field(latitude, distance, 1, g10, a), dipole_field(latitude, distance, moment), rtol=1e-12)
    far, near = zonal_field(0.4, [2 * a, a], 3, 1e-6, a)[:, 1]
    assert np.isclose(far / near, 2.0**-5)


def test_vgp_of_a_tilted_dipole_recovers_its_pole():
    pole_latitude, pole_longitude = np.radians(70), np.radians(300)
    site_latitude, site_longitude = np.radians([10, 45, -30]), np.radians([20, 100, 250])
    cos_p = np.sin(site_latitude) * np.sin(pole_latitude) + np.cos(site_latitude) * np.cos(pole_latitude) * np.cos(
        pole_longitude - site_longitude
    )
    p = np.arccos(cos_p)
    declination = np.arctan2(
        np.cos(pole_latitude) * np.sin(pole_longitude - site_longitude),
        np.cos(site_latitude) * np.sin(pole_latitude)
        - np.sin(site_latitude) * np.cos(pole_latitude) * np.cos(pole_longitude - site_longitude),
    )
    inclination = np.arctan(2 / np.tan(p))
    pole = virtual_geomagnetic_pole(site_latitude, site_longitude, declination, inclination)
    np.testing.assert_allclose(pole, np.column_stack([np.full(3, pole_latitude), np.full(3, pole_longitude)]))


def test_direction_angle_of_a_reversal_is_pi():
    inclination = dipole_inclination(0.6)
    assert np.isclose(direction_angle(0, inclination, np.pi, -inclination), np.pi)
    assert direction_angle(0.3, 0.2, 0.3, 0.2) == 0


def test_block_anomaly_superposes_and_cancels_for_a_layer():
    x = np.linspace(-5e4, 5e4, 11)
    whole = block_magnetic_anomaly(x, -1e4, 1e4, 2500, 3000, 5)
    parts = block_magnetic_anomaly(x, -1e4, 0, 2500, 3000, 5) + block_magnetic_anomaly(x, 0, 1e4, 2500, 3000, 5)
    np.testing.assert_allclose(whole, parts)
    np.testing.assert_allclose(block_magnetic_anomaly(x, -1e12, 1e12, 2500, 3000, 5), 0, atol=1e-6 * np.abs(whole).max())
    np.testing.assert_allclose(block_magnetic_anomaly(x, -1e4, 1e4, 2500, 3000, -5), -whole)


def test_love_numbers_fluid_and_rigid_limits():
    assert love_number_k2(0, 5500, 9.8, 6.371e6) == 1.5
    assert love_number_h2(0, 5500, 9.8, 6.371e6) == 2.5
    stiff = love_number_k2(1e15, 5500, 9.8, 6.371e6)
    assert 0 < stiff < 1e-3
    assert np.isclose(love_number_h2(1e11, 5500, 9.8, 6.371e6), 5 / 3 * love_number_k2(1e11, 5500, 9.8, 6.371e6))


def test_equivalent_water_height_is_linear_and_rejects_noninteger_degree():
    single = equivalent_water_height(1e-10, 10, -0.069, 6.371e6, 5514, 1000)
    assert np.isclose(equivalent_water_height(3e-10, 10, -0.069, 6.371e6, 5514, 1000), 3 * single)
    with pytest.raises(ValueError):
        equivalent_water_height(1e-10, 10.5, -0.069, 6.371e6, 5514, 1000)


def test_rikitake_fixed_point_is_steady_and_flow_reverses():
    mu, k = 1.0, 2.0
    time = np.linspace(0, 20, 41)
    steady = simulate_rikitake_dynamo(mu, k, 1.0, time, initial_state=rikitake_fixed_point(mu, k))
    np.testing.assert_allclose(steady.values, np.broadcast_to(rikitake_fixed_point(mu, k), steady.values.shape), atol=1e-7)
    mirrored = simulate_rikitake_dynamo(mu, k, 1.0, time, initial_state=-rikitake_fixed_point(mu, k) * [1, 1, -1])
    np.testing.assert_allclose(mirrored.values[:, :2], -steady.values[:, :2], atol=1e-7)
    chaotic = simulate_rikitake_dynamo(mu, k, 2.0, np.linspace(0, 400, 2001))
    assert np.count_nonzero(np.diff(np.sign(chaotic.values[:, 0]))) > 4
    assert chaotic.names == ("x", "y", "z") and chaotic.units == ("1", "1", "1")
    for args in [(0, k, 1.0), (mu, k, [1.0, 2.0])]:
        with pytest.raises(ValueError):
            simulate_rikitake_dynamo(*args, time)
    with pytest.raises(ValueError):
        simulate_rikitake_dynamo(mu, k, 1.0, time, initial_state=[1, 0])


def _random_coefficients(degree, seed):
    rng = np.random.default_rng(seed)
    g, h = np.tril(rng.normal(size=(degree + 1, degree + 1))), np.tril(rng.normal(size=(degree + 1, degree + 1)))
    g[0], h[0], h[:, 0] = 0, 0, 0
    return g, h


def test_schmidt_functions_have_the_igrf_normalization():
    from earthsciencekit.geophysics import schmidt_legendre

    for n in range(5):
        for m in range(n + 1):
            # The longitude mean of cos(m lon)**2 is 1/2 for m > 0, which Schmidt's factor 2 compensates.
            mean = quad(lambda lat: schmidt_legendre(n, m, lat) ** 2 * np.cos(lat), -np.pi / 2, np.pi / 2)[0] / 2
            np.testing.assert_allclose(mean * (0.5 if m else 1.0), 1 / (2 * n + 1), rtol=1e-10)


def test_geomagnetic_field_is_minus_the_gradient_of_its_potential():
    from earthsciencekit.geophysics import geomagnetic_field, schmidt_legendre

    g, h = _random_coefficients(4, 7)
    a = 6371.2e3

    def potential(latitude, longitude, radius):
        return a * sum(
            (a / radius) ** (n + 1) * (g[n, m] * np.cos(m * longitude) + h[n, m] * np.sin(m * longitude)) * schmidt_legendre(n, m, latitude)
            for n in range(1, 5)
            for m in range(n + 1)
        )

    lat, lon, r, step = 0.4, -1.1, 1.3 * a, 1e-5
    field = geomagnetic_field(lat, lon, r, g, h, a)
    # North is -B_theta = -(1/r) dV/dlatitude, east -dV/dlon / (r cos(lat)), down dV/dr.
    north = -(potential(lat + step, lon, r) - potential(lat - step, lon, r)) / (2 * step * r)
    east = -(potential(lat, lon + step, r) - potential(lat, lon - step, r)) / (2 * step * r * np.cos(lat))
    down = (potential(lat, lon, r * (1 + step)) - potential(lat, lon, r * (1 - step))) / (2 * step * r)
    np.testing.assert_allclose(field, [north, east, down], rtol=1e-7)


def test_geomagnetic_field_zonal_terms_and_dipole():
    from earthsciencekit.geophysics import dipole_field, gauss_dipole_moment, geomagnetic_field, zonal_field

    g, h = _random_coefficients(3, 3)
    zonal_g = g * (np.arange(4)[None, :] == 0)
    latitude = np.linspace(-1.4, 1.4, 9)
    field = geomagnetic_field(latitude, 0.7, 7e6, zonal_g, np.zeros_like(h), 6.4e6)
    expected = sum(zonal_field(latitude, 7e6, n, zonal_g[n, 0], 6.4e6) for n in range(1, 4))
    np.testing.assert_allclose(field[:, [0, 2]], expected, rtol=1e-12, atol=1e-30)
    np.testing.assert_allclose(field[:, 1], 0, atol=1e-30)
    axial = np.zeros((2, 2))
    axial[1, 0] = -3e-5
    moment = gauss_dipole_moment(-3e-5, 0, 0, 6.4e6)
    dipole = dipole_field(latitude, 7e6, moment)
    np.testing.assert_allclose(geomagnetic_field(latitude, 0.2, 7e6, axial, np.zeros((2, 2)), 6.4e6)[:, [0, 2]], dipole, rtol=1e-12)
