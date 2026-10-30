"""Physical identities, limits and dynamic checks for earthsciencekit.seismology."""

import numpy as np
import pytest

from earthsciencekit.seismology.systems.interfaces import reflection_coefficient, transmission_coefficient
from earthsciencekit.seismology.systems.rays import (
    critical_angle,
    crossover_distance,
    head_wave_travel_time,
    ray_parameter,
    travel_time,
    vertical_travel_time,
)
from earthsciencekit.seismology.systems.sources import local_magnitude, moment_from_magnitude, moment_magnitude, simulate_oscillator
from earthsciencekit.seismology.systems.waves import love_wave_angular_frequency, poisson_ratio, rayleigh_wave_speed, ricker_wavelet


def test_seismic_oscillator_analytic():
    frequency, Q = 1.0, 10.0
    time = np.linspace(0, 2, 101)
    result = simulate_oscillator(frequency, Q, time)
    omega = 2 * np.pi * frequency
    damping = omega / (2 * Q)
    wd = np.sqrt(omega**2 - damping**2)
    expected = np.exp(-damping * time) * (np.cos(wd * time) + damping / wd * np.sin(wd * time))
    np.testing.assert_allclose(result.values[:, 0], expected, atol=3e-8)
    energy = result.values[:, 1] ** 2 + omega**2 * result.values[:, 0] ** 2
    assert np.all(np.diff(energy) <= 1e-8)
    with pytest.raises(ValueError):
        simulate_oscillator([1, 2], Q, time)


def test_wave_interface_energy_balance():
    z1, z2 = 2.0, 6.0
    r = reflection_coefficient(z1, z2)
    t = transmission_coefficient(z1, z2)
    assert np.isclose(r * r + (z1 / z2) * t * t, 1)
    assert np.isclose(moment_magnitude(moment_from_magnitude(7)), 7)
    np.testing.assert_allclose(ricker_wavelet([-1, 0, 1], 2), ricker_wavelet([1, 0, -1], 2))


def test_rayleigh_speed_solves_rayleigh_secular_equation():
    vs = 1.0
    vp = np.array([1.2, np.sqrt(3), 3.0, 10.0])
    c = rayleigh_wave_speed(vp, vs)
    assert np.all(c < vs)
    np.testing.assert_allclose((2 - c**2 / vs**2) ** 2, 4 * np.sqrt(1 - c**2 / vp**2) * np.sqrt(1 - c**2 / vs**2), atol=1e-12)
    # Viktorov's approximation (0.862 + 1.14 nu) / (1 + nu) holds to about 0.3% for nu >= 0.25.
    nu = poisson_ratio(vp[1:], vs)
    np.testing.assert_allclose(c[1:], (0.862 + 1.14 * nu) / (1 + nu), rtol=3e-3)


def test_head_wave_crosses_direct_wave_and_starts_at_reflection():
    h, v1, v2 = 30e3, 6e3, 8e3
    x_cross = crossover_distance(h, v1, v2)
    assert np.isclose(head_wave_travel_time(x_cross, h, v1, v2), travel_time(x_cross, v1))
    # The closed form the model checks against; tan(critical_angle) can round just below it.
    x_crit = 2 * h * v1 / np.sqrt(v2**2 - v1**2)
    assert np.allclose(x_crit, 2 * h * np.tan(critical_angle(v1, v2)))
    assert np.isclose(head_wave_travel_time(x_crit, h, v1, v2), np.hypot(x_crit, 2 * h) / v1)
    x = np.array([2e5, 3e5])
    assert np.isclose(np.diff(head_wave_travel_time(x, h, v1, v2))[0], 1e5 / v2)


def test_love_dispersion_root_and_limits():
    h, b1, b2, m1, m2 = 35e3, 3.5e3, 4.5e3, 3.4e10, 6.7e10
    c = np.linspace(3.6e3, 4.4e3, 9)
    omega = love_wave_angular_frequency(c, h, b1, b2, m1, m2)
    k = omega / c
    s1, s2 = np.sqrt(c**2 / b1**2 - 1), np.sqrt(1 - c**2 / b2**2)
    np.testing.assert_allclose(np.tan(k * h * s1), m2 * s2 / (m1 * s1))
    assert np.all(np.diff(omega) < 0)
    assert love_wave_angular_frequency(b2 * (1 - 1e-12), h, b1, b2, m1, m2) < 1e-5
    assert love_wave_angular_frequency(b1 * (1 + 1e-12), h, b1, b2, m1, m2) > 1e3


def test_ray_parameter_is_snell_on_a_sphere():
    r, v_above, v_below = 3.48e6, 13e3, 8e3
    above = np.deg2rad(50.0)
    p = ray_parameter(r, v_above, above)
    below = np.arcsin(p * v_below / r)
    assert np.isclose(ray_parameter(r, v_below, below), p)
    # A straight ray in a uniform shell keeps r sin(i) fixed, hence the same p at every radius.
    impact = r * np.sin(above)
    for radius in [r, 4e6, 6.371e6]:
        assert np.isclose(ray_parameter(radius, v_above, np.arcsin(impact / radius)), p)


def test_local_magnitude_is_logarithmic_in_amplitude():
    assert np.isclose(local_magnitude(1e-6, 1e5), 0)
    np.testing.assert_allclose(np.diff(local_magnitude([1e-4, 1e-3, 1e-2], 3e5)), 1)
    assert local_magnitude(1e-3, 3e5) > local_magnitude(1e-3, 1e5)


def test_vertical_travel_time_profiles():
    depth = np.linspace(0, 1e5, 11)
    np.testing.assert_allclose(vertical_travel_time(depth, np.full(11, 5e3)), travel_time(depth, 5e3))
    v0, gradient = 5e3, 0.02
    np.testing.assert_allclose(vertical_travel_time(depth, v0 + gradient * depth), np.log1p(gradient * depth / v0) / gradient)
    jump = vertical_travel_time([0, 1e4, 1e4, 2e4], [5e3, 5e3, 8e3, 8e3])
    np.testing.assert_allclose(jump, [0, 2, 2, 2 + 1.25])


def test_layered_rays_in_a_uniform_sphere_are_chords():
    from earthsciencekit.seismology import layered_ray_distance, layered_ray_travel_time

    radius, speed = 6371e3, 8e3
    p = np.linspace(10, 0.999 * radius / speed, 12).reshape(3, 4)
    table = ([0, 0, 3e6, radius - 1e3], [speed] * 4, radius)
    distance = layered_ray_distance(p, *table)
    assert distance.shape == (3, 4)
    np.testing.assert_allclose(distance, 2 * np.arccos(p * speed / radius), rtol=1e-12)
    np.testing.assert_allclose(layered_ray_travel_time(p, *table), 2 * np.sqrt(radius**2 - (p * speed) ** 2) / speed, rtol=1e-12)


def test_layered_rays_match_quadrature_and_benndorf_relation():
    from scipy.integrate import quad
    from scipy.optimize import brentq

    from earthsciencekit.seismology import layered_ray_distance, layered_ray_travel_time

    radius = 6371e3
    depth, speed = np.array([0, 2000e3, 4000e3]), np.array([8e3, 11e3, 13e3])
    b = np.log(speed[1:] / speed[:-1]) / np.log((radius - depth[1:]) / (radius - depth[:-1]))

    def eta(r):
        layer = 0 if r >= radius - depth[1] else 1
        return r / (speed[layer] * (r / (radius - depth[layer])) ** b[layer])

    p = 300.0
    turning = brentq(lambda r: eta(r) - p, radius - depth[2], radius - depth[1])
    width = (radius - depth[1] - turning, depth[1])

    def integral(f, start, span):
        # r = start + u**2 removes the inverse-square-root singularity at the turning radius.
        return quad(lambda u: 2 * u * f(start + u**2) / ((start + u**2) * np.sqrt(eta(start + u**2) ** 2 - p**2)), 0, np.sqrt(span))[0]

    integrals = [
        integral(f, start, span) for f in (lambda r: p, lambda r: eta(r) ** 2) for start, span in zip((turning, radius - depth[1]), width)
    ]
    np.testing.assert_allclose(layered_ray_distance(p, depth, speed, radius), 2 * sum(integrals[:2]), rtol=1e-7)
    np.testing.assert_allclose(layered_ray_travel_time(p, depth, speed, radius), 2 * sum(integrals[2:]), rtol=1e-7)
    near = np.array([p - 0.01, p + 0.01])
    slope = np.diff(layered_ray_travel_time(near, depth, speed, radius)) / np.diff(layered_ray_distance(near, depth, speed, radius))
    np.testing.assert_allclose(slope, p, rtol=1e-6)


def test_layered_rays_reproduce_the_two_shell_core_shadow():
    from earthsciencekit.seismology import layered_ray_distance

    R, rc, vm, vc = 6371e3, 3480e3, 13e3, 8e3
    table = ([0, R - rc, R - rc, R - 1e3], [vm, vm, vc, vc], R)
    p = np.linspace(1, 0.999 * R / vm, 2001)
    expected = np.where(
        p * vm >= rc,
        2 * np.arccos(p * vm / R),
        2 * (np.arccos(p * vm / R) - np.arccos(np.minimum(p * vm / rc, 1))) + np.pi - 2 * np.arcsin(np.minimum(p * vc / rc, 1)),
    )
    np.testing.assert_allclose(layered_ray_distance(p, *table), expected, rtol=1e-10)
