"""Physical identities, limits and dynamic checks for earthsciencekit.atmosphere."""

import numpy as np
import pytest

from earthsciencekit.atmosphere.systems.dynamics import rossby_wave_speed, simulate_barotropic_vorticity
from earthsciencekit.atmosphere.systems.hydrostatics import barometer_column_height, surface_pressure_tendency, two_layer_temperature
from earthsciencekit.atmosphere.systems.thermodynamics import potential_temperature, richardson_number, simulate_dry_ascent


def test_dry_ascent_conserves_potential_temperature():
    time = np.linspace(0, 1000, 51)
    result = simulate_dry_ascent(300, 1e5, 1, time)
    theta = potential_temperature(result.values[:, 0], result.values[:, 1], 1e5)
    np.testing.assert_allclose(theta, 300, rtol=1e-8)
    assert np.all(np.diff(result.values[:, 1]) < 0)
    with pytest.raises(ValueError):
        simulate_dry_ascent([300, 310], 1e5, 1, time)
    with pytest.raises(ValueError):
        simulate_dry_ascent(300, 1e5, 100, time)


def _barotropic_grid(n=16, spacing=6e5):
    x = np.arange(n) * spacing
    return np.meshgrid(x, x), n * spacing


def test_barotropic_single_mode_is_rossby_wave():
    beta, spacing = 1.6e-11, 6e5
    (X, Y), length = _barotropic_grid(spacing=spacing)
    k, l = 2 * np.pi * 2 / length, 2 * np.pi / length
    time = np.linspace(3600, 2 * 86400, 7)
    result = simulate_barotropic_vorticity(1e-5 * np.cos(k * X + l * Y), time, beta, spacing)
    c = rossby_wave_speed(0, beta, np.hypot(k, l))
    expected = 1e-5 * np.cos(k * (X[None] - c * (time - time[0])[:, None, None]) + l * Y[None])
    np.testing.assert_allclose(result.values, expected.reshape(time.size, -1), atol=1e-10)
    assert result.names[17] == "cell_1_1" and result.units[0] == "1/s"


def test_barotropic_conserves_mean_energy_and_enstrophy():
    spacing = 6e5
    (X, Y), length = _barotropic_grid(spacing=spacing)
    rng = np.random.default_rng(1)
    modes = [(i, j) for i in range(1, 4) for j in range(-3, 4)]
    zeta = 2e-6 + sum(1e-5 * rng.normal() * np.cos(2 * np.pi * (i * X + j * Y) / length + rng.uniform(0, 2 * np.pi)) for i, j in modes)
    result = simulate_barotropic_vorticity(zeta, np.linspace(0, 86400, 5), 1.6e-11, spacing)
    fields = result.values.reshape(-1, *X.shape)
    np.testing.assert_allclose(fields.mean(axis=(1, 2)), 0, atol=1e-20)
    spectra = np.abs(np.fft.fft2(fields)) ** 2
    k2 = (2 * np.pi * np.fft.fftfreq(X.shape[0], spacing)[:, None]) ** 2 + (2 * np.pi * np.fft.fftfreq(X.shape[1], spacing)[None, :]) ** 2
    energy = np.sum(np.divide(spectra, k2, out=np.zeros_like(spectra), where=k2 > 0), axis=(1, 2))
    np.testing.assert_allclose(spectra.sum(axis=(1, 2)), spectra[0].sum(), rtol=1e-5)
    np.testing.assert_allclose(energy, energy[0], rtol=1e-5)
    assert not np.allclose(fields[-1], fields[0])


@pytest.mark.parametrize(
    "zeta,beta,spacing",
    [(np.zeros((3, 8)), 1e-11, 1e5), (np.zeros(16), 1e-11, 1e5), (np.zeros((4, 4)), [1e-11, 1e-11], 1e5), (np.zeros((4, 4)), 1e-11, 0)],
)
def test_barotropic_rejects_bad_setup(zeta, beta, spacing):
    with pytest.raises(ValueError):
        simulate_barotropic_vorticity(zeta, [0, 1], beta, spacing)


def test_two_layer_profile_hydrostatic_pressure_matches_standard_atmosphere():
    height = np.linspace(0, 11000, 2001)
    temperature = two_layer_temperature(height, 288.15, 0.0065, 11000)
    integrand = 9.80665 / (287.05 * temperature)
    log_ratio = np.sum((integrand[1:] + integrand[:-1]) / 2 * np.diff(height))
    assert np.isclose(101325 * np.exp(-log_ratio), 101325 * (216.65 / 288.15) ** (9.80665 / (287.05 * 0.0065)), rtol=1e-7)
    assert two_layer_temperature(20000, 288.15, 0.0065, 11000) == two_layer_temperature(11000, 288.15, 0.0065, 11000)


def test_richardson_and_rossby_limits():
    assert richardson_number(-1e-4, 0.01) < 0 < richardson_number(1e-4, -0.01)
    beta, wind = 1.6e-11, 15.0
    assert np.isclose(rossby_wave_speed(wind, beta, np.sqrt(beta / wind)), 0, atol=1e-12)
    assert rossby_wave_speed(wind, 0, 1e-6) == wind
    assert np.isclose(surface_pressure_tendency(-1e-6, 1e5) * 21600, 2160)
    assert np.isclose(barometer_column_height(101325, 13595.1, 9.80665), 0.76, rtol=1e-4)
