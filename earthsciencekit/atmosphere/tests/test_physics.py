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


def test_lcl_is_where_the_lifted_parcel_saturates():
    from earthsciencekit.atmosphere import lcl_height, lcl_pressure, lcl_temperature, mixing_ratio, saturation_vapor_pressure

    temperature, dewpoint = np.array([300.0, 290.0, 280.0]), np.array([290.0, 289.0, 260.0])
    base = lcl_temperature(temperature, dewpoint)
    vapor = mixing_ratio(saturation_vapor_pressure(dewpoint), 1e5)
    at_base = mixing_ratio(saturation_vapor_pressure(base), lcl_pressure(1e5, temperature, base))
    # Bolton's fit is within about 0.1 K of saturation with this package's vapour-pressure formula.
    np.testing.assert_allclose(at_base, vapor, rtol=0.02)
    assert lcl_temperature(300, 300) == 300 and lcl_height(300, 300) == 0
    np.testing.assert_allclose(lcl_height(temperature, base) / (temperature - dewpoint), 125, rtol=0.03)


def test_moist_lapse_rate_tends_to_dry_rate_when_cold_and_is_smaller_when_warm():
    from earthsciencekit.atmosphere import dry_adiabatic_lapse_rate, moist_adiabatic_lapse_rate

    dry = dry_adiabatic_lapse_rate(9.80665, 1004)
    assert np.isclose(moist_adiabatic_lapse_rate(190, 1e5), dry, rtol=1e-3)
    rates = moist_adiabatic_lapse_rate(np.array([300, 280, 260, 240]), 1e5)
    assert np.all(np.diff(rates) > 0) and rates[0] < 0.4 * dry


def test_moist_ascent_follows_dry_adiabat_then_pseudo_adiabat():
    from earthsciencekit.atmosphere import lcl_pressure, lcl_temperature, parcel_temperature_profile, simulate_moist_ascent

    time = np.linspace(0, 2400, 121)
    result = simulate_moist_ascent(300, 293, 1e5, 5, time)
    temperature, pressure, vapor, condensed = result.values.T
    np.testing.assert_allclose(vapor + condensed, vapor[0], rtol=1e-12)
    assert condensed[0] == 0 and np.all(np.diff(condensed) >= -1e-15) and condensed[-1] > 0.01
    base = lcl_pressure(1e5, 300, lcl_temperature(300, 293))
    dry = pressure > base
    np.testing.assert_allclose(potential_temperature(temperature[dry], pressure[dry], 1e5), 300, rtol=1e-8)
    np.testing.assert_allclose(condensed[pressure > base * 1.01], 0)
    profile = parcel_temperature_profile(pressure, 300, 293)
    np.testing.assert_allclose(temperature, profile, atol=0.15)
    assert result.names[2:] == ("mixing_ratio", "condensed_water") and result.units[2] == "kg/kg"


@pytest.mark.parametrize(
    "args",
    [
        (300, 301, 1e5, 1, [0, 1]),
        (300, 290, 1e5, 0, [0, 1]),
        (300, 290, [1e5, 9e4], 1, [0, 1]),
        (400, 290, 1e5, 1, [0, 1]),
    ],
)
def test_moist_ascent_rejects_bad_setup(args):
    from earthsciencekit.atmosphere import simulate_moist_ascent

    with pytest.raises(ValueError):
        simulate_moist_ascent(*args)


def test_cape_and_cin_split_the_buoyancy_integral():
    from earthsciencekit.atmosphere import convective_available_potential_energy, convective_inhibition

    height = np.linspace(0, 12000, 61)
    environment = 300 - 0.0065 * height
    rng = np.random.default_rng(3)
    buoyancy_integral = []
    for offset in (-3.0, 3.0):
        parcel = environment + offset + rng.uniform(-1, 1, height.size)
        buoyancy = 9.80665 * (parcel - environment) / environment
        buoyancy_integral.append(np.sum((buoyancy[1:] + buoyancy[:-1]) / 2 * np.diff(height)))
        cape = convective_available_potential_energy(height, parcel, environment)
        cin = convective_inhibition(height, parcel, environment)
        if offset < 0:
            assert cape == 0 and np.isclose(cin, buoyancy_integral[-1])
        else:
            assert cin == 0 and np.isclose(cape, buoyancy_integral[-1])
    # A stable layer above a buoyant surface parcel does not count once the parcel is free.
    parcel = np.array([301.0, 299.0, 302.0, 299.0])
    assert convective_inhibition([0, 1000, 2000, 3000], parcel, np.full(4, 300.0)) == 0


def test_moist_ascent_saturated_start_and_dry_parcel():
    from earthsciencekit.atmosphere import parcel_temperature_profile, simulate_moist_ascent

    saturated = simulate_moist_ascent(290, 290, 9e4, 2, np.linspace(10, 610, 31))
    assert saturated.values[1, 3] > 0
    np.testing.assert_allclose(saturated.values[:, 0], parcel_temperature_profile(saturated.values[:, 1], 290, 290), atol=1e-6)
    parched = simulate_moist_ascent(300, 180, 1e5, 1, np.linspace(0, 100, 11))
    assert np.all(parched.values[:, 3] == 0)
    np.testing.assert_allclose(parched.values[:, 0], 300 - 9.80665 / 1004 * np.linspace(0, 100, 11))
