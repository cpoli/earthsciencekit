"""Golden checks and invalid-input contracts for every public model."""

from typing import Any

import numpy as np
import pytest

from earthsciencekit import oceanography as domain


def test_hydrostatic_pressure_reference():
    result = domain.hydrostatic_pressure(*[1000, 100, 10])
    np.testing.assert_allclose(result, 1e6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_hydrostatic_pressure_rejects_invalid(bad):
    args = [1000, 100, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hydrostatic_pressure(*args)


def test_linear_density_reference():
    result = domain.linear_density(*[1000, 1, 0, 0.0002, 0.0008])
    np.testing.assert_allclose(result, 999.8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_linear_density_rejects_invalid(bad):
    args = [1000, 1, 0, 0.0002, 0.0008]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.linear_density(*args)


def test_buoyancy_frequency_squared_reference():
    result = domain.buoyancy_frequency_squared(*[1000, -0.01, 10])
    np.testing.assert_allclose(result, 0.0001, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_buoyancy_frequency_squared_rejects_invalid(bad):
    args = [1000, -0.01, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.buoyancy_frequency_squared(*args)


def test_geostrophic_current_reference():
    result = domain.geostrophic_current(*[1e-05, 0.0001, 10])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_geostrophic_current_rejects_invalid(bad):
    args = [1e-05, 0.0001, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.geostrophic_current(*args)


def test_ekman_transport_reference():
    result = domain.ekman_transport(*[0.1, 1000, 0.0001])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ekman_transport_rejects_invalid(bad):
    args = [0.1, 1000, 0.0001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ekman_transport(*args)


def test_inertial_period_reference():
    result = domain.inertial_period(*[0.0001])
    np.testing.assert_allclose(result, 2 * np.pi * 1e4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_inertial_period_rejects_invalid(bad):
    args = [0.0001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.inertial_period(*args)


def test_shallow_wave_speed_reference():
    result = domain.shallow_wave_speed(*[100, 10])
    np.testing.assert_allclose(result, np.sqrt(1000), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_shallow_wave_speed_rejects_invalid(bad):
    args = [100, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shallow_wave_speed(*args)


def test_deep_wave_speed_reference():
    result = domain.deep_wave_speed(*[6.283185307179586, 10])
    np.testing.assert_allclose(result, np.sqrt(10), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_deep_wave_speed_rejects_invalid(bad):
    args = [6.283185307179586, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.deep_wave_speed(*args)


def test_wave_angular_frequency_reference():
    result = domain.wave_angular_frequency(*[1, 1, 10])
    np.testing.assert_allclose(result, np.sqrt(10 * np.tanh(1)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_wave_angular_frequency_rejects_invalid(bad):
    args = [1, 1, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.wave_angular_frequency(*args)


def test_deep_group_speed_reference():
    result = domain.deep_group_speed(*[6.283185307179586, 10])
    np.testing.assert_allclose(result, 0.5 * np.sqrt(10), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_deep_group_speed_rejects_invalid(bad):
    args = [6.283185307179586, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.deep_group_speed(*args)


def test_rossby_radius_reference():
    result = domain.rossby_radius(*[0.01, 1000, 0.0001])
    np.testing.assert_allclose(result, 1e5 / np.pi, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_rossby_radius_rejects_invalid(bad):
    args = [0.01, 1000, 0.0001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rossby_radius(*args)


def test_reduced_gravity_reference():
    result = domain.reduced_gravity(*[2, 1000, 10])
    np.testing.assert_allclose(result, 0.02, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_reduced_gravity_rejects_invalid(bad):
    args = [2, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.reduced_gravity(*args)


def test_internal_wave_speed_reference():
    result = domain.internal_wave_speed(*[0.02, 100, 100])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_internal_wave_speed_rejects_invalid(bad):
    args = [0.02, 100, 100]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.internal_wave_speed(*args)


def test_volume_transport_reference():
    result = domain.volume_transport(*[1, 1000, 100])
    np.testing.assert_allclose(result, 1e5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_volume_transport_rejects_invalid(bad):
    args = [1, 1000, 100]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.volume_transport(*args)


def test_tidal_elevation_reference():
    result = domain.tidal_elevation(*[0, 2, 12, 0])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_tidal_elevation_rejects_invalid(bad):
    args = [0, 2, 12, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.tidal_elevation(*args)


def test_hydrostatic_pressure_broadcasting():
    args: list[Any] = [1000, 100, 10]
    reference = domain.hydrostatic_pressure(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.hydrostatic_pressure(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_hydrostatic_pressure_density_physical_bound():
    args = [1000, 100, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.hydrostatic_pressure(*args)


def test_hydrostatic_pressure_depth_physical_bound():
    args = [1000, 100, 10]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.hydrostatic_pressure(*args)


def test_hydrostatic_pressure_gravity_physical_bound():
    args = [1000, 100, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.hydrostatic_pressure(*args)


def test_linear_density_broadcasting():
    args: list[Any] = [1000, 1, 0, 0.0002, 0.0008]
    reference = domain.linear_density(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.linear_density(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_linear_density_reference_density_physical_bound():
    args = [1000, 1, 0, 0.0002, 0.0008]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.linear_density(*args)


def test_linear_density_thermal_expansion_physical_bound():
    args = [1000, 1, 0, 0.0002, 0.0008]
    args[3] = -1
    with pytest.raises(ValueError):
        domain.linear_density(*args)


def test_linear_density_haline_contraction_physical_bound():
    args = [1000, 1, 0, 0.0002, 0.0008]
    args[4] = -1
    with pytest.raises(ValueError):
        domain.linear_density(*args)


def test_buoyancy_frequency_squared_broadcasting():
    args: list[Any] = [1000, -0.01, 10]
    reference = domain.buoyancy_frequency_squared(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.buoyancy_frequency_squared(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_buoyancy_frequency_squared_density_physical_bound():
    args = [1000, -0.01, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.buoyancy_frequency_squared(*args)


def test_buoyancy_frequency_squared_gravity_physical_bound():
    args = [1000, -0.01, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.buoyancy_frequency_squared(*args)


def test_geostrophic_current_broadcasting():
    args: list[Any] = [1e-05, 0.0001, 10]
    reference = domain.geostrophic_current(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.geostrophic_current(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_geostrophic_current_coriolis_physical_bound():
    args = [1e-05, 0.0001, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.geostrophic_current(*args)


def test_geostrophic_current_gravity_physical_bound():
    args = [1e-05, 0.0001, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.geostrophic_current(*args)


def test_ekman_transport_broadcasting():
    args: list[Any] = [0.1, 1000, 0.0001]
    reference = domain.ekman_transport(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.ekman_transport(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_ekman_transport_density_physical_bound():
    args = [0.1, 1000, 0.0001]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.ekman_transport(*args)


def test_ekman_transport_coriolis_physical_bound():
    args = [0.1, 1000, 0.0001]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.ekman_transport(*args)


def test_inertial_period_broadcasting():
    args: list[Any] = [0.0001]
    reference = domain.inertial_period(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.inertial_period(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_inertial_period_coriolis_physical_bound():
    args = [0.0001]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.inertial_period(*args)


def test_shallow_wave_speed_broadcasting():
    args: list[Any] = [100, 10]
    reference = domain.shallow_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.shallow_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_shallow_wave_speed_depth_physical_bound():
    args = [100, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.shallow_wave_speed(*args)


def test_shallow_wave_speed_gravity_physical_bound():
    args = [100, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.shallow_wave_speed(*args)


def test_deep_wave_speed_broadcasting():
    args: list[Any] = [6.283185307179586, 10]
    reference = domain.deep_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.deep_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_deep_wave_speed_wavelength_physical_bound():
    args = [6.283185307179586, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.deep_wave_speed(*args)


def test_deep_wave_speed_gravity_physical_bound():
    args = [6.283185307179586, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.deep_wave_speed(*args)


def test_wave_angular_frequency_broadcasting():
    args: list[Any] = [1, 1, 10]
    reference = domain.wave_angular_frequency(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.wave_angular_frequency(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_wave_angular_frequency_wavenumber_physical_bound():
    args = [1, 1, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.wave_angular_frequency(*args)


def test_wave_angular_frequency_depth_physical_bound():
    args = [1, 1, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.wave_angular_frequency(*args)


def test_wave_angular_frequency_gravity_physical_bound():
    args = [1, 1, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.wave_angular_frequency(*args)


def test_deep_group_speed_broadcasting():
    args: list[Any] = [6.283185307179586, 10]
    reference = domain.deep_group_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.deep_group_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_deep_group_speed_wavelength_physical_bound():
    args = [6.283185307179586, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.deep_group_speed(*args)


def test_deep_group_speed_gravity_physical_bound():
    args = [6.283185307179586, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.deep_group_speed(*args)


def test_rossby_radius_broadcasting():
    args: list[Any] = [0.01, 1000, 0.0001]
    reference = domain.rossby_radius(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.rossby_radius(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_rossby_radius_buoyancy_frequency_physical_bound():
    args = [0.01, 1000, 0.0001]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.rossby_radius(*args)


def test_rossby_radius_depth_physical_bound():
    args = [0.01, 1000, 0.0001]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.rossby_radius(*args)


def test_rossby_radius_coriolis_physical_bound():
    args = [0.01, 1000, 0.0001]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.rossby_radius(*args)


def test_reduced_gravity_broadcasting():
    args: list[Any] = [2, 1000, 10]
    reference = domain.reduced_gravity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.reduced_gravity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_reduced_gravity_density_difference_physical_bound():
    args = [2, 1000, 10]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.reduced_gravity(*args)


def test_reduced_gravity_reference_density_physical_bound():
    args = [2, 1000, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.reduced_gravity(*args)


def test_reduced_gravity_gravity_physical_bound():
    args = [2, 1000, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.reduced_gravity(*args)


def test_internal_wave_speed_broadcasting():
    args: list[Any] = [0.02, 100, 100]
    reference = domain.internal_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.internal_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_internal_wave_speed_reduced_gravity_physical_bound():
    args = [0.02, 100, 100]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.internal_wave_speed(*args)


def test_internal_wave_speed_upper_depth_physical_bound():
    args = [0.02, 100, 100]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.internal_wave_speed(*args)


def test_internal_wave_speed_lower_depth_physical_bound():
    args = [0.02, 100, 100]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.internal_wave_speed(*args)


def test_volume_transport_broadcasting():
    args: list[Any] = [1, 1000, 100]
    reference = domain.volume_transport(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.volume_transport(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_volume_transport_width_physical_bound():
    args = [1, 1000, 100]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.volume_transport(*args)


def test_volume_transport_depth_physical_bound():
    args = [1, 1000, 100]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.volume_transport(*args)


def test_tidal_elevation_broadcasting():
    args: list[Any] = [0, 2, 12, 0]
    reference = domain.tidal_elevation(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.tidal_elevation(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_tidal_elevation_amplitude_physical_bound():
    args = [0, 2, 12, 0]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.tidal_elevation(*args)


def test_tidal_elevation_period_physical_bound():
    args = [0, 2, 12, 0]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.tidal_elevation(*args)


def test_gerstner_orbit_radius_reference():
    result = domain.gerstner_orbit_radius(*[0.6931471805599453, 0.5, 1])
    np.testing.assert_allclose(result, 0.25, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_gerstner_orbit_radius_rejects_invalid(bad):
    args = [0.6931471805599453, 0.5, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.gerstner_orbit_radius(*args)


def test_gerstner_orbit_radius_broadcasting():
    args: list[Any] = [0.6931471805599453, 0.5, 1]
    reference = domain.gerstner_orbit_radius(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.gerstner_orbit_radius(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_gerstner_orbit_radius_rest_depth_physical_bound():
    args = [0.6931471805599453, 0.5, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.gerstner_orbit_radius(*args)


def test_gerstner_orbit_radius_amplitude_physical_bound():
    args = [0.6931471805599453, 0.5, 1]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.gerstner_orbit_radius(*args)


def test_gerstner_orbit_radius_wavenumber_physical_bound():
    args = [0.6931471805599453, 0.5, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.gerstner_orbit_radius(*args)


def test_sverdrup_transport_reference():
    result = domain.sverdrup_transport(*[6, 2, 3])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_sverdrup_transport_rejects_invalid(bad):
    args = [6, 2, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.sverdrup_transport(*args)


def test_sverdrup_transport_broadcasting():
    args: list[Any] = [6, 2, 3]
    reference = domain.sverdrup_transport(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.sverdrup_transport(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_sverdrup_transport_density_physical_bound():
    args = [6, 2, 3]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.sverdrup_transport(*args)


def test_sverdrup_transport_beta_physical_bound():
    args = [6, 2, 3]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.sverdrup_transport(*args)


def test_stommel_streamfunction_reference():
    result = domain.stommel_streamfunction(*[1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1])
    np.testing.assert_allclose(result, 1 - 1 / np.cosh(1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_stommel_streamfunction_rejects_invalid(bad):
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_broadcasting():
    args: list[Any] = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    reference = domain.stommel_streamfunction(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.stommel_streamfunction(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_stommel_streamfunction_x_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_y_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_basin_width_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_basin_length_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_density_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_beta_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[6] = -1
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_stommel_streamfunction_friction_physical_bound():
    args = [1, 0.5 * np.pi, 1, 2, np.pi, 1, 0, 1]
    args[7] = 0
    with pytest.raises(ValueError):
        domain.stommel_streamfunction(*args)


def test_thermohaline_forcing_reference():
    result = domain.thermohaline_forcing(*[-1, 2, 3])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_thermohaline_forcing_rejects_invalid(bad):
    args = [-1, 2, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.thermohaline_forcing(*args)


def test_thermohaline_forcing_broadcasting():
    args: list[Any] = [-1, 2, 3]
    reference = domain.thermohaline_forcing(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.thermohaline_forcing(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_thermohaline_forcing_thermal_flow_physical_bound():
    args = [-1, 2, 3]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.thermohaline_forcing(*args)


def test_thermohaline_forcing_haline_coefficient_physical_bound():
    args = [-1, 2, 3]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.thermohaline_forcing(*args)


def test_abyssal_temperature_reference():
    result = domain.abyssal_temperature(*[0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_abyssal_temperature_rejects_invalid(bad):
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_broadcasting():
    args: list[Any] = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    reference = domain.abyssal_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.abyssal_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_abyssal_temperature_height_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_bottom_temperature_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_top_temperature_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_column_height_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_diffusivity_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_abyssal_temperature_upwelling_physical_bound():
    args = [0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.abyssal_temperature(*args)


def test_simplified_seawater_density_reference():
    result = domain.simplified_seawater_density(*[293.15, 36, 0])
    np.testing.assert_allclose(result, 1026 - 1.655 * (1 + 0.2976) + 0.76554 * (1 - 0.5 * 7.4914e-4) - 0.024341, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_simplified_seawater_density_rejects_invalid(bad):
    args = [293.15, 36, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.simplified_seawater_density(*args)


def test_simplified_seawater_density_broadcasting():
    args: list[Any] = [293.15, 36, 0]
    reference = domain.simplified_seawater_density(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.simplified_seawater_density(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_simplified_seawater_density_conservative_temperature_physical_bound():
    args = [293.15, 36, 0]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.simplified_seawater_density(*args)


def test_simplified_seawater_density_absolute_salinity_physical_bound():
    args = [293.15, 36, 0]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.simplified_seawater_density(*args)


def test_simplified_seawater_density_depth_physical_bound():
    args = [293.15, 36, 0]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.simplified_seawater_density(*args)


def test_ekman_spiral_velocity_reference():
    result = domain.ekman_spiral_velocity(*[0, 1, 1, 0.5, 2])
    np.testing.assert_allclose(result, [np.sqrt(0.5), np.sqrt(0.5)], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ekman_spiral_velocity_rejects_invalid(bad):
    args = [0, 1, 1, 0.5, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ekman_spiral_velocity(*args)


def test_ekman_spiral_velocity_broadcasting():
    result = domain.ekman_spiral_velocity(np.zeros((2, 1)), 1, 1, 0.5, 2)
    assert result.shape == (2, 1, 2)
    np.testing.assert_allclose(result, np.sqrt(0.5))


@pytest.mark.parametrize("index,value", [(0, -1), (2, 0), (3, 0), (4, 0)])
def test_ekman_spiral_velocity_physical_bound(index, value):
    args = [0, 1, 1, 0.5, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.ekman_spiral_velocity(*args)


@pytest.mark.parametrize(
    "function,args",
    [
        (domain.gerstner_orbit_radius, [0, 2, 1]),
        (domain.stommel_streamfunction, [3, 1, 1, 2, np.pi, 1, 0, 1]),
        (domain.stommel_streamfunction, [1, 4, 1, 2, np.pi, 1, 0, 1]),
        (domain.thermohaline_forcing, [3, 2, 1]),
        (domain.abyssal_temperature, [2, 1, 4, 1, 1, 1]),
    ],
)
def test_new_model_cross_parameter_rejections(function, args):
    with pytest.raises(ValueError):
        function(*args)


def test_teos10_density_reference():
    result = domain.teos10_density(*[283.15, 30, 1e7])
    np.testing.assert_allclose(result, 1 / 9.732819627722660e-4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_teos10_density_rejects_invalid(bad):
    args = [283.15, 30, 1e7]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.teos10_density(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, -1)])
def test_teos10_density_physical_bound(index, value):
    args = [283.15, 30, 1e7]
    args[index] = value
    with pytest.raises(ValueError):
        domain.teos10_density(*args)
