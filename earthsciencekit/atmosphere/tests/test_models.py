"""Golden checks and invalid-input contracts for every public model."""

from typing import Any

import numpy as np
import pytest

from earthsciencekit import atmosphere as domain
from earthsciencekit.constants import EPSILON_WATER


def test_dry_air_density_reference():
    result = domain.dry_air_density(*[86115.0, 300])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_dry_air_density_rejects_invalid(bad):
    args = [86115.0, 300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dry_air_density(*args)


def test_scale_height_reference():
    result = domain.scale_height(*[300, 10])
    np.testing.assert_allclose(result, 8611.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_scale_height_rejects_invalid(bad):
    args = [300, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.scale_height(*args)


def test_barometric_pressure_reference():
    result = domain.barometric_pressure(*[100000.0, 0, 300, 10])
    np.testing.assert_allclose(result, 1e5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_barometric_pressure_rejects_invalid(bad):
    args = [100000.0, 0, 300, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.barometric_pressure(*args)


def test_potential_temperature_reference():
    result = domain.potential_temperature(*[300, 100000.0, 100000.0])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_potential_temperature_rejects_invalid(bad):
    args = [300, 100000.0, 100000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.potential_temperature(*args)


def test_temperature_from_potential_reference():
    result = domain.temperature_from_potential(*[300, 100000.0, 100000.0])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_temperature_from_potential_rejects_invalid(bad):
    args = [300, 100000.0, 100000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.temperature_from_potential(*args)


def test_dry_adiabatic_lapse_rate_reference():
    result = domain.dry_adiabatic_lapse_rate(*[10, 1000])
    np.testing.assert_allclose(result, 0.01, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_dry_adiabatic_lapse_rate_rejects_invalid(bad):
    args = [10, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dry_adiabatic_lapse_rate(*args)


def test_saturation_vapor_pressure_reference():
    result = domain.saturation_vapor_pressure(*[273.15])
    np.testing.assert_allclose(result, 611.2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_saturation_vapor_pressure_rejects_invalid(bad):
    args = [273.15]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.saturation_vapor_pressure(*args)


def test_mixing_ratio_reference():
    result = domain.mixing_ratio(*[1000, 100000])
    np.testing.assert_allclose(result, EPSILON_WATER / 99, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_mixing_ratio_rejects_invalid(bad):
    args = [1000, 100000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.mixing_ratio(*args)


def test_specific_humidity_reference():
    result = domain.specific_humidity(*[1])
    np.testing.assert_allclose(result, 0.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_specific_humidity_rejects_invalid(bad):
    args = [1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.specific_humidity(*args)


def test_virtual_temperature_reference():
    result = domain.virtual_temperature(*[300, 0])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_virtual_temperature_rejects_invalid(bad):
    args = [300, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.virtual_temperature(*args)


def test_relative_humidity_reference():
    result = domain.relative_humidity(*[500, 1000])
    np.testing.assert_allclose(result, 0.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_relative_humidity_rejects_invalid(bad):
    args = [500, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.relative_humidity(*args)


def test_coriolis_parameter_reference():
    result = domain.coriolis_parameter(*[0])
    np.testing.assert_allclose(result, 0, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_coriolis_parameter_rejects_invalid(bad):
    args = [0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.coriolis_parameter(*args)


def test_brunt_vaisala_squared_reference():
    result = domain.brunt_vaisala_squared(*[300, 0.003, 10])
    np.testing.assert_allclose(result, 0.0001, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_brunt_vaisala_squared_rejects_invalid(bad):
    args = [300, 0.003, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.brunt_vaisala_squared(*args)


def test_geostrophic_speed_reference():
    result = domain.geostrophic_speed(*[0.001, 1, 0.0001])
    np.testing.assert_allclose(result, 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_geostrophic_speed_rejects_invalid(bad):
    args = [0.001, 1, 0.0001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.geostrophic_speed(*args)


def test_rossby_number_reference():
    result = domain.rossby_number(*[10, 1000000.0, 0.0001])
    np.testing.assert_allclose(result, 0.1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_rossby_number_rejects_invalid(bad):
    args = [10, 1000000.0, 0.0001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rossby_number(*args)


def test_dry_air_density_broadcasting():
    args: list[Any] = [86115.0, 300]
    reference = domain.dry_air_density(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.dry_air_density(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_dry_air_density_pressure_physical_bound():
    args = [86115.0, 300]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.dry_air_density(*args)


def test_dry_air_density_temperature_physical_bound():
    args = [86115.0, 300]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.dry_air_density(*args)


def test_scale_height_broadcasting():
    args: list[Any] = [300, 10]
    reference = domain.scale_height(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.scale_height(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_scale_height_temperature_physical_bound():
    args = [300, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.scale_height(*args)


def test_scale_height_gravity_physical_bound():
    args = [300, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.scale_height(*args)


def test_barometric_pressure_broadcasting():
    args: list[Any] = [100000.0, 0, 300, 10]
    reference = domain.barometric_pressure(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.barometric_pressure(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_barometric_pressure_surface_pressure_physical_bound():
    args = [100000.0, 0, 300, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.barometric_pressure(*args)


def test_barometric_pressure_temperature_physical_bound():
    args = [100000.0, 0, 300, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.barometric_pressure(*args)


def test_barometric_pressure_gravity_physical_bound():
    args = [100000.0, 0, 300, 10]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.barometric_pressure(*args)


def test_potential_temperature_broadcasting():
    args: list[Any] = [300, 100000.0, 100000.0]
    reference = domain.potential_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.potential_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_potential_temperature_temperature_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.potential_temperature(*args)


def test_potential_temperature_pressure_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.potential_temperature(*args)


def test_potential_temperature_reference_pressure_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.potential_temperature(*args)


def test_temperature_from_potential_broadcasting():
    args: list[Any] = [300, 100000.0, 100000.0]
    reference = domain.temperature_from_potential(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.temperature_from_potential(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_temperature_from_potential_potential_temperature_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.temperature_from_potential(*args)


def test_temperature_from_potential_pressure_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.temperature_from_potential(*args)


def test_temperature_from_potential_reference_pressure_physical_bound():
    args = [300, 100000.0, 100000.0]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.temperature_from_potential(*args)


def test_dry_adiabatic_lapse_rate_broadcasting():
    args: list[Any] = [10, 1000]
    reference = domain.dry_adiabatic_lapse_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.dry_adiabatic_lapse_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_dry_adiabatic_lapse_rate_gravity_physical_bound():
    args = [10, 1000]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.dry_adiabatic_lapse_rate(*args)


def test_dry_adiabatic_lapse_rate_heat_capacity_physical_bound():
    args = [10, 1000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.dry_adiabatic_lapse_rate(*args)


def test_saturation_vapor_pressure_broadcasting():
    args: list[Any] = [273.15]
    reference = domain.saturation_vapor_pressure(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.saturation_vapor_pressure(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_saturation_vapor_pressure_temperature_physical_bound():
    args = [273.15]
    args[0] = 100
    with pytest.raises(ValueError):
        domain.saturation_vapor_pressure(*args)


def test_mixing_ratio_broadcasting():
    args: list[Any] = [1000, 100000]
    reference = domain.mixing_ratio(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.mixing_ratio(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_mixing_ratio_vapor_pressure_physical_bound():
    args = [1000, 100000]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.mixing_ratio(*args)


def test_mixing_ratio_total_pressure_physical_bound():
    args = [1000, 100000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.mixing_ratio(*args)


def test_specific_humidity_broadcasting():
    args: list[Any] = [1]
    reference = domain.specific_humidity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.specific_humidity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_specific_humidity_mixing_ratio_physical_bound():
    args = [1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.specific_humidity(*args)


def test_virtual_temperature_broadcasting():
    args: list[Any] = [300, 0]
    reference = domain.virtual_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.virtual_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_virtual_temperature_temperature_physical_bound():
    args = [300, 0]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.virtual_temperature(*args)


def test_virtual_temperature_specific_humidity_physical_bound():
    args: list[Any] = [300, 0]
    args[1] = 1.1
    with pytest.raises(ValueError):
        domain.virtual_temperature(*args)


def test_relative_humidity_broadcasting():
    args: list[Any] = [500, 1000]
    reference = domain.relative_humidity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.relative_humidity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_relative_humidity_vapor_pressure_physical_bound():
    args = [500, 1000]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.relative_humidity(*args)


def test_relative_humidity_saturation_pressure_physical_bound():
    args = [500, 1000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.relative_humidity(*args)


def test_coriolis_parameter_broadcasting():
    args: list[Any] = [0]
    reference = domain.coriolis_parameter(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.coriolis_parameter(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_coriolis_parameter_latitude_physical_bound():
    args = [0]
    args[0] = 2
    with pytest.raises(ValueError):
        domain.coriolis_parameter(*args)


def test_brunt_vaisala_squared_broadcasting():
    args: list[Any] = [300, 0.003, 10]
    reference = domain.brunt_vaisala_squared(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.brunt_vaisala_squared(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_brunt_vaisala_squared_potential_temperature_physical_bound():
    args = [300, 0.003, 10]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.brunt_vaisala_squared(*args)


def test_brunt_vaisala_squared_gravity_physical_bound():
    args = [300, 0.003, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.brunt_vaisala_squared(*args)


def test_geostrophic_speed_broadcasting():
    args: list[Any] = [0.001, 1, 0.0001]
    reference = domain.geostrophic_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.geostrophic_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_geostrophic_speed_density_physical_bound():
    args = [0.001, 1, 0.0001]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.geostrophic_speed(*args)


def test_geostrophic_speed_coriolis_physical_bound():
    args = [0.001, 1, 0.0001]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.geostrophic_speed(*args)


def test_rossby_number_broadcasting():
    args: list[Any] = [10, 1000000.0, 0.0001]
    reference = domain.rossby_number(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.rossby_number(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_rossby_number_speed_physical_bound():
    args = [10, 1000000.0, 0.0001]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.rossby_number(*args)


def test_rossby_number_length_physical_bound():
    args = [10, 1000000.0, 0.0001]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.rossby_number(*args)


def test_rossby_number_coriolis_physical_bound():
    args = [10, 1000000.0, 0.0001]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.rossby_number(*args)


def test_barometer_column_height_reference():
    result = domain.barometer_column_height(*[12, 2, 3])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_barometer_column_height_rejects_invalid(bad):
    args = [12, 2, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.barometer_column_height(*args)


def test_barometer_column_height_broadcasting():
    args: list[Any] = [12, 2, 3]
    reference = domain.barometer_column_height(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.barometer_column_height(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_barometer_column_height_pressure_physical_bound():
    args = [12, 2, 3]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.barometer_column_height(*args)


def test_barometer_column_height_liquid_density_physical_bound():
    args = [12, 2, 3]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.barometer_column_height(*args)


def test_barometer_column_height_gravity_physical_bound():
    args = [12, 2, 3]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.barometer_column_height(*args)


def test_two_layer_temperature_reference():
    result = domain.two_layer_temperature(*[5000, 288, 0.006, 10000])
    np.testing.assert_allclose(result, 258, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_two_layer_temperature_rejects_invalid(bad):
    args = [5000, 288, 0.006, 10000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.two_layer_temperature(*args)


def test_two_layer_temperature_broadcasting():
    args: list[Any] = [5000, 288, 0.006, 10000]
    reference = domain.two_layer_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.two_layer_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_two_layer_temperature_height_physical_bound():
    args = [5000, 288, 0.006, 10000]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.two_layer_temperature(*args)


def test_two_layer_temperature_surface_temperature_physical_bound():
    args = [5000, 288, 0.006, 10000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.two_layer_temperature(*args)


def test_two_layer_temperature_lapse_rate_physical_bound():
    args = [5000, 288, 0.006, 10000]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.two_layer_temperature(*args)


def test_two_layer_temperature_tropopause_height_physical_bound():
    args = [5000, 288, 0.006, 10000]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.two_layer_temperature(*args)


def test_richardson_number_reference():
    result = domain.richardson_number(*[0.0001, 0.02])
    np.testing.assert_allclose(result, 0.25, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_richardson_number_rejects_invalid(bad):
    args = [0.0001, 0.02]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.richardson_number(*args)


def test_richardson_number_broadcasting():
    args: list[Any] = [0.0001, 0.02]
    reference = domain.richardson_number(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.richardson_number(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_richardson_number_shear_physical_bound():
    args = [0.0001, 0.02]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.richardson_number(*args)


def test_surface_pressure_tendency_reference():
    result = domain.surface_pressure_tendency(*[1e-05, 100000.0])
    np.testing.assert_allclose(result, -1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_surface_pressure_tendency_rejects_invalid(bad):
    args = [1e-05, 100000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.surface_pressure_tendency(*args)


def test_surface_pressure_tendency_broadcasting():
    args: list[Any] = [1e-05, 100000.0]
    reference = domain.surface_pressure_tendency(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.surface_pressure_tendency(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_surface_pressure_tendency_surface_pressure_physical_bound():
    args = [1e-05, 100000.0]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.surface_pressure_tendency(*args)


def test_rossby_wave_speed_reference():
    result = domain.rossby_wave_speed(*[10, 2e-11, 1e-06])
    np.testing.assert_allclose(result, -10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_rossby_wave_speed_rejects_invalid(bad):
    args = [10, 2e-11, 1e-06]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rossby_wave_speed(*args)


def test_rossby_wave_speed_broadcasting():
    args: list[Any] = [10, 2e-11, 1e-06]
    reference = domain.rossby_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.rossby_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_rossby_wave_speed_beta_physical_bound():
    args = [10, 2e-11, 1e-06]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.rossby_wave_speed(*args)


def test_rossby_wave_speed_wavenumber_physical_bound():
    args = [10, 2e-11, 1e-06]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.rossby_wave_speed(*args)


_CONVECTION_CASES: dict[str, tuple[list[Any], Any]] = {
    "lcl_temperature": ([300, 290], 1 / (1 / 234 + np.log(300 / 290) / 800) + 56),
    "lcl_pressure": ([1e5, 300, 270], 1e5 * 0.9 ** (1004 / 287.05)),
    "lcl_height": ([300, 290], 1004 * 10 / 9.80665),
    "moist_adiabatic_lapse_rate": (
        [300, 1e5],
        9.80665
        * (1 + 2.5e6 * (rs := EPSILON_WATER * (es := 611.2 * np.exp(2.5e6 / 461.5 * (1 / 273.15 - 1 / 300))) / (1e5 - es)) / (287.05 * 300))
        / (1004 + 2.5e6**2 * rs * EPSILON_WATER / (287.05 * 300**2)),
    ),
    "parcel_temperature_profile": ([[1e5, 9e4], 300, 250], [300, 300 * 0.9 ** (287.05 / 1004)]),
    "convective_available_potential_energy": ([[0, 1000, 2000], [300, 302, 300], [300, 300, 300]], 9.80665 * 20 / 3),
    "convective_inhibition": ([[0, 1000, 2000], [300, 299, 301], [300, 300, 300]], -2.5 * 9.80665),
}


@pytest.mark.parametrize("name", sorted(_CONVECTION_CASES))
def test_convection_reference(name):
    args, expected = _CONVECTION_CASES[name]
    np.testing.assert_allclose(getattr(domain, name)(*args), expected, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
@pytest.mark.parametrize("name", sorted(_CONVECTION_CASES))
def test_convection_rejects_invalid(name, bad):
    args = list(_CONVECTION_CASES[name][0])
    args[0] = bad
    with pytest.raises(ValueError):
        getattr(domain, name)(*args)


@pytest.mark.parametrize(
    "name,args",
    [
        ("lcl_temperature", [290, 300]),
        ("lcl_temperature", [300, 50]),
        ("lcl_pressure", [1e5, 300, 310]),
        ("lcl_pressure", [0, 300, 290]),
        ("lcl_height", [300, 310]),
        ("moist_adiabatic_lapse_rate", [350, 1e5]),
        ("moist_adiabatic_lapse_rate", [300, 3000]),
        ("parcel_temperature_profile", [[1e5, 1e5], 300, 290]),
        ("parcel_temperature_profile", [[1e5], 300, 290]),
        ("parcel_temperature_profile", [[1e5, 9e4], [300, 301], 290]),
        ("parcel_temperature_profile", [[1e5, 9e4], 300, 310]),
        ("convective_available_potential_energy", [[0, 0, 1], [300, 300, 300], [300, 300, 300]]),
        ("convective_available_potential_energy", [[0, 1], [300, 300, 300], [300, 300, 300]]),
        ("convective_available_potential_energy", [[0], [300], [300]]),
        ("convective_inhibition", [[0, 1, 2], [300, 0, 300], [300, 300, 300]]),
    ],
)
def test_convection_physical_bound(name, args):
    with pytest.raises(ValueError):
        getattr(domain, name)(*args)
