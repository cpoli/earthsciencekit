"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import climate as domain


def test_absorbed_solar_reference():
    result = domain.absorbed_solar(*[1360, 0.3])
    np.testing.assert_allclose(result, 238, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_absorbed_solar_rejects_invalid(bad):
    args = [1360, 0.3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.absorbed_solar(*args)


def test_effective_temperature_reference():
    result = domain.effective_temperature(*[459.300327939, 1])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_effective_temperature_rejects_invalid(bad):
    args = [459.300327939, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.effective_temperature(*args)


def test_outgoing_longwave_reference():
    result = domain.outgoing_longwave(*[300, 1])
    np.testing.assert_allclose(result, 459.300327939, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_outgoing_longwave_rejects_invalid(bad):
    args = [300, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.outgoing_longwave(*args)


def test_planck_feedback_reference():
    result = domain.planck_feedback(*[300, 1])
    np.testing.assert_allclose(result, 6.12400437252, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_planck_feedback_rejects_invalid(bad):
    args = [300, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.planck_feedback(*args)


def test_co2_forcing_reference():
    result = domain.co2_forcing(*[560, 280, 5.35])
    np.testing.assert_allclose(result, 5.35 * np.log(2), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_co2_forcing_rejects_invalid(bad):
    args = [560, 280, 5.35]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.co2_forcing(*args)


def test_equilibrium_warming_reference():
    result = domain.equilibrium_warming(*[4, 2])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_equilibrium_warming_rejects_invalid(bad):
    args = [4, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.equilibrium_warming(*args)


def test_heat_capacity_column_reference():
    result = domain.heat_capacity_column(*[1000, 4000, 50])
    np.testing.assert_allclose(result, 2e8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_heat_capacity_column_rejects_invalid(bad):
    args = [1000, 4000, 50]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.heat_capacity_column(*args)


def test_response_timescale_reference():
    result = domain.response_timescale(*[200000000.0, 2])
    np.testing.assert_allclose(result, 1e8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_response_timescale_rejects_invalid(bad):
    args = [200000000.0, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.response_timescale(*args)


def test_step_temperature_reference():
    result = domain.step_temperature(*[4, 2, 2, 0])
    np.testing.assert_allclose(result, 0, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_step_temperature_rejects_invalid(bad):
    args = [4, 2, 2, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.step_temperature(*args)


def test_heat_uptake_reference():
    result = domain.heat_uptake(*[200000000.0, 2])
    np.testing.assert_allclose(result, 4e8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_heat_uptake_rejects_invalid(bad):
    args = [200000000.0, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.heat_uptake(*args)


def test_ice_albedo_reference():
    result = domain.ice_albedo(*[273, 0.6, 0.2, 273, 5])
    np.testing.assert_allclose(result, 0.4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ice_albedo_rejects_invalid(bad):
    args = [273, 0.6, 0.2, 273, 5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_carbon_airborne_fraction_reference():
    result = domain.carbon_airborne_fraction(*[4, 10])
    np.testing.assert_allclose(result, 0.4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_carbon_airborne_fraction_rejects_invalid(bad):
    args = [4, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.carbon_airborne_fraction(*args)


def test_carbon_decay_reference():
    result = domain.carbon_decay(*[2, 1, np.float64(0.6931471805599453)])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_carbon_decay_rejects_invalid(bad):
    args = [2, 1, np.float64(0.6931471805599453)]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.carbon_decay(*args)


def test_thermal_expansion_sea_level_reference():
    result = domain.thermal_expansion_sea_level(*[0.0002, 1000, 1])
    np.testing.assert_allclose(result, 0.2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_thermal_expansion_sea_level_rejects_invalid(bad):
    args = [0.0002, 1000, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.thermal_expansion_sea_level(*args)


def test_water_equivalent_sea_level_reference():
    result = domain.water_equivalent_sea_level(*[1000000000000.0, 1000000000000.0, 1000])
    np.testing.assert_allclose(result, 0.001, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_water_equivalent_sea_level_rejects_invalid(bad):
    args = [1000000000000.0, 1000000000000.0, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.water_equivalent_sea_level(*args)


def test_absorbed_solar_broadcasting():
    args = [1360, 0.3]
    reference = domain.absorbed_solar(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.absorbed_solar(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_absorbed_solar_solar_constant_physical_bound():
    args = [1360, 0.3]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.absorbed_solar(*args)


def test_absorbed_solar_albedo_physical_bound():
    args = [1360, 0.3]
    args[1] = 1.1
    with pytest.raises(ValueError):
        domain.absorbed_solar(*args)


def test_effective_temperature_broadcasting():
    args = [459.300327939, 1]
    reference = domain.effective_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.effective_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_effective_temperature_absorbed_flux_physical_bound():
    args = [459.300327939, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.effective_temperature(*args)


def test_effective_temperature_emissivity_physical_bound():
    args = [459.300327939, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.effective_temperature(*args)


def test_outgoing_longwave_broadcasting():
    args = [300, 1]
    reference = domain.outgoing_longwave(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.outgoing_longwave(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_outgoing_longwave_temperature_physical_bound():
    args = [300, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.outgoing_longwave(*args)


def test_outgoing_longwave_emissivity_physical_bound():
    args = [300, 1]
    args[1] = 1.1
    with pytest.raises(ValueError):
        domain.outgoing_longwave(*args)


def test_planck_feedback_broadcasting():
    args = [300, 1]
    reference = domain.planck_feedback(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.planck_feedback(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_planck_feedback_temperature_physical_bound():
    args = [300, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.planck_feedback(*args)


def test_planck_feedback_emissivity_physical_bound():
    args = [300, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.planck_feedback(*args)


def test_co2_forcing_broadcasting():
    args = [560, 280, 5.35]
    reference = domain.co2_forcing(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.co2_forcing(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_co2_forcing_concentration_physical_bound():
    args = [560, 280, 5.35]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.co2_forcing(*args)


def test_co2_forcing_reference_concentration_physical_bound():
    args = [560, 280, 5.35]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.co2_forcing(*args)


def test_co2_forcing_coefficient_physical_bound():
    args = [560, 280, 5.35]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.co2_forcing(*args)


def test_equilibrium_warming_broadcasting():
    args = [4, 2]
    reference = domain.equilibrium_warming(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.equilibrium_warming(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_equilibrium_warming_restoring_physical_bound():
    args = [4, 2]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.equilibrium_warming(*args)


def test_heat_capacity_column_broadcasting():
    args = [1000, 4000, 50]
    reference = domain.heat_capacity_column(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.heat_capacity_column(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_heat_capacity_column_density_physical_bound():
    args = [1000, 4000, 50]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.heat_capacity_column(*args)


def test_heat_capacity_column_heat_capacity_physical_bound():
    args = [1000, 4000, 50]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.heat_capacity_column(*args)


def test_heat_capacity_column_depth_physical_bound():
    args = [1000, 4000, 50]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.heat_capacity_column(*args)


def test_response_timescale_broadcasting():
    args = [200000000.0, 2]
    reference = domain.response_timescale(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.response_timescale(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_response_timescale_areal_heat_capacity_physical_bound():
    args = [200000000.0, 2]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.response_timescale(*args)


def test_response_timescale_restoring_physical_bound():
    args = [200000000.0, 2]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.response_timescale(*args)


def test_step_temperature_broadcasting():
    args = [4, 2, 2, 0]
    reference = domain.step_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.step_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_step_temperature_restoring_physical_bound():
    args = [4, 2, 2, 0]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.step_temperature(*args)


def test_step_temperature_areal_heat_capacity_physical_bound():
    args = [4, 2, 2, 0]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.step_temperature(*args)


def test_step_temperature_time_physical_bound():
    args = [4, 2, 2, 0]
    args[3] = -1
    with pytest.raises(ValueError):
        domain.step_temperature(*args)


def test_heat_uptake_broadcasting():
    args = [200000000.0, 2]
    reference = domain.heat_uptake(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.heat_uptake(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_heat_uptake_areal_heat_capacity_physical_bound():
    args = [200000000.0, 2]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.heat_uptake(*args)


def test_ice_albedo_broadcasting():
    args = [273, 0.6, 0.2, 273, 5]
    reference = domain.ice_albedo(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.ice_albedo(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_ice_albedo_temperature_physical_bound():
    args = [273, 0.6, 0.2, 273, 5]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_ice_albedo_cold_albedo_physical_bound():
    args = [273, 0.6, 0.2, 273, 5]
    args[1] = 1.1
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_ice_albedo_warm_albedo_physical_bound():
    args = [273, 0.6, 0.2, 273, 5]
    args[2] = 1.1
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_ice_albedo_transition_temperature_physical_bound():
    args = [273, 0.6, 0.2, 273, 5]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_ice_albedo_width_physical_bound():
    args = [273, 0.6, 0.2, 273, 5]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.ice_albedo(*args)


def test_carbon_airborne_fraction_broadcasting():
    args = [4, 10]
    reference = domain.carbon_airborne_fraction(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.carbon_airborne_fraction(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_carbon_airborne_fraction_atmospheric_increase_physical_bound():
    args = [4, 10]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.carbon_airborne_fraction(*args)


def test_carbon_airborne_fraction_emissions_physical_bound():
    args = [4, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.carbon_airborne_fraction(*args)


def test_carbon_decay_broadcasting():
    args = [2, 1, 0.6931471805599453]
    reference = domain.carbon_decay(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.carbon_decay(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_carbon_decay_initial_excess_physical_bound():
    args = [2, 1, 0.6931471805599453]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.carbon_decay(*args)


def test_carbon_decay_timescale_physical_bound():
    args = [2, 1, 0.6931471805599453]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.carbon_decay(*args)


def test_carbon_decay_time_physical_bound():
    args = [2, 1, 0.6931471805599453]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.carbon_decay(*args)


def test_thermal_expansion_sea_level_broadcasting():
    args = [0.0002, 1000, 1]
    reference = domain.thermal_expansion_sea_level(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.thermal_expansion_sea_level(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_thermal_expansion_sea_level_expansion_coefficient_physical_bound():
    args = [0.0002, 1000, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.thermal_expansion_sea_level(*args)


def test_thermal_expansion_sea_level_depth_physical_bound():
    args = [0.0002, 1000, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.thermal_expansion_sea_level(*args)


def test_water_equivalent_sea_level_broadcasting():
    args = [1000000000000.0, 1000000000000.0, 1000]
    reference = domain.water_equivalent_sea_level(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.water_equivalent_sea_level(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_water_equivalent_sea_level_water_mass_physical_bound():
    args = [1000000000000.0, 1000000000000.0, 1000]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.water_equivalent_sea_level(*args)


def test_water_equivalent_sea_level_ocean_area_physical_bound():
    args = [1000000000000.0, 1000000000000.0, 1000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.water_equivalent_sea_level(*args)


def test_water_equivalent_sea_level_water_density_physical_bound():
    args = [1000000000000.0, 1000000000000.0, 1000]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.water_equivalent_sea_level(*args)


def test_greenhouse_surface_temperature_reference():
    result = domain.greenhouse_surface_temperature(*[459.300327939, 0])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_greenhouse_surface_temperature_rejects_invalid(bad):
    args = [459.300327939, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.greenhouse_surface_temperature(*args)


def test_greenhouse_surface_temperature_broadcasting():
    args = [459.300327939, 0]
    reference = domain.greenhouse_surface_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.greenhouse_surface_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_greenhouse_surface_temperature_absorbed_flux_physical_bound():
    args = [459.300327939, 0]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.greenhouse_surface_temperature(*args)


def test_greenhouse_surface_temperature_layer_emissivity_physical_bound():
    args = [459.300327939, 0]
    args[1] = 1.1
    with pytest.raises(ValueError):
        domain.greenhouse_surface_temperature(*args)


def test_daily_insolation_reference():
    result = domain.daily_insolation(*[1361, 0, 0, 1])
    np.testing.assert_allclose(result, 1361 / np.pi, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_daily_insolation_rejects_invalid(bad):
    args = [1361, 0, 0, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_daily_insolation_broadcasting():
    args = [1361, 0, 0, 1]
    reference = domain.daily_insolation(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.daily_insolation(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_daily_insolation_solar_constant_physical_bound():
    args = [1361, 0, 0, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_daily_insolation_latitude_physical_bound():
    args = [1361, 0, 0, 1]
    args[1] = 2
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_daily_insolation_declination_physical_bound():
    args = [1361, 0, 0, 1]
    args[2] = 2
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_daily_insolation_distance_ratio_physical_bound():
    args = [1361, 0, 0, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_radiative_convective_temperature_reference():
    result = domain.radiative_convective_temperature(*[0, 255, 0.0065, 5000])
    np.testing.assert_allclose(result, 287.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_radiative_convective_temperature_rejects_invalid(bad):
    args = [0, 255, 0.0065, 5000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.radiative_convective_temperature(*args)


def test_radiative_convective_temperature_broadcasting():
    args = [0, 255, 0.0065, 5000]
    reference = domain.radiative_convective_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.radiative_convective_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_radiative_convective_temperature_height_physical_bound():
    args = [0, 255, 0.0065, 5000]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.radiative_convective_temperature(*args)


def test_radiative_convective_temperature_effective_temperature_physical_bound():
    args = [0, 255, 0.0065, 5000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.radiative_convective_temperature(*args)


def test_radiative_convective_temperature_lapse_rate_physical_bound():
    args = [0, 255, 0.0065, 5000]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.radiative_convective_temperature(*args)


def test_radiative_convective_temperature_emission_height_physical_bound():
    args = [0, 255, 0.0065, 5000]
    args[3] = -1
    with pytest.raises(ValueError):
        domain.radiative_convective_temperature(*args)


def test_equilibrium_solar_constant_reference():
    result = domain.equilibrium_solar_constant(*[288, 1, 0.6, 0.2, 288, 5])
    np.testing.assert_allclose(result, 4 * 5.670374419e-8 * 288**4 / 0.6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_equilibrium_solar_constant_rejects_invalid(bad):
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_broadcasting():
    args = [288, 1, 0.6, 0.2, 288, 5]
    reference = domain.equilibrium_solar_constant(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.equilibrium_solar_constant(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_equilibrium_solar_constant_temperature_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_emissivity_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_cold_albedo_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[2] = 1
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_warm_albedo_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[3] = -0.1
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_transition_temperature_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_equilibrium_solar_constant_width_physical_bound():
    args = [288, 1, 0.6, 0.2, 288, 5]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(*args)


def test_red_noise_spectrum_reference():
    result = domain.red_noise_spectrum(*[0, 4, 2, 1])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_red_noise_spectrum_rejects_invalid(bad):
    args = [0, 4, 2, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.red_noise_spectrum(*args)


def test_red_noise_spectrum_broadcasting():
    args = [0, 4, 2, 1]
    reference = domain.red_noise_spectrum(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.red_noise_spectrum(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_red_noise_spectrum_forcing_spectrum_physical_bound():
    args = [0, 4, 2, 1]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.red_noise_spectrum(*args)


def test_red_noise_spectrum_restoring_physical_bound():
    args = [0, 4, 2, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.red_noise_spectrum(*args)


def test_red_noise_spectrum_areal_heat_capacity_physical_bound():
    args = [0, 4, 2, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.red_noise_spectrum(*args)


def test_methane_forcing_reference():
    result = domain.methane_forcing(*[900, 400, 1e-16])
    np.testing.assert_allclose(result, 0.36, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_methane_forcing_rejects_invalid(bad):
    args = [900, 400, 1e-16]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.methane_forcing(*args)


def test_methane_forcing_broadcasting():
    args = [900, 400, 1e-16]
    reference = domain.methane_forcing(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.methane_forcing(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_methane_forcing_concentration_physical_bound():
    args = [900, 400, 1e-16]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.methane_forcing(*args)


def test_methane_forcing_reference_concentration_physical_bound():
    args = [900, 400, 1e-16]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.methane_forcing(*args)


def test_methane_forcing_nitrous_oxide_physical_bound():
    args = [900, 400, 1e-16]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.methane_forcing(*args)


def test_nitrous_oxide_forcing_reference():
    result = domain.nitrous_oxide_forcing(*[400, 100, 1e-16])
    np.testing.assert_allclose(result, 1.2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_nitrous_oxide_forcing_rejects_invalid(bad):
    args = [400, 100, 1e-16]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nitrous_oxide_forcing(*args)


def test_nitrous_oxide_forcing_broadcasting():
    args = [400, 100, 1e-16]
    reference = domain.nitrous_oxide_forcing(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.nitrous_oxide_forcing(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_nitrous_oxide_forcing_concentration_physical_bound():
    args = [400, 100, 1e-16]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.nitrous_oxide_forcing(*args)


def test_nitrous_oxide_forcing_reference_concentration_physical_bound():
    args = [400, 100, 1e-16]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.nitrous_oxide_forcing(*args)


def test_nitrous_oxide_forcing_methane_physical_bound():
    args = [400, 100, 1e-16]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.nitrous_oxide_forcing(*args)


def test_gregory_regression_reference():
    result = domain.gregory_regression([0, 1, 2], [4, 2.8, 1.6])
    np.testing.assert_allclose(result, [4, 1.2], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_gregory_regression_rejects_invalid(bad):
    with pytest.raises(ValueError):
        domain.gregory_regression(bad, [4, 2.8, 1.6])


@pytest.mark.parametrize(
    "temperature, imbalance",
    [(1.0, 1.0), ([1], [1]), ([0, 1], [1, 2, 3]), ([[0, 1]], [[1, 2]]), ([2, 2, 2], [1, 2, 3])],
)
def test_gregory_regression_rejects_bad_shapes(temperature, imbalance):
    with pytest.raises(ValueError):
        domain.gregory_regression(temperature, imbalance)


def test_greenhouse_layer_balances_top_of_atmosphere():
    absorbed, emissivity = 240.0, np.linspace(0, 1, 11)
    surface = domain.greenhouse_surface_temperature(absorbed, emissivity)
    layer = surface / 2**0.25
    outgoing = (1 - emissivity) * domain.outgoing_longwave(surface, 1) + domain.outgoing_longwave(layer, emissivity)
    np.testing.assert_allclose(outgoing, absorbed)
    np.testing.assert_allclose(surface[-1], 2**0.25 * domain.effective_temperature(absorbed, 1))


@pytest.mark.parametrize("declination", [0.0, 0.2, -0.41])
def test_daily_insolation_global_mean_is_quarter_solar_constant(declination):
    latitude = np.linspace(-np.pi / 2, np.pi / 2, 40001)
    flux = domain.daily_insolation(1361, latitude, declination, 1)
    assert np.isclose(np.trapezoid(flux * np.cos(latitude), latitude) / 2, 1361 / 4, rtol=1e-7)


def test_daily_insolation_polar_day_night_and_distance():
    assert np.isclose(domain.daily_insolation(1361, np.pi / 2, 0.4, 1), 1361 * np.sin(0.4))
    assert domain.daily_insolation(1361, -np.pi / 2, 0.4, 1) == 0
    assert np.isclose(domain.daily_insolation(1361, 0.3, 0.1, 1.02), 1.02**2 * domain.daily_insolation(1361, 0.3, 0.1, 1))


@pytest.mark.parametrize("position, value", [(1, 2), (2, -2)])
def test_daily_insolation_angle_upper_bounds(position, value):
    args = [1361, 0, 0, 1]
    args[position] = value
    with pytest.raises(ValueError):
        domain.daily_insolation(*args)


def test_radiative_convective_profile_structure():
    height = np.linspace(0, 30e3, 3001)
    profile = domain.radiative_convective_temperature(height, 255, 6.5e-3, 5e3)
    assert np.isclose(domain.radiative_convective_temperature(5e3, 255, 6.5e-3, 5e3), 255)
    assert np.isclose(profile[0], 255 + 6.5e-3 * 5e3)
    np.testing.assert_allclose(profile[-1], 255 / 2**0.25)
    tropopause = (profile[0] - 255 / 2**0.25) / 6.5e-3
    np.testing.assert_allclose(profile[height < tropopause], profile[0] - 6.5e-3 * height[height < tropopause])


def test_equilibrium_solar_constant_closes_energy_balance():
    temperature = np.linspace(220, 320, 101)
    solar = domain.equilibrium_solar_constant(temperature, 0.62, 0.62, 0.3, 265, 8)
    albedo = domain.ice_albedo(temperature, 0.62, 0.3, 265, 8)
    np.testing.assert_allclose(domain.absorbed_solar(solar, albedo), domain.outgoing_longwave(temperature, 0.62))
    assert np.any(np.diff(solar) < 0), "the ice-albedo feedback must create an unstable branch"


@pytest.mark.parametrize("cold, warm", [(1, 0.3), (0.2, 0.3)])
def test_equilibrium_solar_constant_rejects_albedos(cold, warm):
    with pytest.raises(ValueError):
        domain.equilibrium_solar_constant(288, 1, cold, warm, 288, 5)


def test_red_noise_spectrum_integrates_to_variance():
    from scipy.integrate import quad

    forcing, restoring, capacity = 3.0, 1.2, 2e8
    corner = restoring / (2 * np.pi * capacity)
    integral, _ = quad(lambda f: domain.red_noise_spectrum(f, forcing, restoring, capacity), -np.inf, np.inf, points=None, limit=200)
    assert np.isclose(integral, forcing / (2 * restoring * capacity), rtol=1e-8)
    assert np.isclose(domain.red_noise_spectrum(corner, forcing, restoring, capacity), forcing / (2 * restoring**2))


def test_myhre_forcings_match_ipcc_tar():
    assert np.isclose(domain.methane_forcing(1745, 700, 270), 0.48, atol=0.005)
    assert np.isclose(domain.nitrous_oxide_forcing(314, 270, 700), 0.15, atol=0.005)
    assert domain.methane_forcing(1745, 1745, 314) == 0
    assert domain.nitrous_oxide_forcing(314, 314, 1745) == 0


def test_gregory_regression_recovers_one_box_parameters():
    time = np.linspace(0, 150, 151) * 365.25 * 86400
    temperature = domain.simulate_one_box(7.4, 1.1, 3e8, time).values[:, 0]
    forcing, restoring = domain.gregory_regression(temperature, 7.4 - 1.1 * temperature)
    np.testing.assert_allclose([forcing, restoring], [7.4, 1.1])
