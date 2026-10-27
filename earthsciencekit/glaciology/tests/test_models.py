"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import glaciology as domain
from earthsciencekit.constants import R_GAS, SECONDS_PER_YEAR

BAD = [np.nan, np.inf, -np.inf, [], "bad", True, 1j]
# Herron-Langway stage slopes [1/m] and critical depth [m] at 248.15 K, 1e-9 m/s of ice and 350 kg/m^3 at the surface
HL_UPPER = 0.917 * 11 * np.exp(-10160 / (R_GAS * 248.15))
HL_LOWER = 0.917 * 575 * np.exp(-21400 / (R_GAS * 248.15)) / np.sqrt(1e-09 * SECONDS_PER_YEAR * 0.917)
HL_CRITICAL = (np.log(550 / 367) - np.log(350 / 567)) / HL_UPPER


def test_ela_shift_reference():
    result = domain.ela_shift(*[-6, 0.0065])
    np.testing.assert_allclose(result, -6 / 0.0065, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_ela_shift_rejects_invalid(bad):
    args = [-6, 0.0065]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ela_shift(*args)


@pytest.mark.parametrize("index,value", [(1, 0)])
def test_ela_shift_physical_bound(index, value):
    args = [-6, 0.0065]
    args[index] = value
    with pytest.raises(ValueError):
        domain.ela_shift(*args)


def test_driving_stress_reference():
    result = domain.driving_stress(*[1000, 0.1, 900, 10])
    np.testing.assert_allclose(result, 9e5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_driving_stress_rejects_invalid(bad):
    args = [1000, 0.1, 900, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.driving_stress(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0)])
def test_driving_stress_physical_bound(index, value):
    args = [1000, 0.1, 900, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.driving_stress(*args)


def test_viscous_slab_velocity_reference():
    result = domain.viscous_slab_velocity(*[50, 100, 0.1, 10000000000000.0, 1000, 10])
    np.testing.assert_allclose(result, 3.75e-7, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_viscous_slab_velocity_rejects_invalid(bad):
    args = [50, 100, 0.1, 10000000000000.0, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.viscous_slab_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, -1), (3, 0), (4, 0), (5, 0)])
def test_viscous_slab_velocity_physical_bound(index, value):
    args = [50, 100, 0.1, 10000000000000.0, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.viscous_slab_velocity(*args)


def test_isostatic_depression_reference():
    result = domain.isostatic_depression(*[3000, 900, 3300])
    np.testing.assert_allclose(result, 3000 * 900 / 3300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_isostatic_depression_rejects_invalid(bad):
    args = [3000, 900, 3300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.isostatic_depression(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_isostatic_depression_physical_bound(index, value):
    args = [3000, 900, 3300]
    args[index] = value
    with pytest.raises(ValueError):
        domain.isostatic_depression(*args)


def test_stefan_ice_thickness_reference():
    result = domain.stefan_ice_thickness(*[260000000.0, 2.2, 917, 334000.0])
    np.testing.assert_allclose(result, np.sqrt(2 * 2.2 * 2.6e8 / (917 * 3.34e5)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stefan_ice_thickness_rejects_invalid(bad):
    args = [260000000.0, 2.2, 917, 334000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stefan_ice_thickness(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_stefan_ice_thickness_physical_bound(index, value):
    args = [260000000.0, 2.2, 917, 334000.0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stefan_ice_thickness(*args)


def test_stefan_growth_rate_reference():
    result = domain.stefan_growth_rate(*[1, 20, 2, 1000, 200000.0])
    np.testing.assert_allclose(result, 2e-7, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stefan_growth_rate_rejects_invalid(bad):
    args = [1, 20, 2, 1000, 200000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stefan_growth_rate(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, 0), (3, 0), (4, 0)])
def test_stefan_growth_rate_physical_bound(index, value):
    args = [1, 20, 2, 1000, 200000.0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stefan_growth_rate(*args)


def test_sea_ice_freeboard_reference():
    result = domain.sea_ice_freeboard(*[2, 900, 1000])
    np.testing.assert_allclose(result, 0.2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_sea_ice_freeboard_rejects_invalid(bad):
    args = [2, 900, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.sea_ice_freeboard(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_sea_ice_freeboard_physical_bound(index, value):
    args = [2, 900, 1000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.sea_ice_freeboard(*args)


def test_viscous_relaxation_time_reference():
    result = domain.viscous_relaxation_time(*[1000000.0, 1e21, 3300, 10])
    np.testing.assert_allclose(result, 4 * np.pi * 1e21 / 3.3e10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_viscous_relaxation_time_rejects_invalid(bad):
    args = [1000000.0, 1e21, 3300, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.viscous_relaxation_time(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0)])
def test_viscous_relaxation_time_physical_bound(index, value):
    args = [1000000.0, 1e21, 3300, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.viscous_relaxation_time(*args)


def test_haskell_viscosity_reference():
    result = domain.haskell_viscosity(*[100000000000.0, 1000000.0, 3300, 10])
    np.testing.assert_allclose(result, 3.3e21 / (4 * np.pi), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_haskell_viscosity_rejects_invalid(bad):
    args = [100000000000.0, 1000000.0, 3300, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.haskell_viscosity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0)])
def test_haskell_viscosity_physical_bound(index, value):
    args = [100000000000.0, 1000000.0, 3300, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.haskell_viscosity(*args)


def test_remaining_depression_reference():
    result = domain.remaining_depression(*[1000.0, 300, 1000.0])
    np.testing.assert_allclose(result, 300 / np.e, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_remaining_depression_rejects_invalid(bad):
    args = [1000.0, 300, 1000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.remaining_depression(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0)])
def test_remaining_depression_physical_bound(index, value):
    args = [1000.0, 300, 1000.0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.remaining_depression(*args)


def test_rebound_uplift_rate_reference():
    result = domain.rebound_uplift_rate(*[1000.0, 300, 1000.0])
    np.testing.assert_allclose(result, 0.3 / np.e, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_rebound_uplift_rate_rejects_invalid(bad):
    args = [1000.0, 300, 1000.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rebound_uplift_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0)])
def test_rebound_uplift_rate_physical_bound(index, value):
    args = [1000.0, 300, 1000.0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.rebound_uplift_rate(*args)


def test_nye_plastic_profile_reference():
    result = domain.nye_plastic_profile(*[100000.0, 100000.0, 1000, 10])
    np.testing.assert_allclose(result, np.sqrt(2e6), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_nye_plastic_profile_rejects_invalid(bad):
    args = [100000.0, 100000.0, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nye_plastic_profile(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_nye_plastic_profile_physical_bound(index, value):
    args = [100000.0, 100000.0, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nye_plastic_profile(*args)


def test_glen_strain_rate_reference():
    result = domain.glen_strain_rate(*[100000.0, 2.4e-24, 3])
    np.testing.assert_allclose(result, 2.4e-9, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glen_strain_rate_rejects_invalid(bad):
    args = [100000.0, 2.4e-24, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glen_strain_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5)])
def test_glen_strain_rate_physical_bound(index, value):
    args = [100000.0, 2.4e-24, 3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glen_strain_rate(*args)


def test_glen_stress_reference():
    result = domain.glen_stress(*[2.4e-09, 2.4e-24, 3])
    np.testing.assert_allclose(result, 1e5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glen_stress_rejects_invalid(bad):
    args = [2.4e-09, 2.4e-24, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glen_stress(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5)])
def test_glen_stress_physical_bound(index, value):
    args = [2.4e-09, 2.4e-24, 3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glen_stress(*args)


def test_glen_rate_factor_reference():
    result = domain.glen_rate_factor(*[263.15, 2.4e-24, 60000.0, 273.15])
    np.testing.assert_allclose(result, 2.4e-24 * np.exp(-6e4 / 8.31446261815324 * (1 / 263.15 - 1 / 273.15)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glen_rate_factor_rejects_invalid(bad):
    args = [263.15, 2.4e-24, 60000.0, 273.15]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glen_rate_factor(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1), (3, 0)])
def test_glen_rate_factor_physical_bound(index, value):
    args = [263.15, 2.4e-24, 60000.0, 273.15]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glen_rate_factor(*args)


def test_effective_viscosity_reference():
    result = domain.effective_viscosity(*[100000.0, 2.4e-24, 3])
    np.testing.assert_allclose(result, 1 / 4.8e-14, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_effective_viscosity_rejects_invalid(bad):
    args = [100000.0, 2.4e-24, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.effective_viscosity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0.5)])
def test_effective_viscosity_physical_bound(index, value):
    args = [100000.0, 2.4e-24, 3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.effective_viscosity(*args)


def test_glen_velocity_profile_reference():
    result = domain.glen_velocity_profile(*[50, 100, 0.1, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, 4.6875e-8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glen_velocity_profile_rejects_invalid(bad):
    args = [50, 100, 0.1, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glen_velocity_profile(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, -1), (3, 0), (4, 0.5), (5, 0), (6, 0)])
def test_glen_velocity_profile_physical_bound(index, value):
    args = [50, 100, 0.1, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glen_velocity_profile(*args)


def test_glen_surface_velocity_reference():
    result = domain.glen_surface_velocity(*[100, 0.1, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, 5e-8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glen_surface_velocity_rejects_invalid(bad):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glen_surface_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, 0), (3, 0.5), (4, 0), (5, 0)])
def test_glen_surface_velocity_physical_bound(index, value):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glen_surface_velocity(*args)


def test_nye_crevasse_depth_reference():
    result = domain.nye_crevasse_depth(*[1e-09, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, 20, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_nye_crevasse_depth_rejects_invalid(bad):
    args = [1e-09, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nye_crevasse_depth(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5), (3, 0), (4, 0)])
def test_nye_crevasse_depth_physical_bound(index, value):
    args = [1e-09, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nye_crevasse_depth(*args)


def test_water_filled_crevasse_depth_reference():
    result = domain.water_filled_crevasse_depth(*[1e-09, 1e-24, 3, 10, 900, 1000, 10])
    np.testing.assert_allclose(result, 2e5 / 9000 + 1000 / 900 * 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_water_filled_crevasse_depth_rejects_invalid(bad):
    args = [1e-09, 1e-24, 3, 10, 900, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.water_filled_crevasse_depth(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5), (3, -1), (4, 0), (5, 0), (6, 0)])
def test_water_filled_crevasse_depth_physical_bound(index, value):
    args = [1e-09, 1e-24, 3, 10, 900, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.water_filled_crevasse_depth(*args)


def test_weertman_sliding_velocity_reference():
    result = domain.weertman_sliding_velocity(*[100000.0, 1e-17, 2])
    np.testing.assert_allclose(result, 1e-7, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_weertman_sliding_velocity_rejects_invalid(bad):
    args = [100000.0, 1e-17, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.weertman_sliding_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5)])
def test_weertman_sliding_velocity_physical_bound(index, value):
    args = [100000.0, 1e-17, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.weertman_sliding_velocity(*args)


def test_linear_mass_balance_reference():
    result = domain.linear_mass_balance(*[3000, 2500, 1e-10])
    np.testing.assert_allclose(result, 5e-8, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_linear_mass_balance_rejects_invalid(bad):
    args = [3000, 2500, 1e-10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.linear_mass_balance(*args)


@pytest.mark.parametrize("index,value", [(2, -1)])
def test_linear_mass_balance_physical_bound(index, value):
    args = [3000, 2500, 1e-10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.linear_mass_balance(*args)


def test_accumulation_area_ratio_reference():
    result = domain.accumulation_area_ratio(*[6, 10])
    np.testing.assert_allclose(result, 0.6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_accumulation_area_ratio_rejects_invalid(bad):
    args = [6, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.accumulation_area_ratio(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_accumulation_area_ratio_physical_bound(index, value):
    args = [6, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.accumulation_area_ratio(*args)


def test_dansgaard_delta_reference():
    result = domain.dansgaard_delta(*[253.15, 0.69, -13.6])
    np.testing.assert_allclose(result, 0.69 * -20 - 13.6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dansgaard_delta_rejects_invalid(bad):
    args = [253.15, 0.69, -13.6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dansgaard_delta(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_dansgaard_delta_physical_bound(index, value):
    args = [253.15, 0.69, -13.6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dansgaard_delta(*args)


def test_dansgaard_temperature_reference():
    result = domain.dansgaard_temperature(*[-27.4, 0.69, -13.6])
    np.testing.assert_allclose(result, 253.15, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dansgaard_temperature_rejects_invalid(bad):
    args = [-27.4, 0.69, -13.6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dansgaard_temperature(*args)


@pytest.mark.parametrize("index,value", [(1, 0)])
def test_dansgaard_temperature_physical_bound(index, value):
    args = [-27.4, 0.69, -13.6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dansgaard_temperature(*args)


def test_nye_age_reference():
    result = domain.nye_age(*[2000, 3000, 1e-08])
    np.testing.assert_allclose(result, 3e11 * np.log(3), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_nye_age_rejects_invalid(bad):
    args = [2000, 3000, 1e-08]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nye_age(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_nye_age_physical_bound(index, value):
    args = [2000, 3000, 1e-08]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nye_age(*args)


def test_dansgaard_johnsen_downward_velocity_reference():
    result = domain.dansgaard_johnsen_downward_velocity(*[1000, 3000, 500, 1e-08])
    np.testing.assert_allclose(result, 1e-8 * 3500 / 5500, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dansgaard_johnsen_downward_velocity_rejects_invalid(bad):
    args = [1000, 3000, 500, 1e-08]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dansgaard_johnsen_downward_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_dansgaard_johnsen_downward_velocity_physical_bound(index, value):
    args = [1000, 3000, 500, 1e-08]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dansgaard_johnsen_downward_velocity(*args)


def test_dansgaard_johnsen_age_reference():
    result = domain.dansgaard_johnsen_age(*[2800, 3000, 500, 1e-08])
    np.testing.assert_allclose(result, 5500 / 2e-8 * np.log(11) + 500 * 5500 / 1e-8 * (1 / 200 - 1 / 500), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dansgaard_johnsen_age_rejects_invalid(bad):
    args = [2800, 3000, 500, 1e-08]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dansgaard_johnsen_age(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_dansgaard_johnsen_age_physical_bound(index, value):
    args = [2800, 3000, 500, 1e-08]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dansgaard_johnsen_age(*args)


def test_flotation_thickness_reference():
    result = domain.flotation_thickness(*[900, 900, 1000])
    np.testing.assert_allclose(result, 1000, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_flotation_thickness_rejects_invalid(bad):
    args = [900, 900, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flotation_thickness(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_flotation_thickness_physical_bound(index, value):
    args = [900, 900, 1000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.flotation_thickness(*args)


def test_ice_shelf_strain_rate_reference():
    result = domain.ice_shelf_strain_rate(*[400, 1e-24, 3, 900, 1000, 10])
    np.testing.assert_allclose(result, 7.29e-10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_ice_shelf_strain_rate_rejects_invalid(bad):
    args = [400, 1e-24, 3, 900, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ice_shelf_strain_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5), (3, 0), (4, 0), (5, 0)])
def test_ice_shelf_strain_rate_physical_bound(index, value):
    args = [400, 1e-24, 3, 900, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.ice_shelf_strain_rate(*args)


def test_shallow_ice_flux_reference():
    result = domain.shallow_ice_flux(*[100, 0.1, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, 4e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_shallow_ice_flux_rejects_invalid(bad):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shallow_ice_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0.5), (4, 0), (5, 0)])
def test_shallow_ice_flux_physical_bound(index, value):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.shallow_ice_flux(*args)


def test_shallow_ice_diffusivity_reference():
    result = domain.shallow_ice_diffusivity(*[100, 0.1, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, 4e-5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_shallow_ice_diffusivity_rejects_invalid(bad):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shallow_ice_diffusivity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0.5), (4, 0), (5, 0)])
def test_shallow_ice_diffusivity_physical_bound(index, value):
    args = [100, 0.1, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.shallow_ice_diffusivity(*args)


def test_vialov_dome_height_reference():
    result = domain.vialov_dome_height(*[1000000.0, 1e-08, 1e-24, 3, 1000, 10])
    np.testing.assert_allclose(result, (2 * 25000 ** (1 / 3)) ** (3 / 8) * 1e3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_vialov_dome_height_rejects_invalid(bad):
    args = [1000000.0, 1e-08, 1e-24, 3, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.vialov_dome_height(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0.5), (4, 0), (5, 0)])
def test_vialov_dome_height_physical_bound(index, value):
    args = [1000000.0, 1e-08, 1e-24, 3, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.vialov_dome_height(*args)


def test_vialov_profile_reference():
    result = domain.vialov_profile(*[0.5, 1, 100, 3])
    np.testing.assert_allclose(result, 100 * (1 - 0.5 ** (4 / 3)) ** (3 / 8), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_vialov_profile_rejects_invalid(bad):
    args = [0.5, 1, 100, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.vialov_profile(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, -1), (3, 0.5)])
def test_vialov_profile_physical_bound(index, value):
    args = [0.5, 1, 100, 3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.vialov_profile(*args)


def test_glacier_response_time_reference():
    result = domain.glacier_response_time(*[200, 1e-07])
    np.testing.assert_allclose(result, 2e9, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_glacier_response_time_rejects_invalid(bad):
    args = [200, 1e-07]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.glacier_response_time(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_glacier_response_time_physical_bound(index, value):
    args = [200, 1e-07]
    args[index] = value
    with pytest.raises(ValueError):
        domain.glacier_response_time(*args)


def test_gas_age_difference_reference():
    result = domain.gas_age_difference(*[100, 600, 900, 1e-09])
    np.testing.assert_allclose(result, 6e4 / 9e-7, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_gas_age_difference_rejects_invalid(bad):
    args = [100, 600, 900, 1e-09]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.gas_age_difference(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_gas_age_difference_physical_bound(index, value):
    args = [100, 600, 900, 1e-09]
    args[index] = value
    with pytest.raises(ValueError):
        domain.gas_age_difference(*args)


def test_schoof_grounding_line_flux_reference():
    result = domain.schoof_grounding_line_flux(*[1000, 1e-24, 3, 1000000.0, 1, 900, 1000, 10])
    np.testing.assert_allclose(result, np.sqrt(1e-24 * 9000.0**4 * 0.1**3 / 64e6) * 1000**3.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_schoof_grounding_line_flux_rejects_invalid(bad):
    args = [1000, 1e-24, 3, 1000000.0, 1, 900, 1000, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.schoof_grounding_line_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0.5), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0)])
def test_schoof_grounding_line_flux_physical_bound(index, value):
    args = [1000, 1e-24, 3, 1000000.0, 1, 900, 1000, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.schoof_grounding_line_flux(*args)


def test_herron_langway_density_reference():
    result = domain.herron_langway_density(*[50, 248.15, 1e-09, 350])
    np.testing.assert_allclose(result, 917 / (1 + np.exp(-(np.log(550 / 367) + HL_LOWER * (50 - HL_CRITICAL)))), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_herron_langway_density_rejects_invalid(bad):
    args = [50, 248.15, 1e-09, 350]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.herron_langway_density(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (1, 274), (2, 0), (3, 0), (3, 550)])
def test_herron_langway_density_physical_bound(index, value):
    args = [50, 248.15, 1e-09, 350]
    args[index] = value
    with pytest.raises(ValueError):
        domain.herron_langway_density(*args)


def test_herron_langway_depth_reference():
    result = domain.herron_langway_depth(*[830, 248.15, 1e-09, 350])
    np.testing.assert_allclose(result, HL_CRITICAL + (np.log(830 / 87) - np.log(550 / 367)) / HL_LOWER, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_herron_langway_depth_rejects_invalid(bad):
    args = [830, 248.15, 1e-09, 350]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.herron_langway_depth(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (0, 917), (1, 0), (2, 0), (3, 0)])
def test_herron_langway_depth_physical_bound(index, value):
    args = [830, 248.15, 1e-09, 350]
    args[index] = value
    with pytest.raises(ValueError):
        domain.herron_langway_depth(*args)


def test_herron_langway_mean_density_reference():
    result = domain.herron_langway_mean_density(*[50, 248.15, 1e-09, 350])
    np.testing.assert_allclose(result, 710.4412496420827, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_herron_langway_mean_density_rejects_invalid(bad):
    args = [50, 248.15, 1e-09, 350]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.herron_langway_mean_density(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0)])
def test_herron_langway_mean_density_physical_bound(index, value):
    args = [50, 248.15, 1e-09, 350]
    args[index] = value
    with pytest.raises(ValueError):
        domain.herron_langway_mean_density(*args)


def test_expected_positive_degree_days_reference():
    result = domain.expected_positive_degree_days(*[273.15, 5, 86400])
    np.testing.assert_allclose(result, 86400 * 5 / np.sqrt(2 * np.pi), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_expected_positive_degree_days_rejects_invalid(bad):
    args = [273.15, 5, 86400]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.expected_positive_degree_days(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1)])
def test_expected_positive_degree_days_physical_bound(index, value):
    args = [273.15, 5, 86400]
    args[index] = value
    with pytest.raises(ValueError):
        domain.expected_positive_degree_days(*args)


def test_pdd_surface_melt_reference():
    result = domain.pdd_surface_melt(*[500 * 86400, 0.5, 3e-3 / 86400, 8e-3 / 86400, 0.6])
    np.testing.assert_allclose(result, [0.5, 8 / 3, 8 / 3 + 0.2], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_pdd_surface_melt_rejects_invalid(bad):
    args = [500 * 86400, 0.5, 3e-3 / 86400, 8e-3 / 86400, 0.6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.pdd_surface_melt(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0), (4, -0.1), (4, 1.1)])
def test_pdd_surface_melt_physical_bound(index, value):
    args = [500 * 86400, 0.5, 3e-3 / 86400, 8e-3 / 86400, 0.6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.pdd_surface_melt(*args)


@pytest.mark.parametrize(
    "name,args",
    [
        ("viscous_slab_velocity", [101, 100, 0.1, 10000000000000.0, 1000, 10]),
        ("glen_velocity_profile", [101, 100, 0.1, 1e-24, 3, 1000, 10]),
        ("sea_ice_freeboard", [2, 1000, 1000]),
        ("ice_shelf_strain_rate", [400, 1e-24, 3, 1030, 1000, 10]),
        ("schoof_grounding_line_flux", [1000, 1e-24, 3, 1000000.0, 1, 1030, 1000, 10]),
        ("vialov_profile", [1.5, 1, 100, 3]),
        ("accumulation_area_ratio", [11, 10]),
        ("nye_age", [3000, 3000, 1e-08]),
        ("dansgaard_johnsen_downward_velocity", [3001, 3000, 500, 1e-08]),
        ("dansgaard_johnsen_downward_velocity", [1000, 3000, 3001, 1e-08]),
        ("dansgaard_johnsen_age", [3000, 3000, 500, 1e-08]),
        ("dansgaard_johnsen_age", [1000, 3000, 3001, 1e-08]),
        ("gas_age_difference", [100, 950, 900, 1e-09]),
        ("herron_langway_depth", [300, 248.15, 1e-09, 350]),
    ],
)
def test_cross_parameter_rejections(name, args):
    with pytest.raises(ValueError):
        getattr(domain, name)(*args)
