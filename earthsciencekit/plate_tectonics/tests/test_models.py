"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import plate_tectonics as domain


def test_angular_speed_reference():
    result = domain.angular_speed(*[0.1, 10])
    np.testing.assert_allclose(result, 0.01, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_angular_speed_rejects_invalid(bad):
    args = [0.1, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.angular_speed(*args)


def test_plate_speed_reference():
    result = domain.plate_speed(*[1e-09, 6000000.0, 1.5707963267948966])
    np.testing.assert_allclose(result, 0.006, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_plate_speed_rejects_invalid(bad):
    args = [1e-09, 6000000.0, 1.5707963267948966]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.plate_speed(*args)


def test_euler_velocity_x_reference():
    result = domain.euler_velocity_x(*[1, 2, 3, 4])
    np.testing.assert_allclose(result, -2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_euler_velocity_x_rejects_invalid(bad):
    args = [1, 2, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.euler_velocity_x(*args)


def test_euler_velocity_y_reference():
    result = domain.euler_velocity_y(*[2, 1, 4, 3])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_euler_velocity_y_rejects_invalid(bad):
    args = [2, 1, 4, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.euler_velocity_y(*args)


def test_euler_velocity_z_reference():
    result = domain.euler_velocity_z(*[1, 2, 3, 4])
    np.testing.assert_allclose(result, -2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_euler_velocity_z_rejects_invalid(bad):
    args = [1, 2, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.euler_velocity_z(*args)


def test_plate_displacement_reference():
    result = domain.plate_displacement(*[0.01, 100])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_plate_displacement_rejects_invalid(bad):
    args = [0.01, 100]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.plate_displacement(*args)


def test_spreading_age_reference():
    result = domain.spreading_age(*[100, 2])
    np.testing.assert_allclose(result, 50, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_spreading_age_rejects_invalid(bad):
    args = [100, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.spreading_age(*args)


def test_full_spreading_rate_reference():
    result = domain.full_spreading_rate(*[0.001])
    np.testing.assert_allclose(result, 0.002, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_full_spreading_rate_rejects_invalid(bad):
    args = [0.001]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.full_spreading_rate(*args)


def test_convergence_rate_reference():
    result = domain.convergence_rate(*[2, 1.0471975511965976])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_convergence_rate_rejects_invalid(bad):
    args = [2, 1.0471975511965976]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.convergence_rate(*args)


def test_fault_slip_rate_reference():
    result = domain.fault_slip_rate(*[2, 0.5235987755982988])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_fault_slip_rate_rejects_invalid(bad):
    args = [2, 0.5235987755982988]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.fault_slip_rate(*args)


def test_subduction_travel_time_reference():
    result = domain.subduction_travel_time(*[100, 2, 0.5235987755982988])
    np.testing.assert_allclose(result, 100, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_subduction_travel_time_rejects_invalid(bad):
    args = [100, 2, 0.5235987755982988]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.subduction_travel_time(*args)


def test_slab_pull_force_reference():
    result = domain.slab_pull_force(*[100, 1000000000.0, 10])
    np.testing.assert_allclose(result, 1e12, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_slab_pull_force_rejects_invalid(bad):
    args = [100, 1000000000.0, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.slab_pull_force(*args)


def test_airy_root_reference():
    result = domain.airy_root(*[1000, 2700, 3300])
    np.testing.assert_allclose(result, 4500, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_airy_root_rejects_invalid(bad):
    args = [1000, 2700, 3300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.airy_root(*args)


def test_flexural_rigidity_reference():
    result = domain.flexural_rigidity(*[12, 1, 0])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_flexural_rigidity_rejects_invalid(bad):
    args = [12, 1, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flexural_rigidity(*args)


def test_flexural_parameter_reference():
    result = domain.flexural_parameter(*[4, 1, 1])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_flexural_parameter_rejects_invalid(bad):
    args = [4, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flexural_parameter(*args)


def test_halfspace_temperature_reference():
    result = domain.halfspace_temperature(*[0, 1, 1, 300, 1300])
    np.testing.assert_allclose(result, 300, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_halfspace_temperature_rejects_invalid(bad):
    args = [0, 1, 1, 300, 1300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_heat_flux_reference():
    result = domain.halfspace_heat_flux(*[1, 1, 1, 300, 1300])
    np.testing.assert_allclose(result, 1000 / np.sqrt(np.pi), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_halfspace_heat_flux_rejects_invalid(bad):
    args = [1, 1, 1, 300, 1300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_angular_speed_broadcasting():
    args = [0.1, 10]
    reference = domain.angular_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.angular_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_angular_speed_time_physical_bound():
    args = [0.1, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.angular_speed(*args)


def test_plate_speed_broadcasting():
    args = [1e-09, 6000000.0, 1.5707963267948966]
    reference = domain.plate_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.plate_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_plate_speed_angular_speed_physical_bound():
    args = [1e-09, 6000000.0, 1.5707963267948966]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.plate_speed(*args)


def test_plate_speed_radius_physical_bound():
    args = [1e-09, 6000000.0, 1.5707963267948966]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.plate_speed(*args)


def test_plate_speed_angular_distance_physical_bound():
    args = [1e-09, 6000000.0, 1.5707963267948966]
    args[2] = 7
    with pytest.raises(ValueError):
        domain.plate_speed(*args)


def test_euler_velocity_x_broadcasting():
    args = [1, 2, 3, 4]
    reference = domain.euler_velocity_x(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.euler_velocity_x(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_euler_velocity_y_broadcasting():
    args = [2, 1, 4, 3]
    reference = domain.euler_velocity_y(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.euler_velocity_y(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_euler_velocity_z_broadcasting():
    args = [1, 2, 3, 4]
    reference = domain.euler_velocity_z(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.euler_velocity_z(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_plate_displacement_broadcasting():
    args = [0.01, 100]
    reference = domain.plate_displacement(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.plate_displacement(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_plate_displacement_time_physical_bound():
    args = [0.01, 100]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.plate_displacement(*args)


def test_spreading_age_broadcasting():
    args = [100, 2]
    reference = domain.spreading_age(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.spreading_age(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_spreading_age_distance_physical_bound():
    args = [100, 2]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.spreading_age(*args)


def test_spreading_age_half_spreading_rate_physical_bound():
    args = [100, 2]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.spreading_age(*args)


def test_full_spreading_rate_broadcasting():
    args = [0.001]
    reference = domain.full_spreading_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.full_spreading_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_full_spreading_rate_half_spreading_rate_physical_bound():
    args = [0.001]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.full_spreading_rate(*args)


def test_convergence_rate_broadcasting():
    args = [2, 1.0471975511965976]
    reference = domain.convergence_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.convergence_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_convergence_rate_relative_speed_physical_bound():
    args = [2, 1.0471975511965976]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.convergence_rate(*args)


def test_convergence_rate_obliquity_physical_bound():
    args = [2, 1.0471975511965976]
    args[1] = 7
    with pytest.raises(ValueError):
        domain.convergence_rate(*args)


def test_fault_slip_rate_broadcasting():
    args = [2, 0.5235987755982988]
    reference = domain.fault_slip_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.fault_slip_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_fault_slip_rate_relative_speed_physical_bound():
    args = [2, 0.5235987755982988]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.fault_slip_rate(*args)


def test_fault_slip_rate_obliquity_physical_bound():
    args = [2, 0.5235987755982988]
    args[1] = 7
    with pytest.raises(ValueError):
        domain.fault_slip_rate(*args)


def test_subduction_travel_time_broadcasting():
    args = [100, 2, 0.5235987755982988]
    reference = domain.subduction_travel_time(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.subduction_travel_time(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_subduction_travel_time_vertical_depth_physical_bound():
    args = [100, 2, 0.5235987755982988]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.subduction_travel_time(*args)


def test_subduction_travel_time_slab_speed_physical_bound():
    args = [100, 2, 0.5235987755982988]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.subduction_travel_time(*args)


def test_subduction_travel_time_dip_physical_bound():
    args = [100, 2, 0.5235987755982988]
    args[2] = 2
    with pytest.raises(ValueError):
        domain.subduction_travel_time(*args)


def test_slab_pull_force_broadcasting():
    args = [100, 1000000000.0, 10]
    reference = domain.slab_pull_force(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.slab_pull_force(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_slab_pull_force_density_excess_physical_bound():
    args = [100, 1000000000.0, 10]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.slab_pull_force(*args)


def test_slab_pull_force_slab_volume_physical_bound():
    args = [100, 1000000000.0, 10]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.slab_pull_force(*args)


def test_slab_pull_force_gravity_physical_bound():
    args = [100, 1000000000.0, 10]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.slab_pull_force(*args)


def test_airy_root_broadcasting():
    args = [1000, 2700, 3300]
    reference = domain.airy_root(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.airy_root(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_airy_root_elevation_physical_bound():
    args = [1000, 2700, 3300]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.airy_root(*args)


def test_airy_root_crust_density_physical_bound():
    args = [1000, 2700, 3300]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.airy_root(*args)


def test_airy_root_mantle_density_physical_bound():
    args = [1000, 2700, 3300]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.airy_root(*args)


def test_flexural_rigidity_broadcasting():
    args = [12, 1, 0]
    reference = domain.flexural_rigidity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.flexural_rigidity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_flexural_rigidity_young_modulus_physical_bound():
    args = [12, 1, 0]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.flexural_rigidity(*args)


def test_flexural_rigidity_elastic_thickness_physical_bound():
    args = [12, 1, 0]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.flexural_rigidity(*args)


def test_flexural_rigidity_poisson_ratio_physical_bound():
    args = [12, 1, 0]
    args[2] = 0.5
    with pytest.raises(ValueError):
        domain.flexural_rigidity(*args)


def test_flexural_parameter_broadcasting():
    args = [4, 1, 1]
    reference = domain.flexural_parameter(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.flexural_parameter(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_flexural_parameter_rigidity_physical_bound():
    args = [4, 1, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.flexural_parameter(*args)


def test_flexural_parameter_density_contrast_physical_bound():
    args = [4, 1, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.flexural_parameter(*args)


def test_flexural_parameter_gravity_physical_bound():
    args = [4, 1, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.flexural_parameter(*args)


def test_halfspace_temperature_broadcasting():
    args = [0, 1, 1, 300, 1300]
    reference = domain.halfspace_temperature(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.halfspace_temperature(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_halfspace_temperature_depth_physical_bound():
    args = [0, 1, 1, 300, 1300]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_temperature_age_physical_bound():
    args = [0, 1, 1, 300, 1300]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_temperature_diffusivity_physical_bound():
    args = [0, 1, 1, 300, 1300]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_temperature_surface_temperature_physical_bound():
    args = [0, 1, 1, 300, 1300]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_temperature_mantle_temperature_physical_bound():
    args = [0, 1, 1, 300, 1300]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.halfspace_temperature(*args)


def test_halfspace_heat_flux_broadcasting():
    args = [1, 1, 1, 300, 1300]
    reference = domain.halfspace_heat_flux(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.halfspace_heat_flux(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_halfspace_heat_flux_age_physical_bound():
    args = [1, 1, 1, 300, 1300]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_halfspace_heat_flux_diffusivity_physical_bound():
    args = [1, 1, 1, 300, 1300]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_halfspace_heat_flux_conductivity_physical_bound():
    args = [1, 1, 1, 300, 1300]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_halfspace_heat_flux_surface_temperature_physical_bound():
    args = [1, 1, 1, 300, 1300]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_halfspace_heat_flux_mantle_temperature_physical_bound():
    args = [1, 1, 1, 300, 1300]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.halfspace_heat_flux(*args)


def test_mantle_rayleigh_number_reference():
    result = domain.mantle_rayleigh_number(*[1, 2, 3, 4, 5, 6, 7])
    np.testing.assert_allclose(result, 3000 / 42, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_mantle_rayleigh_number_rejects_invalid(bad):
    args = [1, 2, 3, 4, 5, 6, 7]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_broadcasting():
    args = [1, 2, 3, 4, 5, 6, 7]
    reference = domain.mantle_rayleigh_number(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.mantle_rayleigh_number(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_mantle_rayleigh_number_thermal_expansivity_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_density_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_gravity_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_temperature_difference_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[3] = -1
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_thickness_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_diffusivity_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_mantle_rayleigh_number_viscosity_physical_bound():
    args = [1, 2, 3, 4, 5, 6, 7]
    args[6] = 0
    with pytest.raises(ValueError):
        domain.mantle_rayleigh_number(*args)


def test_flexural_deflection_reference():
    result = domain.flexural_deflection(*[0, 8, 1, 4, 1])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_flexural_deflection_rejects_invalid(bad):
    args = [0, 8, 1, 4, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flexural_deflection(*args)


def test_flexural_deflection_broadcasting():
    args = [0, 8, 1, 4, 1]
    reference = domain.flexural_deflection(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.flexural_deflection(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_flexural_deflection_line_load_physical_bound():
    args = [0, 8, 1, 4, 1]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.flexural_deflection(*args)


def test_flexural_deflection_rigidity_physical_bound():
    args = [0, 8, 1, 4, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.flexural_deflection(*args)


def test_flexural_deflection_density_contrast_physical_bound():
    args = [0, 8, 1, 4, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.flexural_deflection(*args)


def test_flexural_deflection_gravity_physical_bound():
    args = [0, 8, 1, 4, 1]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.flexural_deflection(*args)


def test_hotspot_age_reference():
    result = domain.hotspot_age(*[6, 1, 2, 1.5707963267948966])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_hotspot_age_rejects_invalid(bad):
    args = [6, 1, 2, 1.5707963267948966]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hotspot_age(*args)


def test_hotspot_age_broadcasting():
    args = [6, 1, 2, 1.5707963267948966]
    reference = domain.hotspot_age(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.hotspot_age(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_hotspot_age_track_length_physical_bound():
    args = [6, 1, 2, 1.5707963267948966]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.hotspot_age(*args)


def test_hotspot_age_angular_speed_physical_bound():
    args = [6, 1, 2, 1.5707963267948966]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.hotspot_age(*args)


def test_hotspot_age_radius_physical_bound():
    args = [6, 1, 2, 1.5707963267948966]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.hotspot_age(*args)


def test_hotspot_age_pole_distance_physical_bound():
    args = [6, 1, 2, 1.5707963267948966]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.hotspot_age(*args)


def test_ocean_depth_reference():
    result = domain.ocean_depth(*[0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300])
    np.testing.assert_allclose(result, 3500, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ocean_depth_rejects_invalid(bad):
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_broadcasting():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    reference = domain.ocean_depth(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.ocean_depth(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_ocean_depth_age_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_ridge_depth_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_mantle_density_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_water_density_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_thermal_expansivity_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_diffusivity_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_surface_temperature_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[6] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_ocean_depth_mantle_temperature_physical_bound():
    args = [0.25 * np.pi, 2500, 2, 1, 0.5, 1, 300, 1300]
    args[7] = 0
    with pytest.raises(ValueError):
        domain.ocean_depth(*args)


def test_hotspot_age_pole_distance_upper_bound():
    with pytest.raises(ValueError):
        domain.hotspot_age(6, 1, 2, np.pi)


def test_compose_euler_vectors_reference():
    result = domain.compose_euler_vectors(*[[1, 2, 3], [4, -5, 6]])
    np.testing.assert_allclose(result, [5, -3, 9], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_compose_euler_vectors_rejects_invalid(bad):
    args = [[1, 2, 3], [4, -5, 6]]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.compose_euler_vectors(*args)


def test_compose_euler_vectors_broadcasting():
    result = domain.compose_euler_vectors(np.full((2, 3), 1.0), [4, -5, 6])
    assert result.shape == (2, 3)
    np.testing.assert_allclose(result, [[5, -4, 7]] * 2)


@pytest.mark.parametrize("shape", [(2,), (3, 2)])
def test_compose_euler_vectors_requires_three_components(shape):
    with pytest.raises(ValueError):
        domain.compose_euler_vectors(np.ones(shape), [1, 2, 3])
    with pytest.raises(ValueError):
        domain.compose_euler_vectors([1, 2, 3], np.ones(shape))
