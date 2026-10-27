"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import geomorphology as domain
from earthsciencekit.constants import VON_KARMAN

BAD = [np.nan, np.inf, -np.inf, [], "bad", True, 1j]


def test_stokes_settling_velocity_reference():
    result = domain.stokes_settling_velocity(*[0.0001, 2650, 1000, 0.001, 9.81])
    np.testing.assert_allclose(result, 1650 * 9.81 * 1e-8 / 18e-3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stokes_settling_velocity_rejects_invalid(bad):
    args = [0.0001, 2650, 1000, 0.001, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stokes_settling_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 1000), (2, 0), (3, 0), (4, 0)])
def test_stokes_settling_velocity_physical_bound(index, value):
    args = [0.0001, 2650, 1000, 0.001, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stokes_settling_velocity(*args)


def test_ferguson_church_settling_velocity_reference():
    result = domain.ferguson_church_settling_velocity(*[0.0001, 2650, 1000, 1e-06, 9.81])
    np.testing.assert_allclose(result, 1.65 * 9.81 * 1e-8 / (18e-6 + np.sqrt(0.75 * 1.65 * 9.81 * 1e-12)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_ferguson_church_settling_velocity_rejects_invalid(bad):
    args = [0.0001, 2650, 1000, 1e-06, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ferguson_church_settling_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 900), (2, 0), (3, 0), (4, 0)])
def test_ferguson_church_settling_velocity_physical_bound(index, value):
    args = [0.0001, 2650, 1000, 1e-06, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.ferguson_church_settling_velocity(*args)


def test_bed_shear_stress_reference():
    result = domain.bed_shear_stress(*[2, 0.001, 1000, 9.81])
    np.testing.assert_allclose(result, 19.62, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_bed_shear_stress_rejects_invalid(bad):
    args = [2, 0.001, 1000, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bed_shear_stress(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0)])
def test_bed_shear_stress_physical_bound(index, value):
    args = [2, 0.001, 1000, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bed_shear_stress(*args)


def test_shear_velocity_reference():
    result = domain.shear_velocity(*[19.62, 1000])
    np.testing.assert_allclose(result, np.sqrt(0.01962), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_shear_velocity_rejects_invalid(bad):
    args = [19.62, 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shear_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_shear_velocity_physical_bound(index, value):
    args = [19.62, 1000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.shear_velocity(*args)


def test_shields_number_reference():
    result = domain.shields_number(*[19.62, 0.01, 2650, 1000, 9.81])
    np.testing.assert_allclose(result, 19.62 / (1650 * 9.81 * 0.01), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_shields_number_rejects_invalid(bad):
    args = [19.62, 0.01, 2650, 1000, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shields_number(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 1000), (3, 0), (4, 0)])
def test_shields_number_physical_bound(index, value):
    args = [19.62, 0.01, 2650, 1000, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.shields_number(*args)


def test_dimensionless_grain_size_reference():
    result = domain.dimensionless_grain_size(*[0.001, 2650, 1000, 1e-06, 9.81])
    np.testing.assert_allclose(result, 1e-3 * (1.65 * 9.81 / 1e-12) ** (1 / 3), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dimensionless_grain_size_rejects_invalid(bad):
    args = [0.001, 2650, 1000, 1e-06, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dimensionless_grain_size(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 1000), (2, 0), (3, 0), (4, 0)])
def test_dimensionless_grain_size_physical_bound(index, value):
    args = [0.001, 2650, 1000, 1e-06, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dimensionless_grain_size(*args)


def test_critical_shields_number_reference():
    result = domain.critical_shields_number(*[10])
    np.testing.assert_allclose(result, 0.30 / 13 + 0.055 * (1 - np.exp(-0.2)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_critical_shields_number_rejects_invalid(bad):
    args = [10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.critical_shields_number(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_critical_shields_number_physical_bound(index, value):
    args = [10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.critical_shields_number(*args)


def test_keulegan_mean_velocity_reference():
    result = domain.keulegan_mean_velocity(*[0.1, 2, 0.022])
    np.testing.assert_allclose(result, 0.1 / VON_KARMAN * np.log(1000), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_keulegan_mean_velocity_rejects_invalid(bad):
    args = [0.1, 2, 0.022]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.keulegan_mean_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0.001), (2, 0)])
def test_keulegan_mean_velocity_physical_bound(index, value):
    args = [0.1, 2, 0.022]
    args[index] = value
    with pytest.raises(ValueError):
        domain.keulegan_mean_velocity(*args)


def test_rouse_number_reference():
    result = domain.rouse_number(*[0.041, 0.1])
    np.testing.assert_allclose(result, 0.041 / (VON_KARMAN * 0.1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_rouse_number_rejects_invalid(bad):
    args = [0.041, 0.1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rouse_number(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_rouse_number_physical_bound(index, value):
    args = [0.041, 0.1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.rouse_number(*args)


def test_rouse_concentration_reference():
    result = domain.rouse_concentration(*[0.5, 1, 0.05, 2, 1])
    np.testing.assert_allclose(result, 2 * 0.05 / 0.95, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_rouse_concentration_rejects_invalid(bad):
    args = [0.5, 1, 0.05, 2, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rouse_concentration(*args)


@pytest.mark.parametrize("index,value", [(0, 0.01), (0, 1.5), (1, 0.05), (2, 0), (3, -1), (4, -1)])
def test_rouse_concentration_physical_bound(index, value):
    args = [0.5, 1, 0.05, 2, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.rouse_concentration(*args)


def test_meyer_peter_muller_transport_reference():
    result = domain.meyer_peter_muller_transport(*[0.147, 0.047])
    np.testing.assert_allclose(result, 8 * 0.1**1.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_meyer_peter_muller_transport_rejects_invalid(bad):
    args = [0.147, 0.047]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.meyer_peter_muller_transport(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1)])
def test_meyer_peter_muller_transport_physical_bound(index, value):
    args = [0.147, 0.047]
    args[index] = value
    with pytest.raises(ValueError):
        domain.meyer_peter_muller_transport(*args)


def test_bedload_flux_reference():
    result = domain.bedload_flux(*[1, 0.01, 2650, 1000, 9.81])
    np.testing.assert_allclose(result, np.sqrt(1.65 * 9.81 * 0.01) * 0.01, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_bedload_flux_rejects_invalid(bad):
    args = [1, 0.01, 2650, 1000, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bedload_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 1000), (3, 0), (4, 0)])
def test_bedload_flux_physical_bound(index, value):
    args = [1, 0.01, 2650, 1000, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bedload_flux(*args)


def test_diffusive_hillslope_flux_reference():
    result = domain.diffusive_hillslope_flux(*[0.01, 0.2])
    np.testing.assert_allclose(result, -0.002, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_diffusive_hillslope_flux_rejects_invalid(bad):
    args = [0.01, 0.2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.diffusive_hillslope_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1)])
def test_diffusive_hillslope_flux_physical_bound(index, value):
    args = [0.01, 0.2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.diffusive_hillslope_flux(*args)


def test_nonlinear_hillslope_flux_reference():
    result = domain.nonlinear_hillslope_flux(*[0.01, 0.5, 1.0])
    np.testing.assert_allclose(result, -0.005 / 0.75, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_nonlinear_hillslope_flux_rejects_invalid(bad):
    args = [0.01, 0.5, 1.0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nonlinear_hillslope_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 1.0), (1, -1.5), (2, 0)])
def test_nonlinear_hillslope_flux_physical_bound(index, value):
    args = [0.01, 0.5, 1.0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nonlinear_hillslope_flux(*args)


def test_scarp_elevation_reference():
    result = domain.scarp_elevation(*[1, 2, 0.25, 1, 0.1])
    np.testing.assert_allclose(result, 2 * 0.8427007929497149 + 0.1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_scarp_elevation_rejects_invalid(bad):
    args = [1, 2, 0.25, 1, 0.1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.scarp_elevation(*args)


@pytest.mark.parametrize("index,value", [(1, -1), (2, -1), (3, -1)])
def test_scarp_elevation_physical_bound(index, value):
    args = [1, 2, 0.25, 1, 0.1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.scarp_elevation(*args)


def test_scarp_maximum_gradient_reference():
    result = domain.scarp_maximum_gradient(*[2, 0.25, 4, 0.1])
    np.testing.assert_allclose(result, 2 / np.sqrt(np.pi) + 0.1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_scarp_maximum_gradient_rejects_invalid(bad):
    args = [2, 0.25, 4, 0.1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.scarp_maximum_gradient(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_scarp_maximum_gradient_physical_bound(index, value):
    args = [2, 0.25, 4, 0.1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.scarp_maximum_gradient(*args)


def test_soil_production_rate_reference():
    result = domain.soil_production_rate(*[0.5, 1e-11, 0.5])
    np.testing.assert_allclose(result, 1e-11 * np.exp(-1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_soil_production_rate_rejects_invalid(bad):
    args = [0.5, 1e-11, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.soil_production_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0)])
def test_soil_production_rate_physical_bound(index, value):
    args = [0.5, 1e-11, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.soil_production_rate(*args)


def test_hilltop_curvature_reference():
    result = domain.hilltop_curvature(*[1e-12, 1e-10, 2600, 1300])
    np.testing.assert_allclose(result, -0.02, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_hilltop_curvature_rejects_invalid(bad):
    args = [1e-12, 1e-10, 2600, 1300]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hilltop_curvature(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_hilltop_curvature_physical_bound(index, value):
    args = [1e-12, 1e-10, 2600, 1300]
    args[index] = value
    with pytest.raises(ValueError):
        domain.hilltop_curvature(*args)


def test_stream_power_erosion_reference():
    result = domain.stream_power_erosion(*[1000000.0, 0.01, 1e-10, 0.5, 1])
    np.testing.assert_allclose(result, 1e-9, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stream_power_erosion_rejects_invalid(bad):
    args = [1000000.0, 0.01, 1e-10, 0.5, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stream_power_erosion(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, -1), (3, -1), (4, 0)])
def test_stream_power_erosion_physical_bound(index, value):
    args = [1000000.0, 0.01, 1e-10, 0.5, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stream_power_erosion(*args)


def test_steady_state_slope_reference():
    result = domain.steady_state_slope(*[1e-09, 1000000.0, 1e-10, 0.5, 1])
    np.testing.assert_allclose(result, 0.01, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_steady_state_slope_rejects_invalid(bad):
    args = [1e-09, 1000000.0, 1e-10, 0.5, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.steady_state_slope(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, -1), (4, 0)])
def test_steady_state_slope_physical_bound(index, value):
    args = [1e-09, 1000000.0, 1e-10, 0.5, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.steady_state_slope(*args)


def test_knickpoint_celerity_reference():
    result = domain.knickpoint_celerity(*[1000000.0, 0.01, 1e-10, 0.5, 1])
    np.testing.assert_allclose(result, 1e-7, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_knickpoint_celerity_rejects_invalid(bad):
    args = [1000000.0, 0.01, 1e-10, 0.5, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.knickpoint_celerity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1), (3, -1), (4, 0)])
def test_knickpoint_celerity_physical_bound(index, value):
    args = [1000000.0, 0.01, 1e-10, 0.5, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.knickpoint_celerity(*args)


def test_flint_slope_reference():
    result = domain.flint_slope(*[1000000.0, 10, 0.5])
    np.testing.assert_allclose(result, 0.01, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_flint_slope_rejects_invalid(bad):
    args = [1000000.0, 10, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flint_slope(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, -1)])
def test_flint_slope_physical_bound(index, value):
    args = [1000000.0, 10, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.flint_slope(*args)


def test_chi_coordinate_reference():
    result = domain.chi_coordinate(*[[0, 1000, 2000], [4000000.0, 1000000.0, 250000.0], 1000000.0, 0.5])
    np.testing.assert_allclose(result, [0, 750, 2250], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_chi_coordinate_rejects_invalid(bad):
    args = [[0, 1000, 2000], [4000000.0, 1000000.0, 250000.0], 1000000.0, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.chi_coordinate(*args)


@pytest.mark.parametrize("index,value", [(0, [0, 0, 1]), (0, 5), (1, [1, 1]), (1, [1, 0, 1]), (2, 0), (2, [1, 1]), (3, -1)])
def test_chi_coordinate_physical_bound(index, value):
    args = [[0, 1000, 2000], [4000000.0, 1000000.0, 250000.0], 1000000.0, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.chi_coordinate(*args)


def test_hack_length_reference():
    result = domain.hack_length(*[1000000.0, 0.32, 0.6])
    np.testing.assert_allclose(result, 0.32 * 1e6**0.6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_hack_length_rejects_invalid(bad):
    args = [1000000.0, 0.32, 0.6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hack_length(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_hack_length_physical_bound(index, value):
    args = [1000000.0, 0.32, 0.6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.hack_length(*args)


def test_strahler_order_reference():
    result = domain.strahler_order(*[[2, 2, 4, 4, -1, 3]])
    np.testing.assert_allclose(result, [1, 1, 2, 1, 2, 1], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_strahler_order_rejects_invalid(bad):
    args = [[2, 2, 4, 4, -1, 3]]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.strahler_order(*args)


@pytest.mark.parametrize("index,value", [(0, [1, 0]), (0, [0, -1]), (0, [1.5, -1]), (0, [3, -1]), (0, [-2, -1]), (0, [[1, -1]])])
def test_strahler_order_physical_bound(index, value):
    args = [[2, 2, 4, 4, -1, 3]]
    args[index] = value
    with pytest.raises(ValueError):
        domain.strahler_order(*args)


def test_downstream_hydraulic_geometry_reference():
    result = domain.downstream_hydraulic_geometry(*[100, 4, 0.5, 0.25, 0.4])
    np.testing.assert_allclose(result, [40, 0.25 * 100**0.4, 100 / (40 * 0.25 * 100**0.4)], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_downstream_hydraulic_geometry_rejects_invalid(bad):
    args = [100, 4, 0.5, 0.25, 0.4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.downstream_hydraulic_geometry(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1), (3, 0), (4, -1)])
def test_downstream_hydraulic_geometry_physical_bound(index, value):
    args = [100, 4, 0.5, 0.25, 0.4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.downstream_hydraulic_geometry(*args)


def test_effective_discharge_reference():
    result = domain.effective_discharge(*[50, 1, 2])
    np.testing.assert_allclose(result, 50 * np.exp(1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_effective_discharge_rejects_invalid(bad):
    args = [50, 1, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.effective_discharge(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1), (2, -1)])
def test_effective_discharge_physical_bound(index, value):
    args = [50, 1, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.effective_discharge(*args)


def test_threshold_shear_velocity_reference():
    result = domain.threshold_shear_velocity(*[0.00025, 2650, 1.2, 9.81, 0.1])
    np.testing.assert_allclose(result, 0.1 * np.sqrt(2648.8 / 1.2 * 9.81 * 2.5e-4), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_threshold_shear_velocity_rejects_invalid(bad):
    args = [0.00025, 2650, 1.2, 9.81, 0.1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.threshold_shear_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 1.2), (2, 0), (3, 0), (4, 0)])
def test_threshold_shear_velocity_physical_bound(index, value):
    args = [0.00025, 2650, 1.2, 9.81, 0.1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.threshold_shear_velocity(*args)


def test_bagnold_sand_flux_reference():
    result = domain.bagnold_sand_flux(*[0.5, 0.00025, 1.2, 9.81, 1.8])
    np.testing.assert_allclose(result, 1.8 * 1.2 / 9.81 * 0.125, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_bagnold_sand_flux_rejects_invalid(bad):
    args = [0.5, 0.00025, 1.2, 9.81, 1.8]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bagnold_sand_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0), (4, 0)])
def test_bagnold_sand_flux_physical_bound(index, value):
    args = [0.5, 0.00025, 1.2, 9.81, 1.8]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bagnold_sand_flux(*args)


def test_kawamura_sand_flux_reference():
    result = domain.kawamura_sand_flux(*[0.5, 0.2, 1.2, 9.81, 2.78])
    np.testing.assert_allclose(result, 2.78 * 1.2 / 9.81 * 0.3 * 0.49, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_kawamura_sand_flux_rejects_invalid(bad):
    args = [0.5, 0.2, 1.2, 9.81, 2.78]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.kawamura_sand_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0), (4, 0)])
def test_kawamura_sand_flux_physical_bound(index, value):
    args = [0.5, 0.2, 1.2, 9.81, 2.78]
    args[index] = value
    with pytest.raises(ValueError):
        domain.kawamura_sand_flux(*args)


def test_dune_migration_rate_reference():
    result = domain.dune_migration_rate(*[0.016, 10, 1600])
    np.testing.assert_allclose(result, 1e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dune_migration_rate_rejects_invalid(bad):
    args = [0.016, 10, 1600]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dune_migration_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_dune_migration_rate_physical_bound(index, value):
    args = [0.016, 10, 1600]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dune_migration_rate(*args)


def test_cosmogenic_production_at_depth_reference():
    result = domain.cosmogenic_production_at_depth(*[10, 1, 2000, 1600])
    np.testing.assert_allclose(result, 10 * np.exp(-1.25), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_cosmogenic_production_at_depth_rejects_invalid(bad):
    args = [10, 1, 2000, 1600]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.cosmogenic_production_at_depth(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, 0), (3, 0)])
def test_cosmogenic_production_at_depth_physical_bound(index, value):
    args = [10, 1, 2000, 1600]
    args[index] = value
    with pytest.raises(ValueError):
        domain.cosmogenic_production_at_depth(*args)


def test_cosmogenic_concentration_reference():
    result = domain.cosmogenic_concentration(*[0.1, 10, 0.05, 0, 2000, 1600])
    np.testing.assert_allclose(result, 10 / 0.05 * (1 - np.exp(-0.005)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_cosmogenic_concentration_rejects_invalid(bad):
    args = [0.1, 10, 0.05, 0, 2000, 1600]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.cosmogenic_concentration(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, -1), (3, -1), (4, 0), (5, 0)])
def test_cosmogenic_concentration_physical_bound(index, value):
    args = [0.1, 10, 0.05, 0, 2000, 1600]
    args[index] = value
    with pytest.raises(ValueError):
        domain.cosmogenic_concentration(*args)


def test_steady_state_erosion_rate_reference():
    result = domain.steady_state_erosion_rate(*[100, 10, 0.05, 2000, 1600])
    np.testing.assert_allclose(result, 0.04, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_steady_state_erosion_rate_rejects_invalid(bad):
    args = [100, 10, 0.05, 2000, 1600]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.steady_state_erosion_rate(*args)


@pytest.mark.parametrize("index,value", [(0, 300), (1, 0), (2, -1), (3, 0), (4, 0)])
def test_steady_state_erosion_rate_physical_bound(index, value):
    args = [100, 10, 0.05, 2000, 1600]
    args[index] = value
    with pytest.raises(ValueError):
        domain.steady_state_erosion_rate(*args)


def test_morphologic_age_reference():
    result = domain.morphologic_age(*[2, 0.6, 0.1])
    np.testing.assert_allclose(result, 4 / (np.pi * 0.25), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_morphologic_age_rejects_invalid(bad):
    args = [2, 0.6, 0.1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.morphologic_age(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0.1)])
def test_morphologic_age_physical_bound(index, value):
    args = [2, 0.6, 0.1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.morphologic_age(*args)


def test_luminescence_age_reference():
    result = domain.luminescence_age(*[50, 2e-11])
    np.testing.assert_allclose(result, 2.5e12, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_luminescence_age_rejects_invalid(bad):
    args = [50, 2e-11]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.luminescence_age(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_luminescence_age_physical_bound(index, value):
    args = [50, 2e-11]
    args[index] = value
    with pytest.raises(ValueError):
        domain.luminescence_age(*args)
