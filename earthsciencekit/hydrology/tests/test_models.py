"""Golden checks and invalid-input contracts for every public model."""

from typing import Any

import numpy as np
import pytest

from earthsciencekit import hydrology as domain

SLOPE_0C = 611.2 * 2.5e6 / (461.5 * 273.15**2)  # Pa/K, Clausius-Clapeyron slope at 0 degC
GAMMA_SEA_LEVEL = 1004 * 101325 * 461.5 / (287.05 * 2.5e6)  # Pa/K, psychrometric constant at 101325 Pa
BAD = [np.nan, np.inf, -np.inf, [], "bad", True, 1j]


def test_water_balance_reference():
    result = domain.water_balance(*[3, 1, 0.5])
    np.testing.assert_allclose(result, 1.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_water_balance_rejects_invalid(bad):
    args = [3, 1, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.water_balance(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, -1)])
def test_water_balance_physical_bound(index, value):
    args = [3, 1, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.water_balance(*args)


def test_runoff_coefficient_reference():
    result = domain.runoff_coefficient(*[1, 6])
    np.testing.assert_allclose(result, 1 / 6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_runoff_coefficient_rejects_invalid(bad):
    args = [1, 6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.runoff_coefficient(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_runoff_coefficient_physical_bound(index, value):
    args = [1, 6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.runoff_coefficient(*args)


def test_hydraulic_radius_reference():
    result = domain.hydraulic_radius(*[4, 1])
    np.testing.assert_allclose(result, 2 / 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_hydraulic_radius_rejects_invalid(bad):
    args = [4, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hydraulic_radius(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_hydraulic_radius_physical_bound(index, value):
    args = [4, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.hydraulic_radius(*args)


def test_chezy_velocity_reference():
    result = domain.chezy_velocity(*[50, 2, 0.02])
    np.testing.assert_allclose(result, 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_chezy_velocity_rejects_invalid(bad):
    args = [50, 2, 0.02]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.chezy_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1)])
def test_chezy_velocity_physical_bound(index, value):
    args = [50, 2, 0.02]
    args[index] = value
    with pytest.raises(ValueError):
        domain.chezy_velocity(*args)


def test_dalton_evaporation_reference():
    result = domain.dalton_evaporation(*[3000, 1000, 1e-9])
    np.testing.assert_allclose(result, 2e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dalton_evaporation_rejects_invalid(bad):
    args = [3000, 1000, 1e-9]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dalton_evaporation(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, -1)])
def test_dalton_evaporation_physical_bound(index, value):
    args = [3000, 1000, 1e-9]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dalton_evaporation(*args)


def test_hydraulic_head_reference():
    result = domain.hydraulic_head(*[10, 9810, 1000, 9.81])
    np.testing.assert_allclose(result, 11, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_hydraulic_head_rejects_invalid(bad):
    args = [10, 9810, 1000, 9.81]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.hydraulic_head(*args)


@pytest.mark.parametrize("index,value", [(2, 0), (3, 0)])
def test_hydraulic_head_physical_bound(index, value):
    args = [10, 9810, 1000, 9.81]
    args[index] = value
    with pytest.raises(ValueError):
        domain.hydraulic_head(*args)


def test_darcy_flux_reference():
    result = domain.darcy_flux(*[1e-4, -0.01])
    np.testing.assert_allclose(result, 1e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_darcy_flux_rejects_invalid(bad):
    args = [1e-4, -0.01]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.darcy_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1)])
def test_darcy_flux_physical_bound(index, value):
    args = [1e-4, -0.01]
    args[index] = value
    with pytest.raises(ValueError):
        domain.darcy_flux(*args)


def test_seepage_velocity_reference():
    result = domain.seepage_velocity(*[1e-6, 0.25])
    np.testing.assert_allclose(result, 4e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_seepage_velocity_rejects_invalid(bad):
    args = [1e-6, 0.25]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.seepage_velocity(*args)


@pytest.mark.parametrize("index,value", [(1, 0), (1, 1.5)])
def test_seepage_velocity_physical_bound(index, value):
    args = [1e-6, 0.25]
    args[index] = value
    with pytest.raises(ValueError):
        domain.seepage_velocity(*args)


def test_dupuit_head_reference():
    result = domain.dupuit_head(*[50, 5, 3, 100])
    np.testing.assert_allclose(result, np.sqrt(17), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dupuit_head_rejects_invalid(bad):
    args = [50, 5, 3, 100]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dupuit_head(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (0, 101), (1, -1), (2, -1), (3, 0)])
def test_dupuit_head_physical_bound(index, value):
    args = [50, 5, 3, 100]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dupuit_head(*args)


def test_dupuit_discharge_reference():
    result = domain.dupuit_discharge(*[2, 5, 3, 4])
    np.testing.assert_allclose(result, 4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dupuit_discharge_rejects_invalid(bad):
    args = [2, 5, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dupuit_discharge(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1), (2, -1), (3, 0)])
def test_dupuit_discharge_physical_bound(index, value):
    args = [2, 5, 3, 4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dupuit_discharge(*args)


def test_manning_velocity_reference():
    result = domain.manning_velocity(*[0.04, 8, 0.0004])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_manning_velocity_rejects_invalid(bad):
    args = [0.04, 8, 0.0004]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.manning_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1)])
def test_manning_velocity_physical_bound(index, value):
    args = [0.04, 8, 0.0004]
    args[index] = value
    with pytest.raises(ValueError):
        domain.manning_velocity(*args)


def test_chezy_from_manning_reference():
    result = domain.chezy_from_manning(*[0.05, 64])
    np.testing.assert_allclose(result, 40, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_chezy_from_manning_rejects_invalid(bad):
    args = [0.05, 64]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.chezy_from_manning(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_chezy_from_manning_physical_bound(index, value):
    args = [0.05, 64]
    args[index] = value
    with pytest.raises(ValueError):
        domain.chezy_from_manning(*args)


def test_green_ampt_rate_reference():
    result = domain.green_ampt_rate(*[0.5, 2, 1, 0.5])
    np.testing.assert_allclose(result, 4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_green_ampt_rate_rejects_invalid(bad):
    args = [0.5, 2, 1, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.green_ampt_rate(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1.5)])
def test_green_ampt_rate_physical_bound(index, value):
    args = [0.5, 2, 1, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.green_ampt_rate(*args)


def test_green_ampt_time_reference():
    result = domain.green_ampt_time(*[np.e - 1, 1, 1, 1])
    np.testing.assert_allclose(result, np.e - 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_green_ampt_time_rejects_invalid(bad):
    args = [np.e - 1, 1, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.green_ampt_time(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0), (3, 1.5)])
def test_green_ampt_time_physical_bound(index, value):
    args = [np.e - 1, 1, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.green_ampt_time(*args)


def test_gardner_conductivity_reference():
    result = domain.gardner_conductivity(*[-1, 2, np.log(2)])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_gardner_conductivity_rejects_invalid(bad):
    args = [-1, 2, np.log(2)]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.gardner_conductivity(*args)


@pytest.mark.parametrize("index,value", [(0, 1), (1, 0), (2, 0)])
def test_gardner_conductivity_physical_bound(index, value):
    args = [-1, 2, np.log(2)]
    args[index] = value
    with pytest.raises(ValueError):
        domain.gardner_conductivity(*args)


def test_unsaturated_darcy_flux_reference():
    result = domain.unsaturated_darcy_flux(*[2, -3])
    np.testing.assert_allclose(result, 4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_unsaturated_darcy_flux_rejects_invalid(bad):
    args = [2, -3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.unsaturated_darcy_flux(*args)


@pytest.mark.parametrize("index,value", [(0, -1)])
def test_unsaturated_darcy_flux_physical_bound(index, value):
    args = [2, -3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.unsaturated_darcy_flux(*args)


def test_richards_steady_head_reference():
    result = domain.richards_steady_head(*[np.log(2), -0.5, 1, 1])
    np.testing.assert_allclose(result, np.log(0.75), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_richards_steady_head_rejects_invalid(bad):
    args = [np.log(2), -0.5, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.richards_steady_head(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -2), (1, 1), (2, 0), (3, 0)])
def test_richards_steady_head_physical_bound(index, value):
    args = [np.log(2), -0.5, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.richards_steady_head(*args)


def test_nash_unit_hydrograph_reference():
    result = domain.nash_unit_hydrograph(*[2, 2, 1])
    np.testing.assert_allclose(result, 2 * np.exp(-2), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_nash_unit_hydrograph_rejects_invalid(bad):
    args = [2, 2, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nash_unit_hydrograph(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0.5), (2, 0)])
def test_nash_unit_hydrograph_physical_bound(index, value):
    args = [2, 2, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nash_unit_hydrograph(*args)


def test_unit_hydrograph_convolution_reference():
    result = domain.unit_hydrograph_convolution(*[[1, 2], [3, 4, 5]])
    np.testing.assert_allclose(result, [3, 10, 13, 10], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_unit_hydrograph_convolution_rejects_invalid(bad):
    args = [[1, 2], [3, 4, 5]]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.unit_hydrograph_convolution(*args)


@pytest.mark.parametrize("index,value", [(0, [-1, 2]), (1, [3, -4, 5]), (0, [[1, 2]])])
def test_unit_hydrograph_convolution_physical_bound(index, value):
    args = [[1, 2], [3, 4, 5]]
    args[index] = value
    with pytest.raises(ValueError):
        domain.unit_hydrograph_convolution(*args)


def test_horton_infiltration_rate_reference():
    result = domain.horton_infiltration_rate(*[np.log(2), 5, 1, 1])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_horton_infiltration_rate_rejects_invalid(bad):
    args = [np.log(2), 5, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.horton_infiltration_rate(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0.5), (2, -1), (3, 0)])
def test_horton_infiltration_rate_physical_bound(index, value):
    args = [np.log(2), 5, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.horton_infiltration_rate(*args)


def test_horton_cumulative_infiltration_reference():
    result = domain.horton_cumulative_infiltration(*[np.log(2), 5, 1, 1])
    np.testing.assert_allclose(result, np.log(2) + 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_horton_cumulative_infiltration_rejects_invalid(bad):
    args = [np.log(2), 5, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.horton_cumulative_infiltration(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0.5), (2, -1), (3, 0)])
def test_horton_cumulative_infiltration_physical_bound(index, value):
    args = [np.log(2), 5, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.horton_cumulative_infiltration(*args)


def test_well_function_reference():
    result = domain.well_function(*[1])
    np.testing.assert_allclose(result, 0.21938393439552027, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_well_function_rejects_invalid(bad):
    args = [1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.well_function(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_well_function_physical_bound(index, value):
    args = [1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.well_function(*args)


def test_theis_drawdown_reference():
    result = domain.theis_drawdown(*[4 * np.pi, 1, 4, 1, 1])
    np.testing.assert_allclose(result, 0.21938393439552027, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_theis_drawdown_rejects_invalid(bad):
    args = [4 * np.pi, 1, 4, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.theis_drawdown(*args)


@pytest.mark.parametrize("index,value", [(1, 0), (2, 0), (3, 0), (4, 0)])
def test_theis_drawdown_physical_bound(index, value):
    args = [4 * np.pi, 1, 4, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.theis_drawdown(*args)


def test_muskingum_coefficients_reference():
    result = domain.muskingum_coefficients(*[1, 0.25, 1])
    np.testing.assert_allclose(result, [0.2, 0.6, 0.2], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_muskingum_coefficients_rejects_invalid(bad):
    args = [1, 0.25, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.muskingum_coefficients(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -0.1), (1, 0.6), (2, 0)])
def test_muskingum_coefficients_physical_bound(index, value):
    args = [1, 0.25, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.muskingum_coefficients(*args)


def test_muskingum_routing_reference():
    result = domain.muskingum_routing(*[[1, 3, 2], 1, 0.25, 1])
    np.testing.assert_allclose(result, [1, 1.4, 2.48], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_muskingum_routing_rejects_invalid(bad):
    args = [[1, 3, 2], 1, 0.25, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.muskingum_routing(*args)


@pytest.mark.parametrize("index,value", [(0, [-1, 3, 2]), (0, [1]), (1, 0), (1, [1, 2]), (2, 0.6), (3, 0)])
def test_muskingum_routing_physical_bound(index, value):
    args = [[1, 3, 2], 1, 0.25, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.muskingum_routing(*args)


def test_stream_number_reference():
    result = domain.stream_number(*[1, 3, 4])
    np.testing.assert_allclose(result, 16, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stream_number_rejects_invalid(bad):
    args = [1, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stream_number(*args)


@pytest.mark.parametrize("index,value", [(0, 0.5), (0, 4), (1, 0.5), (2, 0)])
def test_stream_number_physical_bound(index, value):
    args = [1, 3, 4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stream_number(*args)


def test_stream_length_reference():
    result = domain.stream_length(*[3, 100, 2])
    np.testing.assert_allclose(result, 400, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_stream_length_rejects_invalid(bad):
    args = [3, 100, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stream_length(*args)


@pytest.mark.parametrize("index,value", [(0, 0.5), (1, 0), (2, 0)])
def test_stream_length_physical_bound(index, value):
    args = [3, 100, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stream_length(*args)


def test_cooper_jacob_drawdown_reference():
    result = domain.cooper_jacob_drawdown(*[4 * np.pi, 1, 4, 1, 1])
    np.testing.assert_allclose(result, -np.euler_gamma, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_cooper_jacob_drawdown_rejects_invalid(bad):
    args = [4 * np.pi, 1, 4, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.cooper_jacob_drawdown(*args)


@pytest.mark.parametrize("index,value", [(1, 0), (2, 0), (3, 0), (4, 0)])
def test_cooper_jacob_drawdown_physical_bound(index, value):
    args = [4 * np.pi, 1, 4, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.cooper_jacob_drawdown(*args)


def test_cooper_jacob_transmissivity_reference():
    result = domain.cooper_jacob_transmissivity(*[4 * np.pi, np.log(10)])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_cooper_jacob_transmissivity_rejects_invalid(bad):
    args = [4 * np.pi, np.log(10)]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.cooper_jacob_transmissivity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_cooper_jacob_transmissivity_physical_bound(index, value):
    args = [4 * np.pi, np.log(10)]
    args[index] = value
    with pytest.raises(ValueError):
        domain.cooper_jacob_transmissivity(*args)


def test_saturation_slope_reference():
    result = domain.saturation_slope(*[273.15])
    np.testing.assert_allclose(result, SLOPE_0C, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_saturation_slope_rejects_invalid(bad):
    args = [273.15]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.saturation_slope(*args)


@pytest.mark.parametrize("index,value", [(0, 100), (0, 400)])
def test_saturation_slope_physical_bound(index, value):
    args = [273.15]
    args[index] = value
    with pytest.raises(ValueError):
        domain.saturation_slope(*args)


def test_psychrometric_constant_reference():
    result = domain.psychrometric_constant(*[101325])
    np.testing.assert_allclose(result, GAMMA_SEA_LEVEL, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_psychrometric_constant_rejects_invalid(bad):
    args = [101325]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.psychrometric_constant(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_psychrometric_constant_physical_bound(index, value):
    args = [101325]
    args[index] = value
    with pytest.raises(ValueError):
        domain.psychrometric_constant(*args)


def test_penman_evaporation_reference():
    result = domain.penman_evaporation(*[0, 273.15, 1e-5, 101325])
    np.testing.assert_allclose(result, GAMMA_SEA_LEVEL * 1e-5 / (SLOPE_0C + GAMMA_SEA_LEVEL), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_penman_evaporation_rejects_invalid(bad):
    args = [0, 273.15, 1e-5, 101325]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.penman_evaporation(*args)


@pytest.mark.parametrize("index,value", [(1, 100), (3, 0)])
def test_penman_evaporation_physical_bound(index, value):
    args = [0, 273.15, 1e-5, 101325]
    args[index] = value
    with pytest.raises(ValueError):
        domain.penman_evaporation(*args)


def test_aerodynamic_resistance_reference():
    result = domain.aerodynamic_resistance(*[1 / 0.41**2, 2, 0, 2 / np.e, 2 / np.e**2])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_aerodynamic_resistance_rejects_invalid(bad):
    args = [1 / 0.41**2, 2, 0, 2 / np.e, 2 / np.e**2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.aerodynamic_resistance(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1), (2, 1.9), (3, 0), (4, 0)])
def test_aerodynamic_resistance_physical_bound(index, value):
    args = [1 / 0.41**2, 2, 0, 2 / np.e, 2 / np.e**2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.aerodynamic_resistance(*args)


def test_bulk_surface_resistance_reference():
    result = domain.bulk_surface_resistance(*[100, 2])
    np.testing.assert_allclose(result, 50, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_bulk_surface_resistance_rejects_invalid(bad):
    args = [100, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bulk_surface_resistance(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_bulk_surface_resistance_physical_bound(index, value):
    args = [100, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bulk_surface_resistance(*args)


def test_penman_monteith_evaporation_reference():
    result = domain.penman_monteith_evaporation(*[0, 273.15, 1000, 1, 1004, 1004, 101325])
    np.testing.assert_allclose(result, 1000 / (2.5e6 * (SLOPE_0C + 2 * GAMMA_SEA_LEVEL)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_penman_monteith_evaporation_rejects_invalid(bad):
    args = [0, 273.15, 1000, 1, 1004, 1004, 101325]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.penman_monteith_evaporation(*args)


@pytest.mark.parametrize("index,value", [(1, 100), (3, 0), (4, 0), (5, -1), (6, 0)])
def test_penman_monteith_evaporation_physical_bound(index, value):
    args = [0, 273.15, 1000, 1, 1004, 1004, 101325]
    args[index] = value
    with pytest.raises(ValueError):
        domain.penman_monteith_evaporation(*args)


def test_topographic_index_reference():
    result = domain.topographic_index(*[np.e, np.pi / 4])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_topographic_index_rejects_invalid(bad):
    args = [np.e, np.pi / 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.topographic_index(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (1, np.pi / 2)])
def test_topographic_index_physical_bound(index, value):
    args = [np.e, np.pi / 4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.topographic_index(*args)


def test_topmodel_deficit_reference():
    result = domain.topmodel_deficit(*[0.1, 7, 6, 0.05])
    np.testing.assert_allclose(result, 0.05, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_topmodel_deficit_rejects_invalid(bad):
    args = [0.1, 7, 6, 0.05]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.topmodel_deficit(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (3, 0)])
def test_topmodel_deficit_physical_bound(index, value):
    args = [0.1, 7, 6, 0.05]
    args[index] = value
    with pytest.raises(ValueError):
        domain.topmodel_deficit(*args)


def test_topmodel_baseflow_reference():
    result = domain.topmodel_baseflow(*[2, 0, 1, 1])
    np.testing.assert_allclose(result, 2 / np.e, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_topmodel_baseflow_rejects_invalid(bad):
    args = [2, 0, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.topmodel_baseflow(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (2, -1), (3, 0)])
def test_topmodel_baseflow_physical_bound(index, value):
    args = [2, 0, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.topmodel_baseflow(*args)


_DEEPER_CASES: dict[str, tuple[list[Any], Any]] = {
    "exponential_water_content": ([-1, 0.4, 0.05, 2], 0.05 + 0.35 * np.exp(-2)),
    "kinematic_wave_celerity": ([1.5, 5 / 3], 2.5),
    "degree_day_melt": ([278.15, 2, 273.15], 10.0),
    "snowpack_water_equivalent": ([[270, 275], [0.01, 0.002], 1e-3, 273.15, 2], [[0.01, 0], [0.0063, 0.0057]]),
}


@pytest.mark.parametrize("name", sorted(_DEEPER_CASES))
def test_deeper_hydrology_reference(name):
    args, expected = _DEEPER_CASES[name]
    np.testing.assert_allclose(getattr(domain, name)(*args), expected, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
@pytest.mark.parametrize("name", sorted(_DEEPER_CASES))
def test_deeper_hydrology_rejects_invalid(name, bad):
    args = list(_DEEPER_CASES[name][0])
    args[0] = bad
    with pytest.raises(ValueError):
        getattr(domain, name)(*args)


@pytest.mark.parametrize(
    "name,args",
    [
        ("exponential_water_content", [1, 0.4, 0.05, 2]),
        ("exponential_water_content", [-1, 0.04, 0.05, 2]),
        ("exponential_water_content", [-1, 1.2, 0.05, 2]),
        ("exponential_water_content", [-1, 0.4, -0.05, 2]),
        ("exponential_water_content", [-1, 0.4, 0.05, 0]),
        ("kinematic_wave_celerity", [-1, 5 / 3]),
        ("kinematic_wave_celerity", [1, 0]),
        ("degree_day_melt", [0, 2, 273.15]),
        ("degree_day_melt", [278.15, -2, 273.15]),
        ("degree_day_melt", [278.15, 2, 0]),
        ("snowpack_water_equivalent", [[270, 275], [0.01], 1e-3, 273.15, 2]),
        ("snowpack_water_equivalent", [[270, 275], [0.01, -0.01], 1e-3, 273.15, 2]),
        ("snowpack_water_equivalent", [[270, 275], [0.01, 0.002], [1e-3, 1e-3], 273.15, 2]),
        ("snowpack_water_equivalent", [[270, 275], [0.01, 0.002], 1e-3, 273.15, [2, 2]]),
        ("snowpack_water_equivalent", [[270, 275], [0.01, 0.002], 1e-3, 273.15, 0]),
    ],
)
def test_deeper_hydrology_physical_bound(name, args):
    with pytest.raises(ValueError):
        getattr(domain, name)(*args)
