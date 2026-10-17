"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import volcanology as domain
from earthsciencekit.constants import AVOGADRO, LAMBDA_U238, R_GAS

HALF_LIFE_U238 = np.log(2) / LAMBDA_U238
BAD = [np.nan, np.inf, -np.inf, [], "bad", True, 1j]


def test_decay_constant_reference():
    result = domain.decay_constant(*[2])
    np.testing.assert_allclose(result, np.log(2) / 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_decay_constant_rejects_invalid(bad):
    args = [2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.decay_constant(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_decay_constant_physical_bound(index, value):
    args = [2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.decay_constant(*args)


def test_parent_fraction_reference():
    result = domain.parent_fraction(*[np.log(2), 1])
    np.testing.assert_allclose(result, 0.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_parent_fraction_rejects_invalid(bad):
    args = [np.log(2), 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.parent_fraction(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, -1)])
def test_parent_fraction_physical_bound(index, value):
    args = [np.log(2), 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.parent_fraction(*args)


def test_radiometric_age_reference():
    result = domain.radiometric_age(*[np.e - 1, 2])
    np.testing.assert_allclose(result, 0.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_radiometric_age_rejects_invalid(bad):
    args = [np.e - 1, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.radiometric_age(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_radiometric_age_physical_bound(index, value):
    args = [np.e - 1, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.radiometric_age(*args)


def test_radiogenic_daughter_reference():
    result = domain.radiogenic_daughter(*[3, np.log(2), 1])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_radiogenic_daughter_rejects_invalid(bad):
    args = [3, np.log(2), 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.radiogenic_daughter(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, -1)])
def test_radiogenic_daughter_physical_bound(index, value):
    args = [3, np.log(2), 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.radiogenic_daughter(*args)


def test_radiocarbon_age_reference():
    result = domain.radiocarbon_age(*[np.exp(-1)])
    np.testing.assert_allclose(result, 8033 * 365.25 * 86400, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_radiocarbon_age_rejects_invalid(bad):
    args = [np.exp(-1)]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.radiocarbon_age(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_radiocarbon_age_physical_bound(index, value):
    args = [np.exp(-1)]
    args[index] = value
    with pytest.raises(ValueError):
        domain.radiocarbon_age(*args)


def test_radiogenic_lead_ratio_reference():
    result = domain.radiogenic_lead_ratio(*[HALF_LIFE_U238])
    np.testing.assert_allclose(result, (2 ** (9.8485 / 1.55125) - 1) / 137.88, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_radiogenic_lead_ratio_rejects_invalid(bad):
    args = [HALF_LIFE_U238]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.radiogenic_lead_ratio(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_radiogenic_lead_ratio_physical_bound(index, value):
    args = [HALF_LIFE_U238]
    args[index] = value
    with pytest.raises(ValueError):
        domain.radiogenic_lead_ratio(*args)


def test_lead_lead_age_reference():
    result = domain.lead_lead_age(*[(2 ** (9.8485 / 1.55125) - 1) / 137.88])
    np.testing.assert_allclose(result, HALF_LIFE_U238, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_lead_lead_age_rejects_invalid(bad):
    args = [(2 ** (9.8485 / 1.55125) - 1) / 137.88]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.lead_lead_age(*args)


@pytest.mark.parametrize("index,value", [(0, 0.04), (0, 50.0)])
def test_lead_lead_age_physical_bound(index, value):
    args = [(2 ** (9.8485 / 1.55125) - 1) / 137.88]
    args[index] = value
    with pytest.raises(ValueError):
        domain.lead_lead_age(*args)


def test_concordia_ratios_reference():
    result = domain.concordia_ratios(*[HALF_LIFE_U238])
    np.testing.assert_allclose(result, [2 ** (9.8485 / 1.55125) - 1, 1], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_concordia_ratios_rejects_invalid(bad):
    args = [HALF_LIFE_U238]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.concordia_ratios(*args)


@pytest.mark.parametrize("index,value", [(0, -1)])
def test_concordia_ratios_physical_bound(index, value):
    args = [HALF_LIFE_U238]
    args[index] = value
    with pytest.raises(ValueError):
        domain.concordia_ratios(*args)


def test_isochron_fit_reference():
    result = domain.isochron_fit(*[[0, 1, 2], [0.7, 0.7 + np.e - 1, 0.7 + 2 * (np.e - 1)], 1])
    np.testing.assert_allclose(result, [1, 0.7], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_isochron_fit_rejects_invalid(bad):
    args = [[0, 1, 2], [0.7, 0.7 + np.e - 1, 0.7 + 2 * (np.e - 1)], 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.isochron_fit(*args)


@pytest.mark.parametrize("index,value", [(0, [-1, 1, 2]), (1, [-1, 1, 2]), (2, 0)])
def test_isochron_fit_physical_bound(index, value):
    args = [[0, 1, 2], [0.7, 0.7 + np.e - 1, 0.7 + 2 * (np.e - 1)], 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.isochron_fit(*args)


def test_rayleigh_fractionation_reference():
    result = domain.rayleigh_fractionation(*[2, 0.25, 0.5])
    np.testing.assert_allclose(result, 4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_rayleigh_fractionation_rejects_invalid(bad):
    args = [2, 0.25, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rayleigh_fractionation(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (1, 1.5), (2, 0)])
def test_rayleigh_fractionation_physical_bound(index, value):
    args = [2, 0.25, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.rayleigh_fractionation(*args)


def test_delta_value_reference():
    result = domain.delta_value(*[1.01, 1])
    np.testing.assert_allclose(result, 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_delta_value_rejects_invalid(bad):
    args = [1.01, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.delta_value(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_delta_value_physical_bound(index, value):
    args = [1.01, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.delta_value(*args)


def test_fractional_crystallization_melt_reference():
    result = domain.fractional_crystallization_melt(*[10, 0.25, 0.5])
    np.testing.assert_allclose(result, 20, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_fractional_crystallization_melt_rejects_invalid(bad):
    args = [10, 0.25, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.fractional_crystallization_melt(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (1, 1.5), (2, -1)])
def test_fractional_crystallization_melt_physical_bound(index, value):
    args = [10, 0.25, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.fractional_crystallization_melt(*args)


def test_equilibrium_crystallization_melt_reference():
    result = domain.equilibrium_crystallization_melt(*[10, 0.5, 0.5])
    np.testing.assert_allclose(result, 40 / 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_equilibrium_crystallization_melt_rejects_invalid(bad):
    args = [10, 0.5, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.equilibrium_crystallization_melt(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, -0.5), (1, 1.5), (2, 0)])
def test_equilibrium_crystallization_melt_physical_bound(index, value):
    args = [10, 0.5, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.equilibrium_crystallization_melt(*args)


def test_lattice_strain_partition_reference():
    result = domain.lattice_strain_partition(*[2e-10, 1e-10, 3, 1.2 * R_GAS * 1000 / (4 * np.pi * AVOGADRO * 1e-30), 1000])
    np.testing.assert_allclose(result, 3 * np.exp(-1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_lattice_strain_partition_rejects_invalid(bad):
    args = [2e-10, 1e-10, 3, 1.2 * R_GAS * 1000 / (4 * np.pi * AVOGADRO * 1e-30), 1000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.lattice_strain_partition(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
def test_lattice_strain_partition_physical_bound(index, value):
    args = [2e-10, 1e-10, 3, 1.2 * R_GAS * 1000 / (4 * np.pi * AVOGADRO * 1e-30), 1000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.lattice_strain_partition(*args)


def test_isotope_temperature_reference():
    result = domain.isotope_temperature(*[1, 0])
    np.testing.assert_allclose(result, 273.15 + 16.5 - 4.3 + 0.14, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_isotope_temperature_rejects_invalid(bad):
    args = [1, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.isotope_temperature(*args)


def test_arrhenius_viscosity_reference():
    result = domain.arrhenius_viscosity(*[1000, 10, R_GAS * 1000, 500])
    np.testing.assert_allclose(result, 10 * np.exp(-1), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_arrhenius_viscosity_rejects_invalid(bad):
    args = [1000, 10, R_GAS * 1000, 500]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.arrhenius_viscosity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, -1), (3, 0)])
def test_arrhenius_viscosity_physical_bound(index, value):
    args = [1000, 10, R_GAS * 1000, 500]
    args[index] = value
    with pytest.raises(ValueError):
        domain.arrhenius_viscosity(*args)


def test_shaw_viscosity_reference():
    result = domain.shaw_viscosity(*[1e4 / 1.5, 2])
    np.testing.assert_allclose(result, 0.1 * np.exp(-6.4), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_shaw_viscosity_rejects_invalid(bad):
    args = [1e4 / 1.5, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.shaw_viscosity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_shaw_viscosity_physical_bound(index, value):
    args = [1e4 / 1.5, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.shaw_viscosity(*args)


def test_water_solubility_reference():
    result = domain.water_solubility(*[4e8, 4e-6])
    np.testing.assert_allclose(result, 0.08, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_water_solubility_rejects_invalid(bad):
    args = [4e8, 4e-6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.water_solubility(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0)])
def test_water_solubility_physical_bound(index, value):
    args = [4e8, 4e-6]
    args[index] = value
    with pytest.raises(ValueError):
        domain.water_solubility(*args)


def test_dike_ascent_velocity_reference():
    result = domain.dike_ascent_velocity(*[2, 300, 100, 10])
    np.testing.assert_allclose(result, 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_dike_ascent_velocity_rejects_invalid(bad):
    args = [2, 300, 100, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dike_ascent_velocity(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (2, 0), (3, 0)])
def test_dike_ascent_velocity_physical_bound(index, value):
    args = [2, 300, 100, 10]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dike_ascent_velocity(*args)


def test_mogi_displacement_reference():
    result = domain.mogi_displacement(*[0, 1000, np.pi, 0.5])
    np.testing.assert_allclose(result, [0, 5e-7], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_mogi_displacement_rejects_invalid(bad):
    args = [0, 1000, np.pi, 0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.mogi_displacement(*args)


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (3, -0.1), (3, 0.6)])
def test_mogi_displacement_physical_bound(index, value):
    args = [0, 1000, np.pi, 0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.mogi_displacement(*args)


def test_plume_rise_height_reference():
    result = domain.plume_rise_height(*[16, 1 / 16])
    np.testing.assert_allclose(result, 80, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_plume_rise_height_rejects_invalid(bad):
    args = [16, 1 / 16]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.plume_rise_height(*args)


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_plume_rise_height_physical_bound(index, value):
    args = [16, 1 / 16]
    args[index] = value
    with pytest.raises(ValueError):
        domain.plume_rise_height(*args)


def test_column_height_reference():
    result = domain.column_height(*[1e12])
    np.testing.assert_allclose(result, 8200, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_column_height_rejects_invalid(bad):
    args = [1e12]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.column_height(*args)


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_column_height_physical_bound(index, value):
    args = [1e12]
    args[index] = value
    with pytest.raises(ValueError):
        domain.column_height(*args)


def test_volcanic_explosivity_index_reference():
    result = domain.volcanic_explosivity_index(*[5e7])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", BAD)
def test_volcanic_explosivity_index_rejects_invalid(bad):
    args = [5e7]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.volcanic_explosivity_index(*args)


@pytest.mark.parametrize("index,value", [(0, -1)])
def test_volcanic_explosivity_index_physical_bound(index, value):
    args = [5e7]
    args[index] = value
    with pytest.raises(ValueError):
        domain.volcanic_explosivity_index(*args)
