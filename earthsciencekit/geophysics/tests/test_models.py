"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import geophysics as domain


def test_sphere_gravity_reference():
    result = domain.sphere_gravity(*[3, 4, 1, 447111334967.41797])
    np.testing.assert_allclose(result, [-3, 4], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_sphere_gravity_rejects_invalid(bad):
    args = [3, 4, 1, 447111334967.41797]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.sphere_gravity(*args)


def test_sphere_gravity_broadcasting():
    args = [3, 4, 1, 447111334967.41797]
    reference = domain.sphere_gravity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.sphere_gravity(*args)
    assert result.shape == (2, 1, 2)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(2, 0)])
def test_sphere_gravity_physical_bound(index, value):
    args = [3, 4, 1, 447111334967.41797]
    args[index] = value
    with pytest.raises(ValueError):
        domain.sphere_gravity(*args)


def test_bouguer_correction_reference():
    result = domain.bouguer_correction(*[2384593786.492896, 3])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_bouguer_correction_rejects_invalid(bad):
    args = [2384593786.492896, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bouguer_correction(*args)


def test_bouguer_correction_broadcasting():
    args = [2384593786.492896, 3]
    reference = domain.bouguer_correction(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.bouguer_correction(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0)])
def test_bouguer_correction_physical_bound(index, value):
    args = [2384593786.492896, 3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bouguer_correction(*args)


def test_free_air_correction_reference():
    result = domain.free_air_correction(*[100, 9.8, 6400000])
    np.testing.assert_allclose(result, 2 * 9.8 * 100 / 6.4e6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_free_air_correction_rejects_invalid(bad):
    args = [100, 9.8, 6400000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.free_air_correction(*args)


def test_free_air_correction_broadcasting():
    args = [100, 9.8, 6400000]
    reference = domain.free_air_correction(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.free_air_correction(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(1, 0), (2, 0)])
def test_free_air_correction_physical_bound(index, value):
    args = [100, 9.8, 6400000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.free_air_correction(*args)


def test_bouguer_anomaly_reference():
    result = domain.bouguer_anomaly(*[9.81, 9.8, 100, 2384.5937864928956, 6400000])
    np.testing.assert_allclose(result, 0.01 + 3.0625e-4 - 1e-4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_bouguer_anomaly_rejects_invalid(bad):
    args = [9.81, 9.8, 100, 2384.5937864928956, 6400000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.bouguer_anomaly(*args)


def test_bouguer_anomaly_broadcasting():
    args = [9.81, 9.8, 100, 2384.5937864928956, 6400000]
    reference = domain.bouguer_anomaly(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.bouguer_anomaly(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (3, 0), (4, 0)])
def test_bouguer_anomaly_physical_bound(index, value):
    args = [9.81, 9.8, 100, 2384.5937864928956, 6400000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.bouguer_anomaly(*args)


def test_torsion_balance_gravitational_constant_reference():
    result = domain.torsion_balance_gravitational_constant(*[1, 2, 3, 4, 5])
    np.testing.assert_allclose(result, 0.24 * np.pi**2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_torsion_balance_gravitational_constant_rejects_invalid(bad):
    args = [1, 2, 3, 4, 5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.torsion_balance_gravitational_constant(*args)


def test_torsion_balance_gravitational_constant_broadcasting():
    args = [1, 2, 3, 4, 5]
    reference = domain.torsion_balance_gravitational_constant(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.torsion_balance_gravitational_constant(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
def test_torsion_balance_gravitational_constant_physical_bound(index, value):
    args = [1, 2, 3, 4, 5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.torsion_balance_gravitational_constant(*args)


def test_earth_mean_density_reference():
    result = domain.earth_mean_density(*[12.566370614359172, 3, 1])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_earth_mean_density_rejects_invalid(bad):
    args = [12.566370614359172, 3, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.earth_mean_density(*args)


def test_earth_mean_density_broadcasting():
    args = [12.566370614359172, 3, 1]
    reference = domain.earth_mean_density(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.earth_mean_density(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0), (2, 0)])
def test_earth_mean_density_physical_bound(index, value):
    args = [12.566370614359172, 3, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.earth_mean_density(*args)


def test_centrifugal_ratio_reference():
    result = domain.centrifugal_ratio(*[2, 3, 4])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_centrifugal_ratio_rejects_invalid(bad):
    args = [2, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.centrifugal_ratio(*args)


def test_centrifugal_ratio_broadcasting():
    args = [2, 3, 4]
    reference = domain.centrifugal_ratio(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.centrifugal_ratio(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_centrifugal_ratio_physical_bound(index, value):
    args = [2, 3, 4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.centrifugal_ratio(*args)


def test_homogeneous_flattening_reference():
    result = domain.homogeneous_flattening(*[0.08])
    np.testing.assert_allclose(result, 0.1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_homogeneous_flattening_rejects_invalid(bad):
    args = [0.08]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.homogeneous_flattening(*args)


def test_homogeneous_flattening_broadcasting():
    args = [0.08]
    reference = domain.homogeneous_flattening(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.homogeneous_flattening(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -0.01), (0, 0.2)])
def test_homogeneous_flattening_physical_bound(index, value):
    args = [0.08]
    args[index] = value
    with pytest.raises(ValueError):
        domain.homogeneous_flattening(*args)


def test_clairaut_gravity_reference():
    result = domain.clairaut_gravity(*[0.5235987755982988, 10, 0.003, 0.004])
    np.testing.assert_allclose(result, 10.0175, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_clairaut_gravity_rejects_invalid(bad):
    args = [0.5235987755982988, 10, 0.003, 0.004]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.clairaut_gravity(*args)


def test_clairaut_gravity_broadcasting():
    args = [0.5235987755982988, 10, 0.003, 0.004]
    reference = domain.clairaut_gravity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.clairaut_gravity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2), (1, 0), (2, -0.1), (2, 1), (3, -1)])
def test_clairaut_gravity_physical_bound(index, value):
    args = [0.5235987755982988, 10, 0.003, 0.004]
    args[index] = value
    with pytest.raises(ValueError):
        domain.clairaut_gravity(*args)


def test_normal_gravity_reference():
    result = domain.normal_gravity(*[1.5707963267948966])
    np.testing.assert_allclose(result, 9.8321863685, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_normal_gravity_rejects_invalid(bad):
    args = [1.5707963267948966]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.normal_gravity(*args)


def test_normal_gravity_broadcasting():
    args = [1.5707963267948966]
    reference = domain.normal_gravity(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.normal_gravity(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2), (0, -2)])
def test_normal_gravity_physical_bound(index, value):
    args = [1.5707963267948966]
    args[index] = value
    with pytest.raises(ValueError):
        domain.normal_gravity(*args)


def test_flattening_from_j2_reference():
    result = domain.flattening_from_j2(*[0.002, 0.004])
    np.testing.assert_allclose(result, 0.005, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_flattening_from_j2_rejects_invalid(bad):
    args = [0.002, 0.004]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.flattening_from_j2(*args)


def test_flattening_from_j2_broadcasting():
    args = [0.002, 0.004]
    reference = domain.flattening_from_j2(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.flattening_from_j2(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(1, -1)])
def test_flattening_from_j2_physical_bound(index, value):
    args = [0.002, 0.004]
    args[index] = value
    with pytest.raises(ValueError):
        domain.flattening_from_j2(*args)


def test_nodal_precession_rate_reference():
    result = domain.nodal_precession_rate(*[2, 0, 0, 1, 1, 8])
    np.testing.assert_allclose(result, -0.375, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_nodal_precession_rate_rejects_invalid(bad):
    args = [2, 0, 0, 1, 1, 8]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.nodal_precession_rate(*args)


def test_nodal_precession_rate_broadcasting():
    args = [2, 0, 0, 1, 1, 8]
    reference = domain.nodal_precession_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.nodal_precession_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, -0.1), (1, 1), (2, -1), (2, 4), (4, 0), (5, 0)])
def test_nodal_precession_rate_physical_bound(index, value):
    args = [2, 0, 0, 1, 1, 8]
    args[index] = value
    with pytest.raises(ValueError):
        domain.nodal_precession_rate(*args)


def test_stokes_function_reference():
    result = domain.stokes_function(*[3.141592653589793])
    np.testing.assert_allclose(result, 1 + 3 * np.log(2), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_stokes_function_rejects_invalid(bad):
    args = [3.141592653589793]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.stokes_function(*args)


def test_stokes_function_broadcasting():
    args = [3.141592653589793]
    reference = domain.stokes_function(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.stokes_function(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (0, 4)])
def test_stokes_function_physical_bound(index, value):
    args = [3.141592653589793]
    args[index] = value
    with pytest.raises(ValueError):
        domain.stokes_function(*args)


def test_degree_geoid_height_reference():
    result = domain.degree_geoid_height(*[6, 4, 2, 4])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_degree_geoid_height_rejects_invalid(bad):
    args = [6, 4, 2, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.degree_geoid_height(*args)


def test_degree_geoid_height_broadcasting():
    args = [6, 4, 2, 4]
    reference = domain.degree_geoid_height(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.degree_geoid_height(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(1, 1), (1, 2.5), (2, 0), (3, 0)])
def test_degree_geoid_height_physical_bound(index, value):
    args = [6, 4, 2, 4]
    args[index] = value
    with pytest.raises(ValueError):
        domain.degree_geoid_height(*args)


def test_dipole_field_reference():
    result = domain.dipole_field(*[0, 1, 10000000])
    np.testing.assert_allclose(result, [1, 0], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_dipole_field_rejects_invalid(bad):
    args = [0, 1, 10000000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dipole_field(*args)


def test_dipole_field_broadcasting():
    args = [0, 1, 10000000]
    reference = domain.dipole_field(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.dipole_field(*args)
    assert result.shape == (2, 1, 2)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2), (1, 0)])
def test_dipole_field_physical_bound(index, value):
    args = [0, 1, 10000000]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dipole_field(*args)


def test_dipole_inclination_reference():
    result = domain.dipole_inclination(*[0.7853981633974483])
    np.testing.assert_allclose(result, np.arctan(2), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_dipole_inclination_rejects_invalid(bad):
    args = [0.7853981633974483]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.dipole_inclination(*args)


def test_dipole_inclination_broadcasting():
    args = [0.7853981633974483]
    reference = domain.dipole_inclination(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.dipole_inclination(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2)])
def test_dipole_inclination_physical_bound(index, value):
    args = [0.7853981633974483]
    args[index] = value
    with pytest.raises(ValueError):
        domain.dipole_inclination(*args)


def test_paleolatitude_reference():
    result = domain.paleolatitude(*[1.1071487177940906])
    np.testing.assert_allclose(result, np.pi / 4, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_paleolatitude_rejects_invalid(bad):
    args = [1.1071487177940906]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.paleolatitude(*args)


def test_paleolatitude_broadcasting():
    args = [1.1071487177940906]
    reference = domain.paleolatitude(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.paleolatitude(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -2)])
def test_paleolatitude_physical_bound(index, value):
    args = [1.1071487177940906]
    args[index] = value
    with pytest.raises(ValueError):
        domain.paleolatitude(*args)


def test_zonal_field_reference():
    result = domain.zonal_field(*[0.5235987755982988, 1, 2, 1, 1])
    np.testing.assert_allclose(result, [-0.75 * np.sqrt(3), 0.375], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_zonal_field_rejects_invalid(bad):
    args = [0.5235987755982988, 1, 2, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.zonal_field(*args)


def test_zonal_field_broadcasting():
    args = [0.5235987755982988, 1, 2, 1, 1]
    reference = domain.zonal_field(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.zonal_field(*args)
    assert result.shape == (2, 1, 2)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2), (1, 0), (2, 0), (2, 1.5), (4, 0)])
def test_zonal_field_physical_bound(index, value):
    args = [0.5235987755982988, 1, 2, 1, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.zonal_field(*args)


def test_gauss_dipole_moment_reference():
    result = domain.gauss_dipole_moment(*[3, 4, 12, 1])
    np.testing.assert_allclose(result, 4 * np.pi * 13 / 1.2566370614359173e-6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_gauss_dipole_moment_rejects_invalid(bad):
    args = [3, 4, 12, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.gauss_dipole_moment(*args)


def test_gauss_dipole_moment_broadcasting():
    args = [3, 4, 12, 1]
    reference = domain.gauss_dipole_moment(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.gauss_dipole_moment(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(3, 0)])
def test_gauss_dipole_moment_physical_bound(index, value):
    args = [3, 4, 12, 1]
    args[index] = value
    with pytest.raises(ValueError):
        domain.gauss_dipole_moment(*args)


def test_field_elements_reference():
    result = domain.field_elements(*[1, 1, 1.4142135623730951])
    np.testing.assert_allclose(result, [np.pi / 4, np.pi / 4, 2], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_field_elements_rejects_invalid(bad):
    args = [1, 1, 1.4142135623730951]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.field_elements(*args)


def test_field_elements_broadcasting():
    args = [1, 1, 1.4142135623730951]
    reference = domain.field_elements(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.field_elements(*args)
    assert result.shape == (2, 1, 3)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


def test_direction_angle_reference():
    result = domain.direction_angle(*[0, 0.5, 0, -0.5])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_direction_angle_rejects_invalid(bad):
    args = [0, 0.5, 0, -0.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.direction_angle(*args)


def test_direction_angle_broadcasting():
    args = [0, 0.5, 0, -0.5]
    reference = domain.direction_angle(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.direction_angle(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(1, 2), (3, -2)])
def test_direction_angle_physical_bound(index, value):
    args = [0, 0.5, 0, -0.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.direction_angle(*args)


def test_virtual_geomagnetic_pole_reference():
    result = domain.virtual_geomagnetic_pole(*[0, 0, 1.5707963267948966, 0])
    np.testing.assert_allclose(result, [0, np.pi / 2], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_virtual_geomagnetic_pole_rejects_invalid(bad):
    args = [0, 0, 1.5707963267948966, 0]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.virtual_geomagnetic_pole(*args)


def test_virtual_geomagnetic_pole_broadcasting():
    args = [0, 0, 1.5707963267948966, 0]
    reference = domain.virtual_geomagnetic_pole(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.virtual_geomagnetic_pole(*args)
    assert result.shape == (2, 1, 2)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 2), (3, 2)])
def test_virtual_geomagnetic_pole_physical_bound(index, value):
    args = [0, 0, 1.5707963267948966, 0]
    args[index] = value
    with pytest.raises(ValueError):
        domain.virtual_geomagnetic_pole(*args)


def test_block_magnetic_anomaly_reference():
    result = domain.block_magnetic_anomaly(*[0, -1, 1, 1, 1.7320508075688772, 9549296.58551372])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_block_magnetic_anomaly_rejects_invalid(bad):
    args = [0, -1, 1, 1, 1.7320508075688772, 9549296.58551372]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.block_magnetic_anomaly(*args)


def test_block_magnetic_anomaly_broadcasting():
    args = [0, -1, 1, 1, 1.7320508075688772, 9549296.58551372]
    reference = domain.block_magnetic_anomaly(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.block_magnetic_anomaly(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(3, 0), (4, 0)])
def test_block_magnetic_anomaly_physical_bound(index, value):
    args = [0, -1, 1, 1, 1.7320508075688772, 9549296.58551372]
    args[index] = value
    with pytest.raises(ValueError):
        domain.block_magnetic_anomaly(*args)


def test_magnetic_reynolds_number_reference():
    result = domain.magnetic_reynolds_number(*[2, 3, 795774.7154594767])
    np.testing.assert_allclose(result, 6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_magnetic_reynolds_number_rejects_invalid(bad):
    args = [2, 3, 795774.7154594767]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.magnetic_reynolds_number(*args)


def test_magnetic_reynolds_number_broadcasting():
    args = [2, 3, 795774.7154594767]
    reference = domain.magnetic_reynolds_number(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.magnetic_reynolds_number(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0)])
def test_magnetic_reynolds_number_physical_bound(index, value):
    args = [2, 3, 795774.7154594767]
    args[index] = value
    with pytest.raises(ValueError):
        domain.magnetic_reynolds_number(*args)


def test_magnetic_decay_time_reference():
    result = domain.magnetic_decay_time(*[3.141592653589793, 795774.7154594767])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_magnetic_decay_time_rejects_invalid(bad):
    args = [3.141592653589793, 795774.7154594767]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.magnetic_decay_time(*args)


def test_magnetic_decay_time_broadcasting():
    args = [3.141592653589793, 795774.7154594767]
    reference = domain.magnetic_decay_time(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.magnetic_decay_time(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_magnetic_decay_time_physical_bound(index, value):
    args = [3.141592653589793, 795774.7154594767]
    args[index] = value
    with pytest.raises(ValueError):
        domain.magnetic_decay_time(*args)


def test_rikitake_fixed_point_reference():
    result = domain.rikitake_fixed_point(*[3, 2])
    np.testing.assert_allclose(result, [2, 0.5, 12], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_rikitake_fixed_point_rejects_invalid(bad):
    args = [3, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rikitake_fixed_point(*args)


def test_rikitake_fixed_point_broadcasting():
    args = [3, 2]
    reference = domain.rikitake_fixed_point(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.rikitake_fixed_point(*args)
    assert result.shape == (2, 1, 3)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, 0), (1, 0)])
def test_rikitake_fixed_point_physical_bound(index, value):
    args = [3, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.rikitake_fixed_point(*args)


def test_love_number_h2_reference():
    result = domain.love_number_h2(*[1.5, 1, 1, 9.5])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_love_number_h2_rejects_invalid(bad):
    args = [1.5, 1, 1, 9.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.love_number_h2(*args)


def test_love_number_h2_broadcasting():
    args = [1.5, 1, 1, 9.5]
    reference = domain.love_number_h2(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.love_number_h2(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_love_number_h2_physical_bound(index, value):
    args = [1.5, 1, 1, 9.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.love_number_h2(*args)


def test_love_number_k2_reference():
    result = domain.love_number_k2(*[1.5, 1, 1, 9.5])
    np.testing.assert_allclose(result, 0.6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_love_number_k2_rejects_invalid(bad):
    args = [1.5, 1, 1, 9.5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.love_number_k2(*args)


def test_love_number_k2_broadcasting():
    args = [1.5, 1, 1, 9.5]
    reference = domain.love_number_k2(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.love_number_k2(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -1), (1, 0), (2, 0), (3, 0)])
def test_love_number_k2_physical_bound(index, value):
    args = [1.5, 1, 1, 9.5]
    args[index] = value
    with pytest.raises(ValueError):
        domain.love_number_k2(*args)


def test_tidal_gravimetric_factor_reference():
    result = domain.tidal_gravimetric_factor(*[0.6, 0.3])
    np.testing.assert_allclose(result, 1.15, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_tidal_gravimetric_factor_rejects_invalid(bad):
    args = [0.6, 0.3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.tidal_gravimetric_factor(*args)


def test_tidal_gravimetric_factor_broadcasting():
    args = [0.6, 0.3]
    reference = domain.tidal_gravimetric_factor(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.tidal_gravimetric_factor(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(0, -1), (1, -1)])
def test_tidal_gravimetric_factor_physical_bound(index, value):
    args = [0.6, 0.3]
    args[index] = value
    with pytest.raises(ValueError):
        domain.tidal_gravimetric_factor(*args)


def test_equivalent_water_height_reference():
    result = domain.equivalent_water_height(*[2, 1, 0.5, 3, 3, 2])
    np.testing.assert_allclose(result, 6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_equivalent_water_height_rejects_invalid(bad):
    args = [2, 1, 0.5, 3, 3, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.equivalent_water_height(*args)


def test_equivalent_water_height_broadcasting():
    args = [2, 1, 0.5, 3, 3, 2]
    reference = domain.equivalent_water_height(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.equivalent_water_height(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, np.broadcast_to(reference, result.shape))


@pytest.mark.parametrize("index,value", [(1, -1), (1, 0.5), (2, -1), (3, 0), (4, 0), (5, 0)])
def test_equivalent_water_height_physical_bound(index, value):
    args = [2, 1, 0.5, 3, 3, 2]
    args[index] = value
    with pytest.raises(ValueError):
        domain.equivalent_water_height(*args)
