"""Golden checks and invalid-input contracts for every public model."""

import numpy as np
import pytest

from earthsciencekit import seismology as domain


def test_p_wave_speed_reference():
    result = domain.p_wave_speed(*[12, 3, 4])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_p_wave_speed_rejects_invalid(bad):
    args = [12, 3, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.p_wave_speed(*args)


def test_s_wave_speed_reference():
    result = domain.s_wave_speed(*[16, 4])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_s_wave_speed_rejects_invalid(bad):
    args = [16, 4]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.s_wave_speed(*args)


def test_poisson_ratio_reference():
    result = domain.poisson_ratio(*[2, 1])
    np.testing.assert_allclose(result, 1 / 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_poisson_ratio_rejects_invalid(bad):
    args = [2, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.poisson_ratio(*args)


def test_travel_time_reference():
    result = domain.travel_time(*[3000, 1500])
    np.testing.assert_allclose(result, 2, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_travel_time_rejects_invalid(bad):
    args = [3000, 1500]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.travel_time(*args)


def test_seismic_impedance_reference():
    result = domain.seismic_impedance(*[2000, 3000])
    np.testing.assert_allclose(result, 6e6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_seismic_impedance_rejects_invalid(bad):
    args = [2000, 3000]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.seismic_impedance(*args)


def test_reflection_coefficient_reference():
    result = domain.reflection_coefficient(*[2, 6])
    np.testing.assert_allclose(result, 0.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_reflection_coefficient_rejects_invalid(bad):
    args = [2, 6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.reflection_coefficient(*args)


def test_transmission_coefficient_reference():
    result = domain.transmission_coefficient(*[2, 6])
    np.testing.assert_allclose(result, 1.5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_transmission_coefficient_rejects_invalid(bad):
    args = [2, 6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.transmission_coefficient(*args)


def test_critical_angle_reference():
    result = domain.critical_angle(*[1, 2])
    np.testing.assert_allclose(result, np.pi / 6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_critical_angle_rejects_invalid(bad):
    args = [1, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.critical_angle(*args)


def test_seismic_moment_reference():
    result = domain.seismic_moment(*[30000000000.0, 1000000.0, 1])
    np.testing.assert_allclose(result, 3e16, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_seismic_moment_rejects_invalid(bad):
    args = [30000000000.0, 1000000.0, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.seismic_moment(*args)


def test_moment_magnitude_reference():
    result = domain.moment_magnitude(*[1.2589254117941714e18])
    np.testing.assert_allclose(result, 6, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_moment_magnitude_rejects_invalid(bad):
    args = [1.2589254117941714e18]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.moment_magnitude(*args)


def test_moment_from_magnitude_reference():
    result = domain.moment_from_magnitude(*[6])
    np.testing.assert_allclose(result, 10**18.1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_moment_from_magnitude_rejects_invalid(bad):
    args = [6]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.moment_from_magnitude(*args)


def test_gutenberg_richter_reference():
    result = domain.gutenberg_richter(*[4, 1, 3])
    np.testing.assert_allclose(result, 10, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_gutenberg_richter_rejects_invalid(bad):
    args = [4, 1, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.gutenberg_richter(*args)


def test_omori_rate_reference():
    result = domain.omori_rate(*[10, 1, 1, 1])
    np.testing.assert_allclose(result, 5, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_omori_rate_rejects_invalid(bad):
    args = [10, 1, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.omori_rate(*args)


def test_ricker_wavelet_reference():
    result = domain.ricker_wavelet(*[0, 10])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ricker_wavelet_rejects_invalid(bad):
    args = [0, 10]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ricker_wavelet(*args)


def test_attenuated_amplitude_reference():
    result = domain.attenuated_amplitude(*[2, 1, 1, np.float64(4.532360141827194)])
    np.testing.assert_allclose(result, 1, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_attenuated_amplitude_rejects_invalid(bad):
    args = [2, 1, 1, np.float64(4.532360141827194)]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.attenuated_amplitude(*args)


def test_p_wave_speed_broadcasting():
    args = [12, 3, 4]
    reference = domain.p_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.p_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_p_wave_speed_bulk_modulus_physical_bound():
    args = [12, 3, 4]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.p_wave_speed(*args)


def test_p_wave_speed_shear_modulus_physical_bound():
    args = [12, 3, 4]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.p_wave_speed(*args)


def test_p_wave_speed_density_physical_bound():
    args = [12, 3, 4]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.p_wave_speed(*args)


def test_s_wave_speed_broadcasting():
    args = [16, 4]
    reference = domain.s_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.s_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_s_wave_speed_shear_modulus_physical_bound():
    args = [16, 4]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.s_wave_speed(*args)


def test_s_wave_speed_density_physical_bound():
    args = [16, 4]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.s_wave_speed(*args)


def test_poisson_ratio_broadcasting():
    args = [2, 1]
    reference = domain.poisson_ratio(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.poisson_ratio(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_poisson_ratio_p_speed_physical_bound():
    args = [2, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.poisson_ratio(*args)


def test_poisson_ratio_s_speed_physical_bound():
    args = [2, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.poisson_ratio(*args)


def test_travel_time_broadcasting():
    args = [3000, 1500]
    reference = domain.travel_time(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.travel_time(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_travel_time_distance_physical_bound():
    args = [3000, 1500]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.travel_time(*args)


def test_travel_time_speed_physical_bound():
    args = [3000, 1500]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.travel_time(*args)


def test_seismic_impedance_broadcasting():
    args = [2000, 3000]
    reference = domain.seismic_impedance(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.seismic_impedance(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_seismic_impedance_density_physical_bound():
    args = [2000, 3000]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.seismic_impedance(*args)


def test_seismic_impedance_speed_physical_bound():
    args = [2000, 3000]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.seismic_impedance(*args)


def test_reflection_coefficient_broadcasting():
    args = [2, 6]
    reference = domain.reflection_coefficient(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.reflection_coefficient(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_reflection_coefficient_impedance1_physical_bound():
    args = [2, 6]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.reflection_coefficient(*args)


def test_reflection_coefficient_impedance2_physical_bound():
    args = [2, 6]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.reflection_coefficient(*args)


def test_transmission_coefficient_broadcasting():
    args = [2, 6]
    reference = domain.transmission_coefficient(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.transmission_coefficient(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_transmission_coefficient_impedance1_physical_bound():
    args = [2, 6]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.transmission_coefficient(*args)


def test_transmission_coefficient_impedance2_physical_bound():
    args = [2, 6]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.transmission_coefficient(*args)


def test_critical_angle_broadcasting():
    args = [1, 2]
    reference = domain.critical_angle(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.critical_angle(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_critical_angle_speed1_physical_bound():
    args = [1, 2]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.critical_angle(*args)


def test_critical_angle_speed2_physical_bound():
    args = [1, 2]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.critical_angle(*args)


def test_seismic_moment_broadcasting():
    args = [30000000000.0, 1000000.0, 1]
    reference = domain.seismic_moment(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.seismic_moment(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_seismic_moment_rigidity_physical_bound():
    args = [30000000000.0, 1000000.0, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.seismic_moment(*args)


def test_seismic_moment_area_physical_bound():
    args = [30000000000.0, 1000000.0, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.seismic_moment(*args)


def test_seismic_moment_slip_physical_bound():
    args = [30000000000.0, 1000000.0, 1]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.seismic_moment(*args)


def test_moment_magnitude_broadcasting():
    args = [1.2589254117941714e18]
    reference = domain.moment_magnitude(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.moment_magnitude(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_moment_magnitude_moment_physical_bound():
    args = [1.2589254117941714e18]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.moment_magnitude(*args)


def test_moment_from_magnitude_broadcasting():
    args = [6]
    reference = domain.moment_from_magnitude(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.moment_from_magnitude(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_gutenberg_richter_broadcasting():
    args = [4, 1, 3]
    reference = domain.gutenberg_richter(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.gutenberg_richter(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_gutenberg_richter_b_physical_bound():
    args = [4, 1, 3]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.gutenberg_richter(*args)


def test_omori_rate_broadcasting():
    args = [10, 1, 1, 1]
    reference = domain.omori_rate(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.omori_rate(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_omori_rate_productivity_physical_bound():
    args = [10, 1, 1, 1]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.omori_rate(*args)


def test_omori_rate_offset_physical_bound():
    args = [10, 1, 1, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.omori_rate(*args)


def test_omori_rate_exponent_physical_bound():
    args = [10, 1, 1, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.omori_rate(*args)


def test_omori_rate_time_physical_bound():
    args = [10, 1, 1, 1]
    args[3] = -1
    with pytest.raises(ValueError):
        domain.omori_rate(*args)


def test_ricker_wavelet_broadcasting():
    args = [0, 10]
    reference = domain.ricker_wavelet(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.ricker_wavelet(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_ricker_wavelet_frequency_physical_bound():
    args = [0, 10]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.ricker_wavelet(*args)


def test_attenuated_amplitude_broadcasting():
    args = [2, 1, 1, 4.532360141827194]
    reference = domain.attenuated_amplitude(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.attenuated_amplitude(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_attenuated_amplitude_frequency_physical_bound():
    args = [2, 1, 1, 4.532360141827194]
    args[1] = -1
    with pytest.raises(ValueError):
        domain.attenuated_amplitude(*args)


def test_attenuated_amplitude_time_physical_bound():
    args = [2, 1, 1, 4.532360141827194]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.attenuated_amplitude(*args)


def test_attenuated_amplitude_quality_factor_physical_bound():
    args = [2, 1, 1, 4.532360141827194]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.attenuated_amplitude(*args)


def test_rayleigh_wave_speed_reference():
    result = domain.rayleigh_wave_speed(*[1.7320508075688772, 1])
    np.testing.assert_allclose(result, np.sqrt(2 - 2 / np.sqrt(3)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_rayleigh_wave_speed_rejects_invalid(bad):
    args = [1.7320508075688772, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.rayleigh_wave_speed(*args)


def test_rayleigh_wave_speed_broadcasting():
    args = [1.7320508075688772, 1]
    reference = domain.rayleigh_wave_speed(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.rayleigh_wave_speed(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_rayleigh_wave_speed_p_speed_physical_bound():
    args = [1.7320508075688772, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.rayleigh_wave_speed(*args)


def test_rayleigh_wave_speed_s_speed_physical_bound():
    args = [1.7320508075688772, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.rayleigh_wave_speed(*args)


def test_head_wave_travel_time_reference():
    result = domain.head_wave_travel_time(*[4, 1, 1, 2])
    np.testing.assert_allclose(result, 2 + np.sqrt(3), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_head_wave_travel_time_rejects_invalid(bad):
    args = [4, 1, 1, 2]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.head_wave_travel_time(*args)


def test_head_wave_travel_time_broadcasting():
    args = [4, 1, 1, 2]
    reference = domain.head_wave_travel_time(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.head_wave_travel_time(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_head_wave_travel_time_distance_physical_bound():
    args = [4, 1, 1, 2]
    args[0] = -1
    with pytest.raises(ValueError):
        domain.head_wave_travel_time(*args)


def test_head_wave_travel_time_thickness_physical_bound():
    args = [4, 1, 1, 2]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.head_wave_travel_time(*args)


def test_head_wave_travel_time_speed1_physical_bound():
    args = [4, 1, 1, 2]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.head_wave_travel_time(*args)


def test_head_wave_travel_time_speed2_physical_bound():
    args = [4, 1, 1, 2]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.head_wave_travel_time(*args)


def test_crossover_distance_reference():
    result = domain.crossover_distance(*[1, 1, 3])
    np.testing.assert_allclose(result, 2 * np.sqrt(2), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_crossover_distance_rejects_invalid(bad):
    args = [1, 1, 3]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.crossover_distance(*args)


def test_crossover_distance_broadcasting():
    args = [1, 1, 3]
    reference = domain.crossover_distance(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.crossover_distance(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_crossover_distance_thickness_physical_bound():
    args = [1, 1, 3]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.crossover_distance(*args)


def test_crossover_distance_speed1_physical_bound():
    args = [1, 1, 3]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.crossover_distance(*args)


def test_crossover_distance_speed2_physical_bound():
    args = [1, 1, 3]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.crossover_distance(*args)


def test_love_wave_angular_frequency_reference():
    result = domain.love_wave_angular_frequency(*[1.4142135623730951, 1, 1, 2, 1, 1])
    np.testing.assert_allclose(result, np.sqrt(2) * np.arctan(np.sqrt(0.5)), rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_love_wave_angular_frequency_rejects_invalid(bad):
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_broadcasting():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    reference = domain.love_wave_angular_frequency(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.love_wave_angular_frequency(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_love_wave_angular_frequency_phase_speed_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_thickness_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_layer_s_speed_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[2] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_halfspace_s_speed_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[3] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_layer_rigidity_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[4] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_love_wave_angular_frequency_halfspace_rigidity_physical_bound():
    args = [1.4142135623730951, 1, 1, 2, 1, 1]
    args[5] = 0
    with pytest.raises(ValueError):
        domain.love_wave_angular_frequency(*args)


def test_ray_parameter_reference():
    result = domain.ray_parameter(*[2, 4, 0.5235987755982988])
    np.testing.assert_allclose(result, 0.25, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_ray_parameter_rejects_invalid(bad):
    args = [2, 4, 0.5235987755982988]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.ray_parameter(*args)


def test_ray_parameter_broadcasting():
    args = [2, 4, 0.5235987755982988]
    reference = domain.ray_parameter(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.ray_parameter(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_ray_parameter_radius_physical_bound():
    args = [2, 4, 0.5235987755982988]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.ray_parameter(*args)


def test_ray_parameter_speed_physical_bound():
    args = [2, 4, 0.5235987755982988]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.ray_parameter(*args)


def test_ray_parameter_incidence_angle_physical_bound():
    args = [2, 4, 0.5235987755982988]
    args[2] = -1
    with pytest.raises(ValueError):
        domain.ray_parameter(*args)


def test_local_magnitude_reference():
    result = domain.local_magnitude(*[1e-3, 1e5])
    np.testing.assert_allclose(result, 3, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_local_magnitude_rejects_invalid(bad):
    args = [1e-3, 1e5]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.local_magnitude(*args)


def test_local_magnitude_broadcasting():
    args = [1e-3, 1e5]
    reference = domain.local_magnitude(*args)
    args[0] = np.full((2, 1), args[0])
    result = domain.local_magnitude(*args)
    assert result.shape == (2, 1)
    np.testing.assert_allclose(result, reference)


def test_local_magnitude_amplitude_physical_bound():
    args = [1e-3, 1e5]
    args[0] = 0
    with pytest.raises(ValueError):
        domain.local_magnitude(*args)


def test_local_magnitude_distance_physical_bound():
    args = [1e-3, 1e5]
    args[1] = 0
    with pytest.raises(ValueError):
        domain.local_magnitude(*args)


def test_ray_parameter_incidence_angle_upper_bound():
    with pytest.raises(ValueError):
        domain.ray_parameter(2, 4, 2)


def test_vertical_travel_time_reference():
    result = domain.vertical_travel_time(*[[0, 1], [1, np.e]])
    np.testing.assert_allclose(result, [0, 1 / (np.e - 1)], rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, [], "bad", True, 1j])
def test_vertical_travel_time_rejects_invalid(bad):
    args = [[0, 1], [1, np.e]]
    args[0] = bad
    with pytest.raises(ValueError):
        domain.vertical_travel_time(*args)


def test_vertical_travel_time_depth_physical_bound():
    with pytest.raises(ValueError):
        domain.vertical_travel_time([-1, 1], [1, 2])


def test_vertical_travel_time_speed_physical_bound():
    with pytest.raises(ValueError):
        domain.vertical_travel_time([0, 1], [1, 0])


@pytest.mark.parametrize("depth,speed", [([0], [1]), ([[0, 1]], [[1, 2]]), ([0, 1], [1, 2, 3]), ([0, 2, 1], [1, 2, 3])])
def test_vertical_travel_time_requires_profile(depth, speed):
    with pytest.raises(ValueError):
        domain.vertical_travel_time(depth, speed)
