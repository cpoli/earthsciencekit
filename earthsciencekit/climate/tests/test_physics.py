"""Physical identities, limits and dynamic checks for earthsciencekit.climate."""

import numpy as np
import pytest

from earthsciencekit.climate.systems.response import simulate_one_box, step_temperature


def test_one_box_exact_solution_and_nonzero_start():
    time = np.linspace(10, 20, 51)
    result = simulate_one_box(4, 2, 2, time, initial_anomaly=-1)
    np.testing.assert_allclose(result.values[:, 0], 2 - 3 * np.exp(-(time - time[0])), atol=2e-8)
    zero = simulate_one_box(4, 2, 2, time)
    np.testing.assert_allclose(zero.values[:, 0], step_temperature(4, 2, 2, time - time[0]), atol=2e-8)
    with pytest.raises(ValueError):
        simulate_one_box(4, [2, 3], 2, time)


YEAR = 3.15576e7


def test_two_box_matches_matrix_exponential_and_timescales():
    from scipy.linalg import expm

    from earthsciencekit.climate import simulate_two_box, two_box_fast_warming, two_box_timescales

    forcing, feedback, exchange, surface, deep = 3.7, 1.3, 0.7, 8 * YEAR, 100 * YEAR
    time = np.linspace(0, 500, 101) * YEAR
    result = simulate_two_box(forcing, feedback, exchange, surface, deep, time, initial_anomaly=[0.5, 0.2])
    matrix = np.array([[-(feedback + exchange) / surface, exchange / surface], [exchange / deep, -exchange / deep]])
    equilibrium = np.full(2, forcing / feedback)
    expected = [equilibrium + expm(matrix * t) @ (np.array([0.5, 0.2]) - equilibrium) for t in time]
    np.testing.assert_allclose(result.values, expected, atol=1e-6)
    fast, slow = two_box_timescales(feedback, exchange, surface, deep)
    np.testing.assert_allclose(sorted(-1 / np.linalg.eigvals(matrix)), [fast, slow], rtol=1e-10)
    # Held et al.: after the fast adjustment the surface sits near F / (lambda + gamma) while the deep ocean is cold.
    plateau = simulate_two_box(forcing, feedback, exchange, surface, 1e6 * YEAR, [0, 5 * fast]).values[-1, 0]
    assert np.isclose(plateau, two_box_fast_warming(forcing, feedback, exchange), rtol=0.01)
    assert np.isclose(fast, surface / (feedback + exchange), rtol=0.1) and np.isclose(slow, deep * 2 / (1.3 * 0.7), rtol=0.1)


@pytest.mark.parametrize(
    "args,kwargs",
    [(([3.7, 1], 1, 1, 1, 1), {}), ((3.7, 1, 0, 1, 1), {}), ((3.7, 1, 1, 1, 1), {"initial_anomaly": [0.0]})],
)
def test_two_box_rejects_bad_setup(args, kwargs):
    from earthsciencekit.climate import simulate_two_box

    with pytest.raises(ValueError):
        simulate_two_box(*args, [0, 1], **kwargs)


def test_energy_balance_without_ice_feedback_is_norths_two_mode_solution():
    from earthsciencekit.climate import north_temperature_profile, simulate_latitudinal_energy_balance

    kwargs = {"cold_albedo": 0.3, "warm_albedo": 0.3, "warm_albedo_gradient": 0.0, "nodes": 60}
    result = simulate_latitudinal_energy_balance(1361, np.linspace(0, 60 * YEAR, 3), **kwargs)
    x = (np.arange(60) + 0.5) / 60
    expected = north_temperature_profile(x, 1361, 0.3, 210, 2, 0.555, -0.482)
    np.testing.assert_allclose(result.values[-1, 2:], expected, atol=0.01)
    np.testing.assert_allclose(result.values[-1, 0], expected.mean(), atol=0.01)
    assert result.values[-1, 1] == np.pi / 2


def test_energy_balance_closes_globally_and_has_two_stable_states():
    from earthsciencekit.climate import ice_albedo, insolation_distribution, simulate_latitudinal_energy_balance

    time = np.linspace(0, 80 * YEAR, 3)
    warm = simulate_latitudinal_energy_balance(1361, time, nodes=40)
    x = (np.arange(40) + 0.5) / 40
    temperature = warm.values[-1, 2:]
    albedo = ice_albedo(temperature, 0.62, 0.3 + 0.078 * (3 * x**2 - 1) / 2, 263.15, 1.0)
    absorbed = 1361 / 4 * insolation_distribution(x, -0.482) * (1 - albedo)
    # Diffusion only moves heat, so at equilibrium the hemisphere absorbs what it emits.
    assert np.isclose(absorbed.mean(), np.mean(210 + 2 * (temperature - 273.15)), rtol=1e-6)
    edge = warm.values[-1, 1]
    assert np.radians(65) < edge < np.radians(78) and np.all(np.diff(temperature) < 0)
    i = np.searchsorted(-temperature, -263.15)
    assert np.arcsin(x[i - 1]) <= edge <= np.arcsin(x[i])
    snowball = simulate_latitudinal_energy_balance(1361, time, nodes=40, initial_temperature=230)
    assert snowball.values[-1, 1] == 0 and snowball.values[-1, 0] < warm.values[-1, 0] - 40
    assert snowball.names[:3] == ("global_mean_temperature", "ice_edge", "temperature_0") and snowball.units[1] == "rad"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"nodes": 2},
        {"nodes": True},
        {"nodes": 4.0},
        {"initial_temperature": [280.0, 270.0]},
        {"diffusion": [0.5, 0.6]},
        {"warm_albedo_gradient": 0.7},
        {"warm_albedo_gradient": [0.1, 0.1]},
    ],
)
def test_energy_balance_rejects_bad_setup(kwargs):
    from earthsciencekit.climate import simulate_latitudinal_energy_balance

    with pytest.raises(ValueError):
        simulate_latitudinal_energy_balance(1361, [0, 1], **kwargs)
