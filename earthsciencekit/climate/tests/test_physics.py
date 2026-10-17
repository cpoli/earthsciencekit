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
