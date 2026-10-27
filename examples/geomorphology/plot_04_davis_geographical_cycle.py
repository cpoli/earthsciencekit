"""
Davis's geographical cycle (1899)
=================================

William Morris Davis described landscapes as passing through youth,
maturity and old age: after a pulse of uplift, rivers cut down quickly,
then ever more slowly as their slopes decline, until only a low
"peneplain" remains. Here a graded river stops rising and stream power
wears it down. With a slope exponent of 1 the erosion sweeps upstream as a
wave and the relief falls almost linearly; with an exponent of 2, erosion
slows as the slopes flatten and the relief decays ever more slowly, much as
Davis drew it. Airy isostasy (``plate_tectonics.airy_root``) means far more
rock must be removed than the mean elevation loses.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import simulate_river_profile, steady_state_slope
from earthsciencekit.plate_tectonics import airy_root

distance = np.linspace(0, 50e3, 101)
area = ((50e3 - distance + 300) / 0.32) ** (1 / 0.6)
erodibility = 5.4e-6 / SECONDS_PER_YEAR
graded_slope = steady_state_slope(1e-3 / SECONDS_PER_YEAR, area[1:], erodibility, 0.5, 1)
uplifted = np.concatenate([[0.0], np.cumsum(graded_slope * np.diff(distance))])
years = np.linspace(0, 5e6, 101)
linear = simulate_river_profile(distance, area, 0, erodibility, years * SECONDS_PER_YEAR, initial_elevation=uplifted)
# The same initial erosion rate at the outlet with a slope exponent of 2.
steep_k = erodibility / (np.diff(uplifted)[0] / np.diff(distance)[0])
davis = simulate_river_profile(distance, area, 0, steep_k, years * SECONDS_PER_YEAR, slope_exponent=2, initial_elevation=uplifted)
lowering = davis.values[0, 1] - davis.values[-1, 1]
removed = lowering + airy_root(lowering, 2700, 3300)
assert np.all(np.diff(davis.values[:, 0]) < 0) and np.all(np.diff(davis.values[:, 0], 2) > 0)
assert linear.values[-1, 0] < 1 < 20 < davis.values[-1, 0]
assert np.isclose(removed / lowering, 3300 / 600)

fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 8))
for index, stage in [(0, "uplifted"), (5, "youth"), (20, "maturity"), (100, "old age")]:
    top.plot(distance / 1e3, davis.values[index, 2:], label=f"{stage}, {years[index] / 1e6:.1f} Myr")
top.set(xlabel="distance from the outlet [km]", ylabel="elevation [m]", title="Davis's cycle with a slope exponent of 2")
top.legend()
bottom.plot(years / 1e6, linear.values[:, 0], label="slope exponent 1: a wave of erosion")
bottom.plot(years / 1e6, davis.values[:, 0], label="slope exponent 2: ever slower decay")
title = f"{removed:.0f} m of rock removed, mean {lowering:.0f} m lower"
bottom.set(xlabel="time since uplift stopped [Myr]", ylabel="relief [m]", title=title)
bottom.legend()
fig.tight_layout()

plt.show()
