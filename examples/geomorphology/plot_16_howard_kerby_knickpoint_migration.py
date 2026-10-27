"""
Howard and Kerby's stream-power incision (1983)
===============================================

Alan Howard and Gordon Kerby measured channel lowering in the badlands of
Virginia and found it proportional to powers of drainage area and slope,
the stream-power law now used in most models of bedrock rivers. A drop of
base level, or here a doubling of the uplift rate, does not lower the whole
river at once: a knickpoint travels upstream at the speed
``K A**m S**(n - 1)``, fast on the trunk and slowly in the headwaters,
separating the new steep profile below from the old one above.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import knickpoint_celerity, simulate_river_profile, steady_state_slope

distance = np.linspace(0, 50e3, 201)
area = ((50e3 - distance + 300) / 0.32) ** (1 / 0.6)
uplift, erodibility = 1e-3 / SECONDS_PER_YEAR, 5.4e-6 / SECONDS_PER_YEAR
old_slope = steady_state_slope(uplift, area[1:], erodibility, 0.5, 1)
steady = np.concatenate([[0.0], np.cumsum(old_slope * np.diff(distance))])
years = np.linspace(0, 1.5e6, 61)
run = simulate_river_profile(distance, area, 2 * uplift, erodibility, years * SECONDS_PER_YEAR, initial_elevation=steady)
celerity = knickpoint_celerity(area[1:], old_slope, erodibility, 0.5, 1)
travel = np.concatenate([[0.0], np.cumsum(np.diff(distance) / celerity)]) / SECONDS_PER_YEAR
knick = np.array([distance[1:][np.argmax(np.diff(row) / np.diff(distance) / old_slope < 1.5)] for row in run.values[1:, 2:]])
predicted = np.interp(years[1:], travel, distance)
assert np.all(np.abs(knick - predicted) <= 2 * (distance[1] - distance[0]))

fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 8))
for index in (0, 10, 20, 40, 60):
    top.plot(distance / 1e3, run.values[index, 2:] - uplift * years[index] * SECONDS_PER_YEAR, label=f"{years[index] / 1e3:.0f} kyr")
top.set(xlabel="distance from the outlet [km]", ylabel="elevation minus old uplift [m]", title="A knickpoint climbs the river")
top.legend()
bottom.plot(years[1:] / 1e3, knick / 1e3, "o", ms=3, label="simulated slope break")
bottom.plot(years / 1e3, np.interp(years, travel, distance) / 1e3, label="integral of dx / (K A^m)")
bottom.set(xlabel="time since the uplift doubled [kyr]", ylabel="knickpoint position [km]", title="Howard and Kerby: celerity K A^m")
bottom.legend()
fig.tight_layout()

plt.show()
