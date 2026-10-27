"""
Gilbert's graded river (1877)
=============================

In his report on the Henry Mountains of Utah, Grove Karl Gilbert argued
that a river adjusts its slope until it can just do the work demanded of
it: steep in the headwaters, where little water flows, and gentle
downstream, where the discharge is large. The result is the concave
"profile of equilibrium". Here a flat plain rises at 1 mm/year above a
fixed outlet and stream power cuts into it; within about two million years
the river is graded, its profile concave and its incision everywhere equal
to the uplift.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import simulate_river_profile, steady_state_slope, stream_power_erosion

distance = np.linspace(0, 50e3, 201)
area = ((50e3 - distance + 300) / 0.32) ** (1 / 0.6)  # Hack's law from the divide
uplift, erodibility = 1e-3 / SECONDS_PER_YEAR, 5.4e-6 / SECONDS_PER_YEAR
years = np.linspace(0, 3e6, 31)
run = simulate_river_profile(distance, area, uplift, erodibility, years * SECONDS_PER_YEAR)
graded = np.concatenate([[0.0], np.cumsum(steady_state_slope(uplift, area[1:], erodibility, 0.5, 1) * np.diff(distance))])
final = run.values[-1, 2:]
slope = np.diff(final) / np.diff(distance)
assert np.max(np.abs(final - graded)) < 0.5 and np.all(np.diff(slope) > 0)
assert np.allclose(stream_power_erosion(area[1:], slope, erodibility, 0.5, 1), uplift, rtol=1e-3)

fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 8))
for index in (2, 5, 10, 15, 30):
    top.plot(distance / 1e3, run.values[index, 2:], label=f"{years[index] / 1e6:.1f} Myr")
top.plot(distance / 1e3, graded, "k--", label="graded (steady state)")
top.set(xlabel="distance from the outlet [km]", ylabel="elevation [m]", title="Gilbert: a river grades itself to a concave profile")
top.legend()
bottom.plot(years / 1e6, run.values[:, 0])
bottom.axhline(graded[-1], color="k", ls="--")
bottom.set(xlabel="time [Myr]", ylabel="relief [m]", title="Relief rises until incision balances uplift")
fig.tight_layout()

plt.show()
