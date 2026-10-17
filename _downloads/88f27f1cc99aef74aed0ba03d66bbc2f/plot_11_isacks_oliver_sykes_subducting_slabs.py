"""
Isacks, Oliver and Sykes's subducting slabs (1968)
==================================================

Bryan Isacks, Jack Oliver and Lynn Sykes traced deep earthquakes down
inclined zones beneath trenches and read them as cold lithosphere sinking
into the mantle. A shallower dip lengthens the path to a given depth; only the
part of the relative motion normal to the trench is convergence.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.plate_tectonics import convergence_rate, subduction_travel_time
from earthsciencekit.units import convert

relative = convert(9.0, "cm/year", "m/s")
speed = convergence_rate(relative, np.deg2rad(25.0))
depth = np.linspace(0, 660e3, 67)
fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4))
for dip in [30, 45, 60, 90]:
    time = subduction_travel_time(depth, speed, np.deg2rad(dip))
    right.plot(convert(time, "s", "year") / 1e6, depth / 1e3, label=f"dip {dip} degree")
    left.plot(depth / np.tan(np.deg2rad(dip)) / 1e3 if dip < 90 else 0 * depth, depth / 1e3)
assert np.isclose(float(subduction_travel_time(660e3, speed, np.pi / 2)), 660e3 / speed)
assert np.isclose(float(subduction_travel_time(660e3, speed, np.pi / 6)), 2 * 660e3 / speed)
for ax in (left, right):
    ax.invert_yaxis()
left.set(xlabel="distance from trench [km]", ylabel="depth [km]", title="Slab geometry")
right.set(xlabel="time since entering trench [million years]", title="Time to reach depth")
right.legend()
fig.tight_layout()

plt.show()
