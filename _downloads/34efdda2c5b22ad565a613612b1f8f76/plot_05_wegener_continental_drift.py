"""
Wegener's continental drift (1912)
==================================

Alfred Wegener argued from matching coastlines, fossils and rocks that South
America and Africa were once joined. At a few centimetres a year, the
separation since the opening of the South Atlantic is thousands of
kilometres. The rate here is a representative value, not a reconstruction.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.plate_tectonics import plate_displacement
from earthsciencekit.units import convert

speed = convert(3.0, "cm/year", "m/s")
time = convert(np.linspace(0, 120e6, 121), "year", "s")
separation = plate_displacement(speed, time)
assert np.isclose(separation[-1], 3.6e6)

fig, ax = plt.subplots()
ax.plot(convert(time, "s", "year") / 1e6, separation / 1e3)
ax.set(xlabel="time since break-up [million years]", ylabel="South America-Africa separation [km]")
ax.set_title("3 cm/year opens an ocean")

plt.show()
