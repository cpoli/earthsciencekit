"""
Pascal and Périer on the Puy de Dôme (1648)
===========================================

At Blaise Pascal's request, his brother-in-law Florin Périer carried a
mercury barometer from Clermont to the top of the Puy de Dôme, about 500 toises
(975 m) higher. The column fell by 3 inches and 1.5 lines, about 85 mm, while
a second barometer left in the town stayed put. Air has weight, and there is
less of it above a mountain. An isothermal column at 10 °C reproduces the
measured drop to within a few millimetres.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import barometer_column_height, barometric_pressure
from earthsciencekit.constants import GRAVITY

mercury, climb = 13595.1, 975.0
inch, line = 27.07e-3, 27.07e-3 / 12
measured_drop = 3 * inch + 1.5 * line
height = np.linspace(0, 1500, 151)
pressure = barometric_pressure(101325, height, 283.15, GRAVITY)
column = barometer_column_height(pressure, mercury, GRAVITY)
drop = float(column[0] - barometer_column_height(barometric_pressure(101325, climb, 283.15, GRAVITY), mercury, GRAVITY))
assert abs(drop - measured_drop) < 3e-3
print(f"model drop {drop * 1000:.1f} mm, Périer measured {measured_drop * 1000:.1f} mm")

fig, ax = plt.subplots()
ax.plot(column * 1000, height, label="isothermal column, 283 K")
ax.plot([column[0] * 1000, column[0] * 1000 - measured_drop * 1000], [0, climb], "o", color="tab:red", label="Périer, 19 September 1648")
ax.set(xlabel="mercury column [mm]", ylabel="height above Clermont [m]", title="The barometer falls on the mountain")
ax.legend()

plt.show()
