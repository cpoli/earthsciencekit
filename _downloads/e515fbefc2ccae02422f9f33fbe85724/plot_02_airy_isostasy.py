"""
Airy's isostasy (1855)
======================

George Airy explained why the Himalaya deflect a plumb line less than their
visible mass predicts: mountains float on denser rock on light roots. A root
of depth ``r`` under an elevation ``h`` balances the pressure at the base of
the deepest column.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY
from earthsciencekit.plate_tectonics import airy_root

crust, mantle, thickness = 2800.0, 3300.0, 35e3
elevation = np.array([0, 1, 2, 3, 4, 5]) * 1e3
root = airy_root(elevation, crust, mantle)
compensation = thickness + root.max()
pressure = GRAVITY * (crust * (elevation + thickness + root) + mantle * (compensation - thickness - root))
np.testing.assert_allclose(pressure, pressure[0])
assert np.isclose(float(airy_root(4e3, crust, mantle)), 22.4e3)

fig, ax = plt.subplots()
ax.bar(np.arange(elevation.size), (elevation + thickness + root) / 1e3, bottom=-(thickness + root) / 1e3, color="tan", label="crust")
ax.axhline(-thickness / 1e3, color="0.5", ls="--", label="normal Moho")
ax.axhline(-compensation / 1e3, color="k", lw=0.8, label="compensation depth")
ax.set_xticks(np.arange(elevation.size), [f"{h / 1e3:.0f}" for h in elevation])
ax.set(xlabel="elevation [km]", ylabel="height above sea level [km]", title="Light roots under high mountains")
ax.legend(loc="lower left")

plt.show()
