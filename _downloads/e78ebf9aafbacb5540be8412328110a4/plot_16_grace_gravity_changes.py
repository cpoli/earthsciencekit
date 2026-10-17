"""
GRACE gravity changes (2002)
============================

The twin GRACE satellites, launched in 2002, measured the distance between
them to microns. Its changes gave a new gravity field every month, and the
differences between months show where water, ice and groundwater moved. Each
change of a Stokes coefficient becomes a layer of equivalent water, after
allowing for the elastic Earth sinking under the load. A truncated expansion
blurs a river basin's water into a smooth bump: the resolution is limited by
the highest degree, here 60.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import eval_legendre

from earthsciencekit.constants import EARTH_RADIUS
from earthsciencekit.geophysics import equivalent_water_height
from earthsciencekit.units import convert

thickness, cap = 0.1, np.radians(5.0)  # seasonal water in a basin of about 550 km radius
mean_density, water = 5514.0, 1000.0
degrees = np.arange(2, 61)
table_degree = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 40, 50, 70, 100]
table_love = -1e-3 * np.array([303, 194, 132, 104, 89, 81, 76, 72, 69, 64, 58, 51, 40, 33, 27, 20, 14])
load_love = np.interp(degrees, table_degree, table_love)  # approximate load Love numbers k'_n
legendre = thickness / 2 * (eval_legendre(degrees - 1, np.cos(cap)) - eval_legendre(degrees + 1, np.cos(cap)))
stokes = 3 * water * (1 + load_love) / (EARTH_RADIUS * mean_density * (2 * degrees + 1)) * legendre / np.sqrt(2 * degrees + 1)
recovered = equivalent_water_height(stokes, degrees, load_love, EARTH_RADIUS, mean_density, water)
rigid = equivalent_water_height(stokes, degrees, 0, EARTH_RADIUS, mean_density, water)
np.testing.assert_allclose(recovered * np.sqrt(2 * degrees + 1), legendre, rtol=1e-12)
assert np.all(np.abs(rigid) < np.abs(recovered))

angle = np.radians(np.linspace(0, 20, 401))
profile = (recovered * np.sqrt(2 * degrees + 1)) @ eval_legendre(degrees[:, None], np.cos(angle))
profile_rigid = (rigid * np.sqrt(2 * degrees + 1)) @ eval_legendre(degrees[:, None], np.cos(angle))
assert 0.8 * thickness < profile[0] < 1.1 * thickness

fig, ax = plt.subplots()
distance = EARTH_RADIUS * angle / 1e3
ax.plot(distance, convert(np.where(angle <= cap, thickness, 0), "m", "cm"), "k--", label="true water layer")
ax.plot(distance, convert(profile, "m", "cm"), label="from degrees 2-60")
ax.plot(distance, convert(profile_rigid, "m", "cm"), ":", label="ignoring elastic loading")
ax.set(xlabel="distance from basin centre [km]", ylabel="equivalent water [cm]", title="Water mass from monthly gravity")
ax.legend()

plt.show()
