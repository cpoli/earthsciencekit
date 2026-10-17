"""
Newton's equilibrium tide (1687)
================================

Isaac Newton explained the tides as the difference between the Moon's (and the
Sun's) attraction across the Earth. An ocean in equilibrium with that field
bulges by ``(3/4) (M/E) (R/d)**3 R`` in each semidiurnal constituent at the
equator, about 27 cm for the Moon and 12 cm for the Sun. The two constituents
beat into spring and neap tides every 14.8 days.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import EARTH_RADIUS
from earthsciencekit.oceanography import tidal_elevation
from earthsciencekit.units import convert

earth_mass = 5.972e24
lunar = 0.75 * (7.342e22 / earth_mass) * (EARTH_RADIUS / 3.844e8) ** 3 * EARTH_RADIUS
solar = 0.75 * (1.989e30 / earth_mass) * (EARTH_RADIUS / 1.496e11) ** 3 * EARTH_RADIUS
m2, s2 = 12.4206 * 3600, 12.0 * 3600
spring_neap = 1 / (1 / s2 - 1 / m2)
assert 0.45 < solar / lunar < 0.47
assert np.isclose(convert(spring_neap, "s", "day"), 14.77, atol=0.01)

time = np.linspace(0, convert(30, "day", "s"), 20001)
tide = tidal_elevation(time, lunar, m2, 0) + tidal_elevation(time, solar, s2, 0)
assert np.isclose(tide.max(), lunar + solar, rtol=1e-3)

fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(convert(time, "s", "day"), 100 * tide, lw=0.6)
ax.axhline(100 * (lunar + solar), color="tab:red", ls="--", label="spring: Moon + Sun")
ax.axhline(100 * (lunar - solar), color="tab:green", ls="--", label="neap: Moon - Sun")
ax.set(xlabel="time [days]", ylabel="equilibrium tide [cm]", title=f"Spring-neap cycle of {convert(spring_neap, 's', 'day'):.2f} days")
ax.legend(loc="lower right")

plt.show()
