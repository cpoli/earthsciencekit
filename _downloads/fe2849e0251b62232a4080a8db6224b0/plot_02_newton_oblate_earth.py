"""
Newton's oblate Earth (1687)
============================

Newton argued that a rotating, fluid Earth must bulge at the equator. Treating
it as homogeneous, he balanced the weight of two water-filled canals, one to
the pole and one to the equator, and found a flattening of 1/230. The modern
value is 1/298: the Earth's mass is concentrated toward its centre, which makes
it less flattened than a uniform fluid. Jupiter, spinning in under ten hours,
shows the same effect more strongly.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import OMEGA_EARTH
from earthsciencekit.geophysics import centrifugal_ratio, homogeneous_flattening
from earthsciencekit.units import convert

earth = homogeneous_flattening(centrifugal_ratio(OMEGA_EARTH, 6378137.0, 9.7803267715))
assert 225 < 1 / earth < 235
day = np.linspace(6, 48, 200)
flattening = homogeneous_flattening(centrifugal_ratio(2 * np.pi / convert(day / 24, "day", "s"), 6378137.0, 9.7803267715))
jupiter_ratio = centrifugal_ratio(2 * np.pi / convert(9.925 / 24, "day", "s"), 7.1492e7, 24.79)
jupiter = homogeneous_flattening(jupiter_ratio)
assert jupiter > 0.0649

fig, ax = plt.subplots()
ax.loglog(day, flattening, label="homogeneous Earth, 5 m / 4")
ax.plot(convert(2 * np.pi / OMEGA_EARTH, "s", "day") * 24, earth, "o", label=f"Newton: 1/{1 / earth:.0f}")
ax.plot(23.934, 1 / 298.257, "s", label="observed Earth: 1/298")
ax.plot(9.925, jupiter, "^", label=f"homogeneous Jupiter: {jupiter:.3f}")
ax.plot(9.925, 0.0649, "v", label="observed Jupiter: 0.065")
ax.set(xlabel="length of day [hours]", ylabel="flattening (a - c) / a", title="Rotation flattens a fluid planet")
ax.legend()

plt.show()
