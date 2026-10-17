"""
King-Hele's flattening from satellite orbits (1958)
===================================================

The equatorial bulge pulls on an inclined orbit and makes its node drift
westward. Desmond King-Hele and Robert Merson timed the drift of Sputnik 2's
orbit and found a J2 that gave a flattening of 1/298.3, not the 1/297 of the
International Ellipsoid that surveyors had used since 1924. A few months of
tracking one satellite had improved on two centuries of arc measurements.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import OMEGA_EARTH
from earthsciencekit.geophysics import centrifugal_ratio, flattening_from_j2, nodal_precession_rate
from earthsciencekit.units import convert

gm, a, b, j2 = 3.986004418e14, 6378137.0, 6356752.3, 1.08263e-3
perigee, apogee, inclination = 6378137.0 + 212e3, 6378137.0 + 1660e3, np.radians(65.33)  # Sputnik 2
axis, eccentricity = (perigee + apogee) / 2, (apogee - perigee) / (apogee + perigee)
observed = nodal_precession_rate(axis, eccentricity, inclination, j2, a, gm)
recovered = observed / nodal_precession_rate(axis, eccentricity, inclination, 1, a, gm)
flattening = flattening_from_j2(recovered, centrifugal_ratio(OMEGA_EARTH, a, gm / (a * b)))
assert np.isclose(1 / flattening, 298.257, rtol=2e-3) and 1 / flattening - 297 > 1
per_day = convert(observed, "rad", "degree") * convert(1, "day", "s")
assert -3.5 < per_day < -2

angle = np.linspace(0, np.pi, 181)
fig, ax = plt.subplots()
for height in [300e3, 1000e3, 3000e3]:
    rate = nodal_precession_rate(a + height, 0, angle, j2, a, gm)
    ax.plot(
        convert(angle, "rad", "degree"), convert(rate, "rad", "degree") * convert(1, "day", "s"), label=f"circular, {height / 1e3:.0f} km"
    )
ax.plot(convert(inclination, "rad", "degree"), per_day, "ko", label=f"Sputnik 2: f = 1/{1 / flattening:.1f}")
ax.axhline(0, color="k", lw=0.5)
ax.set(xlabel="inclination [degrees]", ylabel="nodal drift [degrees/day]", title="The bulge turns the orbit")
ax.legend()

plt.show()
