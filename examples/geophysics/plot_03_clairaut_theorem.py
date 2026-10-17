"""
Clairaut's theorem (1743)
=========================

Alexis Clairaut showed that a rotating Earth in hydrostatic equilibrium, of any
internal layering, has a gravity that grows from equator to pole by the
fraction ``5 m / 2 - f``. Pendulum gravity measured at a few latitudes therefore
gives the flattening without any arc survey. Applied to today's normal gravity
it returns the flattening to within about 1%, and the first-order formula
matches the exact Somigliana gravity of the reference ellipsoid to 50 parts per
million.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geophysics import clairaut_gravity, normal_gravity
from earthsciencekit.units import convert

ratio = 0.00344978600308  # geodetic m = omega**2 a**2 b / GM of GRS80
equator, pole = normal_gravity(0), normal_gravity(np.pi / 2)
flattening = 5 * ratio / 2 - (pole - equator) / equator
assert 295 < 1 / flattening < 305
latitude = np.linspace(0, np.pi / 2, 91)
clairaut = clairaut_gravity(latitude, equator, 1 / 298.257222101, ratio)
somigliana = normal_gravity(latitude)
np.testing.assert_allclose(clairaut, somigliana, rtol=5e-5)

fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))
degrees = convert(latitude, "rad", "degree")
top.plot(degrees, somigliana, label="Somigliana (GRS80)")
top.plot(degrees, clairaut, "--", label="Clairaut, first order")
top.set(ylabel="gravity [m/s^2]", title=f"Flattening from gravity alone: 1/{1 / flattening:.0f}")
top.legend()
bottom.plot(degrees, convert(clairaut - somigliana, "m/s^2", "mGal"))
bottom.set(xlabel="latitude [degrees]", ylabel="Clairaut - Somigliana [mGal]")
fig.tight_layout()

plt.show()
