"""
Bagnold's physics of blown sand (1941)
======================================

Ralph Bagnold, back from exploring the Libyan desert, built a wind tunnel
and measured how wind moves sand. Grains start to move at a threshold shear
velocity that grows with the square root of their size; above it they hop
in saltation, and the mass they carry grows as the cube of the shear
velocity. Kawamura (1951) added the threshold to the cube law. Bagnold also
showed that a dune keeping its shape advances by the sand flux over its crest
divided by its height, so small dunes overtake large ones.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import bagnold_sand_flux, dune_migration_rate, kawamura_sand_flux, threshold_shear_velocity

air, sand, g, diameter = 1.2, 2650.0, 9.81, 2.5e-4
threshold = threshold_shear_velocity(diameter, sand, air, g, 0.1)
shear = np.linspace(0, 0.8, 401)
bagnold = bagnold_sand_flux(shear, diameter, air, g, 1.8)
kawamura = kawamura_sand_flux(shear, threshold, air, g, 2.78)
assert np.all(kawamura[shear <= threshold] == 0) and kawamura[-1] > bagnold[-1] / 2
height = np.linspace(1, 30, 100)
speed = dune_migration_rate(kawamura_sand_flux(0.45, threshold, air, g, 2.78), height, 1600) * SECONDS_PER_YEAR
assert np.isclose(speed[0] * height[0], speed[-1] * height[-1])

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(shear, bagnold * 1e3, label="Bagnold (1941), cube law")
left.plot(shear, kawamura * 1e3, label="Kawamura (1951), with threshold")
left.axvline(threshold, color="0.6", ls=":", label=f"threshold {threshold:.2f} m/s")
left.set(xlabel="shear velocity [m/s]", ylabel="sand flux [g/(m s)]", title="Bagnold: saltation flux of 0.25 mm sand")
left.legend()
right.plot(height, speed)
right.set(xlabel="dune height [m]", ylabel="migration rate [m/year]", title="Dunes at u* = 0.45 m/s: small ones run")
fig.tight_layout()

plt.show()
