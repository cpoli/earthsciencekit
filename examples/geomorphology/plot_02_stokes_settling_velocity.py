"""
Stokes settling (1851)
======================

George Gabriel Stokes computed the drag on a sphere moving slowly through a
viscous fluid, in a paper on the internal friction of fluids and the motion
of pendulums. Balanced against the submerged weight, it gives a settling
velocity growing as the square of the diameter. The law holds for silt and
clay, but drag on sand and gravel grows faster than viscosity alone, and the
Ferguson-Church (2004) formula bends over to a square-root law. A quartz
grain of 2 mm falls at about 0.2 m/s, not the 3.6 m/s Stokes's law predicts.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import ferguson_church_settling_velocity, stokes_settling_velocity

diameter = np.geomspace(1e-6, 1e-2, 300)
water, viscosity = 1000.0, 1e-3
stokes = stokes_settling_velocity(diameter, 2650, water, viscosity, 9.81)
natural = ferguson_church_settling_velocity(diameter, 2650, water, viscosity / water, 9.81)
reynolds = natural * diameter * water / viscosity
fine = diameter < 2e-5
assert np.allclose(natural[fine], stokes[fine], rtol=0.05)
two_mm = np.argmin(np.abs(diameter - 2e-3))
assert stokes[two_mm] > 10 * natural[two_mm] and 0.1 < natural[two_mm] < 0.3

fig, ax = plt.subplots()
ax.loglog(diameter * 1e3, stokes, "--", label="Stokes (1851)")
ax.loglog(diameter * 1e3, natural, label="Ferguson and Church (2004), natural grains")
ax.axvline(diameter[np.argmin(np.abs(reynolds - 1))] * 1e3, color="0.6", ls=":", label="particle Reynolds number 1")
for name, size in [("clay", 0.002), ("silt", 0.0625), ("sand", 2.0)]:
    ax.axvline(size, color="0.85", lw=0.8)
    ax.text(size * 0.55, 3e-9, name, rotation=90, color="0.4")
ax.set(xlabel="grain diameter [mm]", ylabel="settling velocity [m/s]", ylim=(1e-9, 10), title="Stokes: quartz settling in water")
ax.legend()

plt.show()
