"""
Rouse's suspended-sediment profile (1937)
=========================================

Hunter Rouse balanced the settling of suspended grains against their upward
mixing by turbulence, with the eddy diffusivity of a logarithmic boundary
layer, and found the concentration profile as a power of height whose
exponent, the Rouse number, compares the settling velocity with the shear
velocity. Fine sand fills the whole depth of a river almost uniformly; medium
and coarse sand stay near the bed and travel close to it.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import ferguson_church_settling_velocity, rouse_concentration, rouse_number

depth, reference, shear = 3.0, 0.15, 0.06
height = np.linspace(reference, depth, 300)
fig, ax = plt.subplots()
previous = None
for diameter in (0.0625e-3, 0.125e-3, 0.25e-3, 0.5e-3):
    number = rouse_number(ferguson_church_settling_velocity(diameter, 2650, 1000, 1e-6, 9.81), shear)
    profile = rouse_concentration(height, depth, reference, 1.0, number)
    assert np.all(np.diff(profile) <= 0) and profile[-1] == 0
    if previous is not None:
        assert np.all(profile[1:-1] < previous[1:-1])
    previous = profile
    ax.plot(profile, height, label=f"{diameter * 1e3:g} mm, P = {number:.2f}")
ax.set(xlabel="concentration / concentration at 0.15 m", ylabel="height above the bed [m]", title="Rouse: coarser grains hug the bed")
ax.legend()

plt.show()
