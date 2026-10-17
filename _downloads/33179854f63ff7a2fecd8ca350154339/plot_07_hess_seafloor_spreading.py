"""
Hess's seafloor spreading (1962)
================================

Harry Hess proposed that new ocean floor forms at mid-ocean ridges and moves
away on both sides. Two plates each moving at the half rate open the ridge at
twice that speed, and an ocean basin a few thousand kilometres wide is
renewed in a few hundred million years.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.plate_tectonics import full_spreading_rate
from earthsciencekit.units import convert

half = convert(np.array([1.0, 2.5, 5.0]), "cm/year", "m/s")
full = full_spreading_rate(half)
np.testing.assert_allclose(full, 2 * half)
renewal = 5e6 / full_spreading_rate(convert(1.0, "cm/year", "m/s"))
assert np.isclose(convert(renewal, "s", "year"), 250e6)

time = convert(np.linspace(0, 200e6, 101), "year", "s")
fig, ax = plt.subplots()
for h, f in zip(half, full):
    ax.plot(convert(time, "s", "year") / 1e6, f * time / 1e3, label=f"half rate {convert(h, 'm/s', 'cm/year'):.1f} cm/year")
ax.axhline(5000, color="k", ls="--", lw=0.8)
ax.set(xlabel="time [million years]", ylabel="new ocean floor [km]", title="Both flanks grow from the ridge")
ax.legend()

plt.show()
