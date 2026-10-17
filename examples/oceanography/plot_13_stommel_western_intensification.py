"""
Stommel's western intensification (1948)
========================================

Henry Stommel asked why the strongest currents, like the Gulf Stream and the
Kuroshio, run along the western sides of ocean basins. On a uniformly
rotating Earth his wind-driven gyre is symmetric. Once the Coriolis parameter
increases northward, the return flow crowds into a western boundary layer of
width ``friction / beta`` while the interior obeys Sverdrup's balance.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.oceanography import stommel_streamfunction
from earthsciencekit.units import convert

width, length, stress, density, friction = 6000e3, 3000e3, 0.1, 1025, 1e-6
x, y = np.meshgrid(np.linspace(0, width, 601), np.linspace(0, length, 151))
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, beta, label in zip(axes, [0, 2e-11], ["uniform rotation", "beta plane"]):
    psi = convert(stommel_streamfunction(x, y, stress, width, length, density, beta, friction), "m^3/s", "Sv")
    peak = x[0, np.argmax(psi[75])]
    ax.contour(x / 1e3, y / 1e3, psi, levels=10, cmap="viridis")
    ax.set(xlabel="eastward distance [km]", title=f"{label}: peak at {peak / 1e3:.0f} km")
    if beta == 0:
        assert np.isclose(peak, width / 2)
    else:
        assert peak < 0.1 * width
axes[0].set_ylabel("northward distance [km]")
fig.tight_layout()

plt.show()
