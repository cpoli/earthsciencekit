"""
Bowen's reaction series (1928)
==============================

Norman Bowen's experiments showed that a basaltic melt crystallizes minerals
in a fixed order and that crystals either keep reacting with the melt or are
removed from it. Removing them (fractional crystallization) lets the residual
melt evolve far towards granite; keeping them in equilibrium limits how far it
can go. Trace elements show the difference most clearly: an element the
crystals reject is enriched without limit by fractionation, but at most
``1 / D`` by equilibrium crystallization.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.volcanology import equilibrium_crystallization_melt, fractional_crystallization_melt

melt = np.linspace(1, 0.01, 300)
fig, ax = plt.subplots()
for partition, label in [(0.1, "incompatible, D = 0.1"), (5.0, "compatible (Ni in olivine), D = 5")]:
    ax.semilogy(melt, fractional_crystallization_melt(1, melt, partition), label=f"fractional, {label}")
    ax.semilogy(melt, equilibrium_crystallization_melt(1, melt, partition), "--", label=f"equilibrium, {label}")
assert np.isclose(equilibrium_crystallization_melt(1, 0, 0.1), 10)
assert fractional_crystallization_melt(1, 0.01, 0.1) > 60
ax.invert_xaxis()
ax.set(xlabel="fraction of melt remaining", ylabel="concentration / initial", title="Removing crystals drives the melt further")
ax.legend(fontsize=8)

plt.show()
