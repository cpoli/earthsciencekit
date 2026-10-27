"""
Shields's threshold of motion (1936)
====================================

Albert Shields measured in a flume the stress at which grains of different
sizes and densities begin to move, and collapsed the results with
similarity arguments onto one curve: the critical ratio of bed stress to
submerged grain weight depends only on a grain Reynolds number. The curve
falls through the sand sizes and levels off near 0.055 for gravel. Here it
is drawn with the Soulsby-Whitehouse (1997) fit against the dimensionless
grain size, with quartz in water marked from silt to cobbles.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import critical_shields_number, dimensionless_grain_size

size = np.geomspace(0.1, 1e4, 400)
threshold = critical_shields_number(size)
grains = {"coarse silt": 4e-5, "fine sand": 1.6e-4, "coarse sand": 1e-3, "gravel": 1e-2, "cobble": 1e-1}
marks = {name: dimensionless_grain_size(d, 2650, 1000, 1e-6, 9.81) for name, d in grains.items()}
lowest = size[np.argmin(threshold)]
assert 10 < lowest < 30 and abs(critical_shields_number(marks["cobble"]) - 0.055) < 2e-3
assert critical_shields_number(marks["coarse silt"]) > 2 * threshold.min()

fig, ax = plt.subplots()
ax.loglog(size, threshold, label="Soulsby and Whitehouse (1997)")
for name, value in marks.items():
    ax.plot(value, critical_shields_number(value), "o")
    ax.annotate(name, (value, critical_shields_number(value)), textcoords="offset points", xytext=(4, 6))
ax.axhline(0.055, color="0.6", ls=":", label="0.055 for coarse grains")
ax.set(xlabel="dimensionless grain size D*", ylabel="critical Shields number", title="Shields: the threshold of sediment motion")
ax.legend()

plt.show()
