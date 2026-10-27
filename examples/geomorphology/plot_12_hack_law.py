"""
Hack's law (1957)
=================

John Hack measured streams in the Shenandoah Valley and the Maryland
Piedmont and found that the length of the main stream grows as about the
0.6 power of drainage area, L = 1.4 A^0.6 with miles. An exponent above
one-half means that basins elongate as they grow: large basins are not
scaled-up small ones. The example fits Hack's law to a synthetic sample of
basins whose lengths scatter about it.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import hack_length

rng = np.random.default_rng(11)
area = 10 ** rng.uniform(5, 10, 120)  # m^2, from 0.1 to 10,000 km^2
length = hack_length(area, 0.32, 0.6) * rng.lognormal(0, 0.2, area.size)
exponent, intercept = np.polyfit(np.log(area), np.log(length), 1)
assert abs(exponent - 0.6) < 0.03
# Hack's 1.4 with miles is the SI coefficient 0.32.
mile = 1609.344
assert abs(1.4 * mile ** (1 - 2 * 0.6) - 0.32) < 0.001

fig, ax = plt.subplots()
ax.loglog(area / 1e6, length / 1e3, "o", ms=3, alpha=0.6, label="synthetic basins")
grid = np.geomspace(1e5, 1e10, 50)
ax.loglog(grid / 1e6, np.exp(intercept) * grid**exponent / 1e3, label=f"fit, exponent {exponent:.2f}")
ax.loglog(grid / 1e6, np.sqrt(grid) / 1e3, "--", color="0.5", label="exponent 1/2 (self-similar)")
ax.set(xlabel="drainage area [km^2]", ylabel="main-stream length [km]", title="Hack: basins elongate as they grow")
ax.legend()

plt.show()
