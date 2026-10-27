"""
Flint's slope-area relation (1974)
==================================

J. J. Flint showed that along the streams of a basin the channel slope
falls as a power of the drainage area, with an exponent, now called the
concavity index, of about one-half. It is the quantitative form of the
concave long profile, and the coefficient, the steepness index, has since
become a map of uplift and erosion rates. The example fits both to a
synthetic slope-area plot from a graded stream-power river with scatter.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import flint_slope, steady_state_slope

rng = np.random.default_rng(5)
area = 10 ** rng.uniform(6, 10, 300)
uplift, erodibility, m = 5e-4 / SECONDS_PER_YEAR, 1e-5 / SECONDS_PER_YEAR, 0.45
slope = steady_state_slope(uplift, area, erodibility, m, 1) * rng.lognormal(0, 0.3, area.size)
fit = np.polyfit(np.log(area), np.log(slope), 1)
concavity, steepness = -fit[0], np.exp(fit[1])
assert abs(concavity - m) < 0.03 and abs(np.log(steepness / (uplift / erodibility))) < 0.5

fig, ax = plt.subplots()
ax.loglog(area / 1e6, slope, ".", alpha=0.5, label="synthetic channel reaches")
grid = np.geomspace(1e6, 1e10, 50)
ax.loglog(grid / 1e6, flint_slope(grid, steepness, concavity), label=f"concavity {concavity:.2f}, k_s = {steepness:.0f}")
ax.set(xlabel="drainage area [km^2]", ylabel="channel slope", title="Flint: slope falls as a power of area")
ax.legend()

plt.show()
