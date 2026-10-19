"""
Agassiz's ice age (1840)
========================

Louis Agassiz read the erratic boulders, polished rock and moraines of the
Jura as the work of glaciers far larger than those of the Alps today, and
proposed a great ice age. A modest cooling is enough: the equilibrium line,
above which snow survives the summer, falls by the cooling divided by the
lapse rate, and the area of a mountain range above it grows quickly. The
example uses a synthetic range with a normal distribution of elevations.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from earthsciencekit.glaciology import ela_shift

present_ela, lapse_rate = 2900.0, 0.0065  # m, K/m
mean_elevation, spread = 1400.0, 700.0  # m, of the range
cooling = np.linspace(0, 10, 101)
ela = present_ela + ela_shift(-cooling, lapse_rate)
fraction = norm.sf(ela, mean_elevation, spread)
assert np.isclose(ela_shift(-6.5, lapse_rate), -1000) and fraction[-1] > 20 * fraction[0]

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(cooling, ela, color="tab:blue")
left.set(xlabel="summer cooling [K]", ylabel="equilibrium-line altitude [m]", title="The snow line comes down")
right.semilogy(cooling, 100 * fraction, color="tab:cyan")
right.set(xlabel="summer cooling [K]", ylabel="area above the line [%]", title="Agassiz: a small cooling, a large ice cover")
fig.tight_layout()

plt.show()
