"""
Wolman and Miller's magnitude and frequency (1960)
==================================================

M. Gordon Wolman and John Miller asked which events do most of the work of
shaping landscapes: the rare catastrophe or the ordinary flood. Sediment
transport grows faster than linearly with discharge, but large flows are
rare, so the product of transport rate and frequency peaks at a moderate,
frequent discharge, the "effective discharge", which in many rivers is
close to the bankfull flow that returns every year or two. Here a synthetic
record of daily flows, lognormally distributed, is binned and compared with
the closed form for a lognormal frequency and a power-law rating.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import effective_discharge

rng = np.random.default_rng(1)
median, sigma, exponent = 20.0, 1.0, 2.5
daily = median * rng.lognormal(0, sigma, 200 * 365)  # m^3/s, 200 years of daily flows
edges = np.arange(0, 605, 5.0)  # equal discharge classes, as Wolman and Miller used
counts, _ = np.histogram(daily, edges)
centres = 0.5 * (edges[1:] + edges[:-1])
frequency = counts / counts.sum()
work = frequency * 1e-4 * centres**exponent  # sediment carried in each class
predicted = effective_discharge(median, sigma, exponent)
peak = centres[np.argmax(work)]
assert abs(np.log(peak / predicted)) < 0.25 and np.mean(daily > predicted) > 0.05

fig, ax = plt.subplots()
ax.plot(centres, frequency / frequency.max(), label="frequency of daily flows")
ax.plot(centres, work / work.max(), label="sediment carried (rate x frequency)")
ax.axvline(predicted, color="k", ls="--", label=f"effective discharge {predicted:.0f} m^3/s")
ax.set(xlabel="discharge [m^3/s]", ylabel="normalized", title="Wolman and Miller: moderate floods do the work")
ax.legend()

plt.show()
