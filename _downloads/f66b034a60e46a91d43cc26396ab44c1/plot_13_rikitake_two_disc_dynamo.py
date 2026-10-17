"""
Rikitake's two-disc dynamo (1958)
=================================

Tsuneji Rikitake coupled two disc dynamos, each driven by the same torque and
each powering the coil of the other. The currents never settle: they oscillate
about one of two steady states of opposite sign and, at irregular intervals,
jump to the other. The model was the first to show that a deterministic dynamo
can reverse its field spontaneously, without any outside trigger, with
polarity intervals of very unequal length, as in the geological record.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geophysics import rikitake_fixed_point, simulate_rikitake_dynamo

mu, k = 1.0, 2.0
time = np.linspace(0, 400, 8001)
result = simulate_rikitake_dynamo(mu, k, 1.0, time)
x, y, z = result.values.T
flips = np.flatnonzero(np.diff(np.sign(x)))
assert flips.size > 10
lengths = np.diff(time[flips])
assert lengths.max() > 3 * np.median(lengths)
steady = rikitake_fixed_point(mu, k)

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [2, 1]})
left.plot(time, x, lw=0.8)
left.axhline(0, color="k", lw=0.5)
left.set(xlabel="nondimensional time", ylabel="current x (sign = polarity)", title="Spontaneous reversals")
right.plot(x, y, lw=0.4)
right.plot([steady[0], -steady[0]], [steady[1], -steady[1]], "r*", ms=12, label="steady states")
right.set(xlabel="x", ylabel="y", title="Two unstable centres")
right.legend()
fig.tight_layout()

plt.show()
