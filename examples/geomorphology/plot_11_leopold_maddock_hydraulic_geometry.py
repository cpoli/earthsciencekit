"""
Leopold and Maddock's hydraulic geometry (1953)
===============================================

Luna Leopold and Thomas Maddock plotted width, depth and velocity against
discharge at gauging stations across the American Midwest and Southwest
and found straight lines on log-log axes: downstream, at the mean annual
flow, width grows about as the square root of discharge, depth as the 0.4
power and velocity barely at all. Continuity requires the three exponents
to sum to one. The example fits the power laws to synthetic stations and,
through Manning's equation (``hydrology.manning_velocity``), shows that the
slope those channels need falls downstream: the concave long profile again.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import downstream_hydraulic_geometry
from earthsciencekit.hydrology import hydraulic_radius, manning_velocity

rng = np.random.default_rng(3)
discharge = np.geomspace(1, 3000, 40)
truth = downstream_hydraulic_geometry(discharge, 4.0, 0.5, 0.3, 0.4)
observed = truth * rng.lognormal(0, 0.12, truth.shape)
fits = [np.polyfit(np.log(discharge), np.log(observed[:, i]), 1) for i in range(3)]
exponents = [fit[0] for fit in fits]
assert abs(exponents[0] - 0.5) < 0.05 and abs(exponents[1] - 0.4) < 0.05 and abs(sum(exponents) - 1) < 0.08
width, depth, velocity = truth.T
radius = hydraulic_radius(width, depth)
slope = (velocity * 0.035 / radius ** (2 / 3)) ** 2  # Manning's equation solved for the slope
assert np.allclose(manning_velocity(0.035, radius, slope), velocity) and np.all(np.diff(slope) < 0)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
for i, (name, unit) in enumerate([("width", "m"), ("depth", "m"), ("velocity", "m/s")]):
    left.loglog(discharge, observed[:, i], "o", ms=3)
    left.loglog(discharge, np.exp(np.polyval(fits[i], np.log(discharge))), color=f"C{i}", label=f"{name} [{unit}] ~ Q^{exponents[i]:.2f}")
left.set(xlabel="discharge [m^3/s]", title="Leopold and Maddock: downstream hydraulic geometry")
left.legend()
right.loglog(discharge, slope)
right.set(xlabel="discharge [m^3/s]", ylabel="slope for Manning n = 0.035", title="The slope these channels need")
fig.tight_layout()

plt.show()
