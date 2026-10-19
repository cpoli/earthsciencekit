"""
Dansgaard and Johnsen's age-depth model (1969)
==============================================

Willi Dansgaard and Sigfus Johnsen dated the first deep ice core, drilled to
the bed at Camp Century in Greenland, with a flow model. Annual layers thin
as they sink because the ice spreads sideways; the vertical velocity falls
linearly with depth down to a kink, then quadratically to zero at the frozen
bed. Integrating its reciprocal gives the age, older at depth than in Nye's
(1963) model of uniform thinning.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.glaciology import dansgaard_johnsen_age, dansgaard_johnsen_downward_velocity, nye_age

thickness, kink, accumulation = 1370.0, 400.0, 0.35 / SECONDS_PER_YEAR  # Camp Century
depth = np.linspace(0, 1365, 1000)
age = dansgaard_johnsen_age(depth, thickness, kink, accumulation) / SECONDS_PER_YEAR
layer = dansgaard_johnsen_downward_velocity(depth, thickness, kink, accumulation) * SECONDS_PER_YEAR
assert np.isclose(layer[0], 0.35) and np.all(np.diff(age) > 0)
assert np.all(age[1:] > nye_age(depth[1:], thickness, accumulation) / SECONDS_PER_YEAR)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
left.semilogx(age[1:], depth[1:], label="Dansgaard-Johnsen")
left.semilogx(nye_age(depth[1:], thickness, accumulation) / SECONDS_PER_YEAR, depth[1:], "--", label="Nye")
left.axhline(thickness - kink, color="0.6", ls=":", label="kink")
left.invert_yaxis()
left.set(xlabel="age [years]", ylabel="depth [m of ice]", title="Camp Century time scale")
left.legend()
right.plot(layer * 100, depth)
right.set(xlabel="annual layer thickness [cm]", title="Layers thin with depth")
fig.tight_layout()

plt.show()
