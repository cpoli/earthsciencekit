"""
Weertman's sliding law (1957)
=============================

Johannes Weertman asked how temperate ice gets past the bumps of its bed. It
melts on the upstream side of small obstacles and refreezes behind them, and
it creeps faster around large ones; the obstacle size that limits both sets
a sliding speed proportional to about the square of the basal stress. Sliding
adds to internal deformation, and on soft or wet beds it dominates.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY, ICE_DENSITY, SECONDS_PER_YEAR
from earthsciencekit.glaciology import driving_stress, glen_surface_velocity, weertman_sliding_velocity

thickness = 300.0
slope = np.linspace(0.005, 0.05, 100)
stress = driving_stress(thickness, slope, ICE_DENSITY, GRAVITY)
deformation = glen_surface_velocity(thickness, slope, 2.4e-24, 3, ICE_DENSITY, GRAVITY)
fig, ax = plt.subplots()
ax.plot(stress / 1e3, deformation * SECONDS_PER_YEAR, "k--", label="internal deformation")
for coefficient in (1e-16, 5e-16):
    sliding = weertman_sliding_velocity(stress, coefficient, 2)
    ax.plot(stress / 1e3, (sliding + deformation) * SECONDS_PER_YEAR, label=f"with sliding, C = {coefficient:g}")
assert np.isclose(weertman_sliding_velocity(2e5, 1e-17, 2), 4 * weertman_sliding_velocity(1e5, 1e-17, 2))
ax.set(xlabel="basal stress [kPa]", ylabel="surface velocity [m/year]", title="Weertman: sliding grows as stress squared")
ax.legend()

plt.show()
