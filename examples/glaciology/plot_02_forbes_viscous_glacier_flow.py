"""
Forbes's viscous glacier flow (1843)
====================================

James David Forbes surveyed the Mer de Glace and found that its centre moves
faster than its margins, and its surface faster than its depths, with a speed
that changes with the seasons. He concluded that a glacier is "an imperfect
fluid" flowing under its own weight. For a slab of uniform viscosity the
velocity is a parabola in height, and its depth average is two-thirds of the
surface speed.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from earthsciencekit.constants import GRAVITY, ICE_DENSITY, SECONDS_PER_YEAR
from earthsciencekit.glaciology import driving_stress, viscous_slab_velocity

thickness, slope = 200.0, 0.08
height = np.linspace(0, thickness, 201)
fig, ax = plt.subplots()
for viscosity in (2e13, 5e13, 1e14):
    velocity = viscous_slab_velocity(height, thickness, slope, viscosity, ICE_DENSITY, GRAVITY)
    mean = trapezoid(velocity, height) / thickness
    assert np.isclose(mean, 2 / 3 * velocity[-1], rtol=1e-4)
    ax.plot(velocity * SECONDS_PER_YEAR, height, label=f"viscosity {viscosity:.0e} Pa s")
stress = driving_stress(thickness, slope, ICE_DENSITY, GRAVITY)
ax.set(xlabel="velocity [m/year]", ylabel="height above the bed [m]", title=f"Forbes: a viscous slab, basal stress {stress / 1e3:.0f} kPa")
ax.legend()

plt.show()
