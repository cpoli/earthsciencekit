"""
Richardson's turbulence criterion (1920)
========================================

Lewis Fry Richardson asked when wind shear can stir up a stably stratified
layer. Turbulence must lift heavy air against buoyancy, at a rate set by
``N^2``, using energy drawn from the shear. Their ratio, the Richardson
number, decides: Miles and Howard later proved that a stratified shear flow
can become unstable only where it falls below 1/4. Over a night-time
boundary layer, the critical value marks the height where turbulence dies out.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import brunt_vaisala_squared, richardson_number
from earthsciencekit.constants import GRAVITY

height = np.linspace(10, 500, 200)
theta_gradient = 0.01 * (1 - np.exp(-height / 150))
shear = 8 / 200 * np.exp(-height / 200) + 1e-3
ri = richardson_number(brunt_vaisala_squared(290, theta_gradient, GRAVITY), shear)
assert np.all(np.diff(ri) > 0)
crossing = height[np.argmax(ri > 0.25)]
assert ri[0] < 0.25 < ri[-1]
assert np.isclose(richardson_number(-1e-4, 0.01), -1)
print(f"turbulence can persist below about {crossing:.0f} m")

fig, ax = plt.subplots()
ax.semilogx(ri, height)
ax.axvline(0.25, color="k", ls="--", label="Ri = 1/4")
ax.set(xlabel="gradient Richardson number", ylabel="height [m]", title="A stable night-time layer")
ax.legend()

plt.show()
