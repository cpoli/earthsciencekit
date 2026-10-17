"""
Munk's abyssal recipes (1966)
=============================

Walter Munk noticed that the deep Pacific temperature profile is exponential
with a scale of about 1 km. Cold water sinking near Antarctica must upwell
through the abyss, and the upwelling ``w`` and downward mixing ``kappa``
balance when ``h = kappa / w``. With ``w`` about 1.2 cm per day the profile
needs a diffusivity of about 1e-4 m^2/s, far more than measured in the
interior, a puzzle that took 30 years to resolve.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY
from earthsciencekit.oceanography import abyssal_temperature, buoyancy_frequency_squared, linear_density
from earthsciencekit.units import convert

diffusivity, upwelling, column = 1.3e-4, 1.2e-2 / 86400, 3000.0
height = np.linspace(0, column, 3001)
temperature = abyssal_temperature(height, convert(1.0, "degC", "K"), convert(5.0, "degC", "K"), column, diffusivity, upwelling)
gradient = np.gradient(temperature, height)
scale = 1 / np.polyfit(height[10:-10], np.log(gradient[10:-10]), 1)[0]
assert np.isclose(scale, diffusivity / upwelling, rtol=1e-3)

density = linear_density(1027, temperature - temperature[0], 0, 2e-4, 7.6e-4)
n2 = buoyancy_frequency_squared(density, np.gradient(density, height), GRAVITY)
assert np.all(n2 > 0)

fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
left.plot(convert(temperature, "K", "degC"), (column - height) / 1e3 + 1)
left.set(xlabel="temperature [degC]", ylabel="depth [km]", title=f"Scale height {scale:.0f} m")
left.invert_yaxis()
right.plot(np.sqrt(n2) * 1e3, (column - height) / 1e3 + 1)
right.set(xlabel="buoyancy frequency N [10^-3 1/s]", title="Stratification of the abyss")
fig.tight_layout()

plt.show()
