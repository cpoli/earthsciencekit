"""
Kelvin's cooling Earth (1862)
=============================

William Thomson (Lord Kelvin) treated the Earth as a half-space cooling from a
uniform 7000 degrees Fahrenheit. The surface gradient of the conductive
solution falls as one over the square root of age, so the measured gradient of
about 1 degree Fahrenheit per 50 feet fixes an age near 100 million years.
Kelvin did not know about radioactive heating or mantle convection.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.plate_tectonics import halfspace_heat_flux, halfspace_temperature
from earthsciencekit.units import convert

surface = 273.0
excess = 7000 * 5 / 9
diffusivity = 1.18e-6
measured_gradient = (5 / 9) / (50 * 0.3048)
age = excess**2 / (np.pi * diffusivity * measured_gradient**2)
gradient = halfspace_heat_flux(age, diffusivity, 1.0, surface, surface + excess)
assert np.isclose(float(gradient), measured_gradient)
age_myr = convert(age, "s", "year") / 1e6
assert 90 < age_myr < 105

ages = convert(np.logspace(6, 10, 200), "year", "s")
depth = np.linspace(0, 400e3, 201)
fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4))
left.loglog(convert(ages, "s", "year") / 1e6, halfspace_heat_flux(ages, diffusivity, 1.0, surface, surface + excess) * 1e3)
left.axhline(measured_gradient * 1e3, color="k", ls="--", label="measured gradient")
left.axvline(age_myr, color="tab:red", label=f"{age_myr:.0f} million years")
left.set(xlabel="age [million years]", ylabel="surface gradient [K/km]")
left.legend()
right.plot(halfspace_temperature(depth, age, diffusivity, surface, surface + excess), depth / 1e3)
right.invert_yaxis()
right.set(xlabel="temperature [K]", ylabel="depth [km]", title="Kelvin's geotherm")
fig.tight_layout()

plt.show()
