"""
Culling's hillslope diffusion (1960)
====================================

W. E. H. Culling showed that if soil creeps downhill at a rate proportional
to the slope, as Gilbert had argued in 1909, the land surface obeys the
diffusion equation, and the whole mathematics of heat conduction applies to
hillslopes. A fault scarp is the textbook case: a step that rounds into an
error function whose midpoint slope falls as one over the square root of
age. Hanks et al. (1984) inverted that slope into a "morphologic age", the
diffusivity times the time since the earthquake.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import morphologic_age, scarp_elevation, scarp_maximum_gradient

half, regional = 2.5, 0.05
kappa = 5e-3 / SECONDS_PER_YEAR  # m^2/s, 5 m^2 per kyr
distance = np.linspace(-30, 30, 1201)
ages = np.array([0, 500, 2000, 8000, 32000])
profiles = [scarp_elevation(distance, half, kappa, a * SECONDS_PER_YEAR, regional) for a in ages]
steepest = scarp_maximum_gradient(half, kappa, ages[1:] * SECONDS_PER_YEAR, regional)
recovered = morphologic_age(half, steepest, regional) / kappa / SECONDS_PER_YEAR
measured = [np.gradient(p, distance).max() for p in profiles[1:]]
assert np.allclose(recovered, ages[1:]) and np.allclose(measured, steepest, rtol=1e-3)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
for age, profile in zip(ages, profiles):
    left.plot(distance, profile, label=f"{age} years")
left.set(xlabel="distance from the fault [m]", ylabel="elevation [m]", title="Culling: a scarp degraded by diffusion")
left.legend()
right.loglog(ages[1:], np.degrees(np.arctan(steepest)), "o-")
right.set(xlabel="age [years]", ylabel="steepest slope [degrees]", title="Midpoint slope falls as age^-1/2")
fig.tight_layout()

plt.show()
