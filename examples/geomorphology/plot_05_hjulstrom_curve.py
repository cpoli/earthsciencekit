"""
Hjulström's erosion-deposition diagram (1935)
=============================================

Filip Hjulström gathered flume and river measurements into one diagram of
mean flow velocity against grain size, with a band where grains are eroded
from the bed, a field where they stay in transport and a line below which
they settle out. Its lesson is the wide transport field of fine grains: once
lifted, silt stays in suspension at velocities far below those that eroded
it, while gravel is dropped almost as soon as it stops rolling. The diagram
is rebuilt here for a 1 m deep flow from grain-scale physics: erosion at the
Shields threshold through Keulegan's log-law velocity, suspension where the
shear velocity exceeds the settling velocity, and deposition where the flow
can do neither. Hjulström's erosion curve also rises again for silt and
clay, from the cohesion of mud, which these cohesionless formulas leave out.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import (
    critical_shields_number,
    dimensionless_grain_size,
    ferguson_church_settling_velocity,
    keulegan_mean_velocity,
    shear_velocity,
)

diameter = np.geomspace(1e-6, 0.1, 300)
depth, rho_s, rho, nu, g = 1.0, 2650.0, 1000.0, 1e-6, 9.81
threshold = critical_shields_number(dimensionless_grain_size(diameter, rho_s, rho, nu, g))
roughness = np.maximum(2.5 * diameter, 1e-4)  # a bed of fine grains is no smoother than 0.1 mm
erosion = keulegan_mean_velocity(shear_velocity(threshold * (rho_s - rho) * g * diameter, rho), depth, roughness)
suspension = keulegan_mean_velocity(ferguson_church_settling_velocity(diameter, rho_s, rho, nu, g), depth, roughness)
deposition = np.minimum(erosion, suspension)
silt, gravel = diameter < 2e-5, diameter > 2e-3
assert np.all(deposition[silt] < 0.1 * erosion[silt]) and np.all(deposition[gravel] == erosion[gravel])
assert np.all(np.diff(erosion[diameter > 2e-4]) > 0)

fig, ax = plt.subplots()
ax.loglog(diameter * 1e3, erosion, label="erosion (Shields threshold)")
ax.loglog(diameter * 1e3, deposition, label="deposition (neither suspended nor rolling)")
ax.loglog(diameter * 1e3, suspension, ":", color="0.5", label="suspension (u* = settling velocity)")
ax.fill_between(diameter * 1e3, deposition, erosion, color="0.9")
ax.text(3e-3, 0.03, "transport", color="0.3")
ax.text(3e-3, 1.5, "erosion", color="tab:blue")
ax.text(3, 0.05, "deposition", color="tab:orange")
ax.set(xlabel="grain diameter [mm]", ylabel="mean velocity [m/s]", ylim=(1e-5, 20), title="Hjulström: fine grains stay in transport")
ax.legend(loc="lower right")

plt.show()
