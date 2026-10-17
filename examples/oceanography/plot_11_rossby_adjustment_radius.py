"""
Rossby's adjustment radius (1938)
=================================

Carl-Gustaf Rossby asked what happens when a disturbance in a rotating fluid
is released. Gravity waves carry the imbalance away, but rotation stops them
after a distance equal to the wave speed divided by ``f``: the deformation
radius. In the stratified ocean the first baroclinic radius is a few tens of
kilometres at mid-latitudes, the size of ocean eddies, and grows toward the
equator.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import coriolis_parameter
from earthsciencekit.constants import GRAVITY
from earthsciencekit.oceanography import internal_wave_speed, reduced_gravity, rossby_radius

buoyancy, depth = 3e-3, 4000.0
latitude = np.radians(np.linspace(5, 80, 151))
f = coriolis_parameter(latitude)
radius = rossby_radius(buoyancy, depth, f)
np.testing.assert_allclose(radius * f, buoyancy * depth / np.pi)
assert 20e3 < rossby_radius(buoyancy, depth, coriolis_parameter(np.radians(45))) < 60e3
two_layer = internal_wave_speed(reduced_gravity(2, 1027, GRAVITY), 500, 3500) / f

fig, ax = plt.subplots()
ax.semilogy(np.degrees(latitude), radius / 1e3, label="constant N, first mode")
ax.semilogy(np.degrees(latitude), two_layer / 1e3, "--", label="two layers, 2 kg/m^3 jump")
ax.set(xlabel="latitude [degrees]", ylabel="deformation radius [km]", title="Rossby radius of deformation")
ax.legend()

plt.show()
