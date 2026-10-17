"""
Foote's heated carbon dioxide (1856)
====================================

Eunice Foote set glass cylinders of different gases in sunlight and found
that the one filled with carbon dioxide grew hottest and stayed hot longest.
She suggested that more of it in the air would warm the Earth. The simplest
model of that effect is one atmospheric layer, transparent to sunlight, that
absorbs a fraction of the infrared from the ground and re-emits it both up
and down. A fully absorbing layer warms the surface by a factor of
``2**0.25``.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import absorbed_solar, effective_temperature, greenhouse_surface_temperature, outgoing_longwave

absorbed = absorbed_solar(1361, 0.3)
emissivity = np.linspace(0, 1, 101)
surface = greenhouse_surface_temperature(absorbed, emissivity)
layer = surface / 2**0.25
to_space = (1 - emissivity) * outgoing_longwave(surface, 1) + outgoing_longwave(layer, emissivity)
np.testing.assert_allclose(to_space, absorbed)
assert np.isclose(surface[-1], 2**0.25 * effective_temperature(absorbed, 1))

fig, ax = plt.subplots()
ax.plot(emissivity, surface, label="surface")
ax.plot(emissivity, layer, label="atmospheric layer")
ax.axhline(effective_temperature(absorbed, 1), color="0.5", ls=":", label="no atmosphere")
ax.set(xlabel="infrared emissivity of the layer", ylabel="temperature [K]", title="One absorbing layer warms the ground")
ax.legend()

plt.show()
