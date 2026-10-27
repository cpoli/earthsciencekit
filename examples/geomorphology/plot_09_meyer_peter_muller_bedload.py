"""
Meyer-Peter and Müller bedload (1948)
=====================================

Eugen Meyer-Peter and Robert Müller ran flume experiments at ETH Zürich
with sand and gravel and fitted a bedload law still used in river
engineering: the dimensionless transport rate grows as the 3/2 power of the
excess of the Shields stress over a threshold of 0.047. Applied to a gravel
river in flood, with the bed stress from the depth-slope product and the
hydraulic radius of ``hydrology.hydraulic_radius``, transport starts only
once the water is about a metre deep and then rises steeply. Wong and
Parker (2006) found the original coefficient about twice too high on a
plane bed.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import bed_shear_stress, bedload_flux, meyer_peter_muller_transport, shields_number
from earthsciencekit.hydrology import hydraulic_radius

theta = np.linspace(0, 0.3, 301)
mpm = meyer_peter_muller_transport(theta, 0.047)
wong_parker = 3.97 * np.maximum(theta - 0.0495, 0) ** 1.5
depth = np.linspace(0.1, 3, 300)
stress = bed_shear_stress(hydraulic_radius(30.0, depth), 0.003, 1000, 9.81)
diameter = 0.03
flux = bedload_flux(meyer_peter_muller_transport(shields_number(stress, diameter, 2650, 1000, 9.81), 0.047), diameter, 2650, 1000, 9.81)
onset = depth[np.argmax(flux > 0)]
assert np.all(mpm[theta <= 0.047] == 0) and 0.5 < onset < 1.5 and np.all(np.diff(flux) >= 0)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(theta, mpm, label="Meyer-Peter and Müller (1948)")
left.plot(theta, wong_parker, "--", label="Wong and Parker (2006)")
left.axvline(0.047, color="0.6", ls=":")
left.set(xlabel="Shields number", ylabel="Einstein transport number", title="Bedload against Shields stress")
left.legend()
right.plot(depth, flux * 30 * 2650)
right.axvline(onset, color="0.6", ls=":", label=f"motion starts at {onset:.2f} m")
right.set(xlabel="flow depth [m]", ylabel="bedload [kg/s] across 30 m", title="3 cm gravel, slope 0.003")
right.legend()
fig.tight_layout()

plt.show()
