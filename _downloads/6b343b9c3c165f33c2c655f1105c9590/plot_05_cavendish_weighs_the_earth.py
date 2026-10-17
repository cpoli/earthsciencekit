"""
Cavendish weighs the Earth (1798)
=================================

Henry Cavendish hung a 1.8 m rod carrying two small lead balls from a fine
wire and brought two 158 kg balls close to them. The rod twisted by about a
milliradian, and timing its slow free swing calibrated the wire. Together these
give the gravitational constant, and with surface gravity the mean density of
the Earth: Cavendish found 5.48 times that of water, about twice the density of
surface rock, so the interior must be much denser.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import EARTH_RADIUS, G
from earthsciencekit.geophysics import earth_mean_density, torsion_balance_gravitational_constant

arm, separation, large, period = 1.86, 0.225, 158.0, 420.0
deflection = 1.0e-3
measured = torsion_balance_gravitational_constant(arm, separation, deflection, large, period)
density = earth_mean_density(9.81, EARTH_RADIUS, measured)
assert np.isclose(measured, G, rtol=0.01)
assert 5400 < density < 5600 and density / 2700 > 2

twist = np.linspace(0.5e-3, 1.5e-3, 50)
fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(twist * 1e3, torsion_balance_gravitational_constant(arm, separation, twist, large, period))
left.axhline(G, color="0.5", ls="--", label="CODATA G")
left.plot(deflection * 1e3, measured, "o")
left.set(xlabel="rod deflection [mrad]", ylabel="G [m^3/(kg s^2)]", title="Torsion balance")
left.legend()
constants = np.linspace(5e-11, 9e-11, 50)
right.plot(constants, earth_mean_density(9.81, EARTH_RADIUS, constants))
right.axhline(2700, color="0.5", ls=":", label="granite")
right.plot(measured, density, "o", label=f"{density:.0f} kg/m^3")
right.set(xlabel="G [m^3/(kg s^2)]", ylabel="mean density [kg/m^3]", title="Weighing the Earth")
right.legend()
fig.tight_layout()

plt.show()
