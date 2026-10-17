"""
Dalton's partial pressures (1802)
=================================

John Dalton found that each gas in a mixture exerts the pressure it would
exert alone, and that the total is the sum. Water vapour in air is then an
ideal gas with its own partial pressure ``e``. The mixing ratio is the ratio
of the two partial densities, ``epsilon e / (p - e)``, and the specific
humidity is the vapour fraction of the total mass.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import mixing_ratio, specific_humidity
from earthsciencekit.constants import R_DRY_AIR, R_VAPOR

pressure, temperature = 101325.0, 293.15
vapor = np.linspace(0, 4000, 201)
w = mixing_ratio(vapor, pressure)
q = specific_humidity(w)
vapor_density = vapor / (R_VAPOR * temperature)
dry_density = (pressure - vapor) / (R_DRY_AIR * temperature)
np.testing.assert_allclose(w, vapor_density / dry_density)
np.testing.assert_allclose(q, vapor_density / (vapor_density + dry_density))

fig, ax = plt.subplots()
ax.plot(vapor / 100, w * 1000, label="mixing ratio (per kg dry air)")
ax.plot(vapor / 100, q * 1000, ls="--", label="specific humidity (per kg moist air)")
ax.set(xlabel="vapour partial pressure [hPa]", ylabel="water vapour [g/kg]", title=f"Moisture at {pressure / 100:.0f} hPa total pressure")
ax.legend()

plt.show()
