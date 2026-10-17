"""
Holmes's mantle convection (1931)
=================================

Arthur Holmes proposed that radioactive heat drives slow convection in the
mantle, which could carry continents. A layer heated from below convects when
its Rayleigh number exceeds about 10^3. For any plausible mantle viscosity the
whole mantle is far above that threshold.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY
from earthsciencekit.plate_tectonics import mantle_rayleigh_number

viscosity = np.logspace(19, 26, 200)
rayleigh = mantle_rayleigh_number(3e-5, 4500, GRAVITY, 2500, 2.89e6, 1e-6, viscosity)
critical = 1708.0
assert mantle_rayleigh_number(3e-5, 4500, GRAVITY, 2500, 2.89e6, 1e-6, 1e21) > 1e4 * critical
threshold = viscosity[np.argmin(np.abs(np.log(rayleigh / critical)))]
print(f"convection stops only above a viscosity of about {threshold:.1e} Pa s")

fig, ax = plt.subplots()
ax.loglog(viscosity, rayleigh, label="whole mantle")
ax.axhline(critical, color="k", ls="--", label="onset (rigid walls)")
ax.axvspan(1e21, 1e22, color="tab:orange", alpha=0.2, label="mantle viscosity")
ax.set(xlabel="viscosity [Pa s]", ylabel="Rayleigh number", title="The mantle must convect")
ax.legend()

plt.show()
