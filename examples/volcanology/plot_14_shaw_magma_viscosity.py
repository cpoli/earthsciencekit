"""
Shaw's magma viscosity (1972)
=============================

Herbert Shaw found that the logarithm of the viscosity of silicate melts is a
straight line in ``1/T`` whose slope depends only on composition, and that the
lines for basalt to rhyolite all pass near one point. That makes Shaw's law an
Arrhenius law with activation energy ``s * 1e4 * R``. Rhyolite at eruption
temperature is about a hundred million times stiffer than basalt.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import R_GAS
from earthsciencekit.volcanology import arrhenius_viscosity, shaw_viscosity

temperature = 1e4 / np.linspace(1.4, 10.5, 200)  # out to the pivot, far above any magma temperature
fig, ax = plt.subplots()
for slope, name in [(2.5, "basalt"), (3.2, "andesite"), (4.0, "rhyolite")]:
    viscosity = shaw_viscosity(temperature, slope)
    assert np.allclose(arrhenius_viscosity(temperature, shaw_viscosity(1300, slope), slope * 1e4 * R_GAS, 1300), viscosity)
    ax.semilogy(1e4 / temperature, viscosity, label=f"{name}, s = {slope}")
assert shaw_viscosity(1073, 4.0) / shaw_viscosity(1473, 2.5) > 1e7
ax.axvline(1.5, color="0.5", ls=":", lw=1)
ax.axvspan(1e4 / 1500, 1e4 / 950, color="0.9", label="magmatic temperatures")
ax.set(xlabel="10^4 / T [1/K]", ylabel="viscosity [Pa s]", title="Shaw's lines pivot near 10^4/T = 1.5")
ax.legend()

plt.show()
