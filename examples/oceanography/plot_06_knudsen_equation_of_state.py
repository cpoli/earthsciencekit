"""
Knudsen's equation of state (1901)
==================================

Martin Knudsen, with Carl Forch and Søren Sørensen, measured seawater density
against temperature and salinity and published the Hydrographical Tables used
for the next 70 years. Near a reference state the density is linear in both:
warming makes water lighter and salt makes it heavier, so lines of equal
density in the temperature-salinity plane have slope ``alpha / beta``.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.oceanography import linear_density

alpha, beta = 2e-4, 7.6e-4
temperature, salinity = np.meshgrid(np.linspace(-10, 15, 101), np.linspace(-3, 3, 101))
density = linear_density(1027, temperature, salinity, alpha, beta)
assert np.isclose(linear_density(1027, 1.0, alpha / beta, alpha, beta), 1027)

fig, ax = plt.subplots()
contours = ax.contour(salinity, temperature, density, levels=12, cmap="viridis")
ax.clabel(contours, fmt="%.1f", fontsize=7)
ax.set(xlabel="salinity anomaly [g/kg]", ylabel="temperature anomaly [K]", title="Isopycnals of the linear equation of state [kg/m^3]")

plt.show()
