"""
Hansen's ocean heat uptake (1981)
=================================

James Hansen and colleagues projected the warming from rising CO2 with a
one-dimensional model whose ocean mixed layer stores heat. They showed that
the depth of ocean in contact with the surface sets how fast the warming
appears. Until equilibrium, part of the forcing goes into the ocean rather
than back to space. The heat stored at any time equals the integrated
imbalance between forcing and the radiative restoring.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import heat_capacity_column, heat_uptake, simulate_one_box
from earthsciencekit.units import convert

forcing, restoring = 3.7, 1.2
years = np.linspace(0, 60, 601)
time = convert(years, "year", "s")
fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
for depth in [50, 100, 200]:
    capacity = heat_capacity_column(1025, 3990, depth)
    anomaly = simulate_one_box(forcing, restoring, capacity, time).values[:, 0]
    imbalance = forcing - restoring * anomaly
    stored = np.concatenate([[0], np.cumsum((imbalance[1:] + imbalance[:-1]) / 2 * np.diff(time))])
    np.testing.assert_allclose(stored[-1], heat_uptake(capacity, anomaly[-1]), rtol=1e-4)
    left.plot(years, anomaly, label=f"{depth} m mixed layer")
    right.plot(years, imbalance, label=f"{depth} m mixed layer")
left.set(xlabel="years after step", ylabel="surface warming [K]")
right.set(xlabel="years after step", ylabel="heat flux into the ocean [W/m²]")
left.legend()
fig.tight_layout()

plt.show()
