"""
Lal's cosmogenic-nuclide erosion rates (1991)
=============================================

Devendra Lal worked out how cosmic rays build up nuclides such as
beryllium-10 in quartz near the ground, and how erosion and radioactive
decay remove them. Production falls off exponentially with depth over about
60 cm of rock, so a grain collects nuclides as erosion brings it up, and the
concentration of a steadily eroding surface settles to production divided
by the sum of decay and erosion terms. Measured in river sand, it gives the
erosion rate of a whole catchment averaged over thousands of years.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import LAMBDA_BE10, SECONDS_PER_YEAR
from earthsciencekit.geomorphology import cosmogenic_concentration, cosmogenic_production_at_depth, steady_state_erosion_rate

production = 4.0e3 / SECONDS_PER_YEAR  # atoms/(kg s), 4 atoms per gram per year at sea level
density, attenuation = 2650.0, 1600.0
rates = np.array([1e-6, 1e-5, 1e-4, 1e-3]) / SECONDS_PER_YEAR  # m/s, from 1 m to 1 km per Myr
years = np.geomspace(1e2, 1e7, 200)
growth = cosmogenic_concentration(years[:, None] * SECONDS_PER_YEAR, production, LAMBDA_BE10, rates, density, attenuation)
saturated = growth[-1]
assert np.allclose(steady_state_erosion_rate(saturated, production, LAMBDA_BE10, density, attenuation), rates, rtol=1e-6)
assert np.isclose(cosmogenic_production_at_depth(production, attenuation / density, density, attenuation), production / np.e)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
for rate, curve in zip(rates, growth.T):
    left.loglog(years, curve / 1e3, label=f"{rate * SECONDS_PER_YEAR * 1e6:g} m/Myr")
left.set(xlabel="exposure time [years]", ylabel="10Be [atoms/g]", title="Accumulation until steady state")
left.legend()
grid = np.geomspace(1e-7, 1e-2, 100) / SECONDS_PER_YEAR
steady = cosmogenic_concentration(1e9 * SECONDS_PER_YEAR, production, LAMBDA_BE10, grid, density, attenuation)
right.loglog(grid * SECONDS_PER_YEAR * 1e6, steady / 1e3)
right.set(xlabel="erosion rate [m/Myr]", ylabel="steady 10Be [atoms/g]", title="Lal: concentration measures erosion rate")
fig.tight_layout()

plt.show()
