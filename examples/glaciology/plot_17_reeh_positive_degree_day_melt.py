"""
Reeh's positive-degree-day model (1991)
=======================================

Niels Reeh estimated the melt of the whole Greenland ice sheet from air
temperature alone. A sinusoidal annual cycle with normally distributed daily
departures gives the expected positive degree-days even where the mean stays
below freezing. Those degree-days melt the year's snow at about 3 mm per
degree-day and then the ice beneath at about 8, and up to 60% of the snowfall
can be refilled by refreezing meltwater before any runs off. Along a synthetic
flow line the surface mass balance changes sign at the equilibrium line.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.glaciology import expected_positive_degree_days, pdd_surface_melt

day = 86400.0
days = np.arange(365)
elevation = np.linspace(0, 3000, 61)  # m
mean_annual = 273.15 - 3.0 - 0.0065 * elevation  # K
amplitude = 9.0  # K, July above the annual mean
monthly = mean_annual[:, None] + amplitude * np.cos(2 * np.pi * (days - 196) / 365)[None, :]
pdd = expected_positive_degree_days(monthly, 5.0, day).sum(axis=1)
plain = (np.maximum(monthly - 273.15, 0) * day).sum(axis=1)
assert np.all(pdd >= plain) and pdd[-1] > 0 == plain[-1]  # warm days melt even where no day's mean thaws

snowfall = 0.3 + 0.2 * elevation / 3000  # m water per year
snow, ice, runoff = pdd_surface_melt(pdd, snowfall, 3e-3 / day, 8e-3 / day, 0.6).T
balance = snowfall - runoff
equilibrium = np.interp(0, balance, elevation)
assert np.all(np.diff(balance) > 0) and 500 < equilibrium < 2000
assert np.all(ice[snow < snowfall] == 0)  # ice melts only once the snow is gone
assert np.all(runoff <= snow + ice)
print(f"equilibrium line at {equilibrium:.0f} m; balance {balance[0]:.2f} m/yr at sea level")

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
left.plot(pdd / day, elevation, label="expected, sigma = 5 K")
left.plot(plain / day, elevation, "--", label="from daily means only")
left.set(xlabel="positive degree-days per year [K day]", ylabel="elevation [m]", title="Fluctuations add melt at altitude")
left.legend()
right.plot(snow, elevation, label="snow melt")
right.plot(ice, elevation, label="ice melt")
right.plot(runoff, elevation, label="runoff")
right.plot(balance, elevation, "k", lw=2, label="surface mass balance")
right.axvline(0, color="0.5", lw=0.8)
right.axhline(equilibrium, color="0.5", ls=":", lw=0.8)
right.set(xlabel="m of water per year", title=f"Equilibrium line near {equilibrium:.0f} m")
right.legend(fontsize=8)
fig.tight_layout()

plt.show()
