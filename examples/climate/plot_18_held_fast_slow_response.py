"""
Held's fast and slow responses (2010)
=====================================

Isaac Held and colleagues ran a climate model under quadrupled CO2 and then
returned it abruptly to preindustrial forcing. Most of the warming vanished
within a few years, but a remnant decayed only over centuries: the fast
response of the surface and the slow recalcitrant response of the deep ocean.
A two-box model, a mixed layer coupled to a deep ocean, reproduces both. Under
a step forcing the surface first climbs to a plateau near
``F / (lambda + gamma)`` and then creeps towards ``F / lambda`` as the deep
ocean fills with heat.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import equilibrium_warming, simulate_two_box, two_box_fast_warming, two_box_timescales
from earthsciencekit.units import convert

forcing, feedback, exchange = 7.4, 1.3, 0.7  # W/m^2 for quadrupled CO2, W/(m^2 K)
surface, deep = (convert(c, "year", "s") for c in (8.0, 100.0))  # J/(m^2 K): 8 and 100 W year/(m^2 K)
fast, slow = convert(two_box_timescales(feedback, exchange, surface, deep), "s", "year")
plateau, final = two_box_fast_warming(forcing, feedback, exchange), equilibrium_warming(forcing, feedback)
print(f"fast timescale {fast:.1f} years, slow {slow:.0f} years; plateau {plateau:.2f} K of {final:.2f} K")

years = np.linspace(0, 300, 601)
step = simulate_two_box(forcing, feedback, exchange, surface, deep, convert(years, "year", "s"))
assert abs(step.values[np.searchsorted(years, 5 * fast), 0] - plateau) < 0.15 * final
after = np.linspace(0, 100, 401)
release = simulate_two_box(0.0, feedback, exchange, surface, deep, convert(after, "year", "s"), initial_anomaly=step.values[200])
remnant = release.values[np.searchsorted(after, 5 * fast), 0]
assert 0.1 * step.values[200, 0] < remnant < release.values[0, 0] / 2

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))
left.plot(years, step.values[:, 0], label="surface")
left.plot(years, step.values[:, 1], label="deep ocean")
left.axhline(plateau, color="0.5", ls="--", label="F / (lambda + gamma)")
left.axhline(final, color="0.5", ls=":", label="F / lambda")
left.set(xlabel="years after quadrupling CO2", ylabel="temperature anomaly [K]", title="Step forcing")
left.legend(fontsize="small")
right.plot(after, release.values[:, 0], label="surface")
right.plot(after, release.values[:, 1], label="deep ocean")
right.axvline(5 * fast, color="0.5", ls=":", label="five fast timescales")
right.set(xlabel="years after returning to preindustrial forcing", ylabel="temperature anomaly [K]", title="The slow remnant")
right.legend(fontsize="small")
fig.tight_layout()

plt.show()
