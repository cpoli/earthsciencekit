"""
Budyko and Sellers's ice-albedo feedback (1969)
===============================================

Mikhail Budyko and William Sellers independently built energy-balance models
in which ice reflects more sunlight than ocean or land. Colder climates grow
more ice and absorb less, which cools them further. Both found that a drop of
a few percent in sunlight could freeze the whole planet, and that a frozen
Earth would stay frozen until the Sun grew far brighter: a hysteresis loop.
Inverting the energy balance gives the solar constant that holds each
temperature, which traces every equilibrium at once.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import equilibrium_solar_constant

parameters = (0.62, 0.62, 0.3, 265.0, 8.0)
temperature = np.linspace(220, 310, 2001)
solar = equilibrium_solar_constant(temperature, *parameters)
slope = np.diff(solar)
folds = np.where(np.sign(slope[1:]) != np.sign(slope[:-1]))[0] + 1
escape, collapse = solar[folds]
assert escape > 1361 > collapse
crossings = np.where(np.diff(np.sign(solar - 1361)))[0]
assert crossings.size == 3

stable = np.r_[slope > 0, True]
fig, ax = plt.subplots()
ax.plot(np.where(stable, solar, np.nan), temperature, "C0", label="stable")
ax.plot(np.where(~stable, solar, np.nan), temperature, "C0--", label="unstable")
cold, warm = slice(None, folds[0]), slice(folds[1], None)
freeze = temperature[cold][np.argmin(abs(solar[cold] - collapse))]
thaw = temperature[warm][np.argmin(abs(solar[warm] - escape))]
ax.annotate("", xy=(collapse, freeze), xytext=(collapse, temperature[folds[1]]), arrowprops={"arrowstyle": "->", "color": "tab:blue"})
ax.annotate("", xy=(escape, thaw), xytext=(escape, temperature[folds[0]]), arrowprops={"arrowstyle": "->", "color": "tab:red"})
ax.axvline(1361, color="0.5", ls=":", label="present Sun")
ax.set(xlabel="solar constant [W/m²]", ylabel="equilibrium temperature [K]", title="Snowball hysteresis")
ax.legend(loc="lower right")

plt.show()
