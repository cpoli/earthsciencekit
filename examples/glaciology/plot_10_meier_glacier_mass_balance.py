"""
Meier's glacier mass balance (1962)
===================================

Mark Meier proposed the definitions of glacier mass budget still in use: the
net balance at each point, the equilibrium line where it is zero, and the
accumulation-area ratio. With a balance rising linearly with altitude, a
glacier is in balance with its climate when the ablation below the
equilibrium line cancels the accumulation above it, which on this synthetic
valley glacier leaves about 60% of its area in the accumulation zone.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.glaciology import accumulation_area_ratio, linear_mass_balance

elevation = np.linspace(2000, 3400, 141)  # m, centres of 10 m bands
area = 4e4 * np.exp(-(((elevation - 2900) / 350) ** 2))  # m^2 per band, a wide upper basin
gradient = 0.007 / SECONDS_PER_YEAR


def net_balance(ela):
    return np.sum(linear_mass_balance(elevation, ela, gradient) * area)


balanced = brentq(net_balance, 2000, 3400)
ratio = accumulation_area_ratio(area[elevation > balanced].sum(), area.sum())
assert abs(net_balance(balanced)) < 1e-9 and 0.5 < ratio < 0.7

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
left.barh(elevation, area / 1e3, height=10, color="lightblue")
left.axhline(balanced, color="k", ls="--", label=f"equilibrium line, AAR = {ratio:.2f}")
left.set(xlabel="area per 10 m band [10^3 m^2]", ylabel="elevation [m]", title="Hypsometry")
left.legend()
right.plot(linear_mass_balance(elevation, balanced, gradient) * SECONDS_PER_YEAR, elevation)
right.axvline(0, color="0.5", ls=":")
right.set(xlabel="mass balance [m ice/year]", title="Meier: balance zero at the equilibrium line")
fig.tight_layout()

plt.show()
