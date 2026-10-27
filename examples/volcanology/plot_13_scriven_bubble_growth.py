"""
Scriven's bubble growth (1959)
==============================

L. E. Scriven solved the growth of a bubble that feeds on gas diffusing out of
a supersaturated liquid. The radius grows as ``2 beta sqrt(D t)``, and the
constant ``beta`` depends on the supersaturation and on the flow the
expanding bubble pushes into the liquid. In a magma that has just risen from
200 to 50 MPa, water held at the old pressure diffuses into bubbles. Each
bubble grows until its share of the melt is saturated again, and the magma
ends up 40 % gas by volume.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import R_VAPOR
from earthsciencekit.volcanology import (
    diffusive_supersaturation,
    scriven_bubble_radius,
    scriven_growth_constant,
    simulate_bubble_growth,
    water_solubility,
)

supersaturation = np.geomspace(1e-3, 1e2, 61)
curves = {ratio: scriven_growth_constant(supersaturation[supersaturation * ratio < 0.9], ratio) for ratio in (0.0, 0.01, 0.1)}
assert np.isclose(curves[0.0][0], np.sqrt(supersaturation[0] / 2), rtol=0.03)

temperature, melt_density, diffusivity = 1100.0, 2300.0, 1e-11
initial, final = water_solubility(np.array([2e8, 5e7]), 4.1e-6)
gas_density = 5e7 / (R_VAPOR * temperature)
phi = diffusive_supersaturation(initial, final, melt_density, gas_density)
beta = scriven_growth_constant(phi, gas_density / melt_density)
shell = (3 / (4 * np.pi * 1e12)) ** (1 / 3)  # 10^12 bubbles per cubic metre
time = np.linspace(0, 3000, 301)
cell = simulate_bubble_growth(1e-6, shell, initial, 5e7, temperature, time, diffusivity=diffusivity, melt_density=melt_density)
radius, dissolved, vesicularity = cell.values.T
assert np.isclose(dissolved[-1], final, rtol=1e-4) and 0.35 < vesicularity[-1] < 0.45
print(f"supersaturation {phi:.2f}, beta {beta:.2f}; final vesicularity {vesicularity[-1]:.0%}")

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
for ratio, values in curves.items():
    left.loglog(supersaturation[: values.size], values, label=f"gas / melt density {ratio:g}")
left.loglog(supersaturation, np.sqrt(supersaturation / 2), "k:", label="sqrt(phi / 2), quasi-static")
left.loglog(supersaturation, np.sqrt(3 / np.pi) * supersaturation, "k--", label="sqrt(3 / pi) phi, Plesset-Zwick")
left.set(xlabel="supersaturation phi", ylabel="growth constant beta", ylim=(1e-2, 1e2), title="Scriven's growth constant")
left.legend(fontsize="small")
right.plot(time / 60, scriven_bubble_radius(time, diffusivity, beta) * 1e6, "--", label="Scriven, unbounded melt")
right.plot(time / 60, radius * 1e6, label="one bubble per 10^-12 m^3 of melt")
right.axhline(shell * 1e6, color="0.6", ls=":", label="shell radius")
right.set(xlabel="time after decompression [min]", ylabel="bubble radius [um]", ylim=(0, 80), title="Growth until the melt is saturated")
right.legend(fontsize="small", loc="lower right")
fig.tight_layout()

plt.show()
