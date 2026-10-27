"""
Herron and Langway's firn densification (1980)
==============================================

Michael Herron and Chester Langway fitted the density-depth profiles of
cores from Greenland and Antarctica with a model in two stages: grain
settling and packing up to 550 kg/m^3, then slower sintering to the pore
close-off near 830 kg/m^3. The first stage depends only on temperature, the
second also on accumulation. Warm, snowy sites seal their air near 50 m
within a few centuries; cold, dry plateaus like Vostok seal it near 100 m
after thousands of years, which is the ice-age minus gas-age difference of
an ice core.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import ICE_DENSITY, SECONDS_PER_YEAR
from earthsciencekit.glaciology import gas_age_difference, herron_langway_density, herron_langway_depth, herron_langway_mean_density

sites = {  # mean annual temperature [K] and accumulation [m ice/yr]
    "Dye 3 (Greenland)": (253.15, 0.56),
    "Summit (Greenland)": (241.15, 0.24),
    "Dome C (Antarctica)": (219.15, 0.03),
    "Vostok (Antarctica)": (218.15, 0.024),
}
surface = 350.0
depth = np.linspace(0, 130, 1301)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 5.5))
deltas = {}
for name, (temperature, rate) in sites.items():
    accumulation = rate / SECONDS_PER_YEAR
    density = herron_langway_density(depth, temperature, accumulation, surface)
    critical = herron_langway_depth(550.0, temperature, accumulation, surface)
    close_off = herron_langway_depth(830.0, temperature, accumulation, surface)
    mean = herron_langway_mean_density(close_off, temperature, accumulation, surface)
    deltas[name] = gas_age_difference(close_off, mean, ICE_DENSITY, accumulation) / SECONDS_PER_YEAR
    assert np.all(np.diff(density) > 0) and surface < mean < 830
    line = left.plot(density, depth, label=f"{name}: close-off {close_off:.0f} m")[0]
    left.plot([550, 830], [critical, close_off], "o", color=line.get_color())
    print(f"{name}: 550 kg/m^3 at {critical:.1f} m, close-off at {close_off:.1f} m, ice - gas age {deltas[name]:.0f} yr")

assert deltas["Dye 3 (Greenland)"] < 400 and 2000 < deltas["Vostok (Antarctica)"] < 4000
left.axvline(550, color="0.6", lw=0.8, ls="--")
left.axvline(830, color="0.6", lw=0.8, ls=":")
left.invert_yaxis()
left.set(xlabel="density [kg/m^3]", ylabel="depth [m]", title="Two-stage densification (dots: 550 and 830 kg/m^3)")
left.legend(fontsize=8, loc="lower left")

temperatures = np.linspace(213.15, 258.15, 46)
rates = np.array([0.02, 0.05, 0.1, 0.25, 0.5])
for rate in rates:
    accumulation = rate / SECONDS_PER_YEAR
    close_off = herron_langway_depth(830.0, temperatures, accumulation, surface)
    mean = herron_langway_mean_density(close_off, temperatures, accumulation, surface)
    delta = gas_age_difference(close_off, mean, ICE_DENSITY, accumulation) / SECONDS_PER_YEAR
    assert np.all(np.diff(delta) < 0)
    right.semilogy(temperatures - 273.15, delta, label=f"{rate} m ice/yr")
right.set(xlabel="mean annual temperature [degC]", ylabel="ice age - gas age at close-off [yr]", title="Cold, dry firn seals its air late")
right.legend(title="accumulation", fontsize=8)
fig.tight_layout()

plt.show()
