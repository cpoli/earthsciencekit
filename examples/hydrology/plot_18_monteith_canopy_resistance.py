"""
Monteith's canopy resistance (1965)
===================================

John Monteith treated a vegetated surface as one big leaf whose stomata add a
surface resistance in series with the aerodynamic resistance of the air above
it. Penman's equation becomes the Penman-Monteith equation, which reduces to
Penman's for a wet surface and falls as the stomata close. The example
computes the aerodynamic and surface resistances of the FAO reference grass
and the evaporation it gives, a few millimetres a day.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import CP_DRY_AIR, LATENT_HEAT_VAPORIZATION
from earthsciencekit.hydrology import (
    aerodynamic_resistance,
    bulk_surface_resistance,
    penman_evaporation,
    penman_monteith_evaporation,
    psychrometric_constant,
)

air, pressure, density, deficit, energy = 293.15, 101325.0, 1.2, 1000.0, 150.0
grass = 0.12  # m
aerodynamic = aerodynamic_resistance(2.0, 2.0, 2 / 3 * grass, 0.123 * grass, 0.0123 * grass)
surface = bulk_surface_resistance(100, 0.5 * 24 * grass)  # leaf area index 24 h for clipped grass
reference = penman_monteith_evaporation(energy, air, deficit, density, aerodynamic, surface, pressure)
wet = density * CP_DRY_AIR * deficit / (psychrometric_constant(pressure) * LATENT_HEAT_VAPORIZATION * aerodynamic)
assert np.isclose(aerodynamic, 208 / 2, rtol=0.01) and np.isclose(surface, 70, rtol=0.01)
assert np.isclose(
    penman_monteith_evaporation(energy, air, deficit, density, aerodynamic, 0, pressure), penman_evaporation(energy, air, wet, pressure)
)
assert 3 < reference * 86400 < 6

resistance = np.linspace(0, 1000, 201)
fig, ax = plt.subplots()
for wind in (1.0, 2.0, 5.0):
    ra = aerodynamic_resistance(wind, 2.0, 2 / 3 * grass, 0.123 * grass, 0.0123 * grass)
    ax.plot(
        resistance, penman_monteith_evaporation(energy, air, deficit, density, ra, resistance, pressure) * 86400, label=f"wind {wind:g} m/s"
    )
ax.plot(surface, reference * 86400, "ko", label=f"FAO grass reference, {reference * 86400:.1f} mm/day")
ax.set(xlabel="surface resistance [s/m]", ylabel="evaporation [mm/day]", title="Penman-Monteith: stomata throttle evaporation")
ax.legend()

plt.show()
