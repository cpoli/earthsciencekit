"""
Stommel's two-box thermohaline circulation (1961)
=================================================

Stommel connected a warm, salty equatorial box with a cold polar box. Heat
makes the polar water dense and drives sinking there; salt works against it.
For moderate freshwater forcing the boxes have two stable states, a strong
thermal overturning and a weak reversed haline one, and which one the ocean
settles in depends on where it starts.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.oceanography import simulate_thermohaline_boxes, thermohaline_forcing
from earthsciencekit.units import convert

thermal, haline, volume = 3e7, 6e6, 1e16
forcing = 0.8 * thermal**2 / (4 * haline)
flow = np.linspace(-1.5 * thermal, thermal, 501)
time = convert(np.linspace(0, 1500, 301), "year", "s")
fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(thermohaline_forcing(flow, thermal, haline), convert(flow, "m^3/s", "Sv"))
left.axvline(forcing, color="0.5", ls="--")
for initial in [0.0, 4.0, 8.0]:
    result = simulate_thermohaline_boxes(thermal, haline, forcing, volume, time, initial_salinity_difference=initial)
    final = result.values[-1, 1]
    assert np.isclose(thermohaline_forcing(final, thermal, haline), forcing, rtol=1e-5)
    right.plot(convert(time, "s", "year"), convert(result.values[:, 1], "m^3/s", "Sv"), label=f"initial dS = {initial:g} g/kg")
    left.plot(forcing, convert(final, "m^3/s", "Sv"), "ko")
assert simulate_thermohaline_boxes(thermal, haline, forcing, volume, time).values[-1, 1] > 0
assert simulate_thermohaline_boxes(thermal, haline, forcing, volume, time, initial_salinity_difference=8).values[-1, 1] < 0
left.set(xlabel="salt forcing [g/kg m^3/s]", ylabel="steady overturning [Sv]", title="Equilibria fold at q = q_T / 2")
right.axhline(0, color="k", lw=0.5)
right.set(xlabel="time [years]", ylabel="overturning [Sv]", title="Same forcing, two outcomes")
right.legend()
fig.tight_layout()

plt.show()
