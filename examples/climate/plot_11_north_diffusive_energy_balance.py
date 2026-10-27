"""
North's diffusive energy-balance model (1975)
=============================================

Gerald North solved the Budyko-Sellers model with heat carried poleward by
diffusion, ``D d/dx[(1 - x^2) dT/dx]`` in the sine of latitude ``x``.
Without ice the temperature has just two Legendre modes. With an ice edge
where the temperature falls to -10 degC, the model has several equilibria
for the same Sun. Lowering the solar constant grows the ice cap slowly, until
it passes about 35 degrees and the planet freezes over. Once frozen, it stays
frozen until the Sun is a third brighter, and a small polar cap cannot exist
at all: the hysteresis of the ice-albedo feedback.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import north_temperature_profile, simulate_latitudinal_energy_balance
from earthsciencekit.units import convert

nodes = 45
time = convert(np.array([0.0, 40.0]), "year", "s")
x = (np.arange(nodes) + 0.5) / nodes


def sweep(solar_constants, start):
    state, edges = start, []
    for solar in solar_constants:
        result = simulate_latitudinal_energy_balance(solar, time, nodes=nodes, initial_temperature=state)
        state = result.values[-1, 2:]
        edges.append(result.values[-1, 1])
    return np.array(edges)


cooling = np.arange(1500, 1229, -5.0)
warming = np.arange(1230, 1851, 10.0)
down = sweep(cooling, 300.0)
up = sweep(warming, 230.0)
collapse = cooling[np.argmax(down == 0)]
escape = warming[np.argmax(up > 0)]
assert 1240 < collapse < 1280 and 1700 < escape < 1850
print(f"the ice-covered planet freezes over below {collapse:.0f} W/m^2 and thaws above {escape:.0f} W/m^2")

warm = simulate_latitudinal_energy_balance(1361, time, nodes=nodes).values[-1, 2:]
frozen = simulate_latitudinal_energy_balance(1361, time, nodes=nodes, initial_temperature=230.0).values[-1, 2:]
no_ice = north_temperature_profile(x, 1361, 0.3, 210, 2, 0.555, -0.482)

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
left.plot(cooling, np.sin(down), "o-", ms=3, color="tab:blue", label="Sun dimming from an ice-free start")
left.plot(warming, np.sin(up), "s-", ms=3, color="tab:red", label="Sun brightening from a snowball")
left.axvline(1361, color="0.5", ls=":", label="present Sun")
left.set(xlabel="solar constant [W/m^2]", ylabel="sine of the ice-edge latitude", title="Ice-albedo hysteresis")
left.legend(loc="center right", fontsize="small")
latitude = np.degrees(np.arcsin(x))
right.plot(latitude, warm - 273.15, label="present Sun, warm state")
right.plot(latitude, frozen - 273.15, label="present Sun, snowball")
right.plot(latitude, no_ice - 273.15, "--", label="no ice feedback, albedo 0.3")
right.axhline(-10, color="0.5", ls=":", label="ice edge, -10 degC")
right.set(xlabel="latitude [deg]", ylabel="temperature [degC]", title="Two climates for one Sun")
right.legend(fontsize="small")
fig.tight_layout()

plt.show()
