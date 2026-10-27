"""
Penman's evaporation (1948)
===========================

Howard Penman combined the energy balance of a wet surface with Dalton's law,
and linearized the saturation curve between the air and the surface, so the
unknown surface temperature drops out. Evaporation becomes a weighted mean of
a radiation term and a drying-power term that need only the net radiation and
standard measurements at screen height. The example solves the full energy
balance for the surface temperature and shows that Penman's equation agrees
to within a few per cent.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from earthsciencekit.atmosphere import saturation_vapor_pressure
from earthsciencekit.constants import LATENT_HEAT_VAPORIZATION
from earthsciencekit.hydrology import dalton_evaporation, penman_evaporation, psychrometric_constant, saturation_slope

air, pressure = 293.15, 101325.0
vapor = 0.6 * saturation_vapor_pressure(air)
wind_function = 3.0e-8 * (1 + 0.54 * 2.0)  # kg/(m^2 s Pa) at 2 m/s
gamma = psychrometric_constant(pressure)
drying = dalton_evaporation(saturation_vapor_pressure(air), vapor, wind_function)
radiation = np.linspace(50, 400, 36)


def surface_balance(surface, net):
    latent = LATENT_HEAT_VAPORIZATION * dalton_evaporation(saturation_vapor_pressure(surface), vapor, wind_function)
    return latent + LATENT_HEAT_VAPORIZATION * wind_function * gamma * (surface - air) - net


surface = np.array([brentq(surface_balance, air - 20, air + 20, args=(net,)) for net in radiation])
exact = dalton_evaporation(saturation_vapor_pressure(surface), vapor, wind_function)
penman = penman_evaporation(radiation, air, drying, pressure)
assert np.all(np.abs(penman / exact - 1) < 0.05)
weight = saturation_slope(air) / (saturation_slope(air) + gamma)

fig, ax = plt.subplots()
ax.plot(radiation, exact * 86400, lw=4, alpha=0.4, label="energy balance solved for the surface temperature")
ax.plot(radiation, penman * 86400, "k", label="Penman (1948)")
ax.plot(radiation, weight * radiation / LATENT_HEAT_VAPORIZATION * 86400, ":", label=f"radiation term, weight {weight:.2f}")
ax.set(xlabel="net radiation [W/m^2]", ylabel="evaporation [mm/day]", title="Penman: no surface temperature needed")
ax.legend(fontsize="small")

plt.show()
