"""
Moncrieff and Miller's CAPE (1976)
==================================

Mitchell Moncrieff and Martin Miller modelled tropical cumulonimbus and
squall lines and named the energy that drives them the convective available
potential energy: the work done by buoyancy on a parcel lifted from the
surface, the area between the parcel curve and the sounding. A stable layer
below the level of free convection costs energy instead, the convective
inhibition. If all of CAPE became kinetic energy, the updraught would reach
``sqrt(2 CAPE)``. The sounding is the analytic one of Weisman and Klemp
(1982), with surface moisture capped at 14 g/kg.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import (
    convective_available_potential_energy,
    convective_inhibition,
    mixing_ratio,
    parcel_temperature_profile,
    saturation_vapor_pressure,
    specific_humidity,
    virtual_temperature,
)
from earthsciencekit.constants import CP_DRY_AIR, EPSILON_WATER, GRAVITY, R_DRY_AIR

height = np.linspace(0, 16000, 321)
tropopause, theta_top, t_top = 12000.0, 343.0, 213.0
theta = np.where(
    height <= tropopause,
    300 + (theta_top - 300) * (height / tropopause) ** 1.25,
    theta_top * np.exp(GRAVITY * (height - tropopause) / (CP_DRY_AIR * t_top)),
)
# Hydrostatic balance in the Exner function: d(pi)/dz = -g / (c_p theta).
exner = 1 - np.concatenate([[0], np.cumsum(GRAVITY / CP_DRY_AIR * np.diff(height) * (1 / theta[1:] + 1 / theta[:-1]) / 2)])
pressure = 1e5 * exner ** (CP_DRY_AIR / R_DRY_AIR)
temperature = theta * exner
humidity = np.where(height <= tropopause, 1 - 0.75 * (height / tropopause) ** 1.25, 0.25)
vapor = np.minimum(humidity * mixing_ratio(saturation_vapor_pressure(np.maximum(temperature, 180)), pressure), 0.014)

surface_vapor_pressure = vapor[0] * pressure[0] / (EPSILON_WATER + vapor[0])
dewpoint = 1 / (1 / 273.15 - 461.5 / 2.5e6 * np.log(surface_vapor_pressure / 611.2))
parcel = parcel_temperature_profile(pressure, temperature[0], dewpoint)
parcel_vapor = np.minimum(vapor[0], mixing_ratio(saturation_vapor_pressure(np.maximum(parcel, 180)), pressure))
parcel_tv = virtual_temperature(parcel, specific_humidity(parcel_vapor))
environment_tv = virtual_temperature(temperature, specific_humidity(vapor))

cape = convective_available_potential_energy(height, parcel_tv, environment_tv)
cin = convective_inhibition(height, parcel_tv, environment_tv)
assert 1500 < cape < 3500 and -100 < cin <= 0
print(f"CAPE {cape:.0f} J/kg, CIN {cin:.0f} J/kg, updraught bound {np.sqrt(2 * cape):.0f} m/s")

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(environment_tv, height / 1e3, color="k", label="environment")
ax.plot(parcel_tv, height / 1e3, color="tab:red", label="lifted surface parcel")
ax.fill_betweenx(height / 1e3, environment_tv, parcel_tv, where=parcel_tv > environment_tv, color="tab:red", alpha=0.3, label="CAPE")
low = height < 3000
ax.fill_betweenx(
    height[low] / 1e3,
    environment_tv[low],
    parcel_tv[low],
    where=parcel_tv[low] < environment_tv[low],
    color="tab:blue",
    alpha=0.4,
    label="CIN",
)
ax.set(xlabel="virtual temperature [K]", ylabel="height [km]", title=f"CAPE {cape:.0f} J/kg, CIN {cin:.0f} J/kg")
ax.legend()

plt.show()
