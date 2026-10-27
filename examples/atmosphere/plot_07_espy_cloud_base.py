"""
Espy's cloud base (1841)
========================

James Pollard Espy argued that storms are driven from below: rising air
cools by expansion until its vapour condenses, and the latent heat set free
keeps the cloud rising. Cloud base sits where the cooling parcel reaches its
dew point, so it can be read from the surface dew-point depression, about
125 m per kelvin. Bolton's (1980) condensation temperature makes the rule
exact for this package's thermodynamics: the lifted parcel cools at the dry
rate of 9.8 K/km, its dew point falls at under 2 K/km, and the lines cross at
the lifting condensation level.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import lcl_height, lcl_pressure, lcl_temperature, mixing_ratio, saturation_vapor_pressure
from earthsciencekit.constants import CP_DRY_AIR, EPSILON_WATER, GRAVITY, R_DRY_AIR

depression = np.linspace(0, 25, 51)
bases = {t: lcl_height(t, lcl_temperature(t, t - depression)) for t in (280.0, 300.0, 315.0)}
for base in bases.values():
    np.testing.assert_allclose(base[1:] / depression[1:], 125, rtol=0.04)

temperature, dewpoint = 300.0, 290.0
base_temperature = lcl_temperature(temperature, dewpoint)
height = np.linspace(0, 2000, 201)
parcel = temperature - GRAVITY / CP_DRY_AIR * height
pressure = 1e5 * (parcel / temperature) ** (CP_DRY_AIR / R_DRY_AIR)
vapor = mixing_ratio(saturation_vapor_pressure(dewpoint), 1e5)
partial = vapor * pressure / (EPSILON_WATER + vapor)
parcel_dewpoint = 1 / (1 / 273.15 - 461.5 / 2.5e6 * np.log(partial / 611.2))
crossing = height[np.argmin(np.abs(parcel - parcel_dewpoint))]
assert abs(crossing - lcl_height(temperature, base_temperature)) < 25
print(f"cloud base {lcl_height(temperature, base_temperature):.0f} m at {lcl_pressure(1e5, temperature, base_temperature) / 100:.0f} hPa")

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
for t, base in bases.items():
    left.plot(depression, base / 1e3, label=f"surface air at {t - 273.15:.0f} degC")
left.plot(depression, 0.125 * depression, "k:", label="125 m per kelvin")
left.set(xlabel="dew-point depression [K]", ylabel="cloud base [km]", title="Espy's rule")
left.legend()
right.plot(parcel, height / 1e3, label="parcel temperature (dry adiabat)")
right.plot(parcel_dewpoint, height / 1e3, label="parcel dew point")
right.axhline(lcl_height(temperature, base_temperature) / 1e3, color="0.5", ls="--", label="lifting condensation level")
right.set(xlabel="temperature [K]", ylabel="height [km]", title="Where the lifted parcel saturates")
right.legend()
fig.tight_layout()

plt.show()
