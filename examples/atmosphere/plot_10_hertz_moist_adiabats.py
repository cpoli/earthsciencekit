"""
Hertz's moist adiabats (1884)
=============================

Heinrich Hertz drew the first thermodynamic diagram for moist air: graphs
from which the state of a rising parcel could be read as it passes through
saturation, condensation and freezing. Once saturated, the parcel cools more
slowly than the dry adiabat, because condensation releases latent heat, and
the gap closes in cold air that holds little vapour. Here the condensate
falls out at once (the pseudo-adiabat), and a parcel lifted at 5 m/s shows
the kink at cloud base and the water it rains out.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import (
    dry_adiabatic_lapse_rate,
    moist_adiabatic_lapse_rate,
    parcel_temperature_profile,
    simulate_moist_ascent,
)
from earthsciencekit.constants import CP_DRY_AIR, GRAVITY, R_DRY_AIR

pressure = np.geomspace(1e5, 2.5e4, 161)
starts = [280.0, 290.0, 300.0, 310.0]
moist = [parcel_temperature_profile(pressure, t, t) for t in starts]
dry = [t * (pressure / 1e5) ** (R_DRY_AIR / CP_DRY_AIR) for t in starts]
assert all(np.all(m[1:] > d[1:]) for m, d in zip(moist, dry))
assert moist_adiabatic_lapse_rate(300, 1e5) < 0.4 * dry_adiabatic_lapse_rate(GRAVITY, CP_DRY_AIR)

time = np.linspace(0, 2400, 241)
ascent = simulate_moist_ascent(300, 293, 1e5, 5, time)
temperature, parcel_pressure, vapor, condensed = ascent.values.T
assert np.allclose(vapor + condensed, vapor[0])
print(f"water rained out after {time[-1] * 5 / 1e3:.0f} km: {condensed[-1] * 1e3:.1f} g/kg")

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 5))
for t, m, d in zip(starts, moist, dry):
    left.plot(d, pressure / 100, color="tab:orange", lw=1)
    left.plot(m, pressure / 100, color="tab:blue", lw=1)
left.plot([], [], color="tab:orange", label="dry adiabats")
left.plot([], [], color="tab:blue", label="pseudo-adiabats")
left.set(xlabel="temperature [K]", ylabel="pressure [hPa]", yscale="log", title="Hertz's diagram for moist air")
left.set_yticks([1000, 850, 700, 500, 300, 250], labels=["1000", "850", "700", "500", "300", "250"])
left.invert_yaxis()
left.legend()
height = time * 5 / 1e3
right.plot(temperature, height, label="parcel temperature [K]")
twin = right.twiny()
twin.plot(condensed * 1e3, height, color="tab:green", label="condensed water [g/kg]")
right.set(xlabel="temperature [K]", ylabel="height [km]", title="A parcel lifted through cloud base")
twin.set_xlabel("condensed water [g/kg]", color="tab:green")
fig.tight_layout()

plt.show()
