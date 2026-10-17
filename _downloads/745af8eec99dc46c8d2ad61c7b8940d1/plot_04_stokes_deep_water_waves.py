"""
Stokes's deep-water waves (1847)
================================

George Gabriel Stokes worked out oscillatory waves of small and finite height.
The dispersion relation ``omega**2 = g k tanh(k h)`` makes long waves travel
faster than short ones. In deep water the energy travels at the group speed,
half the speed of the crests, which is why the swell from a distant storm
arrives longest waves first.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY
from earthsciencekit.oceanography import deep_group_speed, deep_wave_speed, shallow_wave_speed, wave_angular_frequency

depth = 50.0
wavelength = np.logspace(0, 4, 2001)
k = 2 * np.pi / wavelength
omega = wave_angular_frequency(k, depth, GRAVITY)
group = np.gradient(omega, k)
deep = wavelength < depth
np.testing.assert_allclose(group[deep][5:-5], deep_group_speed(wavelength[deep][5:-5], GRAVITY), rtol=1e-3)
np.testing.assert_allclose(deep_group_speed(wavelength, GRAVITY), 0.5 * deep_wave_speed(wavelength, GRAVITY))

fig, ax = plt.subplots()
ax.loglog(wavelength, omega / k, label="phase speed")
ax.loglog(wavelength, group, label="group speed")
ax.loglog(wavelength, deep_wave_speed(wavelength, GRAVITY), "k--", lw=0.8, label="deep-water limits")
ax.loglog(wavelength, deep_group_speed(wavelength, GRAVITY), "k--", lw=0.8)
ax.axhline(shallow_wave_speed(depth, GRAVITY), color="0.5", ls=":", label="shallow-water limit")
ax.set(xlabel="wavelength [m]", ylabel="speed [m/s]", title=f"Dispersion in water {depth:.0f} m deep")
ax.legend()

plt.show()
