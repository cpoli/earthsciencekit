"""
Love waves (1911)
=================

Augustus Love showed that horizontal shear waves can be trapped in a slow
surface layer over a faster half-space, which a homogeneous Earth cannot do.
Their speed depends on period: short waves travel at the speed of the layer,
long waves at the speed of the half-space. That dispersion measures the
thickness of the crust. The layer here is a 35 km crust over mantle.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import love_wave_angular_frequency, s_wave_speed

h = 35e3
b1, b2 = 3.5e3, 4.5e3
m1, m2 = 2800 * b1**2, 3300 * b2**2
assert np.isclose(s_wave_speed(m1, 2800), b1)
speed = np.linspace(b1 * 1.001, b2 * 0.999, 400)
period = 2 * np.pi / love_wave_angular_frequency(speed, h, b1, b2, m1, m2)
assert np.all(np.diff(period) > 0)

fig, ax = plt.subplots()
for thickness in [20e3, 35e3, 50e3]:
    omega = love_wave_angular_frequency(speed, thickness, b1, b2, m1, m2)
    ax.semilogx(2 * np.pi / omega, speed / 1e3, label=f"crust {thickness / 1e3:.0f} km")
ax.axhline(b1 / 1e3, color="k", ls=":")
ax.axhline(b2 / 1e3, color="k", ls=":")
ax.set(xlabel="period [s]", ylabel="phase speed [km/s]", title="Long Love waves feel the mantle")
ax.legend()

plt.show()
