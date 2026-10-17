"""
Halley's barometric altitude (1686)
===================================

Edmond Halley combined Boyle's law with the weight of the air and found that
height grows with the logarithm of pressure: equal ratios of pressure mark
equal steps of height. The step for a factor e is the scale height
``R T / g``, about 8 km. Inverting the barometric law turns a barometer into an
altimeter.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import barometric_pressure, scale_height
from earthsciencekit.constants import GRAVITY

height = np.linspace(0, 30e3, 301)
fig, ax = plt.subplots()
for temperature in [230, 260, 290]:
    H = float(scale_height(temperature, GRAVITY))
    pressure = barometric_pressure(101325, height, temperature, GRAVITY)
    np.testing.assert_allclose(H * np.log(101325 / pressure), height, atol=1e-9)
    assert np.isclose(barometric_pressure(101325, H, temperature, GRAVITY), 101325 / np.e)
    ax.semilogx(pressure / 100, height / 1000, label=f"{temperature} K, H = {H / 1000:.2f} km")
ax.set(xlabel="pressure [hPa, log scale]", ylabel="height [km]", title="Height is linear in log pressure")
ax.legend()

plt.show()
