"""
Fourier's greenhouse (1824)
===========================

Joseph Fourier argued that the Earth would be much colder if it simply
radiated away the sunlight it absorbs, and that the atmosphere must hold heat
in. A planet in radiative balance emits what it absorbs, which for the
present albedo of about 0.3 gives an effective temperature near 255 K. The
measured mean surface temperature of about 288 K is 33 K warmer; that gap is
the greenhouse effect.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import absorbed_solar, effective_temperature, outgoing_longwave

albedo = np.linspace(0, 0.8, 201)
flux = absorbed_solar(1361, albedo)
temperature = effective_temperature(flux, 1)
np.testing.assert_allclose(outgoing_longwave(temperature, 1), flux, rtol=1e-14)
earth = float(effective_temperature(absorbed_solar(1361, 0.3), 1))
assert np.isclose(earth, 254.6, atol=0.1)

fig, ax = plt.subplots()
ax.plot(albedo, temperature, label="radiative balance, no atmosphere")
ax.axhline(288, color="tab:red", ls="--", label="observed mean surface")
ax.annotate("", xy=(0.3, 288), xytext=(0.3, earth), arrowprops={"arrowstyle": "<->"})
ax.text(0.31, (288 + earth) / 2, f"greenhouse gap {288 - earth:.0f} K", va="center")
ax.plot(0.3, earth, "ko")
ax.set(xlabel="planetary albedo", ylabel="temperature [K]", title="The Earth is warmer than radiative balance allows")
ax.legend(loc="lower left")

plt.show()
