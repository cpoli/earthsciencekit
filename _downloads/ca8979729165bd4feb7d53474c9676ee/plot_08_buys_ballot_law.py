"""
Buys Ballot's law (1857)
========================

Christophorus Buys Ballot stated the rule from weather maps: in the Northern
Hemisphere, stand with your back to the wind and low pressure is on your
left. The wind is geostrophic, blowing along the isobars rather than down the
pressure gradient, counter-clockwise round a northern low. The geostrophic
components are ``v = (dp/dx) / (rho f)`` and ``u = -(dp/dy) / (rho f)``. The
check is that, at every point, the wind is perpendicular to the gradient and
the low lies to its left.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import coriolis_parameter, geostrophic_speed

x = np.linspace(-1500e3, 1500e3, 31)
X, Y = np.meshgrid(x, x)
pressure = 101325 - 2000 * np.exp(-(X**2 + Y**2) / (2 * (500e3) ** 2))
dpdy, dpdx = np.gradient(pressure, x, x)
f = coriolis_parameter(np.deg2rad(50))
u, v = -geostrophic_speed(dpdy, 1.2, f), geostrophic_speed(dpdx, 1.2, f)
np.testing.assert_allclose(u * dpdx + v * dpdy, 0, atol=1e-12)
to_low = -np.stack([dpdx, dpdy])
assert np.all(u * to_low[1] - v * to_low[0] >= 0)

fig, ax = plt.subplots(figsize=(5.5, 5))
contours = ax.contour(X / 1e3, Y / 1e3, pressure / 100, levels=10, colors="0.5")
ax.clabel(contours, fmt="%.0f")
ax.quiver(X[::3, ::3] / 1e3, Y[::3, ::3] / 1e3, u[::3, ::3], v[::3, ::3])
ax.set(xlabel="east [km]", ylabel="north [km]", title="Geostrophic wind round a low at 50 °N", aspect="equal")

plt.show()
