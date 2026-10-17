"""
Euler's rotation theorem (1776)
===============================

Any motion of a rigid cap on a sphere is a rotation about an axis through the
centre. The velocity of every point is the cross product of the Euler vector
and its position, so it is tangent to the sphere and vanishes at the pole. The
pole and rate here are synthetic.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import EARTH_RADIUS
from earthsciencekit.plate_tectonics import angular_speed, euler_velocity_x, euler_velocity_y, euler_velocity_z
from earthsciencekit.units import convert

rate = angular_speed(convert(1.0, "degree", "rad"), convert(1e6, "year", "s"))
pole_lat, pole_lon = np.deg2rad(50.0), np.deg2rad(-30.0)
omega = rate * np.array([np.cos(pole_lat) * np.cos(pole_lon), np.cos(pole_lat) * np.sin(pole_lon), np.sin(pole_lat)])

lon, lat = np.meshgrid(np.deg2rad(np.arange(-180, 181, 20)), np.deg2rad(np.arange(-75, 76, 15)))
r = EARTH_RADIUS * np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
v = np.stack(
    [
        euler_velocity_x(omega[1], omega[2], r[1], r[2]),
        euler_velocity_y(omega[2], omega[0], r[2], r[0]),
        euler_velocity_z(omega[0], omega[1], r[0], r[1]),
    ]
)
np.testing.assert_allclose(np.sum(v * r, axis=0), 0, atol=1e-9)
pole = EARTH_RADIUS * omega / rate
assert np.allclose(
    [euler_velocity_x(omega[1], omega[2], pole[1], pole[2]), euler_velocity_z(omega[0], omega[1], pole[0], pole[1])], 0, atol=1e-20
)

east = -np.sin(lon) * v[0] + np.cos(lon) * v[1]
north = -np.sin(lat) * np.cos(lon) * v[0] - np.sin(lat) * np.sin(lon) * v[1] + np.cos(lat) * v[2]
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.quiver(np.rad2deg(lon), np.rad2deg(lat), convert(east, "m/s", "cm/year"), convert(north, "m/s", "cm/year"), color="tab:blue")
ax.plot(np.rad2deg(pole_lon), np.rad2deg(pole_lat), "r*", ms=14, label="Euler pole")
ax.set(xlabel="longitude [degree]", ylabel="latitude [degree]", title="Rotation of 1 degree per million years")
ax.legend(loc="lower left")

plt.show()
