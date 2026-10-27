"""
NUVEL-1 relative plate motions (1990)
=====================================

Charles DeMets, Richard Gordon, Donald Argus and Seth Stein fitted spreading
rates, transform azimuths and earthquake slip vectors for all the major
plates at once. The fit requires every circuit of relative Euler vectors to
close: the motion of A relative to C is A-B plus B-C. At a triple junction the
three relative velocities then form a closed triangle. The Euler vectors here
are synthetic.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import EARTH_RADIUS
from earthsciencekit.plate_tectonics import compose_euler_vectors, euler_velocity_x, euler_velocity_y, euler_velocity_z
from earthsciencekit.units import convert


def euler_vector(latitude, longitude, degrees_per_myr):
    rate = convert(degrees_per_myr, "degree", "rad") / convert(1e6, "year", "s")
    lat, lon = np.deg2rad(latitude), np.deg2rad(longitude)
    return rate * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def velocity(omega, r):
    return np.array(
        [
            euler_velocity_x(omega[1], omega[2], r[1], r[2]),
            euler_velocity_y(omega[2], omega[0], r[2], r[0]),
            euler_velocity_z(omega[0], omega[1], r[0], r[1]),
        ]
    )


omega_ab = euler_vector(48.7, -78.2, 0.78)
omega_bc = euler_vector(-10.0, -60.0, 0.5)
omega_ac = compose_euler_vectors(omega_ab, omega_bc)
np.testing.assert_allclose(compose_euler_vectors(omega_ac, -omega_ac), 0)

lat, lon = np.deg2rad(20.0), np.deg2rad(-110.0)
r = EARTH_RADIUS * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
east = np.array([-np.sin(lon), np.cos(lon), 0])
north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
legs = {"A relative to B": velocity(omega_ab, r), "B relative to C": velocity(omega_bc, r), "C relative to A": velocity(-omega_ac, r)}
np.testing.assert_allclose(sum(legs.values()), 0, atol=1e-20)

fig, ax = plt.subplots()
corners = [np.zeros(2)]
for name, v in legs.items():
    step = convert(np.array([east @ v, north @ v]), "m/s", "cm/year")
    ax.annotate("", corners[-1] + step, corners[-1], arrowprops={"arrowstyle": "->", "lw": 2})
    ax.text(*(corners[-1] + step / 2), name, ha="center", va="bottom")
    corners.append(corners[-1] + step)
assert np.allclose(corners[-1], 0, atol=1e-12)
ax.plot(*np.array(corners).T, "o", color="k", ms=4)
ax.set(xlabel="east [cm/year]", ylabel="north [cm/year]", title="The velocity triangle closes", aspect="equal")
ax.margins(0.25)

plt.show()
