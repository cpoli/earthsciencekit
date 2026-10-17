"""
Matuyama's reversed epoch (1929)
================================

Motonori Matuyama measured the magnetization of basalts from Japan, Korea and
Manchuria. The youngest flows were magnetized along today's field, but older,
early Quaternary flows pointed the opposite way. He concluded that the field
had been reversed for a long time before the present one, the interval now
called the Matuyama chron. A flow is reversed when its direction lies more than
90 degrees from the present axial dipole field at the site.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geophysics import dipole_inclination, direction_angle
from earthsciencekit.units import convert

normal_chrons = [(0.0, 0.773), (0.988, 1.072), (1.775, 1.934), (2.595, 3.0)]  # Ma: Brunhes, Jaramillo, Olduvai, Gauss
rng = np.random.default_rng(1929)
age = np.sort(rng.uniform(0.05, 2.9, 40))
normal = np.array([any(start <= t < end for start, end in normal_chrons) for t in age])
site = np.radians(35.0)
present = dipole_inclination(site)
declination = np.where(normal, 0, np.pi) + rng.normal(0, 0.15, age.size)
inclination = np.where(normal, present, -present) + rng.normal(0, 0.1, age.size)
angle = direction_angle(declination, inclination, 0, present)
assert np.array_equal(angle > np.pi / 2, ~normal)
assert np.all(angle[age < 0.773] < np.pi / 2)

fig, ax = plt.subplots(figsize=(8, 4))
for start, end in normal_chrons:
    ax.axvspan(start, end, color="0.85")
ax.plot(age[normal], convert(angle[normal], "rad", "degree"), "ko", mfc="k", label="normal")
ax.plot(age[~normal], convert(angle[~normal], "rad", "degree"), "ko", mfc="w", label="reversed")
ax.axhline(90, color="C3", ls="--")
ax.set(xlabel="age [Ma] (shaded: normal chrons)", ylabel="angle from present field [degrees]", title="Older flows are reversed")
ax.legend()

plt.show()
