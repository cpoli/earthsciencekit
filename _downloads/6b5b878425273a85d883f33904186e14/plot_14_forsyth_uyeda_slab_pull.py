"""
Forsyth and Uyeda's slab pull (1975)
====================================

Donald Forsyth and Seiya Uyeda compared the forces on each plate and found
that plates attached to sinking slabs move fastest. Per metre of trench, the
excess weight of a cold slab exceeds ridge push by about an order of
magnitude. Ridge push uses the half-space formula of Turcotte and Schubert.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import GRAVITY
from earthsciencekit.plate_tectonics import slab_pull_force
from earthsciencekit.units import convert

mantle, water, expansivity, excess_temperature, diffusivity = 3300.0, 1000.0, 3e-5, 1280.0, 1e-6
age = convert(100e6, "year", "s")
buoyancy = 2 * mantle * expansivity * excess_temperature / (np.pi * (mantle - water))
ridge_push = GRAVITY * mantle * expansivity * excess_temperature * diffusivity * age * (1 + buoyancy)
density_excess = mantle * expansivity * excess_temperature / 2
length = np.linspace(0, 700e3, 71)
pull = slab_pull_force(density_excess, 100e3 * length, GRAVITY)
assert float(slab_pull_force(density_excess, 100e3 * 600e3, GRAVITY)) > 5 * ridge_push

fig, ax = plt.subplots()
ax.plot(length / 1e3, pull / 1e12, label="slab pull, 100 km thick slab")
ax.axhline(ridge_push / 1e12, color="tab:red", ls="--", label="ridge push, 100 million year plate")
ax.set(xlabel="slab length [km]", ylabel="force per trench length [10^12 N/m]", title="Slabs pull harder than ridges push")
ax.legend()

plt.show()
