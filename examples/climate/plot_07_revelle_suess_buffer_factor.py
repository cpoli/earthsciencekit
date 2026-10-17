"""
Revelle and Suess's buffer factor (1957)
========================================

Roger Revelle and Hans Suess showed that seawater resists taking up carbon
dioxide. Because of carbonate chemistry, a 10% rise in the CO2 pressure of
the air is balanced by only about a 1% rise in the dissolved carbon of the
surface ocean. This ratio, about 10, is the buffer (Revelle) factor. In
equilibrium the share of added carbon left in the air is ``1 / (1 + M_o /
(R M_a))`` for atmospheric and oceanic inventories ``M_a`` and ``M_o``, so
the ocean holds far less of it than a simple partition by inventory would
suggest. Inventories here are round illustrative values.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.climate import carbon_airborne_fraction

atmosphere, emission = 600e12, 100e12  # kg C
buffer = np.linspace(1, 15, 141)
fig, ax = plt.subplots()
for ocean, label in [(1000e12, "surface ocean only"), (38000e12, "whole ocean, after mixing")]:
    increase = emission / (1 + ocean / (buffer * atmosphere))
    fraction = carbon_airborne_fraction(increase, emission)
    np.testing.assert_allclose(increase / atmosphere, buffer * (emission - increase) / ocean)
    ax.plot(buffer, fraction, label=label)
    assert np.isclose(fraction[0], atmosphere / (atmosphere + ocean))
ax.axvline(10, color="0.5", ls="--", label="Revelle factor of seawater")
ax.set(xlabel="buffer factor R", ylabel="equilibrium airborne fraction", title="A buffered ocean leaves more CO2 in the air")
ax.legend()

plt.show()
