"""
Zoeppritz's reflection coefficients (1919)
==========================================

Karl Zoeppritz derived how a plane wave divides at a boundary between two
elastic media. At normal incidence his equations reduce to the impedance
contrast: the pressure transmission can exceed one, yet the reflected and
transmitted energy fractions always sum to one.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import reflection_coefficient, seismic_impedance, transmission_coefficient

ratio = np.logspace(-2, 2, 201)
r = reflection_coefficient(1, ratio)
t = transmission_coefficient(1, ratio)
reflected, transmitted = r * r, t * t / ratio
np.testing.assert_allclose(reflected + transmitted, 1, atol=1e-14)

# A sediment-on-basement contact.
z1, z2 = seismic_impedance(2200, 2500), seismic_impedance(2700, 5500)
print(f"sediment/basement: R = {reflection_coefficient(z1, z2):.3f}, T = {transmission_coefficient(z1, z2):.3f}")

fig, ax = plt.subplots()
ax.semilogx(ratio, reflected, label="reflected energy")
ax.semilogx(ratio, transmitted, label="transmitted energy")
ax.axvline(z2 / z1, color="k", ls=":", label="sediment on basement")
ax.set(xlabel="impedance ratio Z2/Z1", ylabel="incident energy fraction", title="Energy is shared, never created")
ax.legend()

plt.show()
