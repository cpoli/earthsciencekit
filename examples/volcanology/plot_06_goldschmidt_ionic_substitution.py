"""
Goldschmidt's geochemical rules (1937)
======================================

Victor Goldschmidt classified the elements by where they go in the Earth and
stated the rules for which ions replace one another in a crystal: ions of the
same charge whose radii differ by less than about 15% substitute freely. Blundy
and Wood (1994) made the rule quantitative: the partition coefficient falls off
on both sides of the radius that fits the site, as the strain energy of the
lattice grows.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.volcanology import lattice_strain_partition

optimal, peak, modulus, temperature = 1.0e-10, 2.0, 250e9, 1500  # a clinopyroxene M2-like site
radius = np.linspace(0.6e-10, 1.45e-10, 400)
ions = {"Mg2+": 0.89e-10, "Fe2+": 0.92e-10, "Ca2+": 1.12e-10, "Sr2+": 1.26e-10, "Ba2+": 1.42e-10}  # 8-fold radii
partition = lattice_strain_partition(radius, optimal, peak, modulus, temperature)
assert np.isclose(partition.max(), peak, rtol=1e-4)
assert lattice_strain_partition(ions["Ba2+"], optimal, peak, modulus, temperature) < 1e-3 * peak

fig, ax = plt.subplots()
ax.semilogy(radius * 1e10, partition)
ax.axvspan(0.85, 1.15, color="0.9", label="within 15% of the site radius")
for name, value in ions.items():
    d = lattice_strain_partition(value, optimal, peak, modulus, temperature)
    ax.plot(value * 1e10, d, "ko")
    ax.annotate(name, (value * 1e10, d), textcoords="offset points", xytext=(4, 4), fontsize=8)
ax.set(xlabel="ionic radius [angstrom]", ylabel="mineral/melt partition coefficient", title="Size decides substitution", ylim=(1e-6, 5))
ax.legend()

plt.show()
