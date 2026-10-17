"""
Poisson's two elastic waves (1829)
==================================

Siméon Denis Poisson showed that an elastic solid carries two kinds of wave:
a faster compressional wave and a slower shear wave. When the two Lamé
constants are equal (a Poisson solid) the speed ratio is exactly sqrt(3) and
Poisson's ratio is 1/4. With no rigidity the shear wave vanishes, as in a fluid.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import p_wave_speed, poisson_ratio, s_wave_speed

density, bulk = 3300.0, 1.3e11
shear = np.linspace(0, 1.0e11, 201)
vp = p_wave_speed(bulk, shear, density)
vs = s_wave_speed(shear, density)

# Poisson solid: lambda = mu, so the bulk modulus is 5 mu / 3.
mu = 6e10
vp_poisson, vs_poisson = p_wave_speed(5 * mu / 3, mu, density), s_wave_speed(mu, density)
assert np.isclose(vp_poisson / vs_poisson, np.sqrt(3))
assert np.isclose(poisson_ratio(vp_poisson, vs_poisson), 0.25)
assert vs[0] == 0 and vp[0] > 0

fig, ax = plt.subplots()
ax.plot(shear / 1e9, vp / 1e3, label="P (compressional)")
ax.plot(shear / 1e9, vs / 1e3, label="S (shear)")
ax.axvline(mu / 1e9, color="k", ls=":", lw=1)
ax.set(xlabel="shear modulus [GPa]", ylabel="wave speed [km/s]", title=f"Bulk modulus {bulk / 1e9:.0f} GPa: no rigidity, no S wave")
ax.legend()

plt.show()
