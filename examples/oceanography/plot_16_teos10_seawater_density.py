"""
TEOS-10 (2010)
==============

The Thermodynamic Equation of Seawater 2010 replaced the 1980 equation of
state with one derived from a Gibbs function, in Absolute Salinity and
Conservative Temperature. Its density is curved in temperature: mixing two
waters of equal density gives a denser mixture (cabbeling), which no linear
equation of state can show. The 75-term polynomial of the GSW toolbox
reproduces TEOS-10; the simplified fit of Roquet et al. (2015) keeps its
gradients but not its absolute density, and a linear law keeps neither.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from earthsciencekit.oceanography import linear_density, simplified_seawater_density, teos10_density
from earthsciencekit.units import convert

cold = (convert(2.0, "degC", "K"), 34.6)
target = teos10_density(*cold, 0)
warm_salinity = brentq(lambda s: teos10_density(convert(20.0, "degC", "K"), s, 0) - target, 30, 45, xtol=1e-14)
warm = (convert(20.0, "degC", "K"), warm_salinity)
mixed = teos10_density((cold[0] + warm[0]) / 2, (cold[1] + warm[1]) / 2, 0)
assert mixed > target + 0.1
print(f"mixing two waters of {target:.2f} kg/m^3 gives {mixed:.2f} kg/m^3")

# A linear law fitted at the cold water's state, from the TEOS-10 derivatives.
alpha = -(teos10_density(cold[0] + 0.01, cold[1], 0) - teos10_density(cold[0] - 0.01, cold[1], 0)) / 0.02 / target
beta = (teos10_density(cold[0], cold[1] + 0.01, 0) - teos10_density(cold[0], cold[1] - 0.01, 0)) / 0.02 / target


def linear(kelvin, salinity):
    return linear_density(target, kelvin - cold[0], salinity - cold[1], alpha, beta)


assert np.isclose(linear((cold[0] + warm[0]) / 2, (cold[1] + warm[1]) / 2), (linear(*cold) + linear(*warm)) / 2)

salinity, celsius = np.meshgrid(np.linspace(33, 40, 141), np.linspace(-2, 28, 121))
kelvin = convert(celsius, "degC", "K")
exact = teos10_density(kelvin, salinity, 0)
offset = simplified_seawater_density(kelvin, salinity, 0) - exact
assert np.abs(offset).max() < 3

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
left.contour(salinity, celsius, exact, levels=[target], colors="tab:blue")
left.contour(salinity, celsius, simplified_seawater_density(kelvin, salinity, 0), levels=[target], colors="tab:green", linestyles="-.")
left.contour(salinity, celsius, linear(kelvin, salinity), levels=[target], colors="tab:orange", linestyles="--")
left.plot([cold[1], warm[1]], [2, 20], "k:o")
left.plot((cold[1] + warm[1]) / 2, 11, "r*", markersize=12, label=f"mixture, +{mixed - target:.2f} kg/m^3")
left.plot([], [], color="tab:blue", label="isopycnal, TEOS-10")
left.plot([], [], color="tab:green", ls="-.", label="same density, simplified fit")
left.plot([], [], color="tab:orange", ls="--", label="isopycnal, linear")
left.set(xlabel="Absolute Salinity [g/kg]", ylabel="Conservative Temperature [degC]", title="Cabbeling: mixing makes water denser")
left.legend(loc="lower right", fontsize="small")
image = right.contourf(salinity, celsius, offset, levels=21, cmap="RdBu_r")
fig.colorbar(image, ax=right, label="simplified - TEOS-10 [kg/m^3]")
right.set(xlabel="Absolute Salinity [g/kg]", title="The simplified fit's absolute error at the surface")
fig.tight_layout()

plt.show()
