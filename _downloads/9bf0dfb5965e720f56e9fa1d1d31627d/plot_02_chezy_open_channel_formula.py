"""
Chézy's open-channel formula (1775)
===================================

Antoine de Chézy was asked to size a canal bringing the water of the Yvette to
Paris. He reasoned that in uniform flow the weight of the water down the slope
balances a bed friction that grows as the square of the velocity, so the mean
velocity is ``C sqrt(R S)``, with ``R`` the cross-sectional area over the wetted
perimeter. The example evaluates a rectangular canal: deepening it raises the
hydraulic radius and the velocity, which tends to the wide-channel limit where
``R`` equals the depth.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import chezy_velocity, hydraulic_radius

depth = np.linspace(0.05, 3, 120)
fig, ax = plt.subplots()
for width in (1.0, 4.0, 20.0):
    radius = hydraulic_radius(width, depth)
    ax.plot(depth, chezy_velocity(50, radius, 4e-4), label=f"width {width:g} m")
wide = chezy_velocity(50, depth, 4e-4)
ax.plot(depth, wide, "k--", label="wide channel, R = depth")
assert np.all(chezy_velocity(50, hydraulic_radius(20.0, depth), 4e-4) < wide)
assert np.isclose(chezy_velocity(50, 1.0, 4e-4), 2 * chezy_velocity(50, 0.25, 4e-4))
ax.set(xlabel="flow depth [m]", ylabel="mean velocity [m/s]", title="Chézy: v = C sqrt(R S), C = 50, S = 4e-4")
ax.legend()

plt.show()
