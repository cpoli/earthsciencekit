"""
Manning's formula (1889)
========================

Robert Manning compared the many empirical formulas proposed after Chézy with
the gaugings of rivers and canals, and found that the Chézy coefficient grows
roughly as the sixth root of the hydraulic radius for a given bed. The mean
velocity then becomes ``R**(2/3) S**(1/2) / n``, with one roughness number per
kind of channel. The example shows that Chézy's ``C`` is not a constant, and
that Manning's formula is Chézy's with ``C = R**(1/6) / n``.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import chezy_from_manning, chezy_velocity, manning_velocity

radius = np.logspace(-1, 1, 100)
fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4))
for roughness, bed in [(0.012, "concrete"), (0.03, "clean river"), (0.07, "weedy river")]:
    coefficient = chezy_from_manning(roughness, radius)
    velocity = manning_velocity(roughness, radius, 1e-3)
    assert np.allclose(velocity, chezy_velocity(coefficient, radius, 1e-3))
    left.semilogx(radius, coefficient, label=f"{bed}, n = {roughness}")
    right.loglog(radius, velocity, label=bed)
assert np.isclose(manning_velocity(0.03, 8.0, 1e-3) / manning_velocity(0.03, 1.0, 1e-3), 4)
left.set(xlabel="hydraulic radius [m]", ylabel="Chézy C [m^0.5/s]", title="C = R^(1/6) / n")
right.set(xlabel="hydraulic radius [m]", ylabel="velocity [m/s]", title="v = R^(2/3) S^(1/2) / n, S = 1e-3")
left.legend(fontsize="small")
fig.tight_layout()

plt.show()
