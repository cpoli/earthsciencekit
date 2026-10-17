"""
Larmor's self-exciting dynamo (1919)
====================================

Joseph Larmor asked how a body such as the Sun could be a magnet, and proposed
that conducting fluid moving through a magnetic field could induce the currents
that sustain that field, like a self-exciting dynamo. The Earth needs one: a
field left in a still iron core would decay in tens of thousands of years, yet
rocks record a field for over three billion years. The core's flow of a fraction
of a millimetre per second gives a magnetic Reynolds number in the hundreds,
enough for induction to beat diffusion.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geophysics import magnetic_decay_time, magnetic_reynolds_number
from earthsciencekit.units import convert

core_radius, shell = 3.48e6, 2.26e6
conductivity = np.logspace(5, 6, 100)
decay = convert(magnetic_decay_time(core_radius, conductivity), "s", "kyr")
assert decay.max() < 100 and convert(3.5, "Ga", "kyr") / decay.max() > 1e4
speed = np.logspace(-5, -2, 100)
reynolds = magnetic_reynolds_number(speed, shell, 5e5)
assert magnetic_reynolds_number(5e-4, shell, 5e5) > 100

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.loglog(conductivity, decay)
left.set(xlabel="core conductivity [S/m]", ylabel="free decay time [kyr]", title="A still core forgets its field")
right.loglog(convert(speed, "m/s", "cm/year"), reynolds)
right.axhspan(10, 100, color="0.85", label="typical dynamo threshold")
right.axvline(convert(5e-4, "m/s", "cm/year"), color="C3", ls="--", label="core flow from secular variation")
right.set(xlabel="flow speed [cm/year]", ylabel="magnetic Reynolds number", title="Induction beats diffusion")
right.legend(fontsize="small")
fig.tight_layout()

plt.show()
