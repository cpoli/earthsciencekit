"""
Playfair's law of accordant junctions (1802)
============================================

John Playfair, explaining Hutton's theory, observed that every tributary
joins its main valley at exactly the level of the main stream, neither
hanging above it nor cutting below. Only rivers that had carved their own
valleys could fit together so precisely. In a steady landscape eroded by
stream power the elevation is proportional to the chi coordinate, which
integrates distance scaled by drainage area, so a tributary and its trunk
fall on one straight line in chi and meet at the same height.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import SECONDS_PER_YEAR
from earthsciencekit.geomorphology import chi_coordinate, steady_state_slope

uplift, erodibility = 1e-3 / SECONDS_PER_YEAR, 5.4e-6 / SECONDS_PER_YEAR
junction = 15e3  # m upstream of the outlet on the trunk
trunk = np.linspace(0, 40e3, 401)
tributary = np.linspace(0, 12e3, 121)  # m upstream of the junction
tributary_area = ((12e3 - tributary + 300) / 0.32) ** (1 / 0.6)
upper_trunk_area = ((40e3 - trunk + 300) / 0.32) ** (1 / 0.6)
trunk_area = upper_trunk_area + np.where(trunk <= junction, tributary_area[0], 0.0)


def graded(distance, area, base):
    """Integrate the steady-state slope upstream from a base elevation."""
    slope = steady_state_slope(uplift, 0.5 * (area[1:] + area[:-1]), erodibility, 0.5, 1)
    return base + np.concatenate([[0.0], np.cumsum(slope * np.diff(distance))])


trunk_z = graded(trunk, trunk_area, 0.0)
at_junction = np.interp(junction, trunk, trunk_z)
tributary_z = graded(tributary, tributary_area, at_junction)
trunk_chi = chi_coordinate(trunk, trunk_area, 1e6, 0.5)
tributary_chi = np.interp(junction, trunk, trunk_chi) + chi_coordinate(tributary, tributary_area, 1e6, 0.5)
steepness = uplift / (erodibility * 1e3)
assert tributary_z[0] == at_junction
assert np.allclose(trunk_z, steepness * trunk_chi, rtol=0.03, atol=1)
assert np.allclose(tributary_z, steepness * tributary_chi, rtol=0.03, atol=1)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(trunk / 1e3, trunk_z, label="trunk")
left.plot((junction + tributary) / 1e3, tributary_z, label="tributary")
left.plot(junction / 1e3, at_junction, "ko", label="junction")
left.set(xlabel="distance along the stream from the outlet [km]", ylabel="elevation [m]", title="Long profiles")
left.legend()
right.plot(trunk_chi / 1e3, trunk_z, label="trunk")
right.plot(tributary_chi / 1e3, tributary_z, "--", label="tributary")
right.plot(np.interp(junction, trunk, trunk_chi) / 1e3, at_junction, "ko")
right.set(xlabel="chi [km]", title="Playfair: one line in chi, accordant at the junction")
right.legend()
fig.tight_layout()

plt.show()
