"""
Beven and Kirkby's TOPMODEL (1979)
==================================

Keith Beven and Mike Kirkby argued that runoff comes mostly from saturated
ground that expands and contracts with the wetness of a catchment, and that
where it forms is predictable from the topography. Their topographic index
``ln(a / tan(beta))`` is high where a large area drains through a gentle
slope, at the foot of hillslopes and in hollows. One catchment-mean deficit
then sets the local deficit everywhere, the saturated area and an
exponentially receding baseflow.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import topmodel_baseflow, topmodel_deficit, topographic_index

x = np.linspace(1, 400, 400)  # m from the divide along a hillslope
y = np.linspace(-1, 1, 81)[:, None]  # across the hillslope; convergent toward y = 0
area = x * (1 + 2 * np.exp(-((y / 0.3) ** 2)) * x / 400)  # upslope area per unit contour width [m]
slope = np.arctan(0.02 + 0.25 * np.sin(np.pi * x / 400) ** 2 + 0.01 * np.abs(y))
index = topographic_index(area, slope)
mean_index, scaling = index.mean(), 0.03
mean_deficit = np.linspace(0, 0.15, 61)
saturated = [np.mean(topmodel_deficit(d, index, mean_index, scaling) <= 0) for d in mean_deficit]
assert np.all(np.diff(saturated) <= 0) and saturated[0] > 0.3 and saturated[-1] == 0
baseflow = topmodel_baseflow(5e-3, mean_index, mean_deficit, scaling)
assert np.isclose(baseflow[20] / baseflow[0], np.exp(-mean_deficit[20] / scaling))

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
image = left.pcolormesh(x, y.ravel(), index, shading="auto", cmap="Blues")
left.contour(x, y.ravel(), topmodel_deficit(0.05, index, mean_index, scaling), levels=[0], colors="k")
fig.colorbar(image, ax=left, label="ln(a / tan(beta))")
left.set(xlabel="distance from the divide [m]", ylabel="across slope", title="Saturated at a 50 mm mean deficit (inside line)")
right.plot(mean_deficit * 1e3, saturated, label="saturated fraction")
right.set(xlabel="catchment-mean deficit [mm]", ylabel="saturated fraction", title="Saturated area and baseflow")
twin = right.twinx()
twin.semilogy(mean_deficit * 1e3, baseflow * 3.6e6, "C1--", label="baseflow")
twin.set_ylabel("baseflow [mm/h]")
fig.tight_layout()

plt.show()
