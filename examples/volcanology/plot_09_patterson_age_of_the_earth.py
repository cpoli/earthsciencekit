"""
Patterson's age of the Earth (1956)
===================================

Clair Patterson measured lead in meteorites. The iron meteorite Canyon Diablo
holds almost no uranium, so its lead is primordial; stony meteorites add lead
grown from uranium since they formed. Plotted as 207Pb/204Pb against
206Pb/204Pb they fall on one line whose slope is the radiogenic 207Pb/206Pb
ratio of a 4.55-billion-year age, and terrestrial ocean sediment lies on it
too, so the Earth shares that age.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.units import convert
from earthsciencekit.volcanology import lead_lead_age, radiogenic_lead_ratio

primordial = np.array([9.307, 10.294])  # Canyon Diablo 206Pb/204Pb, 207Pb/204Pb
slope = radiogenic_lead_ratio(convert(4.55, "Ga", "s"))
growth = np.array([0.0, 10.0, 20.0, 30.0, 40.0])  # synthetic radiogenic 206Pb/204Pb added
x = primordial[0] + growth
y = primordial[1] + slope * growth
fit = np.polyfit(x, y, 1)[0]
assert np.isclose(convert(lead_lead_age(fit), "s", "Ga"), 4.55)

fig, ax = plt.subplots()
ax.plot(x, y, "ko")
ax.plot(x[0], y[0], "s", color="C3", label="Canyon Diablo (primordial)")
ax.plot(x, primordial[1] + fit * (x - primordial[0]), label=f"isochron: slope {fit:.3f}, {convert(lead_lead_age(fit), 's', 'Ga'):.2f} Ga")
ax.set(xlabel="206Pb / 204Pb", ylabel="207Pb / 204Pb", title="Meteorite lead isochron")
ax.legend()

plt.show()
