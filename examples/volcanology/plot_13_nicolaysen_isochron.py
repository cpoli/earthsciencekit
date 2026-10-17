"""
Nicolaysen's isochron (1961)
============================

Louis Nicolaysen showed how to date a rock without knowing how much daughter
it started with. Minerals that crystallized together share the initial
87Sr/86Sr but take up different amounts of rubidium. As 87Rb decays they move
along lines of slope one, so at any time they lie on a straight isochron whose
slope gives the age and whose intercept gives the initial ratio.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.units import convert
from earthsciencekit.volcanology import isochron_fit

rate = 1.42e-11 / convert(1, "year", "s")  # 87Rb decay constant, 1/s
age, initial = convert(1.2, "Ga", "s"), 0.7045
rubidium = np.array([0.15, 0.9, 2.1, 4.0, 6.5])
scatter = np.array([2, -3, 1, 3, -2]) * 1e-5  # deterministic analytical noise
strontium = initial + rubidium * np.expm1(rate * age) + scatter
fit_age, fit_initial = isochron_fit(rubidium, strontium, rate)
assert np.isclose(convert(fit_age, "s", "Ga"), 1.2, rtol=2e-3) and np.isclose(fit_initial, initial, atol=1e-4)

fig, ax = plt.subplots()
ax.plot(rubidium, np.full_like(rubidium, initial), "o", mfc="none", color="0.5", label="at crystallization")
ax.plot(rubidium, strontium, "ko", label="today")
line = np.array([0, 7.0])
ax.plot(
    line,
    fit_initial + line * np.expm1(rate * fit_age),
    label=f"isochron: {convert(fit_age, 's', 'Ga'):.3f} Ga, initial {fit_initial:.4f}",
)
ax.set(xlabel="87Rb / 86Sr", ylabel="87Sr / 86Sr", title="Rb-Sr mineral isochron")
ax.legend()

plt.show()
