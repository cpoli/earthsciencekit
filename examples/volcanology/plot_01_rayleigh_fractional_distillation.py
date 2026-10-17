"""
Rayleigh's fractional distillation (1896)
=========================================

Rayleigh asked how far a mixture can be separated when the product is removed
as soon as it forms. Each increment is only slightly richer in one component,
but because the remainder keeps losing it, the remaining reservoir drifts as
the fraction left raised to the power ``alpha - 1``. Evaporating water is the
classic case: the vapour is poorer in oxygen-18, so the last drops of a drying
pool are strongly enriched.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.volcanology import delta_value, rayleigh_fractionation

vsmow, alpha = 2005.2e-6, 1 / 1.0094  # 18O/16O standard; vapour/liquid factor near 25 degC
remaining = np.linspace(1, 0.02, 400)
liquid = rayleigh_fractionation(vsmow, remaining, alpha)
vapour = alpha * liquid
# Mass balance: the vapour collected so far holds what the liquid has lost.
pooled = vsmow * (1 - remaining[1:] ** alpha) / (1 - remaining[1:])
assert np.allclose(remaining[1:] * liquid[1:] + (1 - remaining[1:]) * pooled, vsmow)
assert delta_value(liquid[-1], vsmow) > 30

fig, ax = plt.subplots()
ax.plot(remaining, delta_value(liquid, vsmow), label="remaining water")
ax.plot(remaining, delta_value(vapour, vsmow), label="instantaneous vapour")
ax.plot(remaining[1:], delta_value(pooled, vsmow), label="all vapour collected")
ax.invert_xaxis()
ax.set(xlabel="fraction of water remaining", ylabel="delta 18O [per mil VSMOW]", title="Rayleigh distillation of evaporating water")
ax.legend()

plt.show()
