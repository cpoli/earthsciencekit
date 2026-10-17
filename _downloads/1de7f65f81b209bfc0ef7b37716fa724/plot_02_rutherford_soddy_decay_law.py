"""
Rutherford and Soddy's decay law (1902)
=======================================

Rutherford and Soddy found that the activity of thorium emanation (radon-220)
falls by half every minute or so, whatever its temperature or chemistry, and
concluded that radioactivity is the spontaneous transformation of one element
into another at a fixed fractional rate. Decay is exponential, so a straight
line on a logarithmic plot and its slope gives the half-life.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.volcanology import decay_constant, parent_fraction

half_life = 55.6  # s, radon-220
time = np.linspace(0, 400, 81)
remaining = parent_fraction(decay_constant(half_life), time)
slope = np.polyfit(time, np.log(remaining), 1)[0]
assert np.isclose(np.log(2) / -slope, half_life)
assert np.isclose(parent_fraction(decay_constant(half_life), 2 * half_life), 0.25)

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
left.plot(time, remaining)
for n in range(1, 5):
    left.plot(n * half_life, 0.5**n, "ko", ms=4)
left.set(xlabel="time [s]", ylabel="fraction of radon-220 left", title="Half lost every 55.6 s")
right.semilogy(time, remaining)
right.set(xlabel="time [s]", ylabel="fraction left (log scale)", title=f"Straight line, half-life {np.log(2) / -slope:.1f} s")
fig.tight_layout()

plt.show()
