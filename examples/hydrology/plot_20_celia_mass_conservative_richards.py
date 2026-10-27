"""
Celia's mass-conservative Richards solver (1990)
================================================

Michael Celia, Efthimios Bouloutas and Rebecca Zarba showed why numerical
solutions of Richards's equation lost or gained water: they stepped the
pressure head through a capacity term that cannot be evaluated consistently
over a time step. Writing the storage in terms of water content, the mixed
form, removed the error. Here water content is the state of every cell. Rain
at half the saturated conductivity wets a 3 m column from the top, the
wetting front reaches the water table, and the profile settles on Richards's
steady solution while the water balance closes.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import richards_steady_head, simulate_richards_column

day = 86400.0
flux, conductivity, alpha = -5e-6, 1e-5, 3.0
time = np.linspace(0, 30 * day, 3001)
result = simulate_richards_column(3.0, flux, conductivity, alpha, 0.42, 0.04, time, nodes=120)
storage, bottom, head = result.values[:, 0], result.values[:, 1], result.values[:, 2:]
height = (np.arange(120) + 0.5) * 3.0 / 120

drained = np.concatenate([[0.0], np.cumsum((bottom[1:] + bottom[:-1]) / 2 * np.diff(time))])
error = storage - storage[0] - (drained - flux * time)
steady = richards_steady_head(height, flux, conductivity, alpha)
assert np.abs(error).max() < 1e-5 and np.abs(head[-1] - steady).max() < 1e-4
print(f"largest water-balance error {np.abs(error).max() * 1e3:.4f} mm of {(storage[-1] - storage[0]) * 1e3:.0f} mm stored")

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
for days in (0, 1, 3, 6, 10, 30):
    left.plot(head[np.searchsorted(time, days * day)], height, label=f"day {days}")
left.plot(steady, height, "k:", lw=2, label="steady Gardner profile")
left.set(xlabel="pressure head [m]", ylabel="height above the water table [m]", title="A wetting front reaches the water table")
left.legend(fontsize="small")
right.plot(time / day, (storage - storage[0]) * 1e3, label="change in storage")
right.plot(time / day, (drained - flux * time) * 1e3, "--", label="rain in minus drainage out")
right.plot(time / day, -flux * time * 1e3, ":", label="rain in")
right.set(xlabel="days", ylabel="mm of water", title="The water balance closes")
right.legend(fontsize="small")
fig.tight_layout()

plt.show()
