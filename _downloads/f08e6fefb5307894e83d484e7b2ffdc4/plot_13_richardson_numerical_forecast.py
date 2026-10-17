"""
Richardson's numerical forecast (1922)
======================================

In *Weather Prediction by Numerical Process*, Richardson computed by hand a
six-hour forecast for two points over Germany, from the observations of
20 May 1910. He predicted a surface pressure rise of 145 hPa, about a hundred
times too large. The pressure tendency is the column-mean divergence times the
surface pressure, and the divergence at any one level is much larger than the
column mean, because convergence below nearly cancels divergence above. A
column-mean error of 7e-6 per second, the size of the divergence at a single
level, is enough to reproduce his result.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import surface_pressure_tendency

six_hours = 6 * 3600
surface_pressure = 101325.0
richardson = float(surface_pressure_tendency(-6.6e-6, surface_pressure)) * six_hours
assert abs(richardson - 145e2) < 3e2
divergence = np.logspace(-8, -5, 200)
change = -surface_pressure_tendency(-divergence, surface_pressure) * six_hours / 100
assert np.isclose(-surface_pressure_tendency(1e-7, surface_pressure) * six_hours, 218.8, rtol=1e-3)
print(f"six-hour pressure change: {richardson / 100:.0f} hPa")

fig, ax = plt.subplots()
ax.loglog(divergence, -change)
ax.axhline(145, color="tab:red", ls="--", label="Richardson's 145 hPa")
ax.axvspan(1e-8, 1e-7, color="tab:green", alpha=0.2, label="realistic column mean")
ax.set(xlabel="column-mean convergence [1/s]", ylabel="pressure rise in 6 h [hPa]", title="Small divergence errors, huge tendencies")
ax.legend()

plt.show()
