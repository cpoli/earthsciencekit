"""
Hasselmann's stochastic climate (1976)
======================================

Klaus Hasselmann showed that slow climate variability needs no slow cause.
The ocean integrates fast, random weather forcing, just as a heavy particle
integrates molecular kicks in Brownian motion. The result is red noise: a
spectrum flat below the frequency ``lambda / (2 pi C)`` and falling as
frequency squared above it. Here a seeded white-noise forcing drives a 50 m
mixed layer, sampled every month for 4000 years.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

from earthsciencekit.climate import heat_capacity_column, red_noise_spectrum
from earthsciencekit.units import convert

restoring, capacity = 1.2, heat_capacity_column(1025, 3990, 50)
forcing_spectrum = 4e7  # W^2 m^-4 Hz^-1, two-sided
dt = convert(1 / 12, "year", "s")
decay = np.exp(-restoring * dt / capacity)
variance = forcing_spectrum / (2 * restoring * capacity)
rng = np.random.default_rng(1976)
kicks = rng.normal(0, np.sqrt(variance * (1 - decay**2)), 48000)
temperature = np.zeros(kicks.size)
for i in range(1, kicks.size):
    temperature[i] = decay * temperature[i - 1] + kicks[i]
assert np.isclose(temperature.var(), variance, rtol=0.15)
assert np.isclose(np.corrcoef(temperature[:-1], temperature[1:])[0, 1], decay, atol=0.01)

frequency, density = welch(temperature, fs=1 / dt, nperseg=4096)
per_year = convert(1, "year", "s")
fig, ax = plt.subplots()
ax.loglog(frequency[1:] * per_year, density[1:] / per_year, lw=0.8, label="simulated")
theory = 2 * red_noise_spectrum(frequency[1:], forcing_spectrum, restoring, capacity)
ax.loglog(frequency[1:] * per_year, theory / per_year, label="Hasselmann")
ax.axvline(restoring / (2 * np.pi * capacity) * per_year, color="0.5", ls=":", label="corner frequency")
ax.set(xlabel="frequency [cycles/year]", ylabel="one-sided spectrum [K² year]", title="Weather noise reddened by ocean heat capacity")
ax.legend()

plt.show()
