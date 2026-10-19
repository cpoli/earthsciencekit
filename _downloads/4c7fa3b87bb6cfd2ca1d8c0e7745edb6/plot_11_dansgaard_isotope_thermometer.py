"""
Dansgaard's isotope thermometer (1964)
======================================

Willi Dansgaard collected rain and snow from stations worldwide and found
that their oxygen-18 content falls by about 0.7 permil per kelvin of mean
annual temperature at high latitudes: as moist air cools, its heavy
molecules condense first. Read backwards, the isotope ratio of ice is a
thermometer, and the glacial ice of Greenland, about 7 permil lighter than
Holocene ice, records a cooling of about 10 K.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.glaciology import dansgaard_delta, dansgaard_temperature

slope, offset = 0.69, -13.6
station = np.array([258.0, 248.0, 263.0, 241.0, 268.0, 253.0, 233.0, 271.0])  # K, synthetic stations
scatter = np.array([0.6, -0.8, 0.3, 1.1, -0.4, -0.2, 0.7, -0.9])  # permil
measured = dansgaard_delta(station, slope, offset) + scatter
assert np.isclose(np.polyfit(station, measured, 1)[0], slope, atol=0.05)
holocene, glacial = dansgaard_temperature(np.array([-35.0, -42.0]), slope, offset)
assert np.isclose(holocene - glacial, 7 / slope)

fig, ax = plt.subplots()
line = np.linspace(225, 280, 2)
ax.plot(station - 273.15, measured, "o", label="synthetic stations")
ax.plot(line - 273.15, dansgaard_delta(line, slope, offset), label="0.69 permil/K")
title = f"Dansgaard: glacial ice {holocene - glacial:.0f} K colder"
ax.set(xlabel="mean annual temperature [degC]", ylabel="delta-18O [permil]", title=title)
ax.legend()

plt.show()
