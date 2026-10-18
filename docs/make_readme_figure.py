"""Regenerate the README hero figure: python docs/make_readme_figure.py

Writes docs/source/_static/images/readme_hero.png, one panel per domain,
which the docs landing page and the Open Graph card embed.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import barometric_pressure
from earthsciencekit.climate import simulate_one_box
from earthsciencekit.geophysics import simulate_rikitake_dynamo
from earthsciencekit.glaciology import simulate_flowline_glacier
from earthsciencekit.hydrology import simulate_infiltration
from earthsciencekit.oceanography import simulate_inertial_current
from earthsciencekit.plate_tectonics import halfspace_temperature
from earthsciencekit.seismology import ricker_wavelet
from earthsciencekit.volcanology import mogi_displacement

OUT = Path(__file__).parent / "source" / "_static" / "images" / "readme_hero.png"
YEAR = 365.25 * 86400.0
BACKGROUND = "#f5f8fc"

fig, axes = plt.subplots(2, 5, figsize=(22, 8.4), constrained_layout=True, facecolor=BACKGROUND)
plate, seismic, column, dynamo, climate, current, dome, soil, glacier, label = axes.ravel()

# Half-space cooling of 50 Myr old oceanic lithosphere.
depth = np.linspace(0, 100e3, 201)
plate.plot(halfspace_temperature(depth, 50e6 * YEAR, 1e-6, 273.0, 1600.0), depth / 1e3)
plate.invert_yaxis()
plate.set(title="Plate cooling", xlabel="temperature [K]", ylabel="depth [km]")

# A 10 Hz Ricker source pulse.
time = np.linspace(-0.2, 0.2, 801)
seismic.plot(time, ricker_wavelet(time, 10.0))
seismic.set(title="Seismic source", xlabel="time [s]", ylabel="normalized pulse")

# An isothermal hydrostatic column.
height = np.linspace(0, 15e3, 301)
column.plot(barometric_pressure(101325.0, height, 255.0, 9.80665) / 100, height / 1e3)
column.set(title="Atmospheric column", xlabel="pressure [hPa]", ylabel="height [km]")

# Rikitake's two-disc dynamo reversing its polarity.
time = np.linspace(0, 200, 4001)
dynamo.plot(time, simulate_rikitake_dynamo(1.0, 2.0, 1.0, time).values[:, 0], lw=0.8)
dynamo.axhline(0, color="k", lw=0.5)
dynamo.set(title="Dynamo reversals", xlabel="nondimensional time", ylabel="current x (sign = polarity)")

# A one-box energy balance adjusting to doubled CO2.
time = np.linspace(0, 32 * YEAR, 401)
warming = simulate_one_box(3.7, 1.1, 1.4e8, time)
climate.plot(time / YEAR, warming.values[:, 0])
climate.set(title="Climate adjustment", xlabel="time [year]", ylabel="anomaly [K]")

# A damped inertial current turning clockwise.
time = np.linspace(0, 1.6 * 2 * np.pi / 1e-4, 2001)
east, north = simulate_inertial_current(1e-4, time, damping=1.2e-5).values.T
current.plot(east, north)
current.set_aspect("equal")
current.set(title="Inertial current", xlabel="east [m/s]", ylabel="north [m/s]")

# Mogi uplift above a magma chamber inflating by 1e6 m^3 at 3 km depth.
distance = np.linspace(0, 12e3, 601)
radial, vertical = mogi_displacement(distance, 3000.0, 1e6, 0.25).T
dome.plot(distance / 1e3, vertical * 100, label="uplift")
dome.plot(distance / 1e3, radial * 100, label="radial")
dome.legend(frameon=False)
dome.set(title="Volcano inflation", xlabel="distance from summit [km]", ylabel="displacement [cm]")

# Green-Ampt infiltration of 36 mm/h rain on a loam: the soil ponds and the excess runs off.
time = np.linspace(0, 4 * 3600, 481)
infiltration, runoff = simulate_infiltration(1e-6, 0.11, 0.3, 1e-5, time).values[:, 1:].T
soil.plot(time / 3600, infiltration * 3.6e6, label="infiltration")
soil.plot(time / 3600, runoff * 3.6e6, label="runoff")
soil.legend(frameon=False)
soil.set(title="Storm infiltration", xlabel="time [h]", ylabel="rate [mm/h]")

# A shallow-ice valley glacier grown to steady state on a sloping bed.
x = np.arange(0, 20e3 + 1, 200.0)
bed = 3000 - 0.1 * x
grown = simulate_flowline_glacier(bed, 200.0, 2600.0, 0.007 / YEAR, np.linspace(0, 1500, 76) * YEAR)
glacier.plot(x / 1e3, bed, color="saddlebrown", label="bed")
glacier.plot(x / 1e3, bed + grown.values[-1, 2:], label="ice surface")
glacier.legend(frameon=False)
glacier.set(title="Flowline glacier", xlabel="distance along flow [km]", ylabel="elevation [m]")

for ax in axes.ravel()[:-1]:
    ax.spines[["top", "right"]].set_visible(False)

label.axis("off")
label.text(0.05, 0.62, "EarthScienceKit", fontsize=28, fontweight="bold", transform=label.transAxes)
label.text(
    0.05,
    0.52,
    "Nine domains\nEquations • histories • examples\nInspectable computations",
    fontsize=16,
    linespacing=1.6,
    va="top",
    transform=label.transAxes,
)

fig.savefig(OUT, dpi=100, facecolor=BACKGROUND)
print(f"wrote {OUT}")
