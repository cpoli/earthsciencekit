"""Regenerate the README figures: python docs/make_readme_figure.py

Writes docs/source/_static/images/readme_hero.png, three showcase panels,
which the docs landing page and the Open Graph card embed, and
readme_<domain>.png, three panels per domain, for the landing page's
domain list.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.atmosphere import parcel_temperature_profile, rossby_wave_speed, saturation_vapor_pressure
from earthsciencekit.climate import daily_insolation, equilibrium_solar_constant, simulate_two_box
from earthsciencekit.geomorphology import (
    critical_shields_number,
    dimensionless_grain_size,
    scarp_elevation,
    simulate_river_profile,
)
from earthsciencekit.geophysics import block_magnetic_anomaly, field_elements, geomagnetic_field, simulate_rikitake_dynamo, sphere_gravity
from earthsciencekit.glaciology import glen_velocity_profile, herron_langway_density, vialov_dome_height, vialov_profile
from earthsciencekit.hydrology import muskingum_routing, nash_unit_hydrograph, theis_drawdown
from earthsciencekit.oceanography import ekman_spiral_velocity, stommel_streamfunction, wave_angular_frequency
from earthsciencekit.plate_tectonics import flexural_deflection, flexural_rigidity, halfspace_heat_flux, ocean_depth
from earthsciencekit.seismology import gutenberg_richter, head_wave_travel_time, simulate_oscillator, travel_time
from earthsciencekit.volcanology import (
    concordia_ratios,
    equilibrium_crystallization_melt,
    fractional_crystallization_melt,
    simulate_eruption_column,
)

IMAGES = Path(__file__).parent / "source" / "_static" / "images"
OUT = IMAGES / "readme_hero.png"
YEAR = 365.25 * 86400.0
BACKGROUND = "#f5f8fc"

fig = plt.figure(figsize=(16.5, 5.06), facecolor=BACKGROUND)
gyre = fig.add_axes((0.01, 0.03, 0.31, 0.84))
dynamo = fig.add_axes((0.31, -0.08, 0.37, 1.08), projection="3d")
field = fig.add_axes((0.67, 0.03, 0.32, 0.84))

# Stommel's wind-driven gyre: beta crowds the return flow into a western boundary current.
width, length = 4000e3, 4000e3
x, y = np.linspace(0, width, 401), np.linspace(0, length, 401)
east, north = np.meshgrid(x, y)
psi = stommel_streamfunction(east, north, 0.1, width, length, 1025.0, 2e-11, 3e-6)
v = np.gradient(psi, x, axis=1)
u = -np.gradient(psi, y, axis=0)
speed = np.hypot(u, v) / np.hypot(u, v).max()
gyre.pcolormesh(x / 1e3, y / 1e3, speed**0.5, cmap="viridis", shading="auto")
gyre.streamplot(x / 1e3, y / 1e3, u, v, color="w", density=1.5, linewidth=0.3 + 1.6 * speed**0.5, arrowsize=0.7)
gyre.set_xlim(0, width / 1e3)
gyre.set_ylim(0, length / 1e3)
fig.text(0.165, 0.93, "Stommel gyre: western intensification by $\\beta$", ha="center", fontsize=12)

# Rikitake's two-disc dynamo started twice, 0.001 apart: chaotic polarity reversals.
time = np.linspace(0, 160, 16001)
for start, cmap in (((1.0, 0.0, 0.0), "Blues"), ((1.001, 0.0, 0.0), "Oranges")):
    path = simulate_rikitake_dynamo(1.0, 2.0, 1.0, time, initial_state=start).values
    shades = plt.get_cmap(cmap)(np.linspace(0.35, 1.0, 12))
    for segment, shade in zip(np.array_split(np.arange(time.size), 12), shades):
        stop = min(segment[-1] + 2, time.size)
        dynamo.plot(*path[segment[0] : stop].T, color=shade, lw=0.6, alpha=0.85)
    dynamo.plot(*path[-1], "o", color=shades[-1], ms=6)
dynamo.view_init(elev=24, azim=-58)
dynamo.set_axis_off()
dynamo.set_facecolor(BACKGROUND)
fig.text(0.495, 0.93, "Rikitake dynamo: two starts 0.001 apart, chaotic reversals", ha="center", fontsize=12)

# IGRF-13 (2020) to degree 3: the South Atlantic Anomaly and the magnetic equator.
a = 6371.2e3
g = np.array([[0, 0, 0, 0], [-29404.8, -1450.9, 0, 0], [-2499.6, 2982.0, 1677.0, 0], [1363.2, -2381.2, 1236.2, 525.7]]) * 1e-9
h = np.array([[0, 0, 0, 0], [0, 4652.5, 0, 0], [0, -2991.6, -734.6, 0], [0, -82.1, 241.9, -543.4]]) * 1e-9
latitude, longitude = np.radians(np.linspace(-89, 89, 179)), np.radians(np.linspace(-180, 180, 361))
lat, lon = np.meshgrid(latitude, longitude, indexing="ij")
_, inclination, intensity = np.moveaxis(field_elements(*np.moveaxis(geomagnetic_field(lat, lon, a, g, h, a), -1, 0)), -1, 0)
field.contourf(np.degrees(lon), np.degrees(lat), intensity * 1e9, levels=24, cmap="magma")
field.contour(np.degrees(lon), np.degrees(lat), np.degrees(inclination), levels=np.arange(-80, 81, 20), colors="w", linewidths=0.5)
field.contour(np.degrees(lon), np.degrees(lat), np.degrees(inclination), levels=[0], colors="w", linewidths=2)
weakest = np.unravel_index(np.argmin(intensity), intensity.shape)
field.plot(np.degrees(lon[weakest]), np.degrees(lat[weakest]), "*", color="#7fd4ff", ms=16, mec="k")
fig.text(0.83, 0.93, "IGRF-13 field intensity: the South Atlantic Anomaly", ha="center", fontsize=12)

for ax in (gyre, field):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)

fig.savefig(OUT, dpi=100, facecolor=BACKGROUND)
print(f"wrote {OUT}")


def domain_figure(name):
    """Three panels for one domain, saved as readme_<name>.png."""
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.06), constrained_layout=True, facecolor=BACKGROUND)
    yield axes
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    out = IMAGES / f"readme_{name}.png"
    fig.savefig(out, dpi=100, facecolor=BACKGROUND)
    plt.close(fig)
    print(f"wrote {out}")


MYR = 1e6 * YEAR

for left, middle, right in domain_figure("plate_tectonics"):
    # Sea-floor depth and heat flow from half-space cooling, and a plate bent by a seamount chain.
    age = np.linspace(0, 120, 241) * MYR
    left.plot(age / MYR, ocean_depth(age, 2500.0, 3300.0, 1000.0, 3e-5, 1e-6, 273.0, 1600.0) / 1e3)
    left.invert_yaxis()
    left.set(title="Sea-floor subsidence", xlabel="age [Myr]", ylabel="depth [km]")
    middle.plot(age[1:] / MYR, halfspace_heat_flux(age[1:], 1e-6, 3.3, 273.0, 1600.0) * 1e3)
    middle.set_yscale("log")
    middle.set(title="Oceanic heat flow", xlabel="age [Myr]", ylabel="heat flux [mW/m$^2$]")
    distance = np.linspace(-400e3, 400e3, 801)
    for thickness in (10e3, 25e3, 40e3):
        rigidity = flexural_rigidity(70e9, thickness, 0.25)
        deflection = flexural_deflection(distance, 2e12, rigidity, 2300.0, 9.80665)
        right.plot(distance / 1e3, -deflection / 1e3, label=f"$T_e$ = {thickness / 1e3:.0f} km")
    right.legend(frameon=False)
    right.set(title="Lithospheric flexure", xlabel="distance from load [km]", ylabel="deflection [km]")

for left, middle, right in domain_figure("seismology"):
    # Refraction travel times with the crossover, Gutenberg-Richter, and a damped seismometer.
    distance = np.linspace(0, 300e3, 601)
    critical = 2 * 30e3 * np.tan(np.arcsin(6000.0 / 8000.0))
    beyond = distance[distance >= critical]
    left.plot(distance / 1e3, travel_time(distance, 6000.0), label="direct P")
    left.plot(beyond / 1e3, head_wave_travel_time(beyond, 30e3, 6000.0, 8000.0), label="head wave (Pn)")
    left.legend(frameon=False)
    left.set(title="Refraction travel times", xlabel="distance [km]", ylabel="travel time [s]")
    magnitude = np.linspace(2, 8, 121)
    middle.semilogy(magnitude, gutenberg_richter(8.0, 1.0, magnitude))
    middle.set(title="Gutenberg-Richter law", xlabel="magnitude", ylabel="events with M or larger")
    time = np.linspace(0, 10, 2001)
    for quality in (0.7, 3.0, 20.0):
        right.plot(time, simulate_oscillator(1.0, quality, time).values[:, 0], label=f"Q = {quality:g}")
    right.legend(frameon=False)
    right.set(title="Damped oscillator", xlabel="time [s]", ylabel="displacement [m]")

for left, middle, right in domain_figure("atmosphere"):
    # Clausius-Clapeyron, a lifted parcel's dry and moist adiabats, and Rossby wave speeds.
    temperature = np.linspace(233.15, 313.15, 161)
    left.plot(temperature - 273.15, saturation_vapor_pressure(temperature) / 100)
    left.set(title="Saturation vapour pressure", xlabel="temperature [°C]", ylabel="pressure [hPa]")
    pressure = np.linspace(100000.0, 20000.0, 161)
    middle.plot(parcel_temperature_profile(pressure, 300.0, 300.0) - 273.15, pressure / 100, label="saturated start")
    middle.plot(parcel_temperature_profile(pressure, 300.0, 285.0) - 273.15, pressure / 100, label="dew point 12 °C")
    middle.plot(parcel_temperature_profile(pressure, 300.0, 240.0) - 273.15, pressure / 100, label="dew point -33 °C")
    middle.set_yscale("log")
    middle.invert_yaxis()
    middle.legend(frameon=False)
    middle.set(title="Lifted parcel", xlabel="temperature [°C]", ylabel="pressure [hPa]")
    wavelength = np.linspace(1000e3, 12000e3, 221)
    for wind in (5.0, 15.0, 25.0):
        speed = rossby_wave_speed(wind, 1.6e-11, 2 * np.pi / wavelength)
        right.plot(wavelength / 1e3, speed, label=f"U = {wind:g} m/s")
    right.axhline(0, color="k", lw=0.5)
    right.legend(frameon=False)
    right.set(title="Rossby wave speed", xlabel="wavelength [km]", ylabel="phase speed [m/s]")

for left, middle, right in domain_figure("climate"):
    # Ice-albedo bistability, two-box warming and daily insolation.
    temperature = np.linspace(220, 320, 401)
    left.plot(equilibrium_solar_constant(temperature, 0.61, 0.62, 0.3, 263.15, 8.0), temperature)
    left.axvline(1361, color="k", lw=0.5, ls="--")
    left.set(title="Snowball bistability", xlabel="solar constant [W/m$^2$]", ylabel="equilibrium temperature [K]")
    time = np.linspace(0, 500, 501) * YEAR
    surface, deep = simulate_two_box(3.7, 1.3, 0.7, 2.4e8, 3.0e9, time).values.T
    middle.plot(time / YEAR, surface, label="surface")
    middle.plot(time / YEAR, deep, label="deep ocean")
    middle.legend(frameon=False)
    middle.set(title="Two-box warming", xlabel="time [year]", ylabel="anomaly [K]")
    latitude = np.linspace(-np.pi / 2, np.pi / 2, 361)
    for declination, label in ((23.44, "June solstice"), (0.0, "equinox"), (-23.44, "December solstice")):
        right.plot(np.degrees(latitude), daily_insolation(1361.0, latitude, np.radians(declination), 1.0), label=label)
    right.legend(frameon=False)
    right.set(title="Daily insolation", xlabel="latitude [°]", ylabel="insolation [W/m$^2$]")

for left, middle, right in domain_figure("oceanography"):
    # The Ekman spiral, Stommel's western-intensified gyre and surface-wave dispersion.
    depth = np.linspace(0, 120, 121)
    downwind, crosswind = ekman_spiral_velocity(depth, 0.1, 1025.0, 0.01, 1e-4).T
    left.plot(crosswind, downwind)
    left.plot(crosswind[::10], downwind[::10], "o", ms=3)
    left.set_aspect("equal")
    left.set(title="Ekman spiral", xlabel="crosswind [m/s]", ylabel="downwind [m/s]")
    x, y = np.meshgrid(np.linspace(0, 6000e3, 121), np.linspace(0, 4000e3, 81))
    psi = stommel_streamfunction(x, y, 0.1, 6000e3, 4000e3, 1025.0, 2e-11, 1e-6)
    middle.contourf(x / 1e3, y / 1e3, psi / 1e6, levels=14, cmap="Blues")
    middle.contour(x / 1e3, y / 1e3, psi / 1e6, levels=14, colors="k", linewidths=0.4)
    middle.set_aspect("equal")
    middle.set(title="Stommel gyre", xlabel="east [km]", ylabel="north [km]")
    wavelength = np.logspace(0, 4, 241)
    wavenumber = 2 * np.pi / wavelength
    for water_depth in (5.0, 50.0, 4000.0):
        speed = wave_angular_frequency(wavenumber, water_depth, 9.80665) / wavenumber
        right.semilogx(wavelength, speed, label=f"depth {water_depth:g} m")
    right.legend(frameon=False)
    right.set(title="Wave dispersion", xlabel="wavelength [m]", ylabel="phase speed [m/s]")

for left, middle, right in domain_figure("volcanology"):
    # U-Pb concordia, trace elements in a crystallizing melt, and an eruption column.
    age = np.linspace(0, 4.5e3, 451) * MYR
    ratios = concordia_ratios(age)
    marks = concordia_ratios(np.arange(500, 4501, 500) * MYR)
    left.plot(ratios[:, 0], ratios[:, 1])
    left.plot(marks[:, 0], marks[:, 1], "o", ms=4)
    left.set(title="U-Pb concordia", xlabel="$^{207}$Pb*/$^{235}$U", ylabel="$^{206}$Pb*/$^{238}$U")
    melt = np.linspace(1, 0.05, 191)
    for coefficient in (0.01, 0.5, 3.0):
        line = middle.plot(melt, fractional_crystallization_melt(1.0, melt, coefficient), label=f"D = {coefficient:g}")
        middle.plot(melt, equilibrium_crystallization_melt(1.0, melt, coefficient), ls="--", color=line[0].get_color())
    middle.set_yscale("log")
    middle.invert_xaxis()
    middle.legend(frameon=False, title="solid: fractional, dashed: equilibrium")
    middle.set(title="Crystallizing melt", xlabel="melt fraction", ylabel="concentration / initial")
    time = np.linspace(0, 1200, 1201)
    height = simulate_eruption_column(50.0, 150.0, 5.0, 0.01, time).values[:, 0]
    right.plot(time / 60, height / 1e3)
    right.set(title="Eruption column", xlabel="time [min]", ylabel="parcel height [km]")

for left, middle, right in domain_figure("geophysics"):
    # A buried sphere's gravity, a magnetized block's anomaly and the Rikitake attractor.
    x = np.linspace(-10e3, 10e3, 401)
    for depth in (1e3, 2e3, 4e3):
        left.plot(x / 1e3, sphere_gravity(x, depth, 500.0, 500.0)[:, 1] * 1e5, label=f"depth {depth / 1e3:g} km")
    left.legend(frameon=False)
    left.set(title="Buried sphere", xlabel="distance [km]", ylabel="gravity anomaly [mGal]")
    x = np.linspace(-20e3, 20e3, 801)
    middle.plot(x / 1e3, block_magnetic_anomaly(x, -5e3, 5e3, 1e3, 3e3, 2.0) * 1e9)
    middle.axvspan(-5, 5, color="0.85", zorder=0)
    middle.set(title="Magnetized block", xlabel="distance [km]", ylabel="vertical anomaly [nT]")
    time = np.linspace(0, 300, 12001)
    currents = simulate_rikitake_dynamo(1.0, 2.0, 1.0, time).values
    right.plot(currents[:, 0], currents[:, 1], lw=0.4)
    right.set(title="Rikitake attractor", xlabel="current x", ylabel="current y")

for left, middle, right in domain_figure("hydrology"):
    # Theis drawdown, Nash unit hydrographs and Muskingum flood routing.
    time = np.logspace(1, 6, 251)
    for radius in (10.0, 100.0, 1000.0):
        left.semilogx(time / 3600, theis_drawdown(0.01, 1e-3, 1e-4, radius, time), label=f"r = {radius:g} m")
    left.invert_yaxis()
    left.legend(frameon=False)
    left.set(title="Theis drawdown", xlabel="time [h]", ylabel="drawdown [m]")
    time = np.linspace(0, 48, 481) * 3600
    for reservoirs in (1, 2, 4):
        middle.plot(time / 3600, nash_unit_hydrograph(time, reservoirs, 4 * 3600) * 3600, label=f"n = {reservoirs}")
    middle.legend(frameon=False)
    middle.set(title="Nash unit hydrographs", xlabel="time [h]", ylabel="ordinate [1/h]")
    hours = np.arange(0, 72, 1.0)
    inflow = 20 + 180 * np.exp(-(((hours - 18) / 6) ** 2))
    right.plot(hours, inflow, label="inflow")
    right.plot(hours, muskingum_routing(inflow, 12 * 3600, 0.2, 3600.0), label="outflow")
    right.legend(frameon=False)
    right.set(title="Muskingum routing", xlabel="time [h]", ylabel="discharge [m$^3$/s]")

for left, middle, right in domain_figure("glaciology"):
    # Glen's-law velocity profiles, the Vialov ice-sheet dome and firn densification.
    height = np.linspace(0, 500, 201)
    for exponent in (1.0, 3.0, 5.0):
        rate_factor = 2.4e-24 * (1e5 ** (3 - exponent))
        velocity = glen_velocity_profile(height, 500.0, 0.03, rate_factor, exponent, 917.0, 9.80665)
        left.plot(velocity / velocity[-1], height, label=f"n = {exponent:g}")
    left.legend(frameon=False)
    left.set(title="Glen's flow law", xlabel="velocity / surface velocity", ylabel="height above bed [m]")
    distance = np.linspace(0, 1000e3, 501)
    for accumulation in (0.05, 0.2, 0.5):
        dome = vialov_dome_height(1000e3, accumulation / YEAR, 2.4e-24, 3.0, 917.0, 9.80665)
        middle.plot(distance / 1e3, vialov_profile(distance, 1000e3, dome, 3.0) / 1e3, label=f"{accumulation:g} m/yr")
    middle.legend(frameon=False, title="accumulation")
    middle.set(title="Vialov ice sheet", xlabel="distance from divide [km]", ylabel="thickness [km]")
    depth = np.linspace(0, 120, 241)
    for temperature, accumulation in ((248.15, 0.25), (223.15, 0.03)):
        density = herron_langway_density(depth, temperature, accumulation / YEAR, 350.0)
        right.plot(density, depth, label=f"{temperature - 273.15:.0f} °C, {accumulation:g} m/yr")
    right.invert_yaxis()
    right.legend(frameon=False)
    right.set(title="Firn densification", xlabel="density [kg/m$^3$]", ylabel="depth [m]")

for left, middle, right in domain_figure("geomorphology"):
    # The Shields curve, a diffusing fault scarp and a river profile rising to steady state.
    diameter = np.logspace(-5, -1, 241)
    grain = dimensionless_grain_size(diameter, 2650.0, 1000.0, 1e-6, 9.80665)
    left.loglog(grain, critical_shields_number(grain))
    left.set(title="Shields curve", xlabel="dimensionless grain size $D_*$", ylabel="critical Shields number")
    distance = np.linspace(-30, 30, 241)
    for age in (0.0, 1e3, 10e3, 50e3):
        middle.plot(distance, scarp_elevation(distance, 2.0, 1e-3 / YEAR, age * YEAR, 0.05), label=f"{age / 1e3:g} kyr")
    middle.legend(frameon=False)
    middle.set(title="Degrading fault scarp", xlabel="distance [m]", ylabel="elevation [m]")
    distance = np.linspace(0, 50e3, 51)
    area = 1e6 + 6.7 * (distance[-1] + 1e3 - distance) ** 1.67
    output = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 4.0]) * MYR
    profile = simulate_river_profile(distance, area, 1e-3 / YEAR, 2e-6 / YEAR, output)
    for age, row in zip(output[1:], profile.values[1:]):
        right.plot(distance / 1e3, row[2:], label=f"{age / MYR:g} Myr")
    right.legend(frameon=False)
    right.set(title="River profile", xlabel="distance upstream [km]", ylabel="elevation [m]")
