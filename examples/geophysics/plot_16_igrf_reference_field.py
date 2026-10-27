"""
The International Geomagnetic Reference Field (1968)
====================================================

In 1968 the International Association of Geomagnetism and Aeronomy adopted
one set of Gauss coefficients, for epoch 1965.0, as the reference model of
the main field. Revised every five years since, the IGRF gives the field
everywhere from a short table of ``g`` and ``h`` coefficients. Truncated at
degree 3, the IGRF-13 field for 2020 already shows what a dipole cannot: the
field is weakest over South America and the South Atlantic, and the compass
points up to 20 degrees or more away from geographic north at mid-latitudes.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geophysics import dipole_field, field_elements, gauss_dipole_moment, geomagnetic_field, zonal_field
from earthsciencekit.units import convert

a = 6371.2e3
# IGRF-13 main-field coefficients for 2020.0 [nT], degrees 1 to 3: g[n, m] and h[n, m].
g = np.array([[0, 0, 0, 0], [-29404.8, -1450.9, 0, 0], [-2499.6, 2982.0, 1677.0, 0], [1363.2, -2381.2, 1236.2, 525.7]])
h = np.array([[0, 0, 0, 0], [0, 4652.5, 0, 0], [0, -2991.6, -734.6, 0], [0, -82.1, 241.9, -543.4]])
g, h = convert(g, "nT", "T"), convert(h, "nT", "T")

latitude = np.radians(np.linspace(-85, 85, 69))
longitude = np.radians(np.linspace(-180, 180, 145))
lat, lon = np.meshgrid(latitude, longitude, indexing="ij")
field = geomagnetic_field(lat, lon, a, g, h, a)
declination, inclination, intensity = np.moveaxis(field_elements(*np.moveaxis(field, -1, 0)), -1, 0)

zonal_only = g * (np.arange(4) == 0)  # keep the m = 0 column
dipole_only = np.where((np.arange(4)[:, None] == 1) & (np.arange(4) == 0), g, 0.0)
axial = geomagnetic_field(latitude, 0.0, a, zonal_only, np.zeros_like(h), a)
np.testing.assert_allclose(axial[:, [0, 2]], sum(zonal_field(latitude, a, n, g[n, 0], a) for n in range(1, 4)), rtol=1e-12)
dipole = dipole_field(latitude, a, gauss_dipole_moment(g[1, 0], 0, 0, a))
np.testing.assert_allclose(geomagnetic_field(latitude, 1.0, a, dipole_only, np.zeros_like(h), a)[:, [0, 2]], dipole)
mid = np.abs(lat) < np.radians(60)
assert np.degrees(np.abs(declination[mid])).max() > 20

weakest = np.unravel_index(np.argmin(intensity), intensity.shape)
low_lat, low_lon = np.degrees(lat[weakest]), np.degrees(lon[weakest])
assert -60 < low_lat < 0 and -90 < low_lon < 0
print(f"weakest field {convert(intensity.min(), 'T', 'nT'):.0f} nT at {low_lat:.0f} degN, {low_lon:.0f} degE")

fig, (top, bottom) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
image = top.contourf(np.degrees(lon), np.degrees(lat), convert(intensity, "T", "nT"), levels=20, cmap="viridis")
fig.colorbar(image, ax=top, label="intensity [nT]")
top.plot(low_lon, low_lat, "r*", ms=12)
top.set(ylabel="latitude [deg]", title="IGRF-13 (2020) to degree 3: total intensity")
image = bottom.contourf(
    np.degrees(lon), np.degrees(lat), np.degrees(declination), levels=np.arange(-60, 61, 10), cmap="RdBu_r", extend="both"
)
fig.colorbar(image, ax=bottom, label="declination [deg east]")
bottom.contour(np.degrees(lon), np.degrees(lat), np.degrees(inclination), levels=[0], colors="k", linewidths=1.5)
bottom.set(xlabel="longitude [deg]", ylabel="latitude [deg]", title="Declination, with the magnetic equator (inclination 0)")
fig.tight_layout()

plt.show()
