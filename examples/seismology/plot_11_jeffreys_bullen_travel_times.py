"""
Jeffreys and Bullen's travel-time tables (1940)
===============================================

Harold Jeffreys and Keith Bullen fitted the arrival times of thousands of
earthquakes to a spherically symmetric Earth and published travel times of
P, S and the core phases out to 180 degrees. Their tables located
earthquakes for half a century. Each ray keeps its ray parameter
``r sin(i) / v``, so tracing rays through a layered sphere gives the
travel-time curve, whose slope is the ray parameter. The rounded PREM-like
P speeds below give direct P out to about 98 degrees, where rays graze the
core, and core phases beyond about 113 degrees: the shadow zone between.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import layered_ray_distance, layered_ray_travel_time

# depth [km], P speed [km/s]; repeated depths are discontinuities (Moho, 220, 410, 660 km, core, inner core).
table = np.array(
    [
        (0, 5.8), (15, 5.8), (15, 6.8), (24.4, 6.8), (24.4, 8.11), (220, 7.99), (220, 8.56), (400, 8.91),
        (400, 9.13), (600, 10.16), (670, 10.27), (670, 10.75), (771, 11.07), (971, 11.42), (1171, 11.73),
        (1371, 12.02), (1571, 12.29), (1771, 12.54), (1971, 12.78), (2171, 13.02), (2371, 13.25), (2571, 13.48),
        (2741, 13.68), (2891, 13.72), (2891, 8.06), (3371, 8.65), (3871, 9.2), (4371, 9.65), (4871, 10.0),
        (5150, 10.36), (5150, 11.03), (6370, 11.26),
    ]
)  # fmt: skip
depth, speed, radius = table[:, 0] * 1e3, table[:, 1] * 1e3, 6371e3
core, inner = 3480e3 / 13.72e3, 1221e3 / 10.36e3  # ray parameters that graze the core and the inner core
p = np.geomspace(0.5, 0.9999 * radius / speed[0], 20001)
distance = np.degrees(layered_ray_distance(p, depth, speed, radius))
time = layered_ray_travel_time(p, depth, speed, radius) / 60
distance = np.where(distance > 180, 360 - distance, distance)
phases = {"P": p > core, "PKP": (p < core) & (p > inner), "PKIKP": p < inner}

edge = distance[phases["P"]][p[phases["P"]].argmin()]
core_start = distance[~phases["P"]].min()
p60 = np.interp(60, distance[phases["P"]][::-1], time[phases["P"]][::-1])
assert 95 < edge < 101 and 105 < core_start < 120 and abs(p60 - 10.0) < 0.3
print(f"P to {edge:.1f} deg, core phases from {core_start:.1f} deg; P at 60 deg after {p60:.2f} min")

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
for name, rays in phases.items():
    left.plot(distance[rays], time[rays], ".", ms=1.5, label=name)
    right.plot(p[rays], distance[rays], ".", ms=1.5, label=name)
left.axvspan(edge, core_start, color="0.85", label="shadow zone")
right.axhspan(edge, core_start, color="0.85")
left.set(xlabel="epicentral distance [deg]", ylabel="travel time [min]", title="Travel-time curves")
left.legend(markerscale=8)
right.set(xlabel="ray parameter [s/rad]", ylabel="epicentral distance [deg]", title="Distance against ray parameter")
fig.tight_layout()

plt.show()
