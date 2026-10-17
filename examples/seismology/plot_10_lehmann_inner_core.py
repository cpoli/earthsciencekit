"""
Lehmann's inner core (1936)
===========================

Inge Lehmann noticed weak P arrivals inside the core shadow zone. She
explained them with an inner core faster than the outer core: rays meeting it
beyond the critical angle are totally reflected and emerge in the shadow. In
this straight-ray toy Earth (mantle 13 km/s, outer core 8 km/s, inner core
11 km/s, inner-core radius 1220 km) the post-critical reflections fill the
shadow zone of the 1914 example.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import critical_angle, ray_parameter

R, rc, ri = 6371e3, 3480e3, 1220e3
vm, vc, vi = 13e3, 8e3, 11e3
shadow = (np.degrees(2 * np.arccos(rc / R)), 143.1)
p = ray_parameter(R, vm, np.linspace(1e-4, np.arcsin(rc / R) - 1e-6, 4000))
b_mantle, b_core = p * vm, p * vc
hits = b_core < ri
b_mantle, b_core = b_mantle[hits], b_core[hits]
distance = 2 * (np.arccos(b_mantle / R) - np.arccos(b_mantle / rc)) + 2 * (np.arccos(b_core / rc) - np.arccos(b_core / ri))
incidence = np.arcsin(b_core / ri)
total = incidence >= critical_angle(vc, vi)
lit = np.degrees(distance[total])
assert lit.min() < shadow[0] and lit.max() > shadow[1]
print(f"critical angle {np.degrees(critical_angle(vc, vi)):.1f} deg; total reflections from {lit.min():.0f} to {lit.max():.0f} deg")

fig, ax = plt.subplots()
ax.plot(np.degrees(incidence), np.degrees(distance))
ax.axvline(np.degrees(critical_angle(vc, vi)), color="k", ls=":", label="critical angle")
ax.axhspan(*shadow, color="0.85", label="shadow zone")
ax.set(xlabel="incidence on the inner core [deg]", ylabel="epicentral distance [deg]", title="Inner-core reflections reach the shadow")
ax.legend()

plt.show()
