"""
Lighthill and Whitham's kinematic waves (1955)
==============================================

James Lighthill and Gerald Whitham showed that a flood in a long river is
governed by continuity alone when friction balances the slope. Discharge is
then a function of depth, and a flood crest travels at ``dQ/dA``, faster than
the water: 5/3 of the mean velocity for Manning's law in a wide channel. A
kinematic wave neither spreads nor shrinks, so the attenuation seen in rivers
comes from the pressure-gradient and inertia terms it neglects. Here a flood
is routed down a 50 km reach and compared with Muskingum routing using the
same travel time.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import kinematic_wave_celerity, manning_velocity, muskingum_routing, simulate_kinematic_wave

hour = 3600.0
time = np.linspace(0, 36 * hour, 361)
inflow = 10 + 90 * np.exp(-(((time - 8 * hour) / (3 * hour)) ** 2))
length, width, roughness, slope = 5e4, 50.0, 0.035, 5e-4
result = simulate_kinematic_wave(inflow, length, width, roughness, slope, time, nodes=400)
outflow, storage = result.values[:, 0], result.values[:, 1]

depth = (100 * roughness / (width * np.sqrt(slope))) ** 0.6  # depth of the peak flow in a wide channel
celerity = kinematic_wave_celerity(manning_velocity(roughness, depth, slope), 5 / 3)
travel = length / celerity
delay = time[np.argmax(outflow)] - time[np.argmax(inflow)]
volume = np.sum((inflow - outflow)[1:] + (inflow - outflow)[:-1]) / 2 * (time[1] - time[0])
throughput = np.sum(inflow[1:] + inflow[:-1]) / 2 * (time[1] - time[0])
assert abs(delay / travel - 1) < 0.2 and abs(storage[-1] - storage[0] - volume) < 1e-4 * throughput
print(f"crest celerity {celerity:.2f} m/s, travel time {travel / hour:.1f} h, simulated delay {delay / hour:.1f} h")

muskingum = muskingum_routing(inflow, travel, 0.2, time[1] - time[0])
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(time / hour, inflow, label="inflow")
ax.plot(time / hour, outflow, label="kinematic wave, 50 km downstream")
ax.plot(time / hour, muskingum, "--", label="Muskingum, K = L / c, x = 0.2")
ax.set(xlabel="time [h]", ylabel="discharge [m^3/s]", title="A flood crest travels at the kinematic celerity")
ax.legend()

plt.show()
