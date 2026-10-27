"""
Bullard's fit of the continents (1965)
======================================

Edward Bullard, Jim Everett and Alan Smith fitted the continents around the
Atlantic with a computer. They matched the 500-fathom contours of the
continental slopes, not the coastlines, and found the finite rotation about
one Euler pole that brought them together with the least misfit. South
America closes onto Africa by a rotation of about 57 degrees about a pole near
44 degN, 31 degW. Here a coarse West African margin, and a South American
margin made from it with synthetic irregularities, are closed by that
rotation. Interpolating it shows the Atlantic opening at a steady rate about
the same pole, and composing it with a synthetic motion of Africa gives South
America's rotation in that frame as a single rotation.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.constants import EARTH_RADIUS
from earthsciencekit.plate_tectonics import compose_finite_rotations, interpolate_finite_rotation, rotate_position

africa = np.radians(
    [(35, -6), (28, -13), (21, -17), (15, -17), (10, -14), (5, -7), (5, 2), (4, 9), (0, 9), (-6, 12), (-12, 13.5)]
    + [(-18, 11.8), (-23, 14.5), (-29, 16.5), (-34, 18.5)]
)
fit = np.radians([44.0, -30.6, 57.0])  # pole latitude, longitude and angle carrying South America onto Africa
rng = np.random.default_rng(1965)
opened = rotate_position(africa[:, 0], africa[:, 1], fit[0], fit[1], -fit[2])
america = opened + np.radians(rng.normal(0, 0.25, opened.shape))

closed = rotate_position(america[:, 0], america[:, 1], *fit)


def unit(position):
    return np.stack(
        [np.cos(position[:, 0]) * np.cos(position[:, 1]), np.cos(position[:, 0]) * np.sin(position[:, 1]), np.sin(position[:, 0])]
    )


gap = np.arccos(np.clip(np.sum(unit(closed) * unit(africa), axis=0), -1, 1))
misfit = EARTH_RADIUS * np.sqrt(np.mean(gap**2)) / 1e3
assert misfit < 60
print(f"r.m.s. misfit after closing the Atlantic: {misfit:.0f} km")

fractions = np.linspace(0, 1, 5)
stages = interpolate_finite_rotation(np.pi / 2, 0, 0, *fit, fractions)
assert np.allclose(stages[1:, :2], fit[:2]) and np.allclose(stages[:, 2], fractions * fit[2])

motion = np.radians([10.0, 20.0, 15.0])  # synthetic rotation of Africa in a fixed frame
total = compose_finite_rotations(*fit, *motion)
stepwise = rotate_position(*rotate_position(america[:, 0], america[:, 1], *fit).T, *motion)
assert np.allclose(rotate_position(america[:, 0], america[:, 1], *total), stepwise)
pole = f"{np.degrees(total[0]):.1f} degN, {np.degrees(total[1]):.1f} degE"
print(f"South America in the fixed frame: {np.degrees(total[2]):.1f} deg about {pole}, not {np.degrees(fit[2] + motion[2]):.0f} deg")

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
left.plot(np.degrees(africa[:, 1]), np.degrees(africa[:, 0]), "k-o", ms=3, label="Africa (fixed)")
left.plot(np.degrees(america[:, 1]), np.degrees(america[:, 0]), "-o", ms=3, color="tab:blue", label="South America today")
left.plot(np.degrees(closed[:, 1]), np.degrees(closed[:, 0]), "--o", ms=3, color="tab:red", label="South America, fitted")
left.plot(np.degrees(fit[1]), np.degrees(fit[0]), "k*", ms=12, label="Euler pole of the fit")
left.set(xlabel="longitude [deg]", ylabel="latitude [deg]", title="Closing the Atlantic")
left.legend(fontsize="small", loc="lower left")
right.plot(np.degrees(africa[:, 1]), np.degrees(africa[:, 0]), "k-", label="Africa")
for fraction, rotation in zip(fractions, stages):
    position = rotate_position(america[:, 0], america[:, 1], *rotation)
    right.plot(np.degrees(position[:, 1]), np.degrees(position[:, 0]), color=plt.cm.viridis(fraction), label=f"{1 - fraction:.0%} open")
right.set(xlabel="longitude [deg]", title="Closing at a steady rate about one pole")
right.legend(fontsize="small", loc="lower left")
fig.tight_layout()

plt.show()
