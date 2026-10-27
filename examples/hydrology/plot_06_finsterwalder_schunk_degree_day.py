"""
Finsterwalder and Schunk's degree-day melt (1887)
=================================================

Sebastian Finsterwalder and Hans Schunk, surveying the Suldenferner glacier
in the Ortler Alps, related its melt to the sum of positive air temperatures
rather than to the full energy balance, which nobody could then measure. The
degree-day factor, a few millimetres of water per degree per day, has served
snow hydrology ever since. A synthetic mountain year shows the snowpack
storing winter precipitation and releasing it as a spring melt pulse whose
volume follows the accumulated positive degree-days.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.hydrology import degree_day_melt, snowpack_water_equivalent

day = 86400.0
days = np.arange(365)  # from 1 October
rng = np.random.default_rng(1887)
temperature = 273.15 + 2 - 9 * np.cos(2 * np.pi * (days - 105) / 365) + rng.normal(0, 3, days.size)  # coldest in mid-January
precipitation = np.where(rng.random(days.size) < 0.35, rng.exponential(0.008, days.size), 0.0)  # m per day
factor = 3e-3 / day  # 3 mm per degC per day
pack = snowpack_water_equivalent(temperature, precipitation, factor, 273.15, day)
swe, released = pack.T
assert np.isclose(swe[-1] + released.sum(), precipitation.sum())

degree_days = np.cumsum(np.maximum(temperature - 273.15, 0))
melt = np.cumsum(released - precipitation * (temperature > 273.15))  # released water minus rain
deep = (temperature > 273.15) & (np.r_[0, swe[:-1]] > 0.02)
assert np.allclose(np.diff(np.r_[0, melt])[deep], degree_day_melt(temperature[deep], factor, 273.15) * day)
peak = days[np.argmax(swe)]
gone = days[peak:][np.argmax(swe[peak:] == 0)]
assert 100 < peak < 240 and gone > peak
print(f"peak snow water equivalent {swe.max() * 1e3:.0f} mm on day {peak}; snow gone on day {gone}")

fig, (top, bottom) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
top.plot(days, temperature - 273.15, color="0.5", lw=0.8, label="air temperature [degC]")
top.axhline(0, color="k", lw=0.5)
twin = top.twinx()
twin.plot(days, swe * 1e3, color="tab:blue", label="snow water equivalent [mm]")
top.set(ylabel="air temperature [degC]", title="A degree-day snowpack")
twin.set_ylabel("snow water equivalent [mm]", color="tab:blue")
bottom.plot(days, degree_days * 3, label="3 mm x positive degree-days")
bottom.plot(days, melt * 1e3, label="cumulative melt [mm]")
bottom.set(xlabel="days since 1 October", ylabel="mm of water", title="Melt follows the positive degree-days while snow lasts")
bottom.legend()
fig.tight_layout()

plt.show()
