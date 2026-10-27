"""Snow: degree-day melt and a degree-day snowpack.

SI units unless stated otherwise. Snow, rain and melt are water equivalent.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "degree_day_melt",
    "snowpack_water_equivalent",
]


@model
def degree_day_melt(temperature, degree_day_factor, threshold_temperature):
    """Snowmelt rate from the air temperature by the degree-day method.

    Air temperature alone, the one weather variable measured almost everywhere, predicts melt well over days because it
    tracks most of the energy that reaches the snow.

    Parameters
    ----------
    temperature : float or array_like
        K, mean air temperature over the period; constraint: positive.
    degree_day_factor : float or array_like
        m/(K s) of water equivalent, about 3e-3 / 86400 (3 mm per degC per day) for snow; constraint: nonnegative.
    threshold_temperature : float or array_like
        K, usually 273.15; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s] of water equivalent. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Finsterwalder and Schunk (1887): melt in proportion to the positive degree-days, for periods of a day or longer
    with snow present; the factor varies with season, site and snow age.
    Equation: ``degree_day_factor * np.maximum(temperature - threshold_temperature, 0)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import degree_day_melt
    >>> float(degree_day_melt(278.15, 2, 273.15))
    10.0
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    degree_day_factor = quantity("degree_day_factor", degree_day_factor, lower=0)
    threshold_temperature = quantity("threshold_temperature", threshold_temperature, lower=0, strict_lower=True)
    return degree_day_factor * np.maximum(temperature - threshold_temperature, 0)


@model
def snowpack_water_equivalent(temperature, precipitation, degree_day_factor, threshold_temperature, time_step):
    """Snow water equivalent and water released by a degree-day snowpack.

    Snow stores winter precipitation and releases it as melt in spring, the delay that shapes the hydrograph of every
    mountain river.

    Parameters
    ----------
    temperature : array_like, shape (n,)
        K, mean air temperature of each time step; constraint: positive vector.
    precipitation : array_like, shape (n,)
        m of water equivalent in each time step; constraint: nonnegative vector.
    degree_day_factor : float
        m/(K s); constraint: nonnegative scalar.
    threshold_temperature : float
        K, below or at which precipitation is snow and above which snow melts; constraint: positive scalar.
    time_step : float
        s; constraint: positive scalar.

    Returns
    -------
    ndarray, shape (n, 2)
        Model value [m]: the snow water equivalent at the end of each step and the water released in it (rain plus
        melt) on the last axis.

    Notes
    -----
    A single-layer degree-day snowpack starting bare: snowfall accumulates, melt is :func:`degree_day_melt` times the
    step and cannot exceed the snow present, and rain passes straight through. There is no refreezing, liquid storage
    or sublimation.
    Equation: ``melt = min(swe + snowfall, degree_day_factor * max(T - T0, 0) * time_step)``,
    ``swe = swe + snowfall - melt``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.hydrology import snowpack_water_equivalent
    >>> np.round(snowpack_water_equivalent([270, 275], [0.01, 0.002], 0.001, 273.15, 2), 7).tolist()
    [[0.01, 0.0], [0.0063, 0.0057]]
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    precipitation = quantity("precipitation", precipitation, lower=0)
    degree_day_factor = quantity("degree_day_factor", degree_day_factor, lower=0)
    threshold_temperature = quantity("threshold_temperature", threshold_temperature, lower=0, strict_lower=True)
    time_step = quantity("time_step", time_step, lower=0, strict_lower=True)
    require(temperature.ndim == 1 and precipitation.shape == temperature.shape, "temperature and precipitation must be matching vectors")
    require(degree_day_factor.ndim == threshold_temperature.ndim == time_step.ndim == 0, "snowpack parameters must be scalars")
    potential = degree_day_melt(temperature, degree_day_factor, threshold_temperature) * time_step
    snowing = temperature <= threshold_temperature
    result = np.empty((temperature.size, 2))
    snow = 0.0
    for step in range(temperature.size):
        snow += precipitation[step] if snowing[step] else 0.0
        melt = min(snow, potential[step])
        snow -= melt
        result[step] = snow, melt + (0.0 if snowing[step] else precipitation[step])
    return result
