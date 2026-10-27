"""Surface melt: Reeh's positive-degree-day model of ice-sheet ablation.

SI units unless stated otherwise. Snowfall and melt are water equivalent, as in ``hydrology.degree_day_melt``.
"""

import numpy as np
from scipy.special import erfc

from earthsciencekit._validation import model, quantity

__all__ = [
    "expected_positive_degree_days",
    "pdd_surface_melt",
]


@model
def expected_positive_degree_days(mean_temperature, standard_deviation, duration):
    """Expected positive degree-days of a period with normally distributed temperature.

    A month whose mean is below freezing still has warm days that melt ice, so ice-sheet models count the positive
    part of a temperature that fluctuates about its mean rather than the positive part of the mean.

    Parameters
    ----------
    mean_temperature : float or array_like
        K, mean air temperature of the period; constraint: positive.
    standard_deviation : float or array_like
        K, of the daily temperature about that mean, about 4.5 to 5 on Greenland; constraint: positive.
    duration : float or array_like
        s, length of the period; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K s]; divide by 86400 for degree-days. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Reeh (1991), in the closed form of Calov and Greve (2005): the expected value of ``max(T - 273.15, 0)`` for
    ``T`` normal about ``mean_temperature``, held over the period; day-to-day correlation and diurnal cycles are not
    represented, and the deviation tends to ``max(mean_temperature - 273.15, 0)`` as it shrinks.
    Equation: ``duration * (s / np.sqrt(2 * np.pi) * np.exp(-c**2 / (2 * s**2)) + c / 2 * erfc(-c / (np.sqrt(2) * s)))``
    with ``c = mean_temperature - 273.15`` and ``s = standard_deviation``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import expected_positive_degree_days
    >>> round(float(expected_positive_degree_days(273.15, 5, 86400)), 4)
    172343.0651
    """
    mean_temperature = quantity("mean_temperature", mean_temperature, lower=0, strict_lower=True)
    standard_deviation = quantity("standard_deviation", standard_deviation, lower=0, strict_lower=True)
    duration = quantity("duration", duration, lower=0)
    c, s = mean_temperature - 273.15, standard_deviation
    return duration * (s / np.sqrt(2 * np.pi) * np.exp(-(c**2) / (2 * s**2)) + c / 2 * erfc(-c / (np.sqrt(2) * s)))


@model
def pdd_surface_melt(positive_degree_days, snowfall, snow_factor, ice_factor, refreezing_fraction):
    """Snow melt, ice melt and runoff of an ice-sheet surface from its positive degree-days.

    The positive degree-days first melt the year's snow and only then the darker ice beneath, which melts faster, and
    part of the meltwater refreezes in the cold snow, so runoff is less than melt.

    Parameters
    ----------
    positive_degree_days : float or array_like
        K s, from :func:`expected_positive_degree_days` summed over the year; constraint: nonnegative.
    snowfall : float or array_like
        m of water equivalent over the same period; constraint: nonnegative.
    snow_factor : float or array_like
        m/(K s) of water equivalent, about 3e-3 / 86400 (3 mm per degC per day); constraint: positive.
    ice_factor : float or array_like
        m/(K s) of water equivalent, about 8e-3 / 86400; constraint: positive.
    refreezing_fraction : float or array_like
        1, the largest fraction of the snowfall that meltwater can refreeze into, about 0.6; constraint: between 0
        and 1.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [m] of water equivalent: snow melt, ice melt and runoff on the last axis.

    Notes
    -----
    Reeh (1991), with the Braithwaite factors: snow melts at ``snow_factor`` until the snowfall is gone and the
    remaining degree-days melt ice at ``ice_factor``; refreezing retains melt up to ``refreezing_fraction`` of the
    snowfall. One annual layer, without rain, superimposed ice, albedo feedback or firn warming.
    Equation: ``snow = np.minimum(snow_factor * pdd, snowfall)``, ``ice = ice_factor * np.maximum(pdd - snowfall / snow_factor, 0)``,
    ``runoff = snow + ice - np.minimum(refreezing_fraction * snowfall, snow + ice)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.glaciology import pdd_surface_melt
    >>> np.round(pdd_surface_melt(500 * 86400, 0.5, 3e-3 / 86400, 8e-3 / 86400, 0.6), 4).tolist()
    [0.5, 2.6667, 2.8667]
    """
    positive_degree_days = quantity("positive_degree_days", positive_degree_days, lower=0)
    snowfall = quantity("snowfall", snowfall, lower=0)
    snow_factor = quantity("snow_factor", snow_factor, lower=0, strict_lower=True)
    ice_factor = quantity("ice_factor", ice_factor, lower=0, strict_lower=True)
    refreezing_fraction = quantity("refreezing_fraction", refreezing_fraction, lower=0, upper=1)
    snow = np.minimum(snow_factor * positive_degree_days, snowfall)
    ice = ice_factor * np.maximum(positive_degree_days - snowfall / snow_factor, 0)
    melt = snow + ice
    runoff = melt - np.minimum(refreezing_fraction * snowfall, melt)
    return np.stack(np.broadcast_arrays(snow, ice, runoff), axis=-1)
