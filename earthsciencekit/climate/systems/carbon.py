"""Carbon budget: airborne fraction and single-reservoir decay.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity

__all__ = [
    "carbon_airborne_fraction",
    "carbon_decay",
]


@model
def carbon_airborne_fraction(atmospheric_increase, emissions):
    """Carbon budget fraction.

    A carbon mass budget partitions emissions and atmospheric accumulation.

    Parameters
    ----------
    atmospheric_increase : float or array_like
        kg C; constraint: nonnegative.
    emissions : float or array_like
        kg C; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Same time interval and carbon mass units; feedback-driven values above one are not clipped.
    Equation: ``atmospheric_increase / emissions``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import carbon_airborne_fraction
    >>> float(carbon_airborne_fraction(4, 10))
    0.4
    """
    atmospheric_increase = quantity("atmospheric_increase", atmospheric_increase, lower=0)
    emissions = quantity("emissions", emissions, lower=0, strict_lower=True)
    return atmospheric_increase / emissions


@model
def carbon_decay(initial_excess, timescale, time):
    """Single-reservoir carbon relaxation.

    A reservoir model illustrates first-order removal and residence time.

    Parameters
    ----------
    initial_excess : float or array_like
        kg C; constraint: nonnegative.
    timescale : float or array_like
        s; constraint: positive.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Not a realistic complete atmospheric CO2 impulse response, which has multiple timescales.
    Equation: ``initial_excess * np.exp(-time / timescale)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import carbon_decay
    >>> float(carbon_decay(2, 1, 0.6931471805599453))
    1.0
    """
    initial_excess = quantity("initial_excess", initial_excess, lower=0)
    timescale = quantity("timescale", timescale, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0)
    return initial_excess * np.exp(-time / timescale)
