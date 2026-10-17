"""Moisture: saturation vapour pressure, mixing ratio, specific and relative humidity and virtual temperature.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import EPSILON_WATER, R_VAPOR

__all__ = [
    "saturation_vapor_pressure",
    "mixing_ratio",
    "specific_humidity",
    "virtual_temperature",
    "relative_humidity",
]


@model
def saturation_vapor_pressure(temperature):
    """Integrated Clausius-Clapeyron approximation.

    Phase equilibrium makes saturation vapor pressure strongly temperature dependent.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: temperature.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Liquid-water reference with constant latent heat; restricted to 180-330 K, not an ice saturation formula.
    Equation: ``611.2 * np.exp((2.5e6 / R_VAPOR) * (1 / 273.15 - 1 / temperature))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=180, upper=330)
    return 611.2 * np.exp((2.5e6 / R_VAPOR) * (1 / 273.15 - 1 / temperature))


@model
def mixing_ratio(vapor_pressure, total_pressure):
    """Water-vapor mass mixing ratio.

    Partial pressure determines water mass per unit dry-air mass.

    Parameters
    ----------
    vapor_pressure : float or array_like
        Pa; constraint: nonnegative.
    total_pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Ideal vapor mixture; vapor pressure must be below total pressure.
    Equation: ``EPSILON_WATER * vapor_pressure / (total_pressure - vapor_pressure)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    vapor_pressure = quantity("vapor_pressure", vapor_pressure, lower=0)
    total_pressure = quantity("total_pressure", total_pressure, lower=0, strict_lower=True)
    require(vapor_pressure < total_pressure, "vapor pressure must be below total pressure")
    return EPSILON_WATER * vapor_pressure / (total_pressure - vapor_pressure)


@model
def specific_humidity(mixing_ratio):
    """Specific humidity from mixing ratio.

    Changing the denominator connects two moisture conventions.

    Parameters
    ----------
    mixing_ratio : float or array_like
        kg/kg dry air; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Output is water mass per total moist-air mass.
    Equation: ``mixing_ratio / (1 + mixing_ratio)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    mixing_ratio = quantity("mixing_ratio", mixing_ratio, lower=0)
    return mixing_ratio / (1 + mixing_ratio)


@model
def virtual_temperature(temperature, specific_humidity):
    """Virtual-temperature buoyancy correction.

    Moist air is less dense than dry air at the same pressure and temperature.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    specific_humidity : float or array_like
        kg/kg moist air; constraint: fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Ideal vapor mixture without suspended liquid or ice.
    Equation: ``temperature * (1 + (1 / EPSILON_WATER - 1) * specific_humidity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    specific_humidity = quantity("specific_humidity", specific_humidity, lower=0, upper=1)
    return temperature * (1 + (1 / EPSILON_WATER - 1) * specific_humidity)


@model
def relative_humidity(vapor_pressure, saturation_pressure):
    """Relative humidity fraction.

    Humidity compares actual vapor pressure to a saturation reference.

    Parameters
    ----------
    vapor_pressure : float or array_like
        Pa; constraint: nonnegative.
    saturation_pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Returned fraction may exceed one during supersaturation; no clipping.
    Equation: ``vapor_pressure / saturation_pressure``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    vapor_pressure = quantity("vapor_pressure", vapor_pressure, lower=0)
    saturation_pressure = quantity("saturation_pressure", saturation_pressure, lower=0, strict_lower=True)
    return vapor_pressure / saturation_pressure
