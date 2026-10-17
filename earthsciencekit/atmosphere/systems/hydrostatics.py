"""Hydrostatics: air density, scale height, barometric pressure, the barometer, a layered profile and pressure tendency.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import R_DRY_AIR

__all__ = [
    "dry_air_density",
    "scale_height",
    "barometric_pressure",
    "barometer_column_height",
    "two_layer_temperature",
    "surface_pressure_tendency",
]


@model
def dry_air_density(pressure, temperature):
    """Ideal-gas atmospheric density.

    An equation of state connects atmospheric pressure, temperature and density.

    Parameters
    ----------
    pressure : float or array_like
        Pa; constraint: positive.
    temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dry ideal air with fixed specific gas constant.
    Equation: ``pressure / (R_DRY_AIR * temperature)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    return pressure / (R_DRY_AIR * temperature)


@model
def scale_height(temperature, gravity):
    """Isothermal scale height.

    Hydrostatic balance predicts an exponential pressure scale.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dry isothermal atmosphere with constant gravity.
    Equation: ``R_DRY_AIR * temperature / gravity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return R_DRY_AIR * temperature / gravity


@model
def barometric_pressure(surface_pressure, height, temperature, gravity):
    """Isothermal barometric profile.

    Vertical pressure structure follows from hydrostatic balance and the gas law.

    Parameters
    ----------
    surface_pressure : float or array_like
        Pa; constraint: positive.
    height : float or array_like
        m; constraint: real.
    temperature : float or array_like
        K; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Isothermal column; signed height is relative to the reference level.
    Equation: ``surface_pressure * np.exp(-gravity * height / (R_DRY_AIR * temperature))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    surface_pressure = quantity("surface_pressure", surface_pressure, lower=0, strict_lower=True)
    height = quantity("height", height)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return surface_pressure * np.exp(-gravity * height / (R_DRY_AIR * temperature))


@model
def barometer_column_height(pressure, liquid_density, gravity):
    """Liquid barometer column height.

    The weight of a liquid column balances the pressure of the air above it.

    Parameters
    ----------
    pressure : float or array_like
        Pa; constraint: nonnegative.
    liquid_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Incompressible liquid with a vacuum above the column; vapour pressure and capillarity are neglected.
    Equation: ``pressure / (liquid_density * gravity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    pressure = quantity("pressure", pressure, lower=0)
    liquid_density = quantity("liquid_density", liquid_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return pressure / (liquid_density * gravity)


@model
def two_layer_temperature(height, surface_temperature, lapse_rate, tropopause_height):
    """Troposphere-stratosphere temperature profile.

    Temperature falls with height up to the tropopause and stays constant above it.

    Parameters
    ----------
    height : float or array_like
        m; constraint: nonnegative.
    surface_temperature : float or array_like
        K; constraint: positive.
    lapse_rate : float or array_like
        K/m; constraint: nonnegative.
    tropopause_height : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant-lapse troposphere under an isothermal lower stratosphere; the stratospheric warming above about 20 km is omitted.
    Equation: ``surface_temperature - lapse_rate * np.minimum(height, tropopause_height)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    height = quantity("height", height, lower=0)
    surface_temperature = quantity("surface_temperature", surface_temperature, lower=0, strict_lower=True)
    lapse_rate = quantity("lapse_rate", lapse_rate, lower=0)
    tropopause_height = quantity("tropopause_height", tropopause_height, lower=0, strict_lower=True)
    require(surface_temperature > lapse_rate * tropopause_height, "tropopause temperature must be positive")
    return surface_temperature - lapse_rate * np.minimum(height, tropopause_height)


@model
def surface_pressure_tendency(mean_divergence, surface_pressure):
    """Surface pressure tendency from column divergence.

    Mass leaving the column through net horizontal divergence lowers the surface pressure.

    Parameters
    ----------
    mean_divergence : float or array_like
        1/s; constraint: real. Pressure-weighted column mean of the horizontal divergence.
    surface_pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hydrostatic column with no surface mass flux; divergence at one level is far larger than the column mean.
    Equation: ``-surface_pressure * mean_divergence``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    mean_divergence = quantity("mean_divergence", mean_divergence)
    surface_pressure = quantity("surface_pressure", surface_pressure, lower=0, strict_lower=True)
    return -surface_pressure * mean_divergence
