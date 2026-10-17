"""Radiation: absorbed sunlight, emission, greenhouse layers, insolation and radiative-convective profiles.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity
from earthsciencekit.constants import SIGMA_SB

__all__ = [
    "absorbed_solar",
    "effective_temperature",
    "outgoing_longwave",
    "greenhouse_surface_temperature",
    "daily_insolation",
    "radiative_convective_temperature",
]


@model
def absorbed_solar(solar_constant, albedo):
    """Planetary mean absorbed sunlight.

    Sphere geometry and reflection set the global mean solar input.

    Parameters
    ----------
    solar_constant : float or array_like
        W/m^2; constraint: nonnegative.
    albedo : float or array_like
        1; constraint: fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Global mean rather than local insolation; fixed bolometric albedo.
    Equation: ``solar_constant * (1 - albedo) / 4``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    solar_constant = quantity("solar_constant", solar_constant, lower=0)
    albedo = quantity("albedo", albedo, lower=0, upper=1)
    return solar_constant * (1 - albedo) / 4


@model
def effective_temperature(absorbed_flux, emissivity):
    """Radiative equilibrium temperature.

    Stefan-Boltzmann emission links an energy balance to temperature.

    Parameters
    ----------
    absorbed_flux : float or array_like
        W/m^2; constraint: nonnegative.
    emissivity : float or array_like
        1; constraint: open fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Single emitting level in equilibrium; not the surface temperature of a greenhouse atmosphere.
    Equation: ``(absorbed_flux / (emissivity * SIGMA_SB))**0.25``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    absorbed_flux = quantity("absorbed_flux", absorbed_flux, lower=0)
    emissivity = quantity("emissivity", emissivity, lower=0, upper=1, strict_lower=True)
    return (absorbed_flux / (emissivity * SIGMA_SB)) ** 0.25


@model
def outgoing_longwave(temperature, emissivity):
    """Greybody infrared emission.

    A greybody approximation maps temperature to outgoing radiant flux.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    emissivity : float or array_like
        1; constraint: fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant emissivity; not a spectral atmospheric radiative-transfer solver.
    Equation: ``emissivity * SIGMA_SB * temperature**4``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    emissivity = quantity("emissivity", emissivity, lower=0, upper=1)
    return emissivity * SIGMA_SB * temperature**4


@model
def greenhouse_surface_temperature(absorbed_flux, layer_emissivity):
    """Surface temperature under one absorbing atmospheric layer.

    A layer that absorbs and re-emits infrared returns part of the surface emission and warms the ground.

    Parameters
    ----------
    absorbed_flux : float or array_like
        W/m^2; constraint: nonnegative.
    layer_emissivity : float or array_like
        1; constraint: fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Blackbody surface, a grey isothermal layer transparent to sunlight; the layer temperature is the surface temperature over 2**0.25.
    Equation: ``(2 * absorbed_flux / ((2 - layer_emissivity) * SIGMA_SB))**0.25``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    absorbed_flux = quantity("absorbed_flux", absorbed_flux, lower=0)
    layer_emissivity = quantity("layer_emissivity", layer_emissivity, lower=0, upper=1)
    return (2 * absorbed_flux / ((2 - layer_emissivity) * SIGMA_SB)) ** 0.25


@model
def daily_insolation(solar_constant, latitude, declination, distance_ratio):
    """Daily-mean top-of-atmosphere insolation.

    Latitude, solar declination and Earth-Sun distance set the sunlight a latitude receives over a day.

    Parameters
    ----------
    solar_constant : float or array_like
        W/m^2 at the mean Earth-Sun distance; constraint: nonnegative.
    latitude : float or array_like
        rad; constraint: closed interval [-pi/2, pi/2].
    declination : float or array_like
        rad; constraint: closed interval [-pi/2, pi/2].
    distance_ratio : float or array_like
        1, mean over actual Earth-Sun distance; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Declination fixed over the day, no atmosphere; the sunset hour angle h0 is clipped for polar day and night.
    Equation: ``solar_constant * distance_ratio**2 / pi * (h0 * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(h0))``
    with ``phi = latitude``, ``delta = declination`` and ``h0 = arccos(clip(-tan(phi) * tan(delta), -1, 1))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    solar_constant = quantity("solar_constant", solar_constant, lower=0)
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    declination = quantity("declination", declination, lower=-np.pi / 2, upper=np.pi / 2)
    distance_ratio = quantity("distance_ratio", distance_ratio, lower=0, strict_lower=True)
    hour_angle = np.arccos(np.clip(-np.tan(latitude) * np.tan(declination), -1, 1))
    geometry = hour_angle * np.sin(latitude) * np.sin(declination) + np.cos(latitude) * np.cos(declination) * np.sin(hour_angle)
    return solar_constant * distance_ratio**2 / np.pi * geometry


@model
def radiative_convective_temperature(height, effective_temperature, lapse_rate, emission_height):
    """Radiative-convective equilibrium temperature profile.

    Convection holds the troposphere at a fixed lapse rate below a radiatively balanced stratosphere.

    Parameters
    ----------
    height : float or array_like
        m above the surface; constraint: nonnegative.
    effective_temperature : float or array_like
        K; constraint: positive.
    lapse_rate : float or array_like
        K/m; constraint: positive.
    emission_height : float or array_like
        m; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Grey optically thin stratosphere at the skin temperature effective_temperature / 2**0.25; the troposphere
    passes through effective_temperature at emission_height, so the surface is effective_temperature + lapse_rate * emission_height.
    Equation: ``np.maximum(effective_temperature + lapse_rate * (emission_height - height), effective_temperature / 2**0.25)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    height = quantity("height", height, lower=0)
    effective_temperature = quantity("effective_temperature", effective_temperature, lower=0, strict_lower=True)
    lapse_rate = quantity("lapse_rate", lapse_rate, lower=0, strict_lower=True)
    emission_height = quantity("emission_height", emission_height, lower=0)
    return np.maximum(effective_temperature + lapse_rate * (emission_height - height), effective_temperature / 2**0.25)
