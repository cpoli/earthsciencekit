"""Firn: Herron and Langway's steady densification of snow into ice.

SI units unless stated otherwise. Accumulation is in m/s of ice, as in :func:`gas_age_difference`.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import ICE_DENSITY, R_GAS, SECONDS_PER_YEAR

__all__ = [
    "herron_langway_density",
    "herron_langway_depth",
    "herron_langway_mean_density",
]

CRITICAL_DENSITY = 550.0  # kg/m^3, end of grain settling and packing


def _logit(density):
    return np.log(density / (ICE_DENSITY - density))


def _rates(temperature, accumulation, surface_density):
    """Stage slopes of logit(density) against depth [1/m] and the critical depth [m]."""
    water_equivalent = accumulation * SECONDS_PER_YEAR * ICE_DENSITY / 1000  # m/yr, as fitted
    upper = ICE_DENSITY / 1000 * 11 * np.exp(-10160 / (R_GAS * temperature))
    lower = ICE_DENSITY / 1000 * 575 * np.exp(-21400 / (R_GAS * temperature)) / np.sqrt(water_equivalent)
    critical_depth = (_logit(CRITICAL_DENSITY) - _logit(surface_density)) / upper
    return upper, lower, critical_depth


def _argument(depth, temperature, accumulation, surface_density):
    """logit(density) at each depth, linear in depth within each stage."""
    upper, lower, critical_depth = _rates(temperature, accumulation, surface_density)
    return np.where(
        depth <= critical_depth,
        _logit(surface_density) + upper * depth,
        _logit(CRITICAL_DENSITY) + lower * (depth - critical_depth),
    )


def _firn_inputs(temperature, accumulation, surface_density):
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True, upper=273.15)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    surface_density = quantity("surface_density", surface_density, lower=0, strict_lower=True)
    require(surface_density < CRITICAL_DENSITY, "surface density must be below 550 kg/m^3")
    return temperature, accumulation, surface_density


@model
def herron_langway_density(depth, temperature, accumulation, surface_density):
    """Density of steady-state firn at a depth, by the Herron-Langway model.

    Snow compacts under its own weight into firn and then into ice, so the density profile sets how deep air stays
    open to the atmosphere and how much the ice is younger than the snow it buries.

    Parameters
    ----------
    depth : float or array_like
        m below the surface; constraint: nonnegative.
    temperature : float or array_like
        K, mean annual firn temperature; constraint: positive, at most 273.15.
    accumulation : float or array_like
        m/s of ice; constraint: positive.
    surface_density : float or array_like
        kg/m^3, about 330 to 360; constraint: positive, below 550.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3], tending to ICE_DENSITY at depth. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Herron and Langway (1980): an empirical fit to Greenland and Antarctic cores for dry firn in steady state, with
    ``ln(density / (ICE_DENSITY - density))`` linear in depth in two stages, grain packing above 550 kg/m^3 and
    sintering below. k0 = 11 exp(-10160 / (R T)) and k1 = 575 exp(-21400 / (R T)) take densities in Mg/m^3 and the
    accumulation A in m/yr of water; melt, seasonal layering and transient climate are not represented.
    Equation: ``ICE_DENSITY / (1 + np.exp(-x))``, with ``x = logit(surface_density) + rho_i k0 * depth`` above the
    critical depth and ``x = logit(550) + rho_i k1 / A**0.5 * (depth - critical_depth)`` below it.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import herron_langway_density
    >>> round(float(herron_langway_density(50, 248.15, 1e-09, 350)), 4)
    901.7674
    """
    depth = quantity("depth", depth, lower=0)
    temperature, accumulation, surface_density = _firn_inputs(temperature, accumulation, surface_density)
    return ICE_DENSITY / (1 + np.exp(-_argument(depth, temperature, accumulation, surface_density)))


@model
def herron_langway_depth(density, temperature, accumulation, surface_density):
    """Depth at which steady-state Herron-Langway firn reaches a density.

    The depth of the pore close-off density, about 830 kg/m^3, is where the air in the firn is sealed into bubbles.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: at least surface_density, below ICE_DENSITY.
    temperature : float or array_like
        K, mean annual firn temperature; constraint: positive, at most 273.15.
    accumulation : float or array_like
        m/s of ice; constraint: positive.
    surface_density : float or array_like
        kg/m^3; constraint: positive, below 550.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The inverse of :func:`herron_langway_density`, with the same empirical limits.
    Equation: ``(logit(density) - logit(surface_density)) / (rho_i k0)`` up to 550 kg/m^3, then
    ``critical_depth + (logit(density) - logit(550)) * A**0.5 / (rho_i k1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import herron_langway_depth
    >>> round(float(herron_langway_depth(830, 248.15, 1e-09, 350)), 4)
    31.1813
    """
    density = quantity("density", density, lower=0, strict_lower=True, upper=ICE_DENSITY, strict_upper=True)
    temperature, accumulation, surface_density = _firn_inputs(temperature, accumulation, surface_density)
    require(density >= surface_density, "density must be at least the surface density")
    upper, lower, critical_depth = _rates(temperature, accumulation, surface_density)
    return np.where(
        density <= CRITICAL_DENSITY,
        (_logit(density) - _logit(surface_density)) / upper,
        critical_depth + (_logit(density) - _logit(CRITICAL_DENSITY)) / lower,
    )


@model
def herron_langway_mean_density(depth, temperature, accumulation, surface_density):
    """Mean density of steady-state Herron-Langway firn between the surface and a depth.

    The firn mass above the close-off depth, divided by the yearly accumulation, is the age of the ice there and so the
    ice-age minus gas-age difference of :func:`gas_age_difference`.

    Parameters
    ----------
    depth : float or array_like
        m below the surface; constraint: positive.
    temperature : float or array_like
        K, mean annual firn temperature; constraint: positive, at most 273.15.
    accumulation : float or array_like
        m/s of ice; constraint: positive.
    surface_density : float or array_like
        kg/m^3; constraint: positive, below 550.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The exact integral of :func:`herron_langway_density`, with the same empirical limits: in each stage of slope c the
    mass between densities a and b is ``ICE_DENSITY / c * ln((ICE_DENSITY - a) / (ICE_DENSITY - b))``, evaluated as
    ``ln(1 + exp(x))`` differences so that it stays finite where the density has reached ICE_DENSITY.
    Equation: ``(stage-one mass + stage-two mass) / depth``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import herron_langway_mean_density
    >>> round(float(herron_langway_mean_density(50, 248.15, 1e-09, 350)), 4)
    710.4412
    """
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    temperature, accumulation, surface_density = _firn_inputs(temperature, accumulation, surface_density)
    upper, lower, _ = _rates(temperature, accumulation, surface_density)
    bottom = _argument(depth, temperature, accumulation, surface_density)
    turn = np.minimum(bottom, _logit(CRITICAL_DENSITY))
    mass = ICE_DENSITY / upper * (np.logaddexp(0, turn) - np.logaddexp(0, _logit(surface_density)))
    mass = mass + ICE_DENSITY / lower * (np.logaddexp(0, bottom) - np.logaddexp(0, turn))
    return mass / depth
