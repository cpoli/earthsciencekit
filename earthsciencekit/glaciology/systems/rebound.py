"""Rebound: glacial isostasy and the viscous relaxation of the mantle after deglaciation.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity

__all__ = [
    "isostatic_depression",
    "viscous_relaxation_time",
    "haskell_viscosity",
    "remaining_depression",
    "rebound_uplift_rate",
]


@model
def isostatic_depression(ice_thickness, ice_density, mantle_density):
    """Equilibrium depression of the land under an ice sheet.

    The mantle under a wide ice sheet gives way until the displaced rock weighs as much as the ice, about a third of its
    thickness.

    Parameters
    ----------
    ice_thickness : float or array_like
        m; constraint: nonnegative.
    ice_density : float or array_like
        kg/m^3; constraint: positive.
    mantle_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], positive downward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Jamieson (1865): local isostasy for a load much wider than the flexural wavelength of the lithosphere.
    Equation: ``ice_thickness * ice_density / mantle_density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import isostatic_depression
    >>> round(float(isostatic_depression(3000, 900, 3300)), 3)
    818.182
    """
    ice_thickness = quantity("ice_thickness", ice_thickness, lower=0)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    mantle_density = quantity("mantle_density", mantle_density, lower=0, strict_lower=True)
    return ice_thickness * ice_density / mantle_density


@model
def viscous_relaxation_time(wavelength, viscosity, density, gravity):
    """Relaxation time of a surface depression on a viscous half-space.

    A broad depression in a viscous mantle recovers faster than a narrow one, so the timing of rebound measures the
    viscosity.

    Parameters
    ----------
    wavelength : float or array_like
        m, of the depression; constraint: positive.
    viscosity : float or array_like
        Pa s, of the mantle; constraint: positive.
    density : float or array_like
        kg/m^3, of the mantle; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s], an e-folding time. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Haskell (1935): a uniform Newtonian half-space without an elastic lithosphere, which slows the shortest wavelengths.
    Equation: ``4 * np.pi * viscosity / (density * gravity * wavelength)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import viscous_relaxation_time
    >>> round(float(viscous_relaxation_time(1000000.0, 1e+21, 3300, 10)), 0)
    380799109526.0
    """
    wavelength = quantity("wavelength", wavelength, lower=0, strict_lower=True)
    viscosity = quantity("viscosity", viscosity, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return 4 * np.pi * viscosity / (density * gravity * wavelength)


@model
def haskell_viscosity(relaxation_time, wavelength, density, gravity):
    """Mantle viscosity inferred from a rebound relaxation time.

    Haskell turned the uplift of Fennoscandia into the first estimate of the viscosity of the mantle, about 1e21 Pa s.

    Parameters
    ----------
    relaxation_time : float or array_like
        s; constraint: positive.
    wavelength : float or array_like
        m, of the depression; constraint: positive.
    density : float or array_like
        kg/m^3, of the mantle; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Haskell (1935), the inverse of viscous_relaxation_time for a uniform half-space.
    Equation: ``density * gravity * wavelength * relaxation_time / (4 * np.pi)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import haskell_viscosity
    >>> round(float(haskell_viscosity(100000000000.0, 1000000.0, 3300, 10)), 0)
    2.6260565610162733e+20
    """
    relaxation_time = quantity("relaxation_time", relaxation_time, lower=0, strict_lower=True)
    wavelength = quantity("wavelength", wavelength, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return density * gravity * wavelength * relaxation_time / (4 * np.pi)


@model
def remaining_depression(time, initial_depression, relaxation_time):
    """Depression still to be recovered after the ice has gone.

    After the ice sheet melts the land rises exponentially, and the uplift still to come falls by e every relaxation
    time.

    Parameters
    ----------
    time : float or array_like
        s, since the ice was removed; constraint: nonnegative.
    initial_depression : float or array_like
        m, at removal; constraint: nonnegative.
    relaxation_time : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], positive downward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Haskell (1935): one wavelength on a viscous half-space, with the load removed at once.
    Equation: ``initial_depression * np.exp(-time / relaxation_time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import remaining_depression
    >>> round(float(remaining_depression(1000.0, 300, 1000.0)), 3)
    110.364
    """
    time = quantity("time", time, lower=0)
    initial_depression = quantity("initial_depression", initial_depression, lower=0)
    relaxation_time = quantity("relaxation_time", relaxation_time, lower=0, strict_lower=True)
    return initial_depression * np.exp(-time / relaxation_time)


@model
def rebound_uplift_rate(time, initial_depression, relaxation_time):
    """Rate of postglacial uplift.

    The land rises fastest just after deglaciation; Fennoscandia still rises by up to 9 mm a year.

    Parameters
    ----------
    time : float or array_like
        s, since the ice was removed; constraint: nonnegative.
    initial_depression : float or array_like
        m, at removal; constraint: nonnegative.
    relaxation_time : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], upward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Haskell (1935), the time derivative of remaining_depression.
    Equation: ``initial_depression / relaxation_time * np.exp(-time / relaxation_time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import rebound_uplift_rate
    >>> round(float(rebound_uplift_rate(1000.0, 300, 1000.0)), 6)
    0.110364
    """
    time = quantity("time", time, lower=0)
    initial_depression = quantity("initial_depression", initial_depression, lower=0)
    relaxation_time = quantity("relaxation_time", relaxation_time, lower=0, strict_lower=True)
    return initial_depression / relaxation_time * np.exp(-time / relaxation_time)
