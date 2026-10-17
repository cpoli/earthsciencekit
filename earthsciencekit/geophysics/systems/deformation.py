"""Deformation of the Earth: tidal Love numbers and mass change from time-variable gravity.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "love_number_h2",
    "love_number_k2",
    "tidal_gravimetric_factor",
    "equivalent_water_height",
]


def _rigidity_ratio(rigidity, density, gravity, radius):
    rigidity = quantity("rigidity", rigidity, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    return 19 * rigidity / (2 * density * gravity * radius)


@model
def love_number_h2(rigidity, density, gravity, radius):
    """Degree-2 displacement Love number of a homogeneous elastic Earth.

    The solid Earth rises and falls with the tide by a fraction of the equilibrium tide that measures its rigidity.

    Parameters
    ----------
    rigidity : float or array_like
        Pa; constraint: nonnegative. Shear modulus; zero gives a fluid Earth.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]: radial displacement over the equilibrium tide height.

    Notes
    -----
    Uniform, incompressible, self-gravitating elastic sphere (Kelvin's solution); ``h2 = 5 k2 / 3``.
    Equation: ``2.5 / (1 + 19 * rigidity / (2 * density * gravity * radius))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    return 2.5 / (1 + _rigidity_ratio(rigidity, density, gravity, radius))


@model
def love_number_k2(rigidity, density, gravity, radius):
    """Degree-2 potential Love number of a homogeneous elastic Earth.

    The tidal bulge adds its own gravitational potential, a fraction k2 of the tidal potential.

    Parameters
    ----------
    rigidity : float or array_like
        Pa; constraint: nonnegative. Shear modulus; zero gives a fluid Earth.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uniform, incompressible, self-gravitating elastic sphere (Kelvin's solution); the fluid limit is 3/2.
    Equation: ``1.5 / (1 + 19 * rigidity / (2 * density * gravity * radius))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    return 1.5 / (1 + _rigidity_ratio(rigidity, density, gravity, radius))


@model
def tidal_gravimetric_factor(h2, k2):
    """Gravimetric factor of the body tide.

    A gravimeter on a yielding Earth sees the tide amplified by its uplift and reduced by the bulge's attraction.

    Parameters
    ----------
    h2 : float or array_like
        1; constraint: nonnegative.
    k2 : float or array_like
        1; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]: observed over rigid-Earth tidal gravity change.

    Notes
    -----
    Degree-2 tide on a spherically symmetric elastic Earth, without ocean loading. About 1.16 for the real Earth.
    Equation: ``1 + h2 - 1.5 * k2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    h2 = quantity("h2", h2, lower=0)
    k2 = quantity("k2", k2, lower=0)
    return 1 + h2 - 1.5 * k2


@model
def equivalent_water_height(coefficient_change, degree, load_love_number, radius, mean_density, water_density):
    """Water layer behind a change of one Stokes coefficient.

    Monthly gravity fields from a satellite pair reveal the moving mass of water, ice and groundwater.

    Parameters
    ----------
    coefficient_change : float or array_like
        1; constraint: real. Change of a fully normalized Stokes coefficient of degree ``n``.
    degree : float or array_like
        1; constraint: integer >= 0.
    load_love_number : float or array_like
        1; constraint: greater than -1. Elastic yielding ``k'_n`` of the Earth under the load.
    radius : float or array_like
        m; constraint: positive.
    mean_density : float or array_like
        kg/m^3; constraint: positive. Mean density of the Earth.
    water_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]: coefficient of the same fully normalized harmonic in the water thickness.

    Notes
    -----
    Thin surface layer on an elastic Earth (Wahr et al., 1998); valid while the layer is much thinner than the
    resolution of the truncated expansion.
    Equation: ``radius * mean_density * (2 * degree + 1) * coefficient_change / (3 * water_density * (1 + load_love_number))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    coefficient_change = quantity("coefficient_change", coefficient_change)
    degree = quantity("degree", degree, lower=0)
    require(degree == np.round(degree), "degree must be an integer")
    load_love_number = quantity("load_love_number", load_love_number, lower=-1, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    mean_density = quantity("mean_density", mean_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    return radius * mean_density * (2 * degree + 1) * coefficient_change / (3 * water_density * (1 + load_love_number))
