"""Isostasy and flexure: Airy roots and the bending of an elastic lithosphere under loads.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "airy_root",
    "flexural_rigidity",
    "flexural_parameter",
    "flexural_deflection",
]


@model
def airy_root(elevation, crust_density, mantle_density):
    """Airy crustal compensation.

    Crustal buoyancy provides a baseline for mountain-root support.

    Parameters
    ----------
    elevation : float or array_like
        m; constraint: nonnegative.
    crust_density : float or array_like
        kg/m^3; constraint: positive.
    mantle_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Local compensation without plate flexure; mantle must be denser than crust.
    Equation: ``elevation * crust_density / (mantle_density - crust_density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.plate_tectonics import airy_root
    >>> float(airy_root(1000, 2700, 3300))
    4500.0
    """
    elevation = quantity("elevation", elevation, lower=0)
    crust_density = quantity("crust_density", crust_density, lower=0, strict_lower=True)
    mantle_density = quantity("mantle_density", mantle_density, lower=0, strict_lower=True)
    require(mantle_density > crust_density, "mantle must be denser than crust")
    return elevation * crust_density / (mantle_density - crust_density)


@model
def flexural_rigidity(young_modulus, elastic_thickness, poisson_ratio):
    """Elastic lithosphere rigidity.

    Elastic thickness turns lithospheric strength into a bending resistance.

    Parameters
    ----------
    young_modulus : float or array_like
        Pa; constraint: positive.
    elastic_thickness : float or array_like
        m; constraint: positive.
    poisson_ratio : float or array_like
        1; constraint: poisson.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Thin isotropic elastic plate with stable Poisson ratio, not a yielding viscoelastic lithosphere.
    Equation: ``young_modulus * elastic_thickness**3 / (12 * (1 - poisson_ratio**2))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.plate_tectonics import flexural_rigidity
    >>> float(flexural_rigidity(12, 1, 0))
    1.0
    """
    young_modulus = quantity("young_modulus", young_modulus, lower=0, strict_lower=True)
    elastic_thickness = quantity("elastic_thickness", elastic_thickness, lower=0, strict_lower=True)
    poisson_ratio = quantity("poisson_ratio", poisson_ratio, lower=-1, upper=0.5, strict_lower=True, strict_upper=True)
    return young_modulus * elastic_thickness**3 / (12 * (1 - poisson_ratio**2))


@model
def flexural_parameter(rigidity, density_contrast, gravity):
    """Lithosphere flexural length.

    Bending and buoyancy jointly set the scale of load compensation.

    Parameters
    ----------
    rigidity : float or array_like
        N m; constraint: positive.
    density_contrast : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Infinite beam-like plate convention alpha to the fourth equals 4D over density contrast times g.
    Equation: ``(4 * rigidity / (density_contrast * gravity))**0.25``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.plate_tectonics import flexural_parameter
    >>> float(flexural_parameter(4, 1, 1))
    2.0
    """
    rigidity = quantity("rigidity", rigidity, lower=0, strict_lower=True)
    density_contrast = quantity("density_contrast", density_contrast, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return (4 * rigidity / (density_contrast * gravity)) ** 0.25


@model
def flexural_deflection(distance, line_load, rigidity, density_contrast, gravity):
    """Elastic plate deflection under a line load.

    A load bends the lithosphere into a moat flanked by a peripheral bulge.

    Parameters
    ----------
    distance : float or array_like
        m; constraint: real.
    line_load : float or array_like
        N/m; constraint: nonnegative.
    rigidity : float or array_like
        N m; constraint: positive.
    density_contrast : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Unbroken infinite thin elastic plate, deflection positive downward; density contrast is mantle minus infill.
    Equation: ``line_load * alpha**3 / (8 * rigidity) * np.exp(-x) * (np.cos(x) + np.sin(x))``
    with ``x = abs(distance) / alpha`` and ``alpha = (4 * rigidity / (density_contrast * gravity))**0.25``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.plate_tectonics import flexural_deflection
    >>> float(flexural_deflection(0, 8, 1, 4, 1))
    1.0
    """
    distance = quantity("distance", distance)
    line_load = quantity("line_load", line_load, lower=0)
    rigidity = quantity("rigidity", rigidity, lower=0, strict_lower=True)
    density_contrast = quantity("density_contrast", density_contrast, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    alpha = (4 * rigidity / (density_contrast * gravity)) ** 0.25
    scaled = np.abs(distance) / alpha
    return line_load * alpha**3 / (8 * rigidity) * np.exp(-scaled) * (np.cos(scaled) + np.sin(scaled))
