"""Magma: viscosity, water solubility, ascent in dikes and ground deformation above a magma chamber.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity
from earthsciencekit.constants import R_GAS

__all__ = [
    "arrhenius_viscosity",
    "shaw_viscosity",
    "water_solubility",
    "dike_ascent_velocity",
    "mogi_displacement",
]


@model
def arrhenius_viscosity(temperature, reference_viscosity, activation_energy, reference_temperature):
    """Arrhenius temperature dependence of melt viscosity.

    Silicate melts stiffen by orders of magnitude as they cool, which controls how lava flows and how magma erupts.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    reference_viscosity : float or array_like
        Pa s, viscosity at the reference temperature; constraint: positive.
    activation_energy : float or array_like
        J/mol; constraint: nonnegative.
    reference_temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant activation energy over the range of interest; crystal- and bubble-free Newtonian melt.
    Equation: ``reference_viscosity * np.exp(activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import arrhenius_viscosity
    >>> round(float(arrhenius_viscosity(1000, 10, 8314.46261815324, 500)), 5)
    3.67879
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    reference_viscosity = quantity("reference_viscosity", reference_viscosity, lower=0, strict_lower=True)
    activation_energy = quantity("activation_energy", activation_energy, lower=0)
    reference_temperature = quantity("reference_temperature", reference_temperature, lower=0, strict_lower=True)
    return reference_viscosity * np.exp(activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))


@model
def shaw_viscosity(temperature, slope):
    """Shaw's empirical viscosity of an anhydrous or hydrous silicate melt.

    One composition-dependent slope predicts the viscosity of melts from basalt to rhyolite.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    slope : float or array_like
        1, Shaw's s from the molar composition (about 2-3 for basalt, 3.5-4 for rhyolite); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Shaw (1972), Am. J. Sci. 272: 870-893; all lines pass near ln(eta [poise]) = -6.4 at 1e4/T = 1.5. Crystal- and
    bubble-free Newtonian melt; about 0.2 log units of accuracy.
    Equation: ``0.1 * np.exp(slope * (1e4 / temperature) - 1.5 * slope - 6.4)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import shaw_viscosity
    >>> round(float(shaw_viscosity(6666.666666666667, 2)), 9)
    0.000166156
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0, strict_lower=True)
    return 0.1 * np.exp(slope * (1e4 / temperature) - 1.5 * slope - 6.4)


@model
def water_solubility(pressure, solubility_constant):
    """Dissolved water in a melt saturated at a given pressure.

    Water dissolves as hydroxyl, so its solubility grows as the square root of pressure and ascending magma exsolves
    bubbles.

    Parameters
    ----------
    pressure : float or array_like
        Pa; constraint: nonnegative.
    solubility_constant : float or array_like
        1/sqrt(Pa), mass fraction per root pascal (about 4.1e-6 for rhyolite); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], mass fraction. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Square-root (Henry-type) law for pure water vapour, valid up to a few hundred MPa; ignores CO2 and temperature.
    Equation: ``solubility_constant * np.sqrt(pressure)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import water_solubility
    >>> float(water_solubility(400000000.0, 4e-06))
    0.08
    """
    pressure = quantity("pressure", pressure, lower=0)
    solubility_constant = quantity("solubility_constant", solubility_constant, lower=0, strict_lower=True)
    return solubility_constant * np.sqrt(pressure)


@model
def dike_ascent_velocity(width, density_contrast, viscosity, gravity):
    """Mean laminar velocity of buoyant magma in a dike.

    Buoyancy drives magma up a planar crack, and the viscous drag on its walls sets how fast it rises.

    Parameters
    ----------
    width : float or array_like
        m, full dike width; constraint: positive.
    density_contrast : float or array_like
        kg/m^3, host rock minus magma; constraint: real.
    viscosity : float or array_like
        Pa s; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], positive upward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Wilson and Head (1981): plane Poiseuille flow between parallel walls driven by buoyancy alone; check that the Reynolds
    number ``magma_density * velocity * width / viscosity`` stays below about 2000.
    Equation: ``width**2 * density_contrast * gravity / (12 * viscosity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import dike_ascent_velocity
    >>> float(dike_ascent_velocity(2, 300, 100, 10))
    10.0
    """
    width = quantity("width", width, lower=0, strict_lower=True)
    density_contrast = quantity("density_contrast", density_contrast)
    viscosity = quantity("viscosity", viscosity, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return width**2 * density_contrast * gravity / (12 * viscosity)


@model
def mogi_displacement(radial_distance, depth, volume_change, poisson_ratio):
    """Surface displacement above an inflating point source.

    Magma entering a chamber lifts the ground in a bell-shaped dome whose width reveals the chamber depth.

    Parameters
    ----------
    radial_distance : float or array_like
        m, horizontal distance from the point above the source; constraint: nonnegative.
    depth : float or array_like
        m; constraint: positive.
    volume_change : float or array_like
        m^3, signed change of the source volume; constraint: real.
    poisson_ratio : float or array_like
        1; constraint: between 0 and 0.5.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [m]: ``[radial, vertical]`` displacement, outward and upward positive.

    Notes
    -----
    Mogi (1958) point source in an elastic half-space; the chamber radius must be small compared with its depth. The
    uplift halves at ``radial_distance = 0.766 * depth``.
    Equation: ``(1 - poisson_ratio) * volume_change / np.pi * [radial_distance, depth] / (radial_distance**2 + depth**2)**1.5``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.volcanology import mogi_displacement
    >>> np.round(mogi_displacement(0, 1000, 3.141592653589793, 0.5), 12).tolist()
    [0.0, 5e-07]
    """
    radial_distance = quantity("radial_distance", radial_distance, lower=0)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    volume_change = quantity("volume_change", volume_change)
    poisson_ratio = quantity("poisson_ratio", poisson_ratio, lower=0, upper=0.5)
    scale = (1 - poisson_ratio) * volume_change / np.pi / (radial_distance**2 + depth**2) ** 1.5
    return np.stack(np.broadcast_arrays(scale * radial_distance, scale * depth), axis=-1)
