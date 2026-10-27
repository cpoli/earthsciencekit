"""Hillslopes: diffusive and nonlinear soil creep, fault-scarp degradation and soil production.

SI units unless stated otherwise.
"""

import numpy as np
from scipy.special import erf

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "diffusive_hillslope_flux",
    "nonlinear_hillslope_flux",
    "scarp_elevation",
    "scarp_maximum_gradient",
    "soil_production_rate",
    "hilltop_curvature",
]


@model
def diffusive_hillslope_flux(diffusivity, gradient):
    """Linear diffusive sediment flux of soil creep.

    Creep moves soil downhill at a rate proportional to the slope, so convex hilltops lower and their profiles obey a
    diffusion equation.

    Parameters
    ----------
    diffusivity : float or array_like
        m^2/s, the hillslope transport coefficient; constraint: nonnegative.
    gradient : float or array_like
        1, the elevation gradient dz/dx; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s] of soil volume per unit contour width, positive toward increasing x. Scalar input returns a
        NumPy scalar; arrays broadcast.

    Notes
    -----
    Gilbert (1909) and Culling (1960): gentle soil-mantled slopes; the flux grows faster than linearly near the
    angle of repose (see :func:`nonlinear_hillslope_flux`).
    Equation: ``-diffusivity * gradient``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import diffusive_hillslope_flux
    >>> float(diffusive_hillslope_flux(0.01, 0.2))
    -0.002
    """
    diffusivity = quantity("diffusivity", diffusivity, lower=0)
    gradient = quantity("gradient", gradient)
    return -diffusivity * gradient


@model
def nonlinear_hillslope_flux(diffusivity, gradient, critical_gradient):
    """Nonlinear soil flux that diverges at a critical gradient.

    Steep slopes shed soil far faster than linear creep predicts, which is why hillslopes near their threshold angle
    become planar.

    Parameters
    ----------
    diffusivity : float or array_like
        m^2/s; constraint: nonnegative.
    gradient : float or array_like
        1, the elevation gradient dz/dx; constraint: magnitude less than critical_gradient.
    critical_gradient : float or array_like
        1, about 1.25 on soil-mantled slopes; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s], positive toward increasing x. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Roering, Kirchner and Dietrich (1999): reduces to :func:`diffusive_hillslope_flux` on gentle slopes; landslides
    take over at the critical gradient.
    Equation: ``-diffusivity * gradient / (1 - (gradient / critical_gradient)**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import nonlinear_hillslope_flux
    >>> round(float(nonlinear_hillslope_flux(0.01, 0.5, 1.0)), 14)
    -0.00666666666667
    """
    diffusivity = quantity("diffusivity", diffusivity, lower=0)
    gradient = quantity("gradient", gradient)
    critical_gradient = quantity("critical_gradient", critical_gradient, lower=0, strict_lower=True)
    require(np.abs(gradient) < critical_gradient, "gradient must be gentler than critical_gradient")
    return -diffusivity * gradient / (1 - (gradient / critical_gradient) ** 2)


@model
def scarp_elevation(distance, half_offset, diffusivity, time, far_field_gradient):
    """Elevation across a fault scarp degraded by linear diffusion.

    A fresh step in the land rounds off at a rate set by the hillslope diffusivity, so its shape records its age.

    Parameters
    ----------
    distance : float or array_like
        m from the scarp midpoint, positive toward the upthrown side; constraint: real.
    half_offset : float or array_like
        m, half the vertical offset of the initial step; constraint: nonnegative.
    diffusivity : float or array_like
        m^2/s; constraint: nonnegative.
    time : float or array_like
        s since the scarp formed; constraint: nonnegative.
    far_field_gradient : float or array_like
        1, the regional slope on both sides; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m] relative to the midpoint. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Culling (1960) and Hanks et al. (1984): a vertical step on a uniform slope degraded by linear diffusion without
    base-level change; a fresh step (``diffusivity * time = 0``) is returned exactly.
    Equation: ``half_offset * erf(distance / (2 * np.sqrt(diffusivity * time))) + far_field_gradient * distance``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import scarp_elevation
    >>> round(float(scarp_elevation(1, 2, 0.25, 1, 0.1)), 11)
    1.7854015859
    """
    distance = quantity("distance", distance)
    half_offset = quantity("half_offset", half_offset, lower=0)
    diffusivity = quantity("diffusivity", diffusivity, lower=0)
    time = quantity("time", time, lower=0)
    far_field_gradient = quantity("far_field_gradient", far_field_gradient)
    width = 2 * np.sqrt(diffusivity * time)
    aged = width > 0
    step = np.where(aged, erf(distance / np.where(aged, width, 1.0)), np.sign(distance))
    return half_offset * step + far_field_gradient * distance


@model
def scarp_maximum_gradient(half_offset, diffusivity, time, far_field_gradient):
    """Steepest gradient of a diffusively degraded fault scarp.

    The midpoint slope of a scarp falls as the inverse square root of its age, the measurement behind morphologic
    dating.

    Parameters
    ----------
    half_offset : float or array_like
        m; constraint: nonnegative.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    time : float or array_like
        s; constraint: positive.
    far_field_gradient : float or array_like
        1; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], the gradient at the scarp midpoint. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hanks et al. (1984), the derivative of :func:`scarp_elevation` at its midpoint; same assumptions.
    Equation: ``half_offset / np.sqrt(np.pi * diffusivity * time) + far_field_gradient``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import scarp_maximum_gradient
    >>> round(float(scarp_maximum_gradient(2, 0.25, 4, 0.1)), 11)
    1.2283791671
    """
    half_offset = quantity("half_offset", half_offset, lower=0)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0, strict_lower=True)
    far_field_gradient = quantity("far_field_gradient", far_field_gradient)
    return half_offset / np.sqrt(np.pi * diffusivity * time) + far_field_gradient


@model
def soil_production_rate(soil_depth, maximum_rate, decay_depth):
    """Exponential soil production function.

    Bedrock turns into soil fastest where it is bare and more slowly under a thick soil that shields it, which keeps
    soil-mantled hillslopes in balance.

    Parameters
    ----------
    soil_depth : float or array_like
        m, vertical thickness of mobile soil; constraint: nonnegative.
    maximum_rate : float or array_like
        m/s of bedrock lowering under bare ground; constraint: nonnegative.
    decay_depth : float or array_like
        m, about 0.5 m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s] of bedrock converted to soil. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Heimsath et al. (1997), from cosmogenic nuclides under soils of a few metres at most; some sites show a humped
    function with a maximum under a thin soil.
    Equation: ``maximum_rate * np.exp(-soil_depth / decay_depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import soil_production_rate
    >>> round(float(soil_production_rate(0.5, 1e-11, 0.5)), 23)
    3.67879441171e-12
    """
    soil_depth = quantity("soil_depth", soil_depth, lower=0)
    maximum_rate = quantity("maximum_rate", maximum_rate, lower=0)
    decay_depth = quantity("decay_depth", decay_depth, lower=0, strict_lower=True)
    return maximum_rate * np.exp(-soil_depth / decay_depth)


@model
def hilltop_curvature(erosion_rate, diffusivity, rock_density, soil_density):
    """Hilltop curvature of a diffusive hillslope in steady state.

    Where creep balances uniform lowering, faster erosion makes sharper ridges, so the curvature of a hilltop measures
    its erosion rate.

    Parameters
    ----------
    erosion_rate : float or array_like
        m/s of bedrock lowering; constraint: nonnegative.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    rock_density : float or array_like
        kg/m^3; constraint: positive.
    soil_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1/m], the second derivative of elevation, negative for a convex ridge. Scalar input returns a
        NumPy scalar; arrays broadcast.

    Notes
    -----
    Roering, Perron and Kirchner (2007): steady state with linear creep at the divide, where the gradient is small.
    Equation: ``-rock_density * erosion_rate / (soil_density * diffusivity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import hilltop_curvature
    >>> float(hilltop_curvature(1e-12, 1e-10, 2600, 1300))
    -0.02
    """
    erosion_rate = quantity("erosion_rate", erosion_rate, lower=0)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    rock_density = quantity("rock_density", rock_density, lower=0, strict_lower=True)
    soil_density = quantity("soil_density", soil_density, lower=0, strict_lower=True)
    return -rock_density * erosion_rate / (soil_density * diffusivity)
