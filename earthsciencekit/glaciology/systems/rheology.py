"""Rheology: Forbes's viscous glacier, Glen's flow law, laminar flow profiles and Nye's crevasse depth.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import R_GAS

__all__ = [
    "viscous_slab_velocity",
    "glen_strain_rate",
    "glen_stress",
    "glen_rate_factor",
    "effective_viscosity",
    "glen_velocity_profile",
    "glen_surface_velocity",
    "nye_crevasse_depth",
    "water_filled_crevasse_depth",
]


@model
def viscous_slab_velocity(height, thickness, surface_slope, viscosity, density, gravity):
    """Velocity in a Newtonian viscous slab flowing down a slope.

    Forbes argued that a glacier flows like a very viscous fluid, fastest at the surface and slowest against its bed.

    Parameters
    ----------
    height : float or array_like
        m, above the bed; constraint: between 0 and thickness.
    thickness : float or array_like
        m; constraint: positive.
    surface_slope : float or array_like
        1, tangent of the surface slope; constraint: nonnegative.
    viscosity : float or array_like
        Pa s; constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], downslope, zero at a frozen bed. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Forbes (1843): a wide parallel-sided slab with a linear viscosity, no sliding and small slope.
    Equation: ``density * gravity * surface_slope * (2 * thickness * height - height**2) / (2 * viscosity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import viscous_slab_velocity
    >>> float(viscous_slab_velocity(50, 100, 0.1, 10000000000000.0, 1000, 10))
    3.75e-07
    """
    height = quantity("height", height, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    viscosity = quantity("viscosity", viscosity, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(height <= thickness, "height exceeds the ice thickness")
    return density * gravity * surface_slope * (2 * thickness * height - height**2) / (2 * viscosity)


@model
def glen_strain_rate(stress, rate_factor, exponent):
    """Glen's flow law: strain rate of polycrystalline ice.

    Glen's creep tests showed that ice deforms at a rate proportional to roughly the cube of the stress, so it is far
    from a Newtonian fluid.

    Parameters
    ----------
    stress : float or array_like
        Pa, effective (or shear) stress; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s, the softness A; constraint: positive.
    exponent : float or array_like
        1, Glen's n, about 3; constraint: at least 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1/s], effective strain rate. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Glen (1955): secondary creep of isotropic ice at stresses of about 0.05-1 MPa.
    Equation: ``rate_factor * stress**exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import glen_strain_rate
    >>> float(glen_strain_rate(100000.0, 2.4e-24, 3))
    2.4e-09
    """
    stress = quantity("stress", stress, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    return rate_factor * stress**exponent


@model
def glen_stress(strain_rate, rate_factor, exponent):
    """Stress needed to deform ice at a given strain rate.

    Inverting Glen's law turns a measured stretching rate into the stress that drives it.

    Parameters
    ----------
    strain_rate : float or array_like
        1/s, effective strain rate; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Glen (1955), inverted; secondary creep of isotropic ice.
    Equation: ``(strain_rate / rate_factor) ** (1 / exponent)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import glen_stress
    >>> round(float(glen_stress(2.4e-09, 2.4e-24, 3)), 0)
    100000.0
    """
    strain_rate = quantity("strain_rate", strain_rate, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    return (strain_rate / rate_factor) ** (1 / exponent)


@model
def glen_rate_factor(temperature, reference_factor, activation_energy, reference_temperature):
    """Arrhenius temperature dependence of Glen's rate factor.

    Ice near its melting point is softer by more than an order of magnitude than ice at -30 degC.

    Parameters
    ----------
    temperature : float or array_like
        K, ice temperature (relative to the pressure-melting point); constraint: positive.
    reference_factor : float or array_like
        Pa^-n/s, rate factor at reference_temperature; constraint: positive.
    activation_energy : float or array_like
        J/mol, creep activation energy; constraint: nonnegative.
    reference_temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa^-n/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    An Arrhenius law with one activation energy, about 60 kJ/mol below -10 degC and 115 kJ/mol above (Cuffey and
    Paterson 2010), for temperate-ice water contents neglected.
    Equation: ``reference_factor * np.exp(-activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import glen_rate_factor
    >>> round(float(glen_rate_factor(263.15, 2.4e-24, 60000.0, 273.15)), 30)
    8.79429e-25
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    reference_factor = quantity("reference_factor", reference_factor, lower=0, strict_lower=True)
    activation_energy = quantity("activation_energy", activation_energy, lower=0)
    reference_temperature = quantity("reference_temperature", reference_temperature, lower=0, strict_lower=True)
    return reference_factor * np.exp(-activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))


@model
def effective_viscosity(stress, rate_factor, exponent):
    """Effective viscosity of ice under Glen's law.

    Ice under a large stress behaves as a much less viscous fluid than ice under a small one.

    Parameters
    ----------
    stress : float or array_like
        Pa, effective stress; constraint: positive.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Glen (1955), written as a stress-dependent Newtonian viscosity; with n = 1 it is the viscosity 1 / (2 A).
    Equation: ``1 / (2 * rate_factor * stress ** (exponent - 1))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import effective_viscosity
    >>> round(float(effective_viscosity(100000.0, 2.4e-24, 3)), 0)
    20833333333333.0
    """
    stress = quantity("stress", stress, lower=0, strict_lower=True)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    return 1 / (2 * rate_factor * stress ** (exponent - 1))


@model
def glen_velocity_profile(height, thickness, surface_slope, rate_factor, exponent, density, gravity):
    """Deformation velocity through a glacier obeying Glen's law.

    Most of the shearing in a frozen-bed glacier happens near the bed, so the velocity is nearly uniform in its upper
    half.

    Parameters
    ----------
    height : float or array_like
        m, above the bed; constraint: between 0 and thickness.
    thickness : float or array_like
        m; constraint: positive.
    surface_slope : float or array_like
        1, tangent of the surface slope; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], downslope, zero at a frozen bed. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1952) laminar flow with Glen (1955) rheology: a wide slab, uniform rate factor and no sliding.
    Equation: ``2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * (thickness ** (exponent + 1)
    - (thickness - height) ** (exponent + 1))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import glen_velocity_profile
    >>> round(float(glen_velocity_profile(50, 100, 0.1, 1e-24, 3, 1000, 10)), 13)
    4.6875e-08
    """
    height = quantity("height", height, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(height <= thickness, "height exceeds the ice thickness")
    shape = thickness ** (exponent + 1) - (thickness - height) ** (exponent + 1)
    return 2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * shape


@model
def glen_surface_velocity(thickness, surface_slope, rate_factor, exponent, density, gravity):
    """Surface velocity of a glacier deforming by Glen's law.

    The surface speed rises with the fourth power of thickness and the cube of slope, which is why thick glaciers flow
    so much faster than thin ones.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: positive.
    surface_slope : float or array_like
        1, tangent of the surface slope; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1952) laminar flow with Glen (1955) rheology, no sliding.
    Equation: ``2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import glen_surface_velocity
    >>> round(float(glen_surface_velocity(100, 0.1, 1e-24, 3, 1000, 10)), 13)
    5e-08
    """
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return 2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 1)


@model
def nye_crevasse_depth(strain_rate, rate_factor, exponent, density, gravity):
    """Depth of a dry crevasse in stretching ice.

    A crevasse can only stay open where the tensile stress from stretching exceeds the closing pressure of the ice,
    which caps dry crevasses at a few tens of metres.

    Parameters
    ----------
    strain_rate : float or array_like
        1/s, longitudinal stretching rate; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1955): a field of closely spaced crevasses in uniform ice of density independent of depth; fracture toughness is
    neglected.
    Equation: ``2 * (strain_rate / rate_factor) ** (1 / exponent) / (density * gravity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import nye_crevasse_depth
    >>> round(float(nye_crevasse_depth(1e-09, 1e-24, 3, 1000, 10)), 4)
    20.0
    """
    strain_rate = quantity("strain_rate", strain_rate, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return 2 * (strain_rate / rate_factor) ** (1 / exponent) / (density * gravity)


@model
def water_filled_crevasse_depth(strain_rate, rate_factor, exponent, water_depth, ice_density, water_density, gravity):
    """Depth of a crevasse partly filled with meltwater.

    Water is denser than ice, so a water-filled crevasse is pushed deeper than a dry one and can cut through a whole
    ice shelf.

    Parameters
    ----------
    strain_rate : float or array_like
        1/s, longitudinal stretching rate; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    water_depth : float or array_like
        m, height of water standing in the crevasse; constraint: nonnegative.
    ice_density : float or array_like
        kg/m^3; constraint: positive.
    water_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1955) with the water pressure added (Benn et al. 2007); the result is meaningful only while it is less than the
    ice thickness and at least the water depth.
    Equation: ``2 * (strain_rate / rate_factor) ** (1 / exponent) / (ice_density * gravity) + water_density / ice_density *
    water_depth``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import water_filled_crevasse_depth
    >>> round(float(water_filled_crevasse_depth(1e-09, 1e-24, 3, 10, 900, 1000, 10)), 4)
    33.3333
    """
    strain_rate = quantity("strain_rate", strain_rate, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    water_depth = quantity("water_depth", water_depth, lower=0)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    dry = 2 * (strain_rate / rate_factor) ** (1 / exponent) / (ice_density * gravity)
    return dry + water_density / ice_density * water_depth
