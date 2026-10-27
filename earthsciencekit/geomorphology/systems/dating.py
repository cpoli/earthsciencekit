"""Dating landforms: cosmogenic nuclides, morphologic ages of scarps and luminescence.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "cosmogenic_production_at_depth",
    "cosmogenic_concentration",
    "steady_state_erosion_rate",
    "morphologic_age",
    "luminescence_age",
]


@model
def cosmogenic_production_at_depth(surface_production_rate, depth, density, attenuation_length):
    """Cosmogenic-nuclide production rate beneath the surface.

    Cosmic-ray neutrons are absorbed within the top two metres of rock, so a grain begins to accumulate nuclides only
    as erosion brings it close to the surface.

    Parameters
    ----------
    surface_production_rate : float or array_like
        atoms/(kg s), at the surface; constraint: nonnegative.
    depth : float or array_like
        m; constraint: nonnegative.
    density : float or array_like
        kg/m^3 of the overburden; constraint: positive.
    attenuation_length : float or array_like
        kg/m^2, about 1600 (160 g/cm^2) for spallation; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [atoms/(kg s)]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Lal (1991): neutron spallation only, on a flat surface; muons dominate below a few metres.
    Equation: ``surface_production_rate * np.exp(-density * depth / attenuation_length)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import cosmogenic_production_at_depth
    >>> round(float(cosmogenic_production_at_depth(10, 1, 2000, 1600)), 11)
    2.8650479686
    """
    surface_production_rate = quantity("surface_production_rate", surface_production_rate, lower=0)
    depth = quantity("depth", depth, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    attenuation_length = quantity("attenuation_length", attenuation_length, lower=0, strict_lower=True)
    return surface_production_rate * np.exp(-density * depth / attenuation_length)


@model
def cosmogenic_concentration(time, production_rate, decay_constant, erosion_rate, density, attenuation_length):
    """Concentration of a cosmogenic nuclide in an eroding surface.

    Nuclides build up until decay and erosion remove them as fast as cosmic rays make them, so a surface's
    concentration records either its exposure age or its erosion rate.

    Parameters
    ----------
    time : float or array_like
        s of exposure, starting from zero concentration; constraint: nonnegative.
    production_rate : float or array_like
        atoms/(kg s), at the surface; constraint: nonnegative.
    decay_constant : float or array_like
        1/s, ``LAMBDA_BE10`` for beryllium-10; constraint: nonnegative.
    erosion_rate : float or array_like
        m/s, steady; constraint: nonnegative.
    density : float or array_like
        kg/m^3; constraint: positive.
    attenuation_length : float or array_like
        kg/m^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [atoms/kg]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Lal (1991): steady erosion and production, spallation only, no inheritance; ``production_rate * time`` when
    neither decay nor erosion removes nuclides.
    Equation: ``production_rate / r * (1 - np.exp(-r * time))`` with
    ``r = decay_constant + density * erosion_rate / attenuation_length``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import cosmogenic_concentration
    >>> round(float(cosmogenic_concentration(0.1, 10, 0.05, 0, 2000, 1600)), 12)
    0.997504161464
    """
    time = quantity("time", time, lower=0)
    production_rate = quantity("production_rate", production_rate, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0)
    erosion_rate = quantity("erosion_rate", erosion_rate, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    attenuation_length = quantity("attenuation_length", attenuation_length, lower=0, strict_lower=True)
    rate = decay_constant + density * erosion_rate / attenuation_length
    removing = rate > 0
    exposure = np.where(removing, -np.expm1(-rate * time) / np.where(removing, rate, 1.0), time)
    return production_rate * exposure


@model
def steady_state_erosion_rate(concentration, production_rate, decay_constant, density, attenuation_length):
    """Erosion rate from a steady-state cosmogenic-nuclide concentration.

    A surface that has eroded steadily for long enough carries a concentration inversely proportional to its erosion
    rate, the basis of catchment-wide denudation rates from river sand.

    Parameters
    ----------
    concentration : float or array_like
        atoms/kg; constraint: positive, at most production_rate / decay_constant.
    production_rate : float or array_like
        atoms/(kg s); constraint: positive.
    decay_constant : float or array_like
        1/s; constraint: nonnegative.
    density : float or array_like
        kg/m^3; constraint: positive.
    attenuation_length : float or array_like
        kg/m^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Lal (1991), the inverse of :func:`cosmogenic_concentration` at long times; valid once the surface has been
    exposed for several times ``attenuation_length / (density * erosion_rate)``.
    Equation: ``(production_rate / concentration - decay_constant) * attenuation_length / density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import steady_state_erosion_rate
    >>> float(steady_state_erosion_rate(100, 10, 0.05, 2000, 1600))
    0.04
    """
    concentration = quantity("concentration", concentration, lower=0, strict_lower=True)
    production_rate = quantity("production_rate", production_rate, lower=0, strict_lower=True)
    decay_constant = quantity("decay_constant", decay_constant, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    attenuation_length = quantity("attenuation_length", attenuation_length, lower=0, strict_lower=True)
    require(decay_constant * concentration <= production_rate, "concentration exceeds the saturation value")
    return (production_rate / concentration - decay_constant) * attenuation_length / density


@model
def morphologic_age(half_offset, maximum_gradient, far_field_gradient):
    """Morphologic age of a fault scarp from its steepest gradient.

    Measuring a scarp's height and steepest slope gives the product of diffusivity and age, which dates the
    earthquake once the diffusivity is calibrated.

    Parameters
    ----------
    half_offset : float or array_like
        m; constraint: positive.
    maximum_gradient : float or array_like
        1; constraint: greater than far_field_gradient.
    far_field_gradient : float or array_like
        1; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2], the diffusivity times the age. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hanks et al. (1984), inverting :func:`scarp_maximum_gradient`: one rupture, linear diffusion and no base-level
    change at the foot.
    Equation: ``half_offset**2 / (np.pi * (maximum_gradient - far_field_gradient)**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import morphologic_age
    >>> round(float(morphologic_age(2, 0.6, 0.1)), 11)
    5.09295817894
    """
    half_offset = quantity("half_offset", half_offset, lower=0, strict_lower=True)
    maximum_gradient = quantity("maximum_gradient", maximum_gradient)
    far_field_gradient = quantity("far_field_gradient", far_field_gradient)
    require(maximum_gradient > far_field_gradient, "maximum_gradient must exceed far_field_gradient")
    return half_offset**2 / (np.pi * (maximum_gradient - far_field_gradient) ** 2)


@model
def luminescence_age(equivalent_dose, dose_rate):
    """Luminescence burial age of a sediment.

    Quartz and feldspar grains store the radiation dose received since sunlight last emptied their traps, which
    dates when a dune or flood deposit was buried.

    Parameters
    ----------
    equivalent_dose : float or array_like
        Gy, the laboratory dose that reproduces the natural signal; constraint: nonnegative.
    dose_rate : float or array_like
        Gy/s from the surrounding sediment and cosmic rays; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Aitken (1985): complete bleaching at deposition, a constant dose rate and no fading or saturation of the signal.
    Equation: ``equivalent_dose / dose_rate``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import luminescence_age
    >>> float(luminescence_age(50, 2e-11))
    2500000000000.0
    """
    equivalent_dose = quantity("equivalent_dose", equivalent_dose, lower=0)
    dose_rate = quantity("dose_rate", dose_rate, lower=0, strict_lower=True)
    return equivalent_dose / dose_rate
