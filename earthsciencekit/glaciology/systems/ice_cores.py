"""Ice cores: Dansgaard's isotope thermometer, the Nye and Dansgaard-Johnsen time scales and the gas-age difference.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "dansgaard_delta",
    "dansgaard_temperature",
    "nye_age",
    "dansgaard_johnsen_downward_velocity",
    "dansgaard_johnsen_age",
    "gas_age_difference",
]


@model
def dansgaard_delta(temperature, slope, offset):
    """Oxygen-isotope ratio of precipitation from the air temperature.

    Air loses its heavy isotopes first as it cools and rains out, so the snow at colder sites is lighter.

    Parameters
    ----------
    temperature : float or array_like
        K, mean annual surface air temperature; constraint: positive.
    slope : float or array_like
        permil/K, about 0.69; constraint: positive.
    offset : float or array_like
        permil, the delta at 0 degC, about -13.6; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [permil], delta-18O relative to SMOW. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dansgaard (1964): a spatial relation for mid- and high-latitude coastal and polar stations; the temporal slope at
    one site can differ.
    Equation: ``slope * (temperature - 273.15) + offset``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0, strict_lower=True)
    offset = quantity("offset", offset)
    return slope * (temperature - 273.15) + offset


@model
def dansgaard_temperature(delta, slope, offset):
    """Air temperature recorded by the oxygen-isotope ratio of snow.

    Inverting Dansgaard's relation turns the isotope profile of an ice core into a temperature history.

    Parameters
    ----------
    delta : float or array_like
        permil, delta-18O relative to SMOW; constraint: real.
    slope : float or array_like
        permil/K; constraint: positive.
    offset : float or array_like
        permil; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dansgaard (1964), inverted; the spatial slope is an approximation of the temporal one.
    Equation: ``273.15 + (delta - offset) / slope``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    delta = quantity("delta", delta)
    slope = quantity("slope", slope, lower=0, strict_lower=True)
    offset = quantity("offset", offset)
    return 273.15 + (delta - offset) / slope


@model
def nye_age(depth, thickness, accumulation):
    """Age of ice at depth when layers thin uniformly.

    In a steady ice sheet each layer thins in proportion to its height above the bed, so age grows logarithmically with
    depth.

    Parameters
    ----------
    depth : float or array_like
        m, of ice equivalent below the surface; constraint: between 0 and thickness, exclusive of thickness.
    thickness : float or array_like
        m of ice; constraint: positive.
    accumulation : float or array_like
        m/s of ice; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1963): a steady divide with uniform vertical strain rate, no basal melting or sliding; the age diverges at the
    bed.
    Equation: ``thickness / accumulation * np.log(thickness / (thickness - depth))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    require(depth < thickness, "depth must be above the bed")
    return thickness / accumulation * np.log(thickness / (thickness - depth))


@model
def dansgaard_johnsen_downward_velocity(depth, thickness, kink_height, accumulation):
    """Vertical ice velocity in the Dansgaard-Johnsen flow model.

    Ice frozen to its bed shears only near the bottom, so the vertical velocity falls linearly down to a kink and then
    quadratically to zero at the bed.

    Parameters
    ----------
    depth : float or array_like
        m of ice below the surface; constraint: between 0 and thickness.
    thickness : float or array_like
        m of ice; constraint: positive.
    kink_height : float or array_like
        m above the bed; constraint: between 0 (exclusive) and thickness.
    accumulation : float or array_like
        m/s of ice; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], downward; a year's layer is this speed times one year thick. Scalar input returns a NumPy
        scalar; arrays broadcast.

    Notes
    -----
    Dansgaard and Johnsen (1969): steady state, constant horizontal strain rate above the kink and a linear decrease
    below it, no basal melting. Below, z = thickness - depth and h = kink_height.
    Equation: ``np.where(z >= h, accumulation * (2 * z - h) / (2 * thickness - h), accumulation * z**2 / (h * (2 * thickness - h)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    kink_height = quantity("kink_height", kink_height, lower=0, strict_lower=True)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    require(depth <= thickness, "depth exceeds the ice thickness")
    require(kink_height <= thickness, "kink height exceeds the ice thickness")
    z, h = thickness - depth, kink_height
    return np.where(z >= h, accumulation * (2 * z - h) / (2 * thickness - h), accumulation * z**2 / (h * (2 * thickness - h)))


@model
def dansgaard_johnsen_age(depth, thickness, kink_height, accumulation):
    """Age of ice at depth in the Dansgaard-Johnsen flow model.

    The first time scale for a deep ice core, Camp Century, integrated the thinning of annual layers down to the kink and
    below it.

    Parameters
    ----------
    depth : float or array_like
        m of ice below the surface; constraint: between 0 and thickness, exclusive of thickness.
    thickness : float or array_like
        m of ice; constraint: positive.
    kink_height : float or array_like
        m above the bed; constraint: between 0 (exclusive) and thickness.
    accumulation : float or array_like
        m/s of ice; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dansgaard and Johnsen (1969), the integral of the reciprocal of dansgaard_johnsen_downward_velocity; it tends to
    nye_age as the kink height goes to zero. Below, z = thickness - depth, h = kink_height, a = accumulation and s = 2 * thickness - h.
    Equation: ``s / (2 * a) * np.log(s / (2 * np.maximum(z, h) - h)) + h * s / a * (1 / np.minimum(z, h) - 1 / h)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    kink_height = quantity("kink_height", kink_height, lower=0, strict_lower=True)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    require(depth < thickness, "depth must be above the bed")
    require(kink_height <= thickness, "kink height exceeds the ice thickness")
    z, h = thickness - depth, kink_height
    span = 2 * thickness - h
    upper = span / (2 * accumulation) * np.log(span / (2 * np.maximum(z, h) - h))
    lower = h * span / accumulation * (1 / np.minimum(z, h) - 1 / h)
    return upper + lower


@model
def gas_age_difference(close_off_depth, mean_firn_density, ice_density, accumulation):
    """Age difference between the ice and the air trapped in it.

    Air stays connected to the atmosphere through the porous firn until the bubbles close, so the gas in a core is
    younger than the ice around it.

    Parameters
    ----------
    close_off_depth : float or array_like
        m, depth of bubble close-off; constraint: nonnegative.
    mean_firn_density : float or array_like
        kg/m^3, mean density of the firn above it; constraint: between 0 (exclusive) and ice_density.
    ice_density : float or array_like
        kg/m^3; constraint: positive.
    accumulation : float or array_like
        m/s of ice; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]; about 4 kyr at present-day Vostok. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Schwander and Stauffer (1984), as used for Vostok by Barnola et al. (1987): steady firn densification with no
    thinning above close-off; diffusion through the open pores, which makes the gas slightly older, is neglected.
    Equation: ``close_off_depth * mean_firn_density / (ice_density * accumulation)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    close_off_depth = quantity("close_off_depth", close_off_depth, lower=0)
    mean_firn_density = quantity("mean_firn_density", mean_firn_density, lower=0, strict_lower=True)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    require(mean_firn_density <= ice_density, "firn cannot be denser than ice")
    return close_off_depth * mean_firn_density / (ice_density * accumulation)
