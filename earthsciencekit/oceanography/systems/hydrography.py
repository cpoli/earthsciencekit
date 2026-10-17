"""Hydrography: hydrostatic pressure, equations of state, stratification and the abyssal temperature profile.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "hydrostatic_pressure",
    "linear_density",
    "buoyancy_frequency_squared",
    "abyssal_temperature",
    "simplified_seawater_density",
]


@model
def hydrostatic_pressure(density, depth, gravity):
    """Ocean hydrostatic pressure.

    Water weight connects pressure and depth.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: positive.
    depth : float or array_like
        m; constraint: nonnegative.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Gauge pressure with constant density; atmospheric pressure is excluded.
    Equation: ``density * gravity * depth``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    density = quantity("density", density, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return density * gravity * depth


@model
def linear_density(reference_density, temperature_anomaly, salinity_anomaly, thermal_expansion, haline_contraction):
    """Linear seawater equation of state.

    Temperature and salinity both change buoyancy.

    Parameters
    ----------
    reference_density : float or array_like
        kg/m^3; constraint: positive.
    temperature_anomaly : float or array_like
        K; constraint: real.
    salinity_anomaly : float or array_like
        g/kg; constraint: real.
    thermal_expansion : float or array_like
        1/K; constraint: nonnegative.
    haline_contraction : float or array_like
        kg/g; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Local linearization around a supplied reference; not TEOS-10 and no pressure dependence.
    Equation: ``reference_density * (1 - thermal_expansion * temperature_anomaly + haline_contraction * salinity_anomaly)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    reference_density = quantity("reference_density", reference_density, lower=0, strict_lower=True)
    temperature_anomaly = quantity("temperature_anomaly", temperature_anomaly)
    salinity_anomaly = quantity("salinity_anomaly", salinity_anomaly)
    thermal_expansion = quantity("thermal_expansion", thermal_expansion, lower=0)
    haline_contraction = quantity("haline_contraction", haline_contraction, lower=0)
    require(
        1 - thermal_expansion * temperature_anomaly + haline_contraction * salinity_anomaly > 0, "linearization implies nonpositive density"
    )
    return reference_density * (1 - thermal_expansion * temperature_anomaly + haline_contraction * salinity_anomaly)


@model
def simplified_seawater_density(conservative_temperature, absolute_salinity, depth):
    """Simplified nonlinear seawater density fitted to TEOS-10.

    Cabbeling and thermobaricity make seawater density nonlinear in temperature, salinity and pressure.

    Parameters
    ----------
    conservative_temperature : float or array_like
        K; constraint: positive.
    absolute_salinity : float or array_like
        g/kg; constraint: nonnegative.
    depth : float or array_like
        m; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The S-EOS of Roquet et al. (2015) with NEMO's default coefficients: ``a0 = 0.16550``, ``b0 = 0.76554``,
    ``lambda1 = 5.9520e-2``, ``lambda2 = 7.4914e-4``, ``mu1 = 1.4970e-4``, ``mu2 = 1.1090e-5``, ``nu = 2.4341e-3``
    and ``rho0 = 1026``, with anomalies ``t = conservative_temperature - 283.15`` and ``s = absolute_salinity - 35``.
    It reproduces the TEOS-10 density gradients over typical ocean conditions; absolute density can be off by about
    1 kg/m^3 and it is not TEOS-10 itself.
    Equation: ``1026 - a0 * (1 + lambda1 * t / 2 + mu1 * depth) * t + b0 * (1 - lambda2 * s / 2 - mu2 * depth) * s - nu * t * s``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    conservative_temperature = quantity("conservative_temperature", conservative_temperature, lower=0, strict_lower=True)
    absolute_salinity = quantity("absolute_salinity", absolute_salinity, lower=0)
    depth = quantity("depth", depth, lower=0)
    t = conservative_temperature - 283.15
    s = absolute_salinity - 35
    thermal = 0.16550 * (1 + 0.5 * 5.9520e-2 * t + 1.4970e-4 * depth) * t
    haline = 0.76554 * (1 - 0.5 * 7.4914e-4 * s - 1.1090e-5 * depth) * s
    return 1026 - thermal + haline - 2.4341e-3 * t * s


@model
def buoyancy_frequency_squared(density, upward_density_gradient, gravity):
    """Ocean density stratification.

    Vertical density gradients control internal oscillations and convection.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: positive.
    upward_density_gradient : float or array_like
        kg/m^4; constraint: real.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Coordinate increases upward; negative output denotes unstable stratification.
    Equation: ``-gravity * upward_density_gradient / density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    density = quantity("density", density, lower=0, strict_lower=True)
    upward_density_gradient = quantity("upward_density_gradient", upward_density_gradient)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return -gravity * upward_density_gradient / density


@model
def abyssal_temperature(height, bottom_temperature, top_temperature, column_height, diffusivity, upwelling):
    """Munk's advective-diffusive abyssal temperature profile.

    Slow upwelling balanced by downward mixing of heat sets an exponential thermocline.

    Parameters
    ----------
    height : float or array_like
        m; constraint: 0 <= height <= column_height. Height above the bottom of the column.
    bottom_temperature : float or array_like
        K; constraint: positive.
    top_temperature : float or array_like
        K; constraint: positive.
    column_height : float or array_like
        m; constraint: positive.
    diffusivity : float or array_like
        m^2/s; constraint: positive. Vertical eddy diffusivity.
    upwelling : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Steady one-dimensional balance ``upwelling * dT/dz = diffusivity * d2T/dz2`` with fixed end temperatures;
    the scale height is ``h = diffusivity / upwelling`` (about 1 km in Munk's recipe).
    Equation: ``bottom_temperature + (top_temperature - bottom_temperature)
    * (np.exp((height - column_height) / h) - np.exp(-column_height / h)) / (1 - np.exp(-column_height / h))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    height = quantity("height", height, lower=0)
    bottom_temperature = quantity("bottom_temperature", bottom_temperature, lower=0, strict_lower=True)
    top_temperature = quantity("top_temperature", top_temperature, lower=0, strict_lower=True)
    column_height = quantity("column_height", column_height, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    upwelling = quantity("upwelling", upwelling, lower=0, strict_lower=True)
    require(height <= column_height, "height must lie inside the column")
    scale = diffusivity / upwelling
    shape = np.exp((height - column_height) / scale) * -np.expm1(-height / scale) / -np.expm1(-column_height / scale)
    return bottom_temperature + (top_temperature - bottom_temperature) * shape
