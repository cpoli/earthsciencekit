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
    "teos10_density",
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

    Examples
    --------
    >>> from earthsciencekit.oceanography import hydrostatic_pressure
    >>> float(hydrostatic_pressure(1000, 100, 10))
    1000000.0
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

    Examples
    --------
    >>> from earthsciencekit.oceanography import linear_density
    >>> round(float(linear_density(1000, 1, 0, 0.0002, 0.0008)), 3)
    999.8
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

    Examples
    --------
    >>> from earthsciencekit.oceanography import simplified_seawater_density
    >>> round(float(simplified_seawater_density(293.15, 36, 0)), 2)
    1024.59
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

    Examples
    --------
    >>> from earthsciencekit.oceanography import buoyancy_frequency_squared
    >>> float(buoyancy_frequency_squared(1000, -0.01, 10))
    0.0001
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

    Examples
    --------
    >>> from earthsciencekit.oceanography import abyssal_temperature
    >>> float(abyssal_temperature(0.6931471805599453, 1, 4, 1.3862943611198906, 1, 1))
    2.0
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


# Roquet et al. (2015) 75-term specific volume [m^3/kg]: {(i, j, k): v_ijk} multiplies ys**i * xs**j * z**k, as in GSW-C.
_TEOS10_SPECIFIC_VOLUME = {
    (0, 0, 0): 1.0769995862e-03,
    (0, 0, 1): -6.0799143809e-05,
    (0, 0, 2): 9.9856169219e-06,
    (0, 0, 3): -1.1309361437e-06,
    (0, 0, 4): 1.0531153080e-07,
    (0, 0, 5): -1.2647261286e-08,
    (0, 0, 6): 1.9613503930e-09,
    (0, 1, 0): -3.1038981976e-04,
    (0, 1, 1): 2.4262468747e-05,
    (0, 1, 2): -5.8484432984e-07,
    (0, 1, 3): 3.6310188515e-07,
    (0, 1, 4): -1.1147125423e-07,
    (0, 2, 0): 6.6928067038e-04,
    (0, 2, 1): -3.4792460974e-05,
    (0, 2, 2): -4.8122251597e-06,
    (0, 2, 3): 1.6746303780e-08,
    (0, 3, 0): -8.5047933937e-04,
    (0, 3, 1): 3.7470777305e-05,
    (0, 3, 2): 4.9263106998e-06,
    (0, 4, 0): 5.8086069943e-04,
    (0, 4, 1): -1.7322218612e-05,
    (0, 4, 2): -1.7811974727e-06,
    (0, 5, 0): -2.1092370507e-04,
    (0, 5, 1): 3.0927427253e-06,
    (0, 6, 0): 3.1932457305e-05,
    (1, 0, 0): -1.5649734675e-05,
    (1, 0, 1): 1.8505765429e-05,
    (1, 0, 2): -1.1736386731e-06,
    (1, 0, 3): -3.6527006553e-07,
    (1, 0, 4): 3.1454099902e-07,
    (1, 1, 0): 3.5009599764e-05,
    (1, 1, 1): -9.5677088156e-06,
    (1, 1, 2): -5.5699154557e-06,
    (1, 1, 3): -2.7295696237e-07,
    (1, 2, 0): -4.3592678561e-05,
    (1, 2, 1): 1.1100834765e-05,
    (1, 2, 2): 5.4620748834e-06,
    (1, 3, 0): 3.4532461828e-05,
    (1, 3, 1): -9.8447117844e-06,
    (1, 3, 2): -1.3544185627e-06,
    (1, 4, 0): -1.1959409788e-05,
    (1, 4, 1): 2.5909225260e-06,
    (1, 5, 0): 1.3864594581e-06,
    (2, 0, 0): 2.7762106484e-05,
    (2, 0, 1): -1.1716606853e-05,
    (2, 0, 2): 2.1305028740e-06,
    (2, 0, 3): 2.8695905159e-07,
    (2, 1, 0): -3.7435842344e-05,
    (2, 1, 1): -2.3678308361e-07,
    (2, 1, 2): 3.9137387080e-07,
    (2, 2, 0): 3.5907822760e-05,
    (2, 2, 1): 2.9283346295e-06,
    (2, 2, 2): -6.5731104067e-07,
    (2, 3, 0): -1.8698584187e-05,
    (2, 3, 1): -4.8826139200e-07,
    (2, 4, 0): 3.8595339244e-06,
    (3, 0, 0): -1.6521159259e-05,
    (3, 0, 1): 7.9279656173e-06,
    (3, 0, 2): -4.6132540037e-07,
    (3, 1, 0): 2.4141479483e-05,
    (3, 1, 1): -3.4558773655e-06,
    (3, 1, 2): 7.7618888092e-09,
    (3, 2, 0): -1.4353633048e-05,
    (3, 2, 1): 3.1655306078e-07,
    (3, 3, 0): 2.2863324556e-06,
    (4, 0, 0): 6.9111322702e-06,
    (4, 0, 1): -3.4102187482e-06,
    (4, 0, 2): -6.3352916514e-08,
    (4, 1, 0): -8.7595873154e-06,
    (4, 1, 1): 1.2956717783e-06,
    (4, 2, 0): 4.3703680598e-06,
    (5, 0, 0): -8.0539615540e-07,
    (5, 0, 1): 5.0736766814e-07,
    (5, 1, 0): -3.3052758900e-07,
    (6, 0, 0): 2.0543094268e-07,
}


@model
def teos10_density(conservative_temperature, absolute_salinity, sea_pressure):
    """TEOS-10 in-situ density from the 75-term polynomial.

    One equation of state for all oceanographic work fixes density, and every buoyancy-driven current, to within the
    accuracy of the measurements behind the TEOS-10 Gibbs function.

    Parameters
    ----------
    conservative_temperature : float or array_like
        K; constraint: positive.
    absolute_salinity : float or array_like
        g/kg; constraint: nonnegative.
    sea_pressure : float or array_like
        Pa, absolute pressure minus 101325 Pa; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The 75-term specific-volume polynomial of Roquet et al. (2015) used by the Gibbs-SeaWater (GSW) toolbox for
    ``gsw_rho``; within the oceanographic range (about 0-42 g/kg, -2-40 degC and 0-80 MPa) it matches the full TEOS-10
    Gibbs function to well below the measurement error. Outside that range it is an extrapolation.
    Equation: ``1 / sum(v_ijk * ys**i * xs**j * z**k)`` with ``ys = (conservative_temperature - 273.15) / 40``,
    ``xs = np.sqrt(0.0248826675584615 * absolute_salinity + 0.5971840214030754)`` and ``z = sea_pressure / 1e8``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.oceanography import teos10_density
    >>> round(float(teos10_density(283.15, 30, 10000000.0)), 2)
    1027.45
    """
    conservative_temperature = quantity("conservative_temperature", conservative_temperature, lower=0, strict_lower=True)
    absolute_salinity = quantity("absolute_salinity", absolute_salinity, lower=0)
    sea_pressure = quantity("sea_pressure", sea_pressure, lower=0)
    ys = (conservative_temperature - 273.15) * 0.025
    xs = np.sqrt(0.0248826675584615 * absolute_salinity + 5.971840214030754e-1)
    z = sea_pressure * 1e-8
    specific_volume = sum(v * ys**i * xs**j * z**k for (i, j, k), v in _TEOS10_SPECIFIC_VOLUME.items())
    return 1 / specific_volume
