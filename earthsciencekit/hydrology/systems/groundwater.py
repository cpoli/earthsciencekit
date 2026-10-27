"""Groundwater: hydraulic head, Darcy flow, Dupuit's unconfined aquifer and the Theis and Cooper-Jacob well solutions.

SI units unless stated otherwise.
"""

import numpy as np
from scipy.special import exp1

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "hydraulic_head",
    "darcy_flux",
    "seepage_velocity",
    "dupuit_head",
    "dupuit_discharge",
    "well_function",
    "theis_drawdown",
    "cooper_jacob_drawdown",
    "cooper_jacob_transmissivity",
]


@model
def hydraulic_head(elevation, pressure, density, gravity):
    """Hydraulic head from elevation and pore pressure.

    Groundwater flows from high to low head, the sum of the elevation and pressure heads, not from high to low pressure.

    Parameters
    ----------
    elevation : float or array_like
        m, above a datum; constraint: real.
    pressure : float or array_like
        Pa, gauge pore pressure; constraint: real.
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
    Slow flow with negligible velocity head and uniform fluid density.
    Equation: ``elevation + pressure / (density * gravity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import hydraulic_head
    >>> float(hydraulic_head(10, 9810, 1000, 9.81))
    11.0
    """
    elevation = quantity("elevation", elevation)
    pressure = quantity("pressure", pressure)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return elevation + pressure / (density * gravity)


@model
def darcy_flux(hydraulic_conductivity, head_gradient):
    """Darcy's specific discharge through a porous medium.

    Darcy's filters at Dijon showed that the flow through sand is proportional to the drop in head per unit length.

    Parameters
    ----------
    hydraulic_conductivity : float or array_like
        m/s; constraint: nonnegative.
    head_gradient : float or array_like
        1, change of head per unit length along the flow direction; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], volume per unit cross-sectional area and time, positive down the gradient. Scalar input
        returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Darcy (1856): laminar flow with a pore Reynolds number below about 1-10 in a saturated isotropic medium.
    Equation: ``-hydraulic_conductivity * head_gradient``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import darcy_flux
    >>> round(float(darcy_flux(0.0001, -0.01)), 11)
    1e-06
    """
    hydraulic_conductivity = quantity("hydraulic_conductivity", hydraulic_conductivity, lower=0)
    head_gradient = quantity("head_gradient", head_gradient)
    return -hydraulic_conductivity * head_gradient


@model
def seepage_velocity(darcy_flux, porosity):
    """Mean velocity of water in the pores.

    Water moves only through the connected pores, so a tracer travels faster than the Darcy flux suggests.

    Parameters
    ----------
    darcy_flux : float or array_like
        m/s; constraint: real.
    porosity : float or array_like
        1, effective porosity; constraint: between 0 (exclusive) and 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Average linear velocity; dispersion spreads tracers about it.
    Equation: ``darcy_flux / porosity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import seepage_velocity
    >>> float(seepage_velocity(1e-06, 0.25))
    4e-06
    """
    darcy_flux = quantity("darcy_flux", darcy_flux)
    porosity = quantity("porosity", porosity, lower=0, upper=1, strict_lower=True)
    return darcy_flux / porosity


@model
def dupuit_head(distance, upstream_head, downstream_head, length):
    """Water-table height in an unconfined aquifer between two rivers.

    Treating the flow as horizontal turns the unconfined aquifer into a one-dimensional problem with a parabolic water
    table.

    Parameters
    ----------
    distance : float or array_like
        m, from the upstream boundary; constraint: between 0 and length.
    upstream_head : float or array_like
        m, saturated thickness above the impermeable base; constraint: nonnegative.
    downstream_head : float or array_like
        m; constraint: nonnegative.
    length : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dupuit (1863): no recharge, a horizontal base, and slopes small enough that the head is uniform on each vertical.
    Equation: ``np.sqrt(upstream_head**2 - (upstream_head**2 - downstream_head**2) * distance / length)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import dupuit_head
    >>> round(float(dupuit_head(50, 5, 3, 100)), 5)
    4.12311
    """
    distance = quantity("distance", distance, lower=0)
    upstream_head = quantity("upstream_head", upstream_head, lower=0)
    downstream_head = quantity("downstream_head", downstream_head, lower=0)
    length = quantity("length", length, lower=0, strict_lower=True)
    require(distance <= length, "distance exceeds the aquifer length")
    return np.sqrt(upstream_head**2 - (upstream_head**2 - downstream_head**2) * distance / length)


@model
def dupuit_discharge(hydraulic_conductivity, upstream_head, downstream_head, length):
    """Discharge per unit width through an unconfined aquifer.

    In Dupuit's approximation the flow depends on the difference of the squared heads, not of the heads.

    Parameters
    ----------
    hydraulic_conductivity : float or array_like
        m/s; constraint: nonnegative.
    upstream_head : float or array_like
        m, saturated thickness; constraint: nonnegative.
    downstream_head : float or array_like
        m; constraint: nonnegative.
    length : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s], positive from upstream to downstream. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dupuit (1863), the Dupuit-Forchheimer discharge formula, exact for the discharge even where the head profile is not.
    Equation: ``hydraulic_conductivity * (upstream_head**2 - downstream_head**2) / (2 * length)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import dupuit_discharge
    >>> float(dupuit_discharge(2, 5, 3, 4))
    4.0
    """
    hydraulic_conductivity = quantity("hydraulic_conductivity", hydraulic_conductivity, lower=0)
    upstream_head = quantity("upstream_head", upstream_head, lower=0)
    downstream_head = quantity("downstream_head", downstream_head, lower=0)
    length = quantity("length", length, lower=0, strict_lower=True)
    return hydraulic_conductivity * (upstream_head**2 - downstream_head**2) / (2 * length)


@model
def well_function(u):
    """Theis well function W(u).

    The drawdown around a pumped well is a universal curve of one combination of distance, time and aquifer properties.

    Parameters
    ----------
    u : float or array_like
        1, ``radius**2 * storativity / (4 * transmissivity * time)``; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The exponential integral E1, ``-0.5772 - ln(u) + u - ...`` for small u.
    Equation: ``scipy.special.exp1(u)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import well_function
    >>> round(float(well_function(1)), 6)
    0.219384
    """
    u = quantity("u", u, lower=0, strict_lower=True)
    return exp1(u)


@model
def theis_drawdown(pumping_rate, transmissivity, storativity, radius, time):
    """Drawdown around a well pumping a confined aquifer.

    Theis took the heat-conduction solution for a line sink and found that the cone of depression keeps growing as the
    well draws water from storage.

    Parameters
    ----------
    pumping_rate : float or array_like
        m^3/s; constraint: real.
    transmissivity : float or array_like
        m^2/s; constraint: positive.
    storativity : float or array_like
        1; constraint: positive.
    radius : float or array_like
        m, from the well axis; constraint: positive.
    time : float or array_like
        s, since pumping began; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], positive downward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Theis (1935): constant rate, a fully penetrating well of negligible radius in an infinite, homogeneous, isotropic
    confined aquifer with water released instantly from storage.
    Equation: ``pumping_rate / (4 * np.pi * transmissivity) * scipy.special.exp1(radius**2 * storativity / (4 *
    transmissivity * time))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import theis_drawdown
    >>> round(float(theis_drawdown(12.566370614359172, 1, 4, 1, 1)), 6)
    0.219384
    """
    pumping_rate = quantity("pumping_rate", pumping_rate)
    transmissivity = quantity("transmissivity", transmissivity, lower=0, strict_lower=True)
    storativity = quantity("storativity", storativity, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0, strict_lower=True)
    return pumping_rate / (4 * np.pi * transmissivity) * exp1(radius**2 * storativity / (4 * transmissivity * time))


@model
def cooper_jacob_drawdown(pumping_rate, transmissivity, storativity, radius, time):
    """Cooper and Jacob's late-time approximation to the Theis drawdown.

    Late in a pumping test the drawdown grows linearly with the logarithm of time, so a straight line on a semilog plot
    gives the aquifer properties.

    Parameters
    ----------
    pumping_rate : float or array_like
        m^3/s; constraint: real.
    transmissivity : float or array_like
        m^2/s; constraint: positive.
    storativity : float or array_like
        1; constraint: positive.
    radius : float or array_like
        m; constraint: positive.
    time : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Cooper and Jacob (1946): the first two terms of the Theis series, within 1% for ``u < 0.01`` and meaningless (it goes
    negative) for u above about 0.56.
    Equation: ``pumping_rate / (4 * np.pi * transmissivity) * (-np.euler_gamma - np.log(radius**2 * storativity / (4 *
    transmissivity * time)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import cooper_jacob_drawdown
    >>> round(float(cooper_jacob_drawdown(12.566370614359172, 1, 4, 1, 1)), 6)
    -0.577216
    """
    pumping_rate = quantity("pumping_rate", pumping_rate)
    transmissivity = quantity("transmissivity", transmissivity, lower=0, strict_lower=True)
    storativity = quantity("storativity", storativity, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0, strict_lower=True)
    u = radius**2 * storativity / (4 * transmissivity * time)
    return pumping_rate / (4 * np.pi * transmissivity) * (-np.euler_gamma - np.log(u))


@model
def cooper_jacob_transmissivity(pumping_rate, drawdown_per_log_cycle):
    """Transmissivity from the slope of a Cooper-Jacob straight line.

    The drawdown gained per tenfold increase in time depends only on the pumping rate and the transmissivity.

    Parameters
    ----------
    pumping_rate : float or array_like
        m^3/s; constraint: positive.
    drawdown_per_log_cycle : float or array_like
        m, drawdown increase per factor of ten in time; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Cooper and Jacob (1946), from late-time data where ``u < 0.01``; the storativity follows from the zero-drawdown
    intercept ``t0`` as ``2.25 * transmissivity * t0 / radius**2``.
    Equation: ``np.log(10) * pumping_rate / (4 * np.pi * drawdown_per_log_cycle)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import cooper_jacob_transmissivity
    >>> float(cooper_jacob_transmissivity(12.566370614359172, 2.302585092994046))
    1.0
    """
    pumping_rate = quantity("pumping_rate", pumping_rate, lower=0, strict_lower=True)
    drawdown_per_log_cycle = quantity("drawdown_per_log_cycle", drawdown_per_log_cycle, lower=0, strict_lower=True)
    return np.log(10) * pumping_rate / (4 * np.pi * drawdown_per_log_cycle)
