"""Rivers: stream power, steady and transient long profiles, network scaling and hydraulic geometry.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.results import TimeSeries

__all__ = [
    "stream_power_erosion",
    "steady_state_slope",
    "knickpoint_celerity",
    "flint_slope",
    "chi_coordinate",
    "hack_length",
    "strahler_order",
    "downstream_hydraulic_geometry",
    "effective_discharge",
    "simulate_river_profile",
]


@model
def stream_power_erosion(drainage_area, slope, erodibility, area_exponent, slope_exponent):
    """Bedrock incision rate of the stream-power law.

    A river cuts its bed faster where it carries more water down a steeper slope, the rule behind most models of
    mountain landscapes.

    Parameters
    ----------
    drainage_area : float or array_like
        m^2, the area upstream, a proxy for discharge; constraint: positive.
    slope : float or array_like
        1, tangent of the channel slope; constraint: nonnegative.
    erodibility : float or array_like
        m^(1 - 2 area_exponent)/s for a slope exponent of 1; constraint: nonnegative.
    area_exponent : float or array_like
        1, ``m``, about 0.5; constraint: nonnegative.
    slope_exponent : float or array_like
        1, ``n``, about 1; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s] of bed lowering. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Howard and Kerby (1983): detachment-limited bedrock channels, with discharge and width folded into powers of the
    area; no threshold, sediment cover or tools effect.
    Equation: ``erodibility * drainage_area**area_exponent * slope**slope_exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import stream_power_erosion
    >>> round(float(stream_power_erosion(1000000.0, 0.01, 1e-10, 0.5, 1)), 20)
    1e-09
    """
    drainage_area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0)
    erodibility = quantity("erodibility", erodibility, lower=0)
    area_exponent = quantity("area_exponent", area_exponent, lower=0)
    slope_exponent = quantity("slope_exponent", slope_exponent, lower=0, strict_lower=True)
    return erodibility * drainage_area**area_exponent * slope**slope_exponent


@model
def steady_state_slope(uplift_rate, drainage_area, erodibility, area_exponent, slope_exponent):
    """Channel slope at which stream-power incision balances uplift.

    A graded bedrock river steepens until it cuts down as fast as the rock rises, so its slope falls with drainage
    area along a concave profile.

    Parameters
    ----------
    uplift_rate : float or array_like
        m/s relative to base level; constraint: nonnegative.
    drainage_area : float or array_like
        m^2; constraint: positive.
    erodibility : float or array_like
        m^(1 - 2 area_exponent)/s for a slope exponent of 1; constraint: positive.
    area_exponent : float or array_like
        1; constraint: nonnegative.
    slope_exponent : float or array_like
        1; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Steady state of :func:`stream_power_erosion` under uniform uplift; the slope follows Flint's law with concavity
    ``area_exponent / slope_exponent``.
    Equation: ``(uplift_rate / (erodibility * drainage_area**area_exponent))**(1 / slope_exponent)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import steady_state_slope
    >>> float(steady_state_slope(1e-09, 1000000.0, 1e-10, 0.5, 1))
    0.01
    """
    uplift_rate = quantity("uplift_rate", uplift_rate, lower=0)
    drainage_area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    erodibility = quantity("erodibility", erodibility, lower=0, strict_lower=True)
    area_exponent = quantity("area_exponent", area_exponent, lower=0)
    slope_exponent = quantity("slope_exponent", slope_exponent, lower=0, strict_lower=True)
    return (uplift_rate / (erodibility * drainage_area**area_exponent)) ** (1 / slope_exponent)


@model
def knickpoint_celerity(drainage_area, slope, erodibility, area_exponent, slope_exponent):
    """Upstream migration speed of a knickpoint under stream power.

    A drop in base level travels up a river as a step that moves fastest on the large trunk streams and slows as it
    climbs into the headwaters.

    Parameters
    ----------
    drainage_area : float or array_like
        m^2; constraint: positive.
    slope : float or array_like
        1, the channel slope at the knickpoint; constraint: positive.
    erodibility : float or array_like
        m^(1 - 2 area_exponent)/s for a slope exponent of 1; constraint: nonnegative.
    area_exponent : float or array_like
        1; constraint: nonnegative.
    slope_exponent : float or array_like
        1; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], horizontal and upstream. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Whipple and Tucker (1999), the characteristic speed of :func:`stream_power_erosion`; independent of the slope
    only when ``slope_exponent = 1``, where knickpoints keep their shape.
    Equation: ``erodibility * drainage_area**area_exponent * slope**(slope_exponent - 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import knickpoint_celerity
    >>> round(float(knickpoint_celerity(1000000.0, 0.01, 1e-10, 0.5, 1)), 18)
    1e-07
    """
    drainage_area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0, strict_lower=True)
    erodibility = quantity("erodibility", erodibility, lower=0)
    area_exponent = quantity("area_exponent", area_exponent, lower=0)
    slope_exponent = quantity("slope_exponent", slope_exponent, lower=0, strict_lower=True)
    return erodibility * drainage_area**area_exponent * slope ** (slope_exponent - 1)


@model
def flint_slope(drainage_area, steepness_index, concavity):
    """Flint's power law of channel slope against drainage area.

    Along almost every river the slope falls as a power of the area drained, so a log-log plot of slope against area
    is a straight line whose slope is the concavity.

    Parameters
    ----------
    drainage_area : float or array_like
        m^2; constraint: positive.
    steepness_index : float or array_like
        m^(2 concavity), the channel steepness; constraint: nonnegative.
    concavity : float or array_like
        1, about 0.4-0.6; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Flint (1974): below the colluvial-fluvial transition at about 1 km^2, in rivers without knickpoints or changes of
    lithology.
    Equation: ``steepness_index * drainage_area**(-concavity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import flint_slope
    >>> float(flint_slope(1000000.0, 10, 0.5))
    0.01
    """
    drainage_area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    steepness_index = quantity("steepness_index", steepness_index, lower=0)
    concavity = quantity("concavity", concavity, lower=0)
    return steepness_index * drainage_area ** (-concavity)


@model
def chi_coordinate(distance, drainage_area, reference_area, concavity):
    """The chi coordinate of a river long profile.

    Rescaling distance by drainage area turns a graded profile into a straight line, so tributaries that share a base
    level and erosion rate plot on top of each other.

    Parameters
    ----------
    distance : array_like, shape (n,)
        m upstream from the outlet; constraint: strictly increasing vector.
    drainage_area : array_like, shape (n,)
        m^2 at each node; constraint: positive vector.
    reference_area : float
        m^2; constraint: positive scalar.
    concavity : float
        1; constraint: nonnegative scalar.

    Returns
    -------
    ndarray, shape (n,)
        Model value [m], zero at the first node.

    Notes
    -----
    Perron and Royden (2013), integrated by the trapezoidal rule between nodes; in steady state with
    :func:`steady_state_slope` the elevation is linear in chi.
    Equation: ``chi = integral of (reference_area / drainage_area)**concavity d distance``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import chi_coordinate
    >>> import numpy as np
    >>> np.round(chi_coordinate([0, 1000, 2000], [4000000.0, 1000000.0, 250000.0], 1000000.0, 0.5), 6).tolist()
    [0.0, 750.0, 2250.0]
    """
    distance = quantity("distance", distance)
    drainage_area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    reference_area = quantity("reference_area", reference_area, lower=0, strict_lower=True)
    concavity = quantity("concavity", concavity, lower=0)
    require(distance.ndim == 1 and distance.size >= 2 and np.all(np.diff(distance) > 0), "distance must strictly increase")
    require(drainage_area.shape == distance.shape, "drainage_area must match distance")
    require(reference_area.ndim == concavity.ndim == 0, "chi parameters must be scalars")
    integrand = (reference_area / drainage_area) ** concavity
    return np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(distance))])


@model
def hack_length(drainage_area, coefficient, exponent):
    """Hack's law of main-stream length against drainage area.

    Larger basins are longer than they are wide, because the main stream grows faster than the square root of the
    area it drains.

    Parameters
    ----------
    drainage_area : float or array_like
        m^2; constraint: nonnegative.
    coefficient : float or array_like
        m^(1 - 2 exponent), about 0.32 for an exponent of 0.6 in SI (1.4 in miles); constraint: positive.
    exponent : float or array_like
        1, about 0.6; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hack (1957), from basins in Virginia and Maryland; the exponent drifts toward 0.5 in the largest basins.
    Equation: ``coefficient * drainage_area**exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import hack_length
    >>> round(float(hack_length(1000000.0, 0.32, 0.6)), 9)
    1273.942945771
    """
    drainage_area = quantity("drainage_area", drainage_area, lower=0)
    coefficient = quantity("coefficient", coefficient, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=0, strict_lower=True)
    return coefficient * drainage_area**exponent


@model
def strahler_order(downstream):
    """Strahler order of every link in a river network.

    Ordering streams from the headwaters down, rising by one only where two streams of equal order meet, gives the
    network a hierarchy whose counts and lengths follow Horton's laws.

    Parameters
    ----------
    downstream : array_like, shape (n,)
        Index of the link each link drains into, -1 at an outlet; constraint: integer vector forming a tree.

    Returns
    -------
    ndarray, shape (n,)
        Model value [1], the order of each link, 1 for a source.

    Notes
    -----
    Strahler (1952, 1957): a link takes the highest order among those draining into it, plus one when two or more
    share that highest order; several outlets are allowed.
    Equation: ``order = max(upstream orders) + (count of that maximum >= 2)``, 1 at sources.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import strahler_order
    >>> import numpy as np
    >>> np.round(strahler_order([2, 2, 4, 4, -1, 3]), 6).tolist()
    [1.0, 1.0, 2.0, 1.0, 2.0, 1.0]
    """
    links = quantity("downstream", downstream, lower=-1)
    require(links.ndim == 1 and np.all(links == np.round(links)), "downstream must be a vector of link indices")
    links = links.astype(int)
    size = links.size
    require(np.all(links < size) and np.all(links != np.arange(size)), "downstream must point to other links")
    inflow = np.bincount(links[links >= 0], minlength=size)
    remaining = inflow.copy()
    highest = np.zeros(size)
    ties = np.zeros(size, dtype=int)
    order = np.zeros(size)
    ready = list(np.flatnonzero(inflow == 0))
    done = 0
    while ready:
        link = ready.pop()
        done += 1
        order[link] = 1.0 if inflow[link] == 0 else highest[link] + (ties[link] >= 2)
        target = links[link]
        if target >= 0:
            if order[link] > highest[target]:
                highest[target], ties[target] = order[link], 1
            elif order[link] == highest[target]:
                ties[target] += 1
            remaining[target] -= 1
            if remaining[target] == 0:
                ready.append(target)
    require(done == size, "downstream links must form a tree without cycles")
    return order


@model
def downstream_hydraulic_geometry(discharge, width_coefficient, width_exponent, depth_coefficient, depth_exponent):
    """Width, depth and velocity of a channel as powers of discharge.

    Rivers adjust their shape to the water they carry, widening fastest, deepening more slowly and speeding up least.

    Parameters
    ----------
    discharge : float or array_like
        m^3/s, the bankfull or mean discharge; constraint: positive.
    width_coefficient : float or array_like
        m^(1 - 3 width_exponent) s^width_exponent; constraint: positive.
    width_exponent : float or array_like
        1, about 0.5 downstream; constraint: nonnegative.
    depth_coefficient : float or array_like
        m^(1 - 3 depth_exponent) s^depth_exponent; constraint: positive.
    depth_exponent : float or array_like
        1, about 0.4 downstream; constraint: nonnegative.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value: width [m], depth [m] and velocity [m/s] on the last axis.

    Notes
    -----
    Leopold and Maddock (1953): power laws fitted to gauging stations, with the velocity from continuity so that its
    coefficient and exponent are ``1 / (a c)`` and ``1 - b - f``.
    Equation: ``width = a * Q**b``, ``depth = c * Q**f``, ``velocity = Q / (width * depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import downstream_hydraulic_geometry
    >>> import numpy as np
    >>> np.round(downstream_hydraulic_geometry(100, 4, 0.5, 0.25, 0.4), 6).tolist()
    [40.0, 1.577393, 1.584893]
    """
    discharge = quantity("discharge", discharge, lower=0, strict_lower=True)
    width_coefficient = quantity("width_coefficient", width_coefficient, lower=0, strict_lower=True)
    width_exponent = quantity("width_exponent", width_exponent, lower=0)
    depth_coefficient = quantity("depth_coefficient", depth_coefficient, lower=0, strict_lower=True)
    depth_exponent = quantity("depth_exponent", depth_exponent, lower=0)
    width = width_coefficient * discharge**width_exponent
    depth = depth_coefficient * discharge**depth_exponent
    return np.stack(np.broadcast_arrays(width, depth, discharge / (width * depth)), axis=-1)


@model
def effective_discharge(median_discharge, log_standard_deviation, rating_exponent):
    """Effective discharge that carries the most sediment over time.

    Rare floods move the most sediment per event but ordinary flows happen far more often; their product peaks at a
    moderate, frequent flood that does most of the geomorphic work.

    Parameters
    ----------
    median_discharge : float or array_like
        m^3/s; constraint: positive.
    log_standard_deviation : float or array_like
        1, standard deviation of the natural logarithm of discharge; constraint: nonnegative.
    rating_exponent : float or array_like
        1, ``b`` in the sediment rating ``Qs = a Q**b``, usually 1.5-3; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^3/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Wolman and Miller (1960), with a lognormal frequency of discharge and a power-law sediment rating; the maximum of
    ``Q**b`` times the lognormal density.
    Equation: ``median_discharge * np.exp((rating_exponent - 1) * log_standard_deviation**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import effective_discharge
    >>> round(float(effective_discharge(50, 1, 2)), 9)
    135.914091423
    """
    median_discharge = quantity("median_discharge", median_discharge, lower=0, strict_lower=True)
    log_standard_deviation = quantity("log_standard_deviation", log_standard_deviation, lower=0)
    rating_exponent = quantity("rating_exponent", rating_exponent, lower=0)
    return median_discharge * np.exp((rating_exponent - 1) * log_standard_deviation**2)


@safe_arithmetic
def simulate_river_profile(
    distance,
    drainage_area,
    uplift_rate,
    erodibility,
    time,
    *,
    area_exponent=0.5,
    slope_exponent=1.0,
    initial_elevation=None,
):
    """Detachment-limited river long profile under uplift and stream-power incision.

    Parameters
    ----------
    distance : array_like, shape (m,), m >= 2
        Strictly increasing distance upstream of the outlet, outlet first [m].
    drainage_area : array_like, shape (m,)
        Positive drainage area at each node [m^2].
    uplift_rate : float
        Nonnegative scalar rock uplift rate relative to the fixed outlet [m/s].
    erodibility : float
        Positive scalar stream-power erodibility [m^(1 - 2 area_exponent)/s for a slope exponent of 1].
    time : array_like
        Nonnegative strictly increasing output times [s]; the initial state belongs to the first sample.
    area_exponent : float, optional
        Nonnegative scalar ``m`` [1].
    slope_exponent : float, optional
        Positive scalar ``n`` [1].
    initial_elevation : array_like, shape (m,), optional
        Elevation at the first sample [m], flat at zero by default; the outlet keeps its value throughout. Pass the
        steady profile for one uplift rate with another to send a knickpoint upstream.

    Returns
    -------
    TimeSeries
        Relief (highest node above the outlet) [m], mean elevation of the nodes [m] and the elevation of every node
        [m], labeled elevation_<node> from the outlet.

    Notes
    -----
    Solves ``dz/dt = U - K A**m max(dz/dx, 0)**n`` by the method of lines, with the slope of each node taken toward
    its downstream neighbour (first-order upwind, since knickpoints travel upstream) and ``integrators.integrate_ode``
    (BDF) in time. The outlet is a fixed base level; reversed slopes do not erode, and there is no sediment flux,
    erosion threshold, hillslope coupling, channel-width change or isostasy. The upwind scheme smears a knickpoint over
    a few nodes.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    distance = quantity("distance", distance)
    require(distance.ndim == 1 and distance.size >= 2 and np.all(np.diff(distance) > 0), "distance must strictly increase")
    area = quantity("drainage_area", drainage_area, lower=0, strict_lower=True)
    require(area.shape == distance.shape, "drainage_area must match distance")
    uplift = quantity("uplift_rate", uplift_rate, lower=0)
    erodibility = quantity("erodibility", erodibility, lower=0, strict_lower=True)
    m = quantity("area_exponent", area_exponent, lower=0)
    n = quantity("slope_exponent", slope_exponent, lower=0, strict_lower=True)
    require(all(v.ndim == 0 for v in (uplift, erodibility, m, n)), "profile parameters must be scalars")
    time = time_grid(time)
    initial = np.zeros_like(distance) if initial_elevation is None else quantity("initial_elevation", initial_elevation)
    require(initial.shape == distance.shape, "initial_elevation must match distance")
    coefficient = erodibility * area[1:] ** m
    spacing = np.diff(distance)

    def rhs(t, elevation):
        slope = np.maximum(np.diff(elevation) / spacing, 0.0)
        return np.concatenate([[0.0], uplift - coefficient * slope**n])

    elevation = integrate_ode(rhs, initial, time, method="BDF", rtol=1e-7, atol=1e-6).values
    relief = elevation.max(axis=1) - elevation[:, 0]
    names = ("relief", "mean_elevation") + tuple(f"elevation_{i}" for i in range(distance.size))
    units = ("m",) * (distance.size + 2)
    return TimeSeries(time, np.column_stack([relief, elevation.mean(axis=1), elevation]), names, units)
