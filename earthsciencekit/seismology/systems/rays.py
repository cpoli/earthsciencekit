"""Rays and travel times: direct, head and vertical waves, crossover distance, critical angle and ray parameter.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "travel_time",
    "critical_angle",
    "head_wave_travel_time",
    "crossover_distance",
    "ray_parameter",
    "vertical_travel_time",
]


@model
def travel_time(distance, speed):
    """Homogeneous travel time.

    Arrival times connect source distance and propagation speed.

    Parameters
    ----------
    distance : float or array_like
        m; constraint: nonnegative.
    speed : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Straight path in a homogeneous medium; not a ray tracer.
    Equation: ``distance / speed``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    distance = quantity("distance", distance, lower=0)
    speed = quantity("speed", speed, lower=0, strict_lower=True)
    return distance / speed


@model
def critical_angle(speed1, speed2):
    """Snell refraction critical angle.

    Snell refraction explains head waves above a faster layer.

    Parameters
    ----------
    speed1 : float or array_like
        m/s; constraint: positive.
    speed2 : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Incidence from the slower medium; second speed must exceed first.
    Equation: ``np.arcsin(speed1 / speed2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    speed1 = quantity("speed1", speed1, lower=0, strict_lower=True)
    speed2 = quantity("speed2", speed2, lower=0, strict_lower=True)
    require(speed2 > speed1, "second medium must be faster")
    return np.arcsin(speed1 / speed2)


@model
def head_wave_travel_time(distance, thickness, speed1, speed2):
    """Head-wave travel time below a single layer.

    A wave critically refracted along a faster layer overtakes the direct wave at long range.

    Parameters
    ----------
    distance : float or array_like
        m; constraint: nonnegative, at least the critical distance.
    thickness : float or array_like
        m; constraint: positive.
    speed1 : float or array_like
        m/s; constraint: positive.
    speed2 : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Flat layer of speed1 over a faster half-space of speed2, surface source and receiver.
    The head wave exists only beyond the critical distance ``2 * thickness * tan(critical_angle)``.
    Equation: ``distance / speed2 + 2 * thickness * np.sqrt(speed2**2 - speed1**2) / (speed1 * speed2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    distance = quantity("distance", distance, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    speed1 = quantity("speed1", speed1, lower=0, strict_lower=True)
    speed2 = quantity("speed2", speed2, lower=0, strict_lower=True)
    require(speed2 > speed1, "second medium must be faster")
    require(distance >= 2 * thickness * speed1 / np.sqrt(speed2**2 - speed1**2), "distance is inside the critical distance")
    return distance / speed2 + 2 * thickness * np.sqrt(speed2**2 - speed1**2) / (speed1 * speed2)


@model
def crossover_distance(thickness, speed1, speed2):
    """Head-wave crossover distance.

    The range at which the refracted wave arrives first measures the depth of the layer.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: positive.
    speed1 : float or array_like
        m/s; constraint: positive.
    speed2 : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Flat layer of speed1 over a faster half-space of speed2; direct and head-wave times are equal there.
    Equation: ``2 * thickness * np.sqrt((speed2 + speed1) / (speed2 - speed1))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    speed1 = quantity("speed1", speed1, lower=0, strict_lower=True)
    speed2 = quantity("speed2", speed2, lower=0, strict_lower=True)
    require(speed2 > speed1, "second medium must be faster")
    return 2 * thickness * np.sqrt((speed2 + speed1) / (speed2 - speed1))


@model
def ray_parameter(radius, speed, incidence_angle):
    """Ray parameter in a spherically layered Earth.

    Snell's law on a sphere keeps one quantity fixed along each ray, which fixes where it emerges.

    Parameters
    ----------
    radius : float or array_like
        m; constraint: positive.
    speed : float or array_like
        m/s; constraint: positive.
    incidence_angle : float or array_like
        rad; constraint: between 0 and pi/2, measured from the radial direction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s/rad]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Spherically symmetric isotropic Earth; the value is the same at every point of one ray.
    Equation: ``radius * np.sin(incidence_angle) / speed``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    speed = quantity("speed", speed, lower=0, strict_lower=True)
    incidence_angle = quantity("incidence_angle", incidence_angle, lower=0, upper=np.pi / 2)
    return radius * np.sin(incidence_angle) / speed


@model
def vertical_travel_time(depth, speed):
    """Vertical travel time through a tabulated speed profile.

    Integrating slowness through a reference Earth model predicts arrival times.

    Parameters
    ----------
    depth : array_like, shape (n,), n >= 2
        m; constraint: nonnegative and nondecreasing. A repeated depth marks a discontinuity.
    speed : array_like, shape (n,)
        m/s; constraint: positive.

    Returns
    -------
    ndarray, shape (n,)
        Model value [s]: one-way vertical time from the first depth to each depth.

    Notes
    -----
    Speed varies linearly in depth between tabulated points; the integral of each segment is exact.
    Equation: ``cumsum(dz * np.log(v2 / v1) / (v2 - v1))``, which is ``dz / v1`` where ``v2 == v1``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    speed = quantity("speed", speed, lower=0, strict_lower=True)
    require(depth.ndim == 1 and depth.size >= 2, "depth must be a vector with at least two samples")
    require(speed.shape == depth.shape, "speed must have the same shape as depth")
    require(np.all(np.diff(depth) >= 0), "depth must not decrease")
    change = np.diff(speed) / speed[:-1]
    factor = np.ones_like(change)
    varying = change != 0
    factor[varying] = np.log1p(change[varying]) / change[varying]
    return np.concatenate([[0.0], np.cumsum(np.diff(depth) / speed[:-1] * factor)])
