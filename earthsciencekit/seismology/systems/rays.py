"""Rays and travel times: direct, head and vertical waves, crossover distance, critical angle, ray parameter and rays in a layered sphere.

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
    "layered_ray_distance",
    "layered_ray_travel_time",
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

    Examples
    --------
    >>> from earthsciencekit.seismology import travel_time
    >>> float(travel_time(3000, 1500))
    2.0
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

    Examples
    --------
    >>> from earthsciencekit.seismology import critical_angle
    >>> round(float(critical_angle(1, 2)), 6)
    0.523599
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

    Examples
    --------
    >>> from earthsciencekit.seismology import head_wave_travel_time
    >>> round(float(head_wave_travel_time(4, 1, 1, 2)), 5)
    3.73205
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

    Examples
    --------
    >>> from earthsciencekit.seismology import crossover_distance
    >>> round(float(crossover_distance(1, 1, 3)), 5)
    2.82843
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

    Examples
    --------
    >>> from earthsciencekit.seismology import ray_parameter
    >>> round(float(ray_parameter(2, 4, 0.5235987755982988)), 6)
    0.25
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

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.seismology import vertical_travel_time
    >>> np.round(vertical_travel_time([0, 1], [1, 2.718281828459045]), 6).tolist()
    [0.0, 0.581977]
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


def _layered_ray(ray_parameter, depth, speed, planet_radius):
    """Epicentral distance and travel time of surface-to-surface rays through Mohorovicic-law layers."""
    p = quantity("ray_parameter", ray_parameter, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0)
    speed = quantity("speed", speed, lower=0, strict_lower=True)
    planet_radius = quantity("planet_radius", planet_radius, lower=0, strict_lower=True)
    require(depth.ndim == 1 and depth.size >= 2, "depth must be a vector with at least two samples")
    require(speed.shape == depth.shape, "speed must have the same shape as depth")
    require(planet_radius.ndim == 0, "planet_radius must be a scalar")
    require(depth[0] == 0, "depth must start at the surface")
    require(np.all(np.diff(depth) >= 0), "depth must not decrease")
    require(depth[-1] < planet_radius, "depth must stay above the centre")
    radius = planet_radius - depth
    eta = radius / speed
    layer = np.flatnonzero(np.diff(depth) > 0)
    top, bottom = layer, layer + 1
    require(np.all(eta[top] > eta[bottom]), "radius / speed must decrease with depth within every layer")
    require(p < eta[0], "ray_parameter must be below the surface radius / speed")
    distance, time = np.zeros_like(p), np.zeros_like(p)
    active = np.ones(p.shape, dtype=bool)
    for upper, lower in zip(top, bottom):
        # In a layer with speed = a * r**b, theta and T are exact in eta = r / speed (Mohorovicic 1910).
        scale = np.log(radius[upper] / radius[lower]) / np.log(eta[upper] / eta[lower])
        enters = active & (p < eta[upper])
        turns = enters & (p > eta[lower])
        through = enters & ~turns
        cosine_top = np.where(enters, p / eta[upper], 1.0)
        cosine_bottom = np.where(through, p / eta[lower], 1.0)
        distance += scale * (np.arccos(cosine_top) - np.arccos(cosine_bottom))
        time += (
            scale
            * eta[upper]
            * (np.sqrt(1 - cosine_top**2) - np.where(through, eta[lower] / eta[upper], 0.0) * np.sqrt(1 - cosine_bottom**2))
        )
        active = through
    require(~active, "ray passes below the deepest tabulated level")
    return 2 * distance, 2 * time


@model
def layered_ray_distance(ray_parameter, depth, speed, planet_radius):
    """Epicentral distance of a ray through a radially layered Earth.

    Rays bend back to the surface because speed grows with depth, and the distance each ray travels maps the
    velocity structure, including the shadow cast by a slower core.

    Parameters
    ----------
    ray_parameter : float or array_like
        s/rad, ``r sin(i) / v`` conserved along the ray; constraint: positive and below the surface ``radius / speed``.
    depth : array_like, shape (n,), n >= 2
        m; constraint: starts at 0, nondecreasing, below planet_radius. A repeated depth marks a discontinuity.
    speed : array_like, shape (n,)
        m/s; constraint: positive, with ``(planet_radius - depth) / speed`` decreasing with depth within every layer.
    planet_radius : float
        m; constraint: positive scalar.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad], surface source to surface receiver; it can exceed pi for rays that pass the antipode.
        Scalar ray parameters return a NumPy scalar; arrays keep their shape.

    Notes
    -----
    Speed follows the Mohorovicic law ``v = a r**b`` between tabulated points, so the ray integrals are exact. A ray
    turns where ``r / v`` falls to the ray parameter, or at the top of a layer it cannot enter, and must turn above the
    deepest tabulated level.
    Equation: ``2 * sum(np.arccos(p / eta_top) - np.arccos(p / eta_bottom)) / (1 - b)`` with ``eta = r / v``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import layered_ray_distance
    >>> round(float(layered_ray_distance(1, [0, 1.5], [1, 1], 2)), 5)
    2.0944
    """
    return _layered_ray(ray_parameter, depth, speed, planet_radius)[0]


@model
def layered_ray_travel_time(ray_parameter, depth, speed, planet_radius):
    """Travel time of a ray through a radially layered Earth.

    Paired with the epicentral distance it gives the travel-time curve, the table every earthquake location rests on.

    Parameters
    ----------
    ray_parameter : float or array_like
        s/rad; constraint: positive and below the surface ``radius / speed``.
    depth : array_like, shape (n,), n >= 2
        m; constraint: starts at 0, nondecreasing, below planet_radius. A repeated depth marks a discontinuity.
    speed : array_like, shape (n,)
        m/s; constraint: positive, with ``(planet_radius - depth) / speed`` decreasing with depth within every layer.
    planet_radius : float
        m; constraint: positive scalar.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s], surface source to surface receiver. Scalar ray parameters return a NumPy scalar.

    Notes
    -----
    Same rays as :func:`layered_ray_distance`; the slope of the travel-time curve ``dT/d(distance)`` is the ray
    parameter (Benndorf's relation).
    Equation: ``2 * sum(np.sqrt(eta_top**2 - p**2) - np.sqrt(eta_bottom**2 - p**2)) / (1 - b)`` with ``eta = r / v``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import layered_ray_travel_time
    >>> round(float(layered_ray_travel_time(1, [0, 1.5], [1, 1], 2)), 5)
    3.4641
    """
    return _layered_ray(ray_parameter, depth, speed, planet_radius)[1]
