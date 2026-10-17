"""Gravity: buried masses, the Bouguer and free-air reductions and the weighing of the Earth.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import G

__all__ = [
    "sphere_gravity",
    "bouguer_correction",
    "free_air_correction",
    "bouguer_anomaly",
    "torsion_balance_gravitational_constant",
    "earth_mean_density",
]


@model
def sphere_gravity(x, depth, radius, density_contrast):
    """Gravitational attraction of a uniform sphere.

    A compact body with excess density pulls a plumb line sideways and adds to gravity above it.

    Parameters
    ----------
    x : float or array_like
        m; constraint: real. Horizontal distance of the observer from the point above the centre.
    depth : float or array_like
        m; constraint: real. Depth of the centre below the observer, positive downward.
    radius : float or array_like
        m; constraint: positive.
    density_contrast : float or array_like
        kg/m^3; constraint: real. Density minus that of the surroundings.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [m/s^2]: the horizontal component, positive toward increasing ``x``, and the downward component.

    Notes
    -----
    The observer must lie outside the sphere, where it attracts as a point mass ``M = 4/3 pi radius**3 density_contrast``.
    The horizontal component always points toward the centre for positive contrast.
    Equation: ``G * M * [-x, depth] / (x**2 + depth**2) ** 1.5``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    x = quantity("x", x)
    depth = quantity("depth", depth)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    density_contrast = quantity("density_contrast", density_contrast)
    distance_squared = x**2 + depth**2
    require(distance_squared >= radius**2, "the observer must lie outside the sphere")
    mass = 4 / 3 * np.pi * radius**3 * density_contrast
    scale = G * mass / distance_squared**1.5
    return np.stack(np.broadcast_arrays(-scale * x, scale * depth), axis=-1)


@model
def bouguer_correction(density, thickness):
    """Attraction of an infinite horizontal slab.

    Rock between the station and the reference level adds gravity that does not depend on the slab's extent.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: positive.
    thickness : float or array_like
        m; constraint: real. Signed station height above the reference level.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Flat, uniform, laterally infinite slab; no terrain or curvature correction. About 0.112 mGal per metre at
    2670 kg/m^3.
    Equation: ``2 * np.pi * G * density * thickness``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    density = quantity("density", density, lower=0, strict_lower=True)
    thickness = quantity("thickness", thickness)
    return 2 * np.pi * G * density * thickness


@model
def free_air_correction(height, gravity, radius):
    """Free-air gravity correction.

    Gravity falls off with distance from the Earth's centre, about 0.3086 mGal per metre at the surface.

    Parameters
    ----------
    height : float or array_like
        m; constraint: real. Station height above the reference surface.
    gravity : float or array_like
        m/s^2; constraint: positive. Gravity on the reference surface.
    radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s^2]: the amount to add to gravity observed at ``height``.

    Notes
    -----
    First order in height / radius for a spherical point-mass Earth; no rock between station and reference.
    Equation: ``2 * gravity * height / radius``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    height = quantity("height", height)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    return 2 * gravity * height / radius


@model
def bouguer_anomaly(observed_gravity, normal_gravity, height, density, radius):
    """Simple Bouguer gravity anomaly.

    Removing latitude, height and the rock slab under a station leaves the attraction of buried density contrasts.

    Parameters
    ----------
    observed_gravity : float or array_like
        m/s^2; constraint: positive.
    normal_gravity : float or array_like
        m/s^2; constraint: positive. Normal gravity on the reference surface at the station's latitude.
    height : float or array_like
        m; constraint: real. Station height above the reference surface.
    density : float or array_like
        kg/m^3; constraint: positive. Reduction density of the slab.
    radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Free-air and infinite-slab reductions only, without terrain correction. The free-air anomaly is the same
    expression without the slab term.
    Equation: ``observed_gravity - normal_gravity + 2 * normal_gravity * height / radius - 2 * np.pi * G * density * height``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    observed_gravity = quantity("observed_gravity", observed_gravity, lower=0, strict_lower=True)
    normal_gravity = quantity("normal_gravity", normal_gravity, lower=0, strict_lower=True)
    height = quantity("height", height)
    density = quantity("density", density, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    return observed_gravity - normal_gravity + 2 * normal_gravity * height / radius - 2 * np.pi * G * density * height


@model
def torsion_balance_gravitational_constant(arm_length, separation, deflection, attracting_mass, period):
    """Gravitational constant from a torsion balance.

    The twist of a suspended rod between large masses, timed by its free swing, measures G.

    Parameters
    ----------
    arm_length : float or array_like
        m; constraint: positive. Distance between the two small balls.
    separation : float or array_like
        m; constraint: positive. Centre distance between each small ball and its large ball.
    deflection : float or array_like
        rad; constraint: positive. Equilibrium twist of the rod.
    attracting_mass : float or array_like
        kg; constraint: positive. Mass of each large ball.
    period : float or array_like
        s; constraint: positive. Free oscillation period of the rod.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^3/(kg s^2)]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Point masses, a massless rod and attraction by the near large ball only; the small mass cancels.
    Equation: ``2 * np.pi**2 * arm_length * separation**2 * deflection / (attracting_mass * period**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    arm_length = quantity("arm_length", arm_length, lower=0, strict_lower=True)
    separation = quantity("separation", separation, lower=0, strict_lower=True)
    deflection = quantity("deflection", deflection, lower=0, strict_lower=True)
    attracting_mass = quantity("attracting_mass", attracting_mass, lower=0, strict_lower=True)
    period = quantity("period", period, lower=0, strict_lower=True)
    return 2 * np.pi**2 * arm_length * separation**2 * deflection / (attracting_mass * period**2)


@model
def earth_mean_density(gravity, radius, gravitational_constant=G):
    """Mean density of the Earth.

    Surface gravity and the gravitational constant together weigh the Earth.

    Parameters
    ----------
    gravity : float or array_like
        m/s^2; constraint: positive.
    radius : float or array_like
        m; constraint: positive.
    gravitational_constant : float or array_like
        m^3/(kg s^2); constraint: positive. Defaults to :data:`earthsciencekit.constants.G`.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Spherically symmetric, nonrotating Earth, so ``gravity = G M / radius**2``.
    Equation: ``3 * gravity / (4 * np.pi * gravitational_constant * radius)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    gravitational_constant = quantity("gravitational_constant", gravitational_constant, lower=0, strict_lower=True)
    return 3 * gravity / (4 * np.pi * gravitational_constant * radius)
