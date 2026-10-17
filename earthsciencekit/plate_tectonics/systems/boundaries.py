"""Plate boundaries: spreading ridges, transform and convergent margins, subducting slabs and slab pull.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "spreading_age",
    "full_spreading_rate",
    "convergence_rate",
    "fault_slip_rate",
    "subduction_travel_time",
    "slab_pull_force",
]


@model
def spreading_age(distance, half_spreading_rate):
    """Magnetic-stripe spreading clock.

    Symmetric ocean-floor ages link ridge spreading to plate creation.

    Parameters
    ----------
    distance : float or array_like
        m; constraint: nonnegative.
    half_spreading_rate : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Distance on one side and half-spreading rate; changes in rate are omitted.
    Equation: ``distance / half_spreading_rate``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    distance = quantity("distance", distance, lower=0)
    half_spreading_rate = quantity("half_spreading_rate", half_spreading_rate, lower=0, strict_lower=True)
    return distance / half_spreading_rate


@model
def full_spreading_rate(half_spreading_rate):
    """Symmetric full spreading rate.

    Opposite plate velocities combine into ridge opening speed.

    Parameters
    ----------
    half_spreading_rate : float or array_like
        m/s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Symmetric spreading; asymmetric velocities require their own relative-vector calculation.
    Equation: ``2 * half_spreading_rate``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    half_spreading_rate = quantity("half_spreading_rate", half_spreading_rate, lower=0)
    return 2 * half_spreading_rate


@model
def convergence_rate(relative_speed, obliquity):
    """Signed plate-boundary convergence.

    Resolving relative motion normal to a boundary separates convergence and shear.

    Parameters
    ----------
    relative_speed : float or array_like
        m/s; constraint: nonnegative.
    obliquity : float or array_like
        rad; constraint: width.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Positive normal component denotes convergence; negative denotes divergence.
    Equation: ``relative_speed * np.cos(obliquity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    relative_speed = quantity("relative_speed", relative_speed, lower=0)
    obliquity = quantity("obliquity", obliquity, lower=0, upper=2 * np.pi)
    return relative_speed * np.cos(obliquity)


@model
def fault_slip_rate(relative_speed, obliquity):
    """Tangential boundary slip rate.

    The tangential motion component connects oblique convergence to strike slip.

    Parameters
    ----------
    relative_speed : float or array_like
        m/s; constraint: nonnegative.
    obliquity : float or array_like
        rad; constraint: width.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Supplied obliquity is oriented from the positive convergence normal.
    Equation: ``relative_speed * np.sin(obliquity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    relative_speed = quantity("relative_speed", relative_speed, lower=0)
    obliquity = quantity("obliquity", obliquity, lower=0, upper=2 * np.pi)
    return relative_speed * np.sin(obliquity)


@model
def subduction_travel_time(vertical_depth, slab_speed, dip):
    """Slab path travel time.

    Slab dip changes the path length required to reach a specified depth.

    Parameters
    ----------
    vertical_depth : float or array_like
        m; constraint: nonnegative.
    slab_speed : float or array_like
        m/s; constraint: positive.
    dip : float or array_like
        rad; constraint: inclination.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Straight slab at constant speed; zero dip cannot reach depth and is rejected.
    Equation: ``vertical_depth / (slab_speed * np.sin(dip))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    vertical_depth = quantity("vertical_depth", vertical_depth, lower=0)
    slab_speed = quantity("slab_speed", slab_speed, lower=0, strict_lower=True)
    dip = quantity("dip", dip, lower=0, upper=np.pi / 2)
    require(dip > 0, "slab dip must be positive")
    return vertical_depth / (slab_speed * np.sin(dip))


@model
def slab_pull_force(density_excess, slab_volume, gravity):
    """Slab negative-buoyancy force.

    A cold slab can drive motion through its excess weight.

    Parameters
    ----------
    density_excess : float or array_like
        kg/m^3; constraint: nonnegative.
    slab_volume : float or array_like
        m^3; constraint: nonnegative.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Simple buoyancy force magnitude; drag, phase changes and mantle flow are excluded.
    Equation: ``density_excess * slab_volume * gravity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    density_excess = quantity("density_excess", density_excess, lower=0)
    slab_volume = quantity("slab_volume", slab_volume, lower=0)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return density_excess * slab_volume * gravity
