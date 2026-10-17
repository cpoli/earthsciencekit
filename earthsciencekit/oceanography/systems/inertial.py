"""Inertial motion: the inertial period and damped inertial currents.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "inertial_period",
    "simulate_inertial_current",
]


@model
def inertial_period(coriolis):
    """Inertial oscillation period.

    Rotation supplies a natural timescale for ocean motion.

    Parameters
    ----------
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    No pressure gradient or friction; invalid at zero Coriolis parameter.
    Equation: ``2 * np.pi / np.abs(coriolis)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    return 2 * np.pi / np.abs(coriolis)


@safe_arithmetic
def simulate_inertial_current(coriolis, time, *, east_speed=1.0, north_speed=0.0, damping=0.0):
    """Integrate rotating, optionally damped horizontal ocean currents.

    Parameters
    ----------
    coriolis : float
        Finite nonzero scalar Coriolis parameter [1/s], signed by hemisphere.
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    east_speed, north_speed : float
        Initial signed horizontal velocity components [m/s].
    damping : float
        Nonnegative scalar linear decay rate [1/s].

    Returns
    -------
    TimeSeries
        Eastward and northward velocity [m/s].

    Notes
    -----
    du/dt = f v - r u; dv/dt = -f u - r v. With no damping, speed is
    conserved and northern-hemisphere motion rotates clockwise. Pressure
    gradients and wind forcing are omitted.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    coriolis = quantity("coriolis", coriolis)
    damping = quantity("damping", damping, lower=0)
    east_speed = quantity("east_speed", east_speed)
    north_speed = quantity("north_speed", north_speed)
    require(all(v.ndim == 0 for v in (coriolis, damping, east_speed, north_speed)), "current parameters must be scalars")
    require(coriolis != 0, "coriolis must be nonzero")
    time = time_grid(time)
    return integrate_ode(
        lambda t, y: [coriolis * y[1] - damping * y[0], -coriolis * y[0] - damping * y[1]],
        [east_speed, north_speed],
        time,
        names=("east_velocity", "north_velocity"),
        units=("m/s", "m/s"),
    )
