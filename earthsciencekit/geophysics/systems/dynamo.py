"""Geodynamo: magnetic diffusion, the magnetic Reynolds number and Rikitake's two-disc dynamo.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.constants import MU0

__all__ = [
    "magnetic_reynolds_number",
    "magnetic_decay_time",
    "rikitake_fixed_point",
    "simulate_rikitake_dynamo",
]


@model
def magnetic_reynolds_number(speed, length, conductivity):
    """Magnetic Reynolds number.

    A conducting fluid can stretch and regenerate a magnetic field only if advection beats diffusion.

    Parameters
    ----------
    speed : float or array_like
        m/s; constraint: nonnegative.
    length : float or array_like
        m; constraint: positive.
    conductivity : float or array_like
        S/m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Order-of-magnitude scale; dynamos need values above a geometry-dependent threshold of order 10 to 100.
    Equation: ``MU0 * conductivity * speed * length``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import magnetic_reynolds_number
    >>> float(magnetic_reynolds_number(2, 3, 795774.7154594767))
    6.0
    """
    speed = quantity("speed", speed, lower=0)
    length = quantity("length", length, lower=0, strict_lower=True)
    conductivity = quantity("conductivity", conductivity, lower=0, strict_lower=True)
    return MU0 * conductivity * speed * length


@model
def magnetic_decay_time(radius, conductivity):
    """Free decay time of the slowest magnetic mode of a conducting sphere.

    Without a dynamo the Earth's field would die in tens of thousands of years, far less than its age.

    Parameters
    ----------
    radius : float or array_like
        m; constraint: positive.
    conductivity : float or array_like
        S/m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]: e-folding time of the field of the slowest dipole mode.

    Notes
    -----
    Uniform stationary conductor surrounded by an insulator.
    Equation: ``MU0 * conductivity * radius**2 / np.pi**2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import magnetic_decay_time
    >>> float(magnetic_decay_time(3.141592653589793, 795774.7154594767))
    1.0
    """
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    conductivity = quantity("conductivity", conductivity, lower=0, strict_lower=True)
    return MU0 * conductivity * radius**2 / np.pi**2


@model
def rikitake_fixed_point(mu, k):
    """Steady state of Rikitake's two-disc dynamo.

    The two coupled discs have a pair of steady states of opposite polarity, which the chaotic flow visits in turn.

    Parameters
    ----------
    mu : float or array_like
        1; constraint: positive. Scaled resistance.
    k : float or array_like
        1; constraint: positive. Sets ``A = mu * (k**2 - 1 / k**2)``.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [1]: the scaled currents ``x``, ``y`` and the scaled angular speed ``z`` of the positive state; the
        other state is ``(-x, -y, z)``.

    Notes
    -----
    Nondimensional system of :func:`simulate_rikitake_dynamo`. Both states are unstable for the usual parameters.
    Equation: ``[k, 1 / k, mu * k**2]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import rikitake_fixed_point
    >>> np.round(rikitake_fixed_point(3, 2), 4).tolist()
    [2.0, 0.5, 12.0]
    """
    mu = quantity("mu", mu, lower=0, strict_lower=True)
    k = quantity("k", k, lower=0, strict_lower=True)
    return np.stack(np.broadcast_arrays(k, 1 / k, mu * k**2), axis=-1)


@safe_arithmetic
def simulate_rikitake_dynamo(mu, k, time_scale, time, *, initial_state=(1.0, 0.0, 0.0)):
    """Integrate Rikitake's coupled two-disc dynamo.

    Parameters
    ----------
    mu : float
        Positive scalar scaled resistance [1].
    k : float
        Positive scalar setting the torque asymmetry ``A = mu * (k**2 - 1 / k**2)`` [1].
    time_scale : float
        Positive scalar seconds per nondimensional time unit [s].
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    initial_state : sequence of float
        Scaled currents ``x``, ``y`` and angular speed ``z`` at the first time [1].

    Returns
    -------
    TimeSeries
        The scaled currents ``x`` and ``y`` (whose sign is the polarity) and the scaled angular speed ``z`` [1].

    Notes
    -----
    ``dx/dt = -mu * x + y * z``, ``dy/dt = -mu * y + (z - A) * x``, ``dz/dt = 1 - x * y`` in units of ``time_scale``.
    The steady states are :func:`rikitake_fixed_point`. Trajectories oscillate about one and jump irregularly to the
    other, reversing both currents: a chaotic analogue of geomagnetic reversals, not a model of the core.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.geophysics.core.base import TimeSeries
    from earthsciencekit.integrators import integrate_ode

    mu = quantity("mu", mu, lower=0, strict_lower=True)
    k = quantity("k", k, lower=0, strict_lower=True)
    time_scale = quantity("time_scale", time_scale, lower=0, strict_lower=True)
    initial_state = quantity("initial_state", initial_state)
    require(all(v.ndim == 0 for v in (mu, k, time_scale)), "dynamo parameters must be scalars")
    require(initial_state.shape == (3,), "initial_state must hold x, y and z")
    time = time_grid(time)
    a = mu * (k**2 - 1 / k**2)

    def rhs(t, state):
        x, y, z = state
        return [(-mu * x + y * z) / time_scale, (-mu * y + (z - a) * x) / time_scale, (1 - x * y) / time_scale]

    result = integrate_ode(rhs, initial_state, time, names=("x", "y", "z"), units=("1", "1", "1"), method="DOP853", rtol=1e-10)
    return TimeSeries(result.time, result.values, result.names, result.units)
