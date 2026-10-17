"""Validated access to SciPy initial-value integration."""

import numpy as np
from scipy.integrate import solve_ivp

from earthsciencekit._validation import quantity, require
from earthsciencekit.results import TimeSeries

__all__ = ["integrate_ode"]


def integrate_ode(rhs, initial, time, *, names=None, units=None, rtol=1e-8, atol=1e-10, method="RK45"):
    """Integrate a real initial-value system and retain the requested samples.

    Parameters
    ----------
    rhs : callable
        ``rhs(t, state)`` returns derivatives with the same vector shape as state.
    initial : array_like, shape (m,)
        State at the first supplied time, in the units specified by ``units``.
    time : array_like, shape (n,)
        Strictly increasing output times in seconds, with at least two samples.
    names, units : sequence of str or None
        State labels; defaults are x0, x1, ... and dimensionless unit labels.
    rtol, atol : float
        Positive finite scalar solver tolerances; atol uses each state's units.
    method : str
        One of RK23, RK45, DOP853, Radau, BDF or LSODA.

    Returns
    -------
    TimeSeries
        Immutable sampled states. Numerical solver errors raise ValueError.

    Notes
    -----
    SciPy adapts its internal steps; the output grid is not the integration step.
    Solver tolerances do not validate a physical model or guarantee conservation.
    """
    require(callable(rhs), "rhs must be callable")
    initial = quantity("initial", initial)
    time = quantity("time", time)
    require(initial.ndim == 1, "initial must be a nonempty state vector")
    require(time.ndim == 1 and time.size >= 2, "time must have at least two samples")
    require(np.all(np.diff(time) > 0), "time must strictly increase")
    for name, value in [("rtol", rtol), ("atol", atol)]:
        tolerance = quantity(name, value, lower=0, strict_lower=True)
        require(tolerance.ndim == 0, "tolerances must be scalars")
    require(method in ("RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"), "unsupported integration method")
    names = tuple(f"x{i}" for i in range(initial.size)) if names is None else tuple(names)
    units = tuple("1" for _ in initial) if units is None else tuple(units)
    # Validate labels before invoking potentially expensive or stateful callbacks.
    TimeSeries(time, np.broadcast_to(initial, (time.size, initial.size)), names, units)

    def checked_rhs(t, state):
        derivative = quantity("rhs derivative", rhs(t, state))
        require(derivative.shape == initial.shape, "rhs must return a matching state vector")
        return derivative

    solution = solve_ivp(checked_rhs, (time[0], time[-1]), initial, t_eval=time, rtol=rtol, atol=atol, method=method)
    require(solution.success, f"integration failed: {solution.message}")
    return TimeSeries(solution.t, solution.y.T, names, units)
