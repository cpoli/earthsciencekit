"""Central contracts for real, finite, nonempty numerical quantities."""

from functools import wraps

import numpy as np


def quantity(name, value, *, lower=None, upper=None, strict_lower=False, strict_upper=False):
    """Return a floating NumPy quantity after validating physical bounds."""
    array = np.asarray(value)
    if array.dtype.kind not in "iuf" or not array.size or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite real numbers, excluding booleans")
    array = array.astype(float)
    if lower is not None:
        require(array > lower if strict_lower else array >= lower, f"{name} violates its lower bound")
    if upper is not None:
        require(array < upper if strict_upper else array <= upper, f"{name} violates its upper bound")
    return array


def require(condition, message):
    """Reject any violated relation, including array-valued physical constraints."""
    if not np.all(condition):
        raise ValueError(message)


def model(function):
    """Preserve signatures and turn arithmetic singularities into ValueError."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                result = function(*args, **kwargs)
        except FloatingPointError as error:
            raise ValueError(f"{function.__name__}: singular or overflowing arithmetic") from error
        require(np.all(np.isfinite(result)), f"{function.__name__}: nonfinite result")
        return result

    return wrapped


def time_grid(value):
    """Validate a nonnegative strictly increasing vector of sampled seconds."""
    time = quantity("time", value, lower=0)
    require(time.ndim == 1 and time.size >= 2, "time must be a vector with at least two samples")
    require(np.all(np.diff(time) > 0), "time must strictly increase")
    return time


def safe_arithmetic(function):
    """Translate floating arithmetic failures in structured numerical models."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                return function(*args, **kwargs)
        except FloatingPointError as error:
            raise ValueError(f"{function.__name__}: singular or overflowing arithmetic") from error

    return wrapped
