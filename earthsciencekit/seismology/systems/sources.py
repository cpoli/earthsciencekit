"""Earthquake sources: seismic moment, magnitudes, frequency-magnitude and aftershock laws, and a damped oscillator.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "seismic_moment",
    "moment_magnitude",
    "moment_from_magnitude",
    "gutenberg_richter",
    "omori_rate",
    "local_magnitude",
    "simulate_oscillator",
]


@model
def seismic_moment(rigidity, area, slip):
    """Fault seismic moment.

    Moment connects an earthquake to fault area and displacement.

    Parameters
    ----------
    rigidity : float or array_like
        Pa; constraint: positive.
    area : float or array_like
        m^2; constraint: positive.
    slip : float or array_like
        m; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Scalar mean-slip estimate; no spatial rupture distribution.
    Equation: ``rigidity * area * slip``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    rigidity = quantity("rigidity", rigidity, lower=0, strict_lower=True)
    area = quantity("area", area, lower=0, strict_lower=True)
    slip = quantity("slip", slip, lower=0)
    return rigidity * area * slip


@model
def moment_magnitude(moment):
    """Moment magnitude in SI.

    Moment magnitude relates logarithmic event size to scalar moment.

    Parameters
    ----------
    moment : float or array_like
        N m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uses the 9.1 intercept with moment in N m, not dyne cm.
    Equation: ``(2 / 3) * (np.log10(moment) - 9.1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    moment = quantity("moment", moment, lower=0, strict_lower=True)
    return (2 / 3) * (np.log10(moment) - 9.1)


@model
def moment_from_magnitude(magnitude):
    """Inverse moment magnitude.

    Inverting magnitude enables comparisons on a physical moment scale.

    Parameters
    ----------
    magnitude : float or array_like
        1; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Same SI convention as moment_magnitude; very large inputs can overflow.
    Equation: ``10 ** (1.5 * magnitude + 9.1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    magnitude = quantity("magnitude", magnitude)
    return 10 ** (1.5 * magnitude + 9.1)


@model
def local_magnitude(amplitude, distance):
    """Richter local magnitude.

    A logarithm of recorded amplitude, corrected for distance, ranks earthquakes by size.

    Parameters
    ----------
    amplitude : float or array_like
        m; constraint: positive. Peak Wood-Anderson displacement amplitude.
    distance : float or array_like
        m; constraint: positive. Hypocentral distance.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Richter's anchor (1 mm at 100 km is magnitude 3) with the Hutton and Boore (1987)
    southern California distance correction, fitted from about 10 to 700 km.
    Equation: ``np.log10(amplitude / 1e-3) + 1.11 * np.log10(distance / 1e5) + 1.89e-6 * (distance - 1e5) + 3``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    amplitude = quantity("amplitude", amplitude, lower=0, strict_lower=True)
    distance = quantity("distance", distance, lower=0, strict_lower=True)
    return np.log10(amplitude / 1e-3) + 1.11 * np.log10(distance / 1e5) + 1.89e-6 * (distance - 1e5) + 3


@model
def gutenberg_richter(a, b, magnitude):
    """Earthquake frequency-size relation.

    A logarithmic frequency-size law summarizes an earthquake catalogue.

    Parameters
    ----------
    a : float or array_like
        1; constraint: real.
    b : float or array_like
        1; constraint: positive.
    magnitude : float or array_like
        1; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Expected cumulative count for the fitted catalogue interval; completeness matters.
    Equation: ``10 ** (a - b * magnitude)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    a = quantity("a", a)
    b = quantity("b", b, lower=0, strict_lower=True)
    magnitude = quantity("magnitude", magnitude)
    return 10 ** (a - b * magnitude)


@model
def omori_rate(productivity, offset, exponent, time):
    """Modified Omori aftershock rate.

    Aftershock rates often decay with time after a main shock.

    Parameters
    ----------
    productivity : float or array_like
        count s^(p-1); constraint: nonnegative.
    offset : float or array_like
        s; constraint: positive.
    exponent : float or array_like
        1; constraint: positive.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Empirical fitted rate; no earthquake timing or location prediction.
    Equation: ``productivity / (offset + time)**exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    productivity = quantity("productivity", productivity, lower=0)
    offset = quantity("offset", offset, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0)
    return productivity / (offset + time) ** exponent


@safe_arithmetic
def simulate_oscillator(frequency, quality_factor, time, *, displacement=1.0, velocity=0.0):
    """Integrate an unforced damped seismic oscillator.

    Parameters
    ----------
    frequency : float
        Positive scalar undamped natural frequency [Hz].
    quality_factor : float
        Positive scalar oscillator Q; damping ratio is 1/(2Q).
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    displacement, velocity : float
        Finite scalar initial displacement [m] and velocity [m/s].

    Returns
    -------
    TimeSeries
        Displacement [m] and velocity [m/s].

    Notes
    -----
    Solves x'' + omega/Q x' + omega^2 x = 0. This single oscillator can
    illustrate a seismometer response; it is not an earthquake rupture model.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    frequency = quantity("frequency", frequency, lower=0, strict_lower=True)
    quality_factor = quantity("quality_factor", quality_factor, lower=0, strict_lower=True)
    displacement = quantity("displacement", displacement)
    velocity = quantity("velocity", velocity)
    require(all(v.ndim == 0 for v in (frequency, quality_factor, displacement, velocity)), "parameters must be scalars")
    time = time_grid(time)
    omega = 2 * np.pi * frequency
    return integrate_ode(
        lambda t, y: [y[1], -(omega**2) * y[0] - omega / quality_factor * y[1]],
        [displacement, velocity],
        time,
        names=("displacement", "velocity"),
        units=("m", "m/s"),
    )
