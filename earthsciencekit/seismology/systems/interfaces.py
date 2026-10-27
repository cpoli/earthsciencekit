"""Interfaces: seismic impedance and normal-incidence reflection and transmission.

SI units unless stated otherwise.
"""

from earthsciencekit._validation import model, quantity

__all__ = [
    "seismic_impedance",
    "reflection_coefficient",
    "transmission_coefficient",
]


@model
def seismic_impedance(density, speed):
    """Acoustic impedance.

    Impedance contrasts control normal-incidence reflection.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: positive.
    speed : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Scalar acoustic impedance; elastic mode conversion is excluded.
    Equation: ``density * speed``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import seismic_impedance
    >>> float(seismic_impedance(2000, 3000))
    6000000.0
    """
    density = quantity("density", density, lower=0, strict_lower=True)
    speed = quantity("speed", speed, lower=0, strict_lower=True)
    return density * speed


@model
def reflection_coefficient(impedance1, impedance2):
    """Normal-incidence pressure reflection.

    An interface divides incident energy into reflected and transmitted waves.

    Parameters
    ----------
    impedance1 : float or array_like
        kg/(m^2 s); constraint: positive.
    impedance2 : float or array_like
        kg/(m^2 s); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Pressure amplitude convention at normal incidence; coefficient can be negative.
    Equation: ``(impedance2 - impedance1) / (impedance2 + impedance1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import reflection_coefficient
    >>> float(reflection_coefficient(2, 6))
    0.5
    """
    impedance1 = quantity("impedance1", impedance1, lower=0, strict_lower=True)
    impedance2 = quantity("impedance2", impedance2, lower=0, strict_lower=True)
    return (impedance2 - impedance1) / (impedance2 + impedance1)


@model
def transmission_coefficient(impedance1, impedance2):
    """Normal-incidence pressure transmission.

    Pressure transmission differs from transmitted energy fraction.

    Parameters
    ----------
    impedance1 : float or array_like
        kg/(m^2 s); constraint: positive.
    impedance2 : float or array_like
        kg/(m^2 s); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Normal-incidence acoustic pressure amplitudes, not displacement amplitudes.
    Equation: ``2 * impedance2 / (impedance2 + impedance1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import transmission_coefficient
    >>> float(transmission_coefficient(2, 6))
    1.5
    """
    impedance1 = quantity("impedance1", impedance1, lower=0, strict_lower=True)
    impedance2 = quantity("impedance2", impedance2, lower=0, strict_lower=True)
    return 2 * impedance2 / (impedance2 + impedance1)
