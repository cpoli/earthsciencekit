"""Elastic waves: body and surface wave speeds, the Ricker wavelet and anelastic attenuation.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "p_wave_speed",
    "s_wave_speed",
    "poisson_ratio",
    "ricker_wavelet",
    "attenuated_amplitude",
    "rayleigh_wave_speed",
    "love_wave_angular_frequency",
]


@model
def p_wave_speed(bulk_modulus, shear_modulus, density):
    """Compressional-wave elasticity.

    Elastic moduli and density determine compressional wave speed.

    Parameters
    ----------
    bulk_modulus : float or array_like
        Pa; constraint: positive.
    shear_modulus : float or array_like
        Pa; constraint: nonnegative.
    density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Homogeneous isotropic small-strain medium; shear modulus may vanish in a fluid.
    Equation: ``np.sqrt((bulk_modulus + 4 * shear_modulus / 3) / density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import p_wave_speed
    >>> float(p_wave_speed(12, 3, 4))
    2.0
    """
    bulk_modulus = quantity("bulk_modulus", bulk_modulus, lower=0, strict_lower=True)
    shear_modulus = quantity("shear_modulus", shear_modulus, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    return np.sqrt((bulk_modulus + 4 * shear_modulus / 3) / density)


@model
def s_wave_speed(shear_modulus, density):
    """Shear-wave elasticity.

    Shear waves reveal rigidity and disappear in an ideal fluid.

    Parameters
    ----------
    shear_modulus : float or array_like
        Pa; constraint: nonnegative.
    density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Homogeneous isotropic medium; zero shear modulus returns zero speed.
    Equation: ``np.sqrt(shear_modulus / density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import s_wave_speed
    >>> float(s_wave_speed(16, 4))
    2.0
    """
    shear_modulus = quantity("shear_modulus", shear_modulus, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    return np.sqrt(shear_modulus / density)


@model
def poisson_ratio(p_speed, s_speed):
    """Elastic velocity ratio.

    Wave-speed ratios constrain isotropic elastic response.

    Parameters
    ----------
    p_speed : float or array_like
        m/s; constraint: positive.
    s_speed : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Requires positive bulk modulus: vp squared must exceed 4 vs squared over 3.
    Equation: ``(p_speed**2 - 2 * s_speed**2) / (2 * (p_speed**2 - s_speed**2))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import poisson_ratio
    >>> round(float(poisson_ratio(2, 1)), 6)
    0.333333
    """
    p_speed = quantity("p_speed", p_speed, lower=0, strict_lower=True)
    s_speed = quantity("s_speed", s_speed, lower=0, strict_lower=True)
    require(p_speed**2 > 4 * s_speed**2 / 3, "wave speeds imply nonpositive bulk modulus")
    return (p_speed**2 - 2 * s_speed**2) / (2 * (p_speed**2 - s_speed**2))


@model
def rayleigh_wave_speed(p_speed, s_speed):
    """Rayleigh surface-wave speed.

    A free surface guides a slower wave that combines compression and shear.

    Parameters
    ----------
    p_speed : float or array_like
        m/s; constraint: positive.
    s_speed : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Homogeneous isotropic half-space with a traction-free surface; nondispersive.
    Requires positive bulk modulus: vp squared must exceed 4 vs squared over 3.
    Equation: ``s_speed * np.sqrt(xi)``, where ``xi`` is the root in (0, 1) of
    ``xi**3 - 8 * xi**2 + (24 - 16 * gamma) * xi - 16 * (1 - gamma) = 0`` and
    ``gamma = (s_speed / p_speed)**2``, found by bisection to machine precision.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import rayleigh_wave_speed
    >>> round(float(rayleigh_wave_speed(1.7320508075688772, 1)), 6)
    0.919402
    """
    p_speed = quantity("p_speed", p_speed, lower=0, strict_lower=True)
    s_speed = quantity("s_speed", s_speed, lower=0, strict_lower=True)
    require(p_speed**2 > 4 * s_speed**2 / 3, "wave speeds imply nonpositive bulk modulus")
    gamma = (s_speed / p_speed) ** 2
    # The cubic is -16(1 - gamma) < 0 at xi = 0 and 1 at xi = 1, with exactly one root between.
    lower, upper = np.zeros_like(gamma), np.ones_like(gamma)
    for _ in range(64):
        middle = (lower + upper) / 2
        below = middle**3 - 8 * middle**2 + (24 - 16 * gamma) * middle - 16 * (1 - gamma) < 0
        lower, upper = np.where(below, middle, lower), np.where(below, upper, middle)
    return s_speed * np.sqrt((lower + upper) / 2)


@model
def love_wave_angular_frequency(phase_speed, thickness, layer_s_speed, halfspace_s_speed, layer_rigidity, halfspace_rigidity):
    """Fundamental Love-wave dispersion in a layer over a half-space.

    Shear waves trapped in a slow surface layer travel at a speed that depends on their period.

    Parameters
    ----------
    phase_speed : float or array_like
        m/s; constraint: strictly between the layer and half-space shear speeds.
    thickness : float or array_like
        m; constraint: positive.
    layer_s_speed : float or array_like
        m/s; constraint: positive.
    halfspace_s_speed : float or array_like
        m/s; constraint: positive.
    layer_rigidity : float or array_like
        Pa; constraint: positive.
    halfspace_rigidity : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Fundamental mode of a uniform elastic layer over a faster uniform half-space.
    Solves ``tan(k * thickness * s1) = halfspace_rigidity * s2 / (layer_rigidity * s1)`` for the
    wavenumber k at the given phase speed, with ``s1 = np.sqrt(phase_speed**2 / layer_s_speed**2 - 1)``
    and ``s2 = np.sqrt(1 - phase_speed**2 / halfspace_s_speed**2)``.
    Equation: ``phase_speed * np.arctan(halfspace_rigidity * s2 / (layer_rigidity * s1)) / (thickness * s1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import love_wave_angular_frequency
    >>> round(float(love_wave_angular_frequency(1.4142135623730951, 1, 1, 2, 1, 1)), 6)
    0.87042
    """
    phase_speed = quantity("phase_speed", phase_speed, lower=0, strict_lower=True)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    layer_s_speed = quantity("layer_s_speed", layer_s_speed, lower=0, strict_lower=True)
    halfspace_s_speed = quantity("halfspace_s_speed", halfspace_s_speed, lower=0, strict_lower=True)
    layer_rigidity = quantity("layer_rigidity", layer_rigidity, lower=0, strict_lower=True)
    halfspace_rigidity = quantity("halfspace_rigidity", halfspace_rigidity, lower=0, strict_lower=True)
    require(phase_speed > layer_s_speed, "phase speed must exceed the layer shear speed")
    require(phase_speed < halfspace_s_speed, "phase speed must be below the half-space shear speed")
    s1 = np.sqrt(phase_speed**2 / layer_s_speed**2 - 1)
    s2 = np.sqrt(1 - phase_speed**2 / halfspace_s_speed**2)
    return phase_speed * np.arctan(halfspace_rigidity * s2 / (layer_rigidity * s1)) / (thickness * s1)


@model
def ricker_wavelet(time, frequency):
    """Ricker seismic source pulse.

    A compact pulse makes propagation and filtering examples inspectable.

    Parameters
    ----------
    time : float or array_like
        s; constraint: real.
    frequency : float or array_like
        Hz; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Unit peak at zero time; frequency is the peak-spectrum convention.
    Equation: ``(1 - 2 * (np.pi * frequency * time)**2) * np.exp(-(np.pi * frequency * time)**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import ricker_wavelet
    >>> float(ricker_wavelet(0, 10))
    1.0
    """
    time = quantity("time", time)
    frequency = quantity("frequency", frequency, lower=0, strict_lower=True)
    return (1 - 2 * (np.pi * frequency * time) ** 2) * np.exp(-((np.pi * frequency * time) ** 2))


@model
def attenuated_amplitude(amplitude, frequency, time, quality_factor):
    """Constant-Q wave attenuation.

    Finite quality factor models intrinsic amplitude decay during propagation.

    Parameters
    ----------
    amplitude : float or array_like
        arbitrary amplitude; constraint: real.
    frequency : float or array_like
        Hz; constraint: nonnegative.
    time : float or array_like
        s; constraint: nonnegative.
    quality_factor : float or array_like
        1; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant-Q narrowband attenuation; geometric spreading is separate.
    Equation: ``amplitude * np.exp(-np.pi * frequency * time / quality_factor)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.seismology import attenuated_amplitude
    >>> round(float(attenuated_amplitude(2, 1, 1, 4.532360141827194)), 5)
    1.0
    """
    amplitude = quantity("amplitude", amplitude)
    frequency = quantity("frequency", frequency, lower=0)
    time = quantity("time", time, lower=0)
    quality_factor = quantity("quality_factor", quality_factor, lower=0, strict_lower=True)
    return amplitude * np.exp(-np.pi * frequency * time / quality_factor)
