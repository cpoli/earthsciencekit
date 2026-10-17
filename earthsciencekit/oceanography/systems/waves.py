"""Waves and tides: surface and internal wave dispersion, trochoidal orbits and the equilibrium tide.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "shallow_wave_speed",
    "deep_wave_speed",
    "wave_angular_frequency",
    "deep_group_speed",
    "reduced_gravity",
    "internal_wave_speed",
    "tidal_elevation",
    "gerstner_orbit_radius",
]


@model
def shallow_wave_speed(depth, gravity):
    """Shallow-water gravity wave.

    Long surface waves respond to water depth through gravity.

    Parameters
    ----------
    depth : float or array_like
        m; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Wavelength much greater than depth, small amplitude, no dispersion.
    Equation: ``np.sqrt(gravity * depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return np.sqrt(gravity * depth)


@model
def deep_wave_speed(wavelength, gravity):
    """Deep-water phase speed.

    Short waves relative to depth propagate dispersively.

    Parameters
    ----------
    wavelength : float or array_like
        m; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Deep-water small-amplitude gravity waves without capillarity.
    Equation: ``np.sqrt(gravity * wavelength / (2 * np.pi))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wavelength = quantity("wavelength", wavelength, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return np.sqrt(gravity * wavelength / (2 * np.pi))


@model
def wave_angular_frequency(wavenumber, depth, gravity):
    """Finite-depth dispersion relation.

    The dispersion relation bridges shallow and deep gravity waves.

    Parameters
    ----------
    wavenumber : float or array_like
        1/m; constraint: positive.
    depth : float or array_like
        m; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Linear inviscid surface gravity wave, no current or surface tension.
    Equation: ``np.sqrt(gravity * wavenumber * np.tanh(wavenumber * depth))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wavenumber = quantity("wavenumber", wavenumber, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return np.sqrt(gravity * wavenumber * np.tanh(wavenumber * depth))


@model
def deep_group_speed(wavelength, gravity):
    """Deep-water group speed.

    Wave energy travels at group speed rather than crest phase speed.

    Parameters
    ----------
    wavelength : float or array_like
        m; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Deep-water gravity waves; group speed is half phase speed.
    Equation: ``0.5 * np.sqrt(gravity * wavelength / (2 * np.pi))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wavelength = quantity("wavelength", wavelength, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return 0.5 * np.sqrt(gravity * wavelength / (2 * np.pi))


@model
def gerstner_orbit_radius(rest_depth, amplitude, wavenumber):
    """Particle-orbit radius in a Gerstner trochoidal wave.

    Water particles in a deep-water wave move on circles that shrink with depth.

    Parameters
    ----------
    rest_depth : float or array_like
        m; constraint: nonnegative. Depth of the particle's orbit centre below the mean surface.
    amplitude : float or array_like
        m; constraint: nonnegative. Orbit radius at the surface.
    wavenumber : float or array_like
        1/m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Exact nonlinear deep-water solution with closed orbits and no Stokes drift; the flow has vorticity.
    ``amplitude * wavenumber`` must not exceed 1, where the surface trochoid forms a cusp.
    A particle at rest position ``(a, b)`` moves to ``(a - r sin(k a - omega t), -b + r cos(k a - omega t))``
    with ``omega = sqrt(g k)``.
    Equation: ``amplitude * np.exp(-wavenumber * rest_depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    rest_depth = quantity("rest_depth", rest_depth, lower=0)
    amplitude = quantity("amplitude", amplitude, lower=0)
    wavenumber = quantity("wavenumber", wavenumber, lower=0, strict_lower=True)
    require(amplitude * wavenumber <= 1, "a trochoid steeper than amplitude * wavenumber = 1 overturns")
    return amplitude * np.exp(-wavenumber * rest_depth)


@model
def reduced_gravity(density_difference, reference_density, gravity):
    """Two-layer reduced gravity.

    A small density contrast slows internal waves relative to surface waves.

    Parameters
    ----------
    density_difference : float or array_like
        kg/m^3; constraint: nonnegative.
    reference_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Boussinesq small density contrast; output is acceleration.
    Equation: ``gravity * density_difference / reference_density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    density_difference = quantity("density_difference", density_difference, lower=0)
    reference_density = quantity("reference_density", reference_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return gravity * density_difference / reference_density


@model
def internal_wave_speed(reduced_gravity, upper_depth, lower_depth):
    """Two-layer internal long wave.

    Layer depths and buoyancy jointly set internal long-wave speed.

    Parameters
    ----------
    reduced_gravity : float or array_like
        m/s^2; constraint: positive.
    upper_depth : float or array_like
        m; constraint: positive.
    lower_depth : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hydrostatic two-layer fluid under a rigid lid.
    Equation: ``np.sqrt(reduced_gravity * upper_depth * lower_depth / (upper_depth + lower_depth))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    reduced_gravity = quantity("reduced_gravity", reduced_gravity, lower=0, strict_lower=True)
    upper_depth = quantity("upper_depth", upper_depth, lower=0, strict_lower=True)
    lower_depth = quantity("lower_depth", lower_depth, lower=0, strict_lower=True)
    return np.sqrt(reduced_gravity * upper_depth * lower_depth / (upper_depth + lower_depth))


@model
def tidal_elevation(time, amplitude, period, phase):
    """Single harmonic tide.

    Harmonic constituents provide an inspectable tide representation.

    Parameters
    ----------
    time : float or array_like
        s; constraint: real.
    amplitude : float or array_like
        m; constraint: nonnegative.
    period : float or array_like
        s; constraint: positive.
    phase : float or array_like
        rad; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    One supplied constituent; not a site-specific tide prediction or navigational product.
    Equation: ``amplitude * np.cos(2 * np.pi * time / period + phase)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    time = quantity("time", time)
    amplitude = quantity("amplitude", amplitude, lower=0)
    period = quantity("period", period, lower=0, strict_lower=True)
    phase = quantity("phase", phase)
    return amplitude * np.cos(2 * np.pi * time / period + phase)
