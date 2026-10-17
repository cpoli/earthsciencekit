"""Cooling and convection: half-space cooling, heat flux, ocean depth, mantle Rayleigh number and periodic diffusion.

SI units unless stated otherwise.
"""

import numpy as np
from scipy.special import erf

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "halfspace_temperature",
    "halfspace_heat_flux",
    "mantle_rayleigh_number",
    "ocean_depth",
    "diffuse_temperature",
]


@model
def halfspace_temperature(depth, age, diffusivity, surface_temperature, mantle_temperature):
    """Cooling oceanic lithosphere.

    Oceanic plates cool as their age and distance from a ridge increase.

    Parameters
    ----------
    depth : float or array_like
        m; constraint: nonnegative.
    age : float or array_like
        s; constraint: positive.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    surface_temperature : float or array_like
        K; constraint: positive.
    mantle_temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Semi-infinite conduction with fixed surface temperature; no hydrothermal or advective heat transfer.
    Equation: ``surface_temperature + (mantle_temperature - surface_temperature) * erf(depth / (2 * np.sqrt(diffusivity * age)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    age = quantity("age", age, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    surface_temperature = quantity("surface_temperature", surface_temperature, lower=0, strict_lower=True)
    mantle_temperature = quantity("mantle_temperature", mantle_temperature, lower=0, strict_lower=True)
    require(mantle_temperature >= surface_temperature, "mantle must be at least as warm as surface")
    return surface_temperature + (mantle_temperature - surface_temperature) * erf(depth / (2 * np.sqrt(diffusivity * age)))


@model
def halfspace_heat_flux(age, diffusivity, conductivity, surface_temperature, mantle_temperature):
    """Cooling-plate surface heat flux.

    Surface heat flux decreases as an oceanic plate cools.

    Parameters
    ----------
    age : float or array_like
        s; constraint: positive.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    conductivity : float or array_like
        W/(m K); constraint: positive.
    surface_temperature : float or array_like
        K; constraint: positive.
    mantle_temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Positive upward flux; semi-infinite conductive model excludes the zero-age singularity.
    Equation: ``conductivity * (mantle_temperature - surface_temperature) / np.sqrt(np.pi * diffusivity * age)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    age = quantity("age", age, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    conductivity = quantity("conductivity", conductivity, lower=0, strict_lower=True)
    surface_temperature = quantity("surface_temperature", surface_temperature, lower=0, strict_lower=True)
    mantle_temperature = quantity("mantle_temperature", mantle_temperature, lower=0, strict_lower=True)
    require(mantle_temperature >= surface_temperature, "mantle must be at least as warm as surface")
    return conductivity * (mantle_temperature - surface_temperature) / np.sqrt(np.pi * diffusivity * age)


@model
def ocean_depth(age, ridge_depth, mantle_density, water_density, thermal_expansivity, diffusivity, surface_temperature, mantle_temperature):
    """Ocean depth of a cooling half-space.

    Thermal contraction of the cooling plate makes the sea floor deepen with the square root of its age.

    Parameters
    ----------
    age : float or array_like
        s; constraint: nonnegative.
    ridge_depth : float or array_like
        m; constraint: nonnegative.
    mantle_density : float or array_like
        kg/m^3; constraint: positive.
    water_density : float or array_like
        kg/m^3; constraint: positive.
    thermal_expansivity : float or array_like
        1/K; constraint: positive.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    surface_temperature : float or array_like
        K; constraint: positive.
    mantle_temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Isostatic half-space cooling, depth positive downward; it overestimates depth on sea floor older than about 70 million years.
    Equation: ``ridge_depth + 2 * mantle_density * thermal_expansivity * (mantle_temperature - surface_temperature)
    / (mantle_density - water_density) * np.sqrt(diffusivity * age / np.pi)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    age = quantity("age", age, lower=0)
    ridge_depth = quantity("ridge_depth", ridge_depth, lower=0)
    mantle_density = quantity("mantle_density", mantle_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    thermal_expansivity = quantity("thermal_expansivity", thermal_expansivity, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    surface_temperature = quantity("surface_temperature", surface_temperature, lower=0, strict_lower=True)
    mantle_temperature = quantity("mantle_temperature", mantle_temperature, lower=0, strict_lower=True)
    require(mantle_density > water_density, "mantle must be denser than water")
    require(mantle_temperature >= surface_temperature, "mantle must be at least as warm as surface")
    subsidence = 2 * mantle_density * thermal_expansivity * (mantle_temperature - surface_temperature) / (mantle_density - water_density)
    return ridge_depth + subsidence * np.sqrt(diffusivity * age / np.pi)


@model
def mantle_rayleigh_number(thermal_expansivity, density, gravity, temperature_difference, thickness, diffusivity, viscosity):
    """Rayleigh number of a convecting mantle layer.

    Buoyancy against viscous and thermal diffusion decides whether the mantle convects.

    Parameters
    ----------
    thermal_expansivity : float or array_like
        1/K; constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    temperature_difference : float or array_like
        K; constraint: nonnegative.
    thickness : float or array_like
        m; constraint: positive.
    diffusivity : float or array_like
        m^2/s; constraint: positive.
    viscosity : float or array_like
        Pa s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Layer heated from below with uniform properties; convection starts above about 657.5 (free slip) or 1708 (rigid walls).
    Equation: ``thermal_expansivity * density * gravity * temperature_difference * thickness**3 / (diffusivity * viscosity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thermal_expansivity = quantity("thermal_expansivity", thermal_expansivity, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    temperature_difference = quantity("temperature_difference", temperature_difference, lower=0)
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    viscosity = quantity("viscosity", viscosity, lower=0, strict_lower=True)
    return thermal_expansivity * density * gravity * temperature_difference * thickness**3 / (diffusivity * viscosity)


@safe_arithmetic
def diffuse_temperature(initial, time, diffusivity, spacing):
    """Evolve a periodic one-dimensional conductive temperature grid.

    Parameters
    ----------
    initial : array_like, shape (m,), m >= 3
        Positive temperatures [K] on an equally spaced periodic grid.
    time : array_like
        Strictly increasing nonnegative output times [s]. Initial state belongs
        to the first supplied time; values need not start at zero.
    diffusivity : float
        Positive scalar thermal diffusivity [m^2/s].
    spacing : float
        Positive scalar grid spacing [m].

    Returns
    -------
    TimeSeries
        Temperature at each grid point [K]. Columns are labeled cell_0, ... .

    Notes
    -----
    The centered discrete Laplacian is diagonalized with a real FFT. Time
    evolution of that discrete system is exact to floating-point accuracy.
    Periodicity conserves mean temperature. Spatial discretization error remains;
    this is not the fixed-surface half-space cooling boundary condition.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.plate_tectonics.core.base import TimeSeries

    initial = quantity("initial", initial, lower=0, strict_lower=True)
    require(initial.ndim == 1 and initial.size >= 3, "at least three temperature cells are required")
    time = time_grid(time)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    spacing = quantity("spacing", spacing, lower=0, strict_lower=True)
    require(diffusivity.ndim == spacing.ndim == 0, "diffusivity and spacing must be scalars")
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        eigenvalues = -4 * diffusivity * np.sin(np.pi * np.fft.rfftfreq(initial.size)) ** 2 / spacing**2
        spectra = np.fft.rfft(initial)[None, :] * np.exp((time - time[0])[:, None] * eigenvalues)
        values = np.fft.irfft(spectra, n=initial.size, axis=1)
    return TimeSeries(time, values, tuple(f"cell_{i}" for i in range(initial.size)), ("K",) * initial.size)
