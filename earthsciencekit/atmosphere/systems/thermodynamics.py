"""Dry thermodynamics: potential temperature, the dry adiabat, static stability, the Richardson number and parcel ascent.

SI units unless stated otherwise.
"""

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.constants import CP_DRY_AIR, R_DRY_AIR

__all__ = [
    "potential_temperature",
    "temperature_from_potential",
    "dry_adiabatic_lapse_rate",
    "brunt_vaisala_squared",
    "richardson_number",
    "simulate_dry_ascent",
]


@model
def potential_temperature(temperature, pressure, reference_pressure):
    """Poisson potential temperature.

    Adiabatic temperature adjustment permits comparisons across pressure levels.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    pressure : float or array_like
        Pa; constraint: positive.
    reference_pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dry ideal gas with constant heat capacity; not equivalent potential temperature.
    Equation: ``temperature * (reference_pressure / pressure)**(R_DRY_AIR / CP_DRY_AIR)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    reference_pressure = quantity("reference_pressure", reference_pressure, lower=0, strict_lower=True)
    return temperature * (reference_pressure / pressure) ** (R_DRY_AIR / CP_DRY_AIR)


@model
def temperature_from_potential(potential_temperature, pressure, reference_pressure):
    """Inverse potential temperature.

    Potential-temperature inversion recovers the local dry-air temperature.

    Parameters
    ----------
    potential_temperature : float or array_like
        K; constraint: positive.
    pressure : float or array_like
        Pa; constraint: positive.
    reference_pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uses the same dry-air heat capacity and reference pressure convention.
    Equation: ``potential_temperature * (pressure / reference_pressure)**(R_DRY_AIR / CP_DRY_AIR)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    potential_temperature = quantity("potential_temperature", potential_temperature, lower=0, strict_lower=True)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    reference_pressure = quantity("reference_pressure", reference_pressure, lower=0, strict_lower=True)
    return potential_temperature * (pressure / reference_pressure) ** (R_DRY_AIR / CP_DRY_AIR)


@model
def dry_adiabatic_lapse_rate(gravity, heat_capacity):
    """Dry adiabatic lapse rate.

    Rising dry parcels cool as they perform expansion work.

    Parameters
    ----------
    gravity : float or array_like
        m/s^2; constraint: positive.
    heat_capacity : float or array_like
        J/(kg K); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Returned positive lapse magnitude in K per metre; temperature decreases upward.
    Equation: ``gravity / heat_capacity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    heat_capacity = quantity("heat_capacity", heat_capacity, lower=0, strict_lower=True)
    return gravity / heat_capacity


@model
def brunt_vaisala_squared(potential_temperature, vertical_gradient, gravity):
    """Dry static stability.

    A potential-temperature gradient determines dry stratification stability.

    Parameters
    ----------
    potential_temperature : float or array_like
        K; constraint: positive.
    vertical_gradient : float or array_like
        K/m; constraint: real.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Returns signed N squared; negative values represent instability rather than an imaginary frequency.
    Equation: ``gravity * vertical_gradient / potential_temperature``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    potential_temperature = quantity("potential_temperature", potential_temperature, lower=0, strict_lower=True)
    vertical_gradient = quantity("vertical_gradient", vertical_gradient)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return gravity * vertical_gradient / potential_temperature


@model
def richardson_number(brunt_vaisala_squared, shear):
    """Gradient Richardson number.

    Stratification suppresses the turbulence that vertical wind shear generates.

    Parameters
    ----------
    brunt_vaisala_squared : float or array_like
        1/s^2; constraint: real.
    shear : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Signed; stratified shear flow can become unstable only where it falls below 1/4 (Miles-Howard).
    Equation: ``brunt_vaisala_squared / shear**2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    brunt_vaisala_squared = quantity("brunt_vaisala_squared", brunt_vaisala_squared)
    shear = quantity("shear", shear)
    require(shear != 0, "shear must be nonzero")
    return brunt_vaisala_squared / shear**2


@safe_arithmetic
def simulate_dry_ascent(initial_temperature, initial_pressure, vertical_speed, time):
    """Follow an ideal dry parcel along a hydrostatic constant-speed ascent.

    Parameters
    ----------
    initial_temperature : float
        Positive scalar initial temperature [K].
    initial_pressure : float
        Positive scalar initial pressure [Pa].
    vertical_speed : float
        Finite scalar vertical velocity [m/s], positive upward.
    time : array_like
        Nonnegative strictly increasing sampled times [s].

    Returns
    -------
    TimeSeries
        Parcel temperature [K] and pressure [Pa].

    Notes
    -----
    Dry lapse and hydrostatic pressure equations use constant standard gravity
    and dry-air heat capacity. Potential temperature should remain constant.
    The parcel is unsaturated; condensation and environmental entrainment are
    excluded. Reject trajectories whose predicted temperature reaches zero.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.constants import GRAVITY
    from earthsciencekit.integrators import integrate_ode

    temperature = quantity("initial_temperature", initial_temperature, lower=0, strict_lower=True)
    pressure = quantity("initial_pressure", initial_pressure, lower=0, strict_lower=True)
    speed = quantity("vertical_speed", vertical_speed)
    require(temperature.ndim == pressure.ndim == speed.ndim == 0, "parcel parameters must be scalars")
    time = time_grid(time)
    require(temperature - GRAVITY / CP_DRY_AIR * speed * (time[-1] - time[0]) > 0, "trajectory reaches nonpositive temperature")
    return integrate_ode(
        lambda t, y: [-GRAVITY / CP_DRY_AIR * speed, -GRAVITY * y[1] * speed / (R_DRY_AIR * y[0])],
        [temperature, pressure],
        time,
        names=("temperature", "pressure"),
        units=("K", "Pa"),
    )
