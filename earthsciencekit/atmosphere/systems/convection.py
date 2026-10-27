"""Moist convection: the lifting condensation level, the moist pseudo-adiabat, CAPE, CIN and moist parcel ascent.

SI units unless stated otherwise. Moisture follows :mod:`earthsciencekit.atmosphere.systems.moisture`: a liquid-water
saturation vapour pressure with the constant latent heat ``LATENT_HEAT_VAPORIZATION``.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.atmosphere.systems.moisture import mixing_ratio, saturation_vapor_pressure
from earthsciencekit.constants import CP_DRY_AIR, EPSILON_WATER, GRAVITY, LATENT_HEAT_VAPORIZATION, R_DRY_AIR

__all__ = [
    "lcl_temperature",
    "lcl_pressure",
    "lcl_height",
    "moist_adiabatic_lapse_rate",
    "parcel_temperature_profile",
    "convective_available_potential_energy",
    "convective_inhibition",
    "simulate_moist_ascent",
]


@model
def lcl_temperature(temperature, dewpoint):
    """Temperature at the lifting condensation level.

    Air lifted dry-adiabatically cools faster than its dew point falls, and condenses where the two meet: the base of
    a cumulus cloud.

    Parameters
    ----------
    temperature : float or array_like
        K, of the air before lifting; constraint: at least dewpoint.
    dewpoint : float or array_like
        K; constraint: above 56 K.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Bolton (1980, his equation 15), accurate to about 0.1 K for terrestrial air; it does not depend on pressure.
    Equation: ``1 / (1 / (dewpoint - 56) + np.log(temperature / dewpoint) / 800) + 56``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import lcl_temperature
    >>> round(float(lcl_temperature(300, 290)), 3)
    287.702
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    dewpoint = quantity("dewpoint", dewpoint, lower=56, strict_lower=True)
    require(dewpoint <= temperature, "dewpoint must not exceed the temperature")
    return 1 / (1 / (dewpoint - 56) + np.log(temperature / dewpoint) / 800) + 56


@model
def lcl_pressure(pressure, temperature, lcl_temperature):
    """Pressure at the lifting condensation level.

    Below cloud base the lifted parcel keeps its potential temperature, so the condensation temperature fixes the
    condensation pressure.

    Parameters
    ----------
    pressure : float or array_like
        Pa, of the air before lifting; constraint: positive.
    temperature : float or array_like
        K, of the air before lifting; constraint: positive.
    lcl_temperature : float or array_like
        K; constraint: positive and at most temperature.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Poisson's dry adiabat with constant heat capacity, the inverse of ``potential_temperature``.
    Equation: ``pressure * (lcl_temperature / temperature)**(CP_DRY_AIR / R_DRY_AIR)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import lcl_pressure
    >>> round(float(lcl_pressure(100000.0, 300, 270)), 1)
    69176.1
    """
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    lcl_temperature = quantity("lcl_temperature", lcl_temperature, lower=0, strict_lower=True)
    require(lcl_temperature <= temperature, "lcl_temperature must not exceed the temperature")
    return pressure * (lcl_temperature / temperature) ** (CP_DRY_AIR / R_DRY_AIR)


@model
def lcl_height(temperature, lcl_temperature):
    """Height of the lifting condensation level above the starting level.

    The dry adiabatic lapse rate turns the cooling needed for condensation into a cloud-base height, about 125 m per
    kelvin of dew-point depression.

    Parameters
    ----------
    temperature : float or array_like
        K, of the air before lifting; constraint: positive.
    lcl_temperature : float or array_like
        K; constraint: positive and at most temperature.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dry adiabatic ascent at ``GRAVITY / CP_DRY_AIR`` with no mixing; pass ``lcl_temperature`` from
    :func:`lcl_temperature`.
    Equation: ``CP_DRY_AIR * (temperature - lcl_temperature) / GRAVITY``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import lcl_height
    >>> round(float(lcl_height(300, 290)), 2)
    1023.8
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    lcl_temperature = quantity("lcl_temperature", lcl_temperature, lower=0, strict_lower=True)
    require(lcl_temperature <= temperature, "lcl_temperature must not exceed the temperature")
    return CP_DRY_AIR * (temperature - lcl_temperature) / GRAVITY


@model
def moist_adiabatic_lapse_rate(temperature, pressure):
    """Pseudo-adiabatic lapse rate of saturated rising air.

    Condensation releases latent heat inside a rising cloud, so saturated air cools more slowly with height than dry
    air: about 4 K/km in warm, moist air and close to the dry rate in cold air.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: 180 to 330 K, the range of ``saturation_vapor_pressure``.
    pressure : float or array_like
        Pa; constraint: above the saturation vapour pressure.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K/m], a positive magnitude; temperature decreases upward. Scalar input returns a NumPy scalar;
        arrays broadcast.

    Notes
    -----
    Pseudo-adiabatic: condensate falls out at once, and the heat capacities of vapour and water are neglected (AMS
    Glossary form). ``r_s`` is the saturation mixing ratio from ``saturation_vapor_pressure`` and ``L`` is
    ``LATENT_HEAT_VAPORIZATION``.
    Equation: ``GRAVITY * (1 + L * r_s / (R_DRY_AIR * temperature)) / (CP_DRY_AIR + L**2 * r_s * EPSILON_WATER /
    (R_DRY_AIR * temperature**2))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import moist_adiabatic_lapse_rate
    >>> round(float(moist_adiabatic_lapse_rate(300, 100000.0)), 8)
    0.003647
    """
    temperature = quantity("temperature", temperature, lower=180, upper=330)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    return _moist_lapse(temperature, pressure)


def _moist_lapse(temperature, pressure):
    saturation = mixing_ratio(saturation_vapor_pressure(temperature), pressure)
    latent = LATENT_HEAT_VAPORIZATION * saturation / (R_DRY_AIR * temperature)
    return GRAVITY * (1 + latent) / (CP_DRY_AIR + latent * LATENT_HEAT_VAPORIZATION * EPSILON_WATER / temperature)


@model
def parcel_temperature_profile(pressure, temperature, dewpoint):
    """Temperature of a lifted surface parcel: dry below cloud base, moist pseudo-adiabatic above.

    The parcel curve on a sounding is the reference against which buoyancy, and hence convective energy, is
    measured.

    Parameters
    ----------
    pressure : array_like, shape (n,), n >= 2
        Pa, strictly decreasing, starting at the parcel's level; constraint: positive.
    temperature : float
        K, scalar starting temperature; constraint: positive.
    dewpoint : float
        K, scalar starting dew point; constraint: above 56 K and at most temperature.

    Returns
    -------
    ndarray, shape (n,)
        Model value [K]: the parcel temperature at each pressure.

    Notes
    -----
    Poisson's dry adiabat up to the lifting condensation level of :func:`lcl_temperature` (Bolton 1980), then
    ``dT/d(ln p) = moist_adiabatic_lapse_rate * R_DRY_AIR * T / GRAVITY`` integrated with
    ``integrators.integrate_ode``, using the dry-air hydrostatic thickness. The moist part must stay within 180-330 K.
    Equation: ``T (p / p0)**(R_DRY_AIR / CP_DRY_AIR)`` below cloud base, the pseudo-adiabat above.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.atmosphere import parcel_temperature_profile
    >>> np.round(parcel_temperature_profile([100000.0, 90000.0], 300, 250), 3).tolist()
    [300.0, 291.098]
    """
    from earthsciencekit.integrators import integrate_ode

    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    dewpoint = quantity("dewpoint", dewpoint, lower=56, strict_lower=True)
    require(pressure.ndim == 1 and pressure.size >= 2, "pressure must be a vector with at least two levels")
    require(np.all(np.diff(pressure) < 0), "pressure must strictly decrease")
    require(temperature.ndim == dewpoint.ndim == 0, "temperature and dewpoint must be scalars")
    base_temperature = lcl_temperature(temperature, dewpoint)
    base_pressure = lcl_pressure(pressure[0], temperature, base_temperature)
    profile = temperature * (pressure / pressure[0]) ** (R_DRY_AIR / CP_DRY_AIR)
    moist = pressure < base_pressure
    if np.any(moist):
        log_pressure = -np.log(np.concatenate([[base_pressure], pressure[moist]]))
        result = integrate_ode(
            lambda s, state: _moist_lapse(state, np.exp(-s)) * R_DRY_AIR * state / GRAVITY * -1,
            [base_temperature],
            log_pressure,
            rtol=1e-10,
            atol=1e-8,
        )
        profile[moist] = result.values[1:, 0]
    return profile


def _buoyancy(height, parcel_virtual_temperature, environment_virtual_temperature):
    height = quantity("height", height)
    parcel = quantity("parcel_virtual_temperature", parcel_virtual_temperature, lower=0, strict_lower=True)
    environment = quantity("environment_virtual_temperature", environment_virtual_temperature, lower=0, strict_lower=True)
    require(height.ndim == 1 and height.size >= 2, "height must be a vector with at least two levels")
    require(parcel.shape == environment.shape == height.shape, "temperatures must have the same shape as height")
    require(np.all(np.diff(height) > 0), "height must strictly increase")
    return GRAVITY * (parcel - environment) / environment, np.diff(height)


def _positive_area(buoyancy, thickness):
    """Exact integral of the positive part of a piecewise-linear profile over each layer."""
    lower, upper = buoyancy[:-1], buoyancy[1:]
    same_sign = lower * upper >= 0
    span = np.where(same_sign, 1.0, np.abs(lower) + np.abs(upper))
    crossing = np.maximum(lower, upper) ** 2 / (2 * span)
    return thickness * np.where(same_sign, (np.maximum(lower, 0) + np.maximum(upper, 0)) / 2, crossing)


@model
def convective_available_potential_energy(height, parcel_virtual_temperature, environment_virtual_temperature):
    """Convective available potential energy (CAPE) of a sounding.

    The work buoyancy can do on a lifted parcel bounds the updraught speed of a thunderstorm, about
    ``np.sqrt(2 * CAPE)``.

    Parameters
    ----------
    height : array_like, shape (n,), n >= 2
        m, strictly increasing; constraint: real.
    parcel_virtual_temperature : array_like, shape (n,)
        K, of the lifted parcel at each height; constraint: positive.
    environment_virtual_temperature : array_like, shape (n,)
        K, of the environment at each height; constraint: positive.

    Returns
    -------
    numpy scalar
        Model value [J/kg], nonnegative.

    Notes
    -----
    Moncrieff and Miller (1976): the integral of positive buoyancy ``GRAVITY * (Tv_parcel - Tv_env) / Tv_env``, taken
    exactly for buoyancy linear between levels, over every positively buoyant layer. Pass virtual temperatures to
    include the buoyancy of vapour; condensate loading and entrainment are neglected.
    Equation: ``sum(integral of max(GRAVITY * (Tv_parcel - Tv_env) / Tv_env, 0) dz)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import convective_available_potential_energy
    >>> round(float(convective_available_potential_energy([0, 1000, 2000], [300, 302, 300], [300, 300, 300])), 4)
    65.3777
    """
    buoyancy, thickness = _buoyancy(height, parcel_virtual_temperature, environment_virtual_temperature)
    return np.sum(_positive_area(buoyancy, thickness))


@model
def convective_inhibition(height, parcel_virtual_temperature, environment_virtual_temperature):
    """Convective inhibition (CIN) below the level of free convection.

    A stable layer above the surface can hold back a very unstable sounding until lifting supplies this much energy.

    Parameters
    ----------
    height : array_like, shape (n,), n >= 2
        m, strictly increasing from the parcel's starting level; constraint: real.
    parcel_virtual_temperature : array_like, shape (n,)
        K; constraint: positive.
    environment_virtual_temperature : array_like, shape (n,)
        K; constraint: positive.

    Returns
    -------
    numpy scalar
        Model value [J/kg], nonpositive.

    Notes
    -----
    The negative buoyancy integrated from the first level up to the level of free convection, where the parcel first
    becomes buoyant (over the whole profile if it never does), exactly for buoyancy linear between levels.
    Equation: ``-sum(integral of max(-GRAVITY * (Tv_parcel - Tv_env) / Tv_env, 0) dz)`` below the level of free
    convection.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import convective_inhibition
    >>> float(convective_inhibition([0, 1000, 2000], [300, 299, 301], [300, 300, 300]))
    -24.516625
    """
    buoyancy, thickness = _buoyancy(height, parcel_virtual_temperature, environment_virtual_temperature)
    buoyant = np.flatnonzero(buoyancy > 0)
    layers = buoyant[0] if buoyant.size else thickness.size
    return -np.sum(_positive_area(-buoyancy, thickness)[:layers])


@safe_arithmetic
def simulate_moist_ascent(initial_temperature, initial_dewpoint, initial_pressure, vertical_speed, time):
    """Follow a rising parcel that condenses its vapour above the lifting condensation level.

    Parameters
    ----------
    initial_temperature : float
        Scalar starting temperature [K].
    initial_dewpoint : float
        Scalar starting dew point, at most the temperature [K].
    initial_pressure : float
        Positive scalar starting pressure [Pa].
    vertical_speed : float
        Positive scalar ascent rate [m/s].
    time : array_like
        Nonnegative strictly increasing sampled times [s].

    Returns
    -------
    TimeSeries
        Parcel temperature [K], pressure [Pa], vapour mixing ratio [kg/kg] and the water condensed so far [kg/kg].

    Notes
    -----
    The moist successor of ``simulate_dry_ascent``: hydrostatic pressure with the dry-air density and constant
    standard gravity, the dry lapse rate while the vapour mixing ratio is below saturation and
    :func:`moist_adiabatic_lapse_rate` once it reaches it. The dry ascent is evaluated in closed form and the moment of
    saturation found by root finding, so the integrator never steps across the kink at cloud base. Condensate is
    removed at once (pseudo-adiabatic), and there is no entrainment or ice. The saturation formula restricts the
    parcel to 180-330 K.
    """
    from scipy.optimize import brentq

    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode
    from earthsciencekit.results import TimeSeries

    temperature = quantity("initial_temperature", initial_temperature, lower=180, upper=330)
    dewpoint = quantity("initial_dewpoint", initial_dewpoint, lower=180, upper=330)
    pressure = quantity("initial_pressure", initial_pressure, lower=0, strict_lower=True)
    speed = quantity("vertical_speed", vertical_speed, lower=0, strict_lower=True)
    require(all(v.ndim == 0 for v in (temperature, dewpoint, pressure, speed)), "parcel parameters must be scalars")
    require(dewpoint <= temperature, "initial_dewpoint must not exceed the initial temperature")
    time = time_grid(time)
    vapor = mixing_ratio(saturation_vapor_pressure(dewpoint), pressure)
    dry_rate = GRAVITY / CP_DRY_AIR * speed
    elapsed = time - time[0]

    def dry(tau):
        cooled = temperature - dry_rate * tau
        return cooled, pressure * (cooled / temperature) ** (CP_DRY_AIR / R_DRY_AIR)

    def excess(tau):
        return mixing_ratio(saturation_vapor_pressure(dry(tau)[0]), dry(tau)[1]) - vapor

    # The dry ascent is analytic; the parcel saturates where its saturation mixing ratio falls to its own.
    coldest = (temperature - 180) / dry_rate
    if excess(0) <= 0:
        saturation = 0.0
    elif excess(coldest) > 0:
        saturation = np.inf
    else:
        saturation = brentq(excess, 0, coldest, xtol=1e-12, rtol=1e-14)
    state = np.column_stack(dry(np.minimum(elapsed, saturation)))
    moist = elapsed > saturation
    if np.any(moist):
        result = integrate_ode(
            lambda t, y: [-_moist_lapse(y[0], y[1]) * speed, -GRAVITY * y[1] * speed / (R_DRY_AIR * y[0])],
            np.array(dry(saturation)),
            np.concatenate([[saturation], elapsed[moist]]),
            names=("temperature", "pressure"),
            units=("K", "Pa"),
        )
        state[moist] = result.values[1:]
    remaining = np.minimum(vapor, mixing_ratio(saturation_vapor_pressure(state[:, 0]), state[:, 1]))
    values = np.column_stack([state, remaining, vapor - remaining])
    names = ("temperature", "pressure", "mixing_ratio", "condensed_water")
    return TimeSeries(time, values, names, ("K", "Pa", "kg/kg", "kg/kg"))
