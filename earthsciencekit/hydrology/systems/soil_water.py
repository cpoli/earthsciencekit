"""Soil water: unsaturated Darcy flow, steady Richards profiles and the Green-Ampt and Horton infiltration laws.

SI units unless stated otherwise. Heights and fluxes are positive upward in the unsaturated zone.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.results import TimeSeries

__all__ = [
    "gardner_conductivity",
    "unsaturated_darcy_flux",
    "richards_steady_head",
    "green_ampt_rate",
    "green_ampt_time",
    "horton_infiltration_rate",
    "horton_cumulative_infiltration",
    "simulate_infiltration",
]


@model
def gardner_conductivity(pressure_head, saturated_conductivity, sorptive_number):
    """Gardner's exponential unsaturated hydraulic conductivity.

    Drying soil conducts water through ever thinner films and smaller pores, so its conductivity falls by orders of
    magnitude as suction rises.

    Parameters
    ----------
    pressure_head : float or array_like
        m, negative under suction; constraint: nonpositive.
    saturated_conductivity : float or array_like
        m/s; constraint: positive.
    sorptive_number : float or array_like
        1/m, Gardner's alpha (about 1 for clay to 10 for sand); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Gardner (1958), a one-parameter fit without hysteresis; it makes the steady Richards equation solvable in closed form.
    Equation: ``saturated_conductivity * np.exp(sorptive_number * pressure_head)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    pressure_head = quantity("pressure_head", pressure_head, upper=0)
    saturated_conductivity = quantity("saturated_conductivity", saturated_conductivity, lower=0, strict_lower=True)
    sorptive_number = quantity("sorptive_number", sorptive_number, lower=0, strict_lower=True)
    return saturated_conductivity * np.exp(sorptive_number * pressure_head)


@model
def unsaturated_darcy_flux(conductivity, pressure_head_gradient):
    """Vertical Darcy-Buckingham flux in an unsaturated soil.

    Richards extended Darcy's law to unsaturated soil, where capillary suction and gravity together drive the water.

    Parameters
    ----------
    conductivity : float or array_like
        m/s, at the local pressure head; constraint: nonnegative.
    pressure_head_gradient : float or array_like
        1, d(pressure_head)/dz with z upward; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], positive upward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Richards (1931): the total head is ``pressure_head + z`` and air stays at atmospheric pressure.
    Equation: ``-conductivity * (pressure_head_gradient + 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    conductivity = quantity("conductivity", conductivity, lower=0)
    pressure_head_gradient = quantity("pressure_head_gradient", pressure_head_gradient)
    return -conductivity * (pressure_head_gradient + 1)


@model
def richards_steady_head(height, flux, saturated_conductivity, sorptive_number):
    """Steady pressure-head profile above a water table in a Gardner soil.

    With a constant flux through the unsaturated zone, Richards's equation fixes the suction at every height: capillary
    rise for evaporation, a unit-gradient tail for infiltration.

    Parameters
    ----------
    height : float or array_like
        m, above the water table; constraint: nonnegative.
    flux : float or array_like
        m/s, positive upward (evaporation), negative downward (infiltration); constraint: greater than
        ``-saturated_conductivity``.
    saturated_conductivity : float or array_like
        m/s; constraint: positive.
    sorptive_number : float or array_like
        1/m, Gardner's alpha; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], zero at the water table and negative above it. Scalar input returns a NumPy scalar; arrays
        broadcast.

    Notes
    -----
    Richards (1931) steady equation ``flux = -K(psi) (dpsi/dz + 1)`` with Gardner's (1958) conductivity; an upward flux
    can only be sustained below ``np.log(1 + saturated_conductivity / flux) / sorptive_number``.
    Equation: ``np.log(((flux + saturated_conductivity) * np.exp(-sorptive_number * height) - flux) /
    saturated_conductivity) / sorptive_number``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    height = quantity("height", height, lower=0)
    flux = quantity("flux", flux)
    saturated_conductivity = quantity("saturated_conductivity", saturated_conductivity, lower=0, strict_lower=True)
    sorptive_number = quantity("sorptive_number", sorptive_number, lower=0, strict_lower=True)
    require(flux > -saturated_conductivity, "downward flux must be below the saturated conductivity")
    conductivity = (flux + saturated_conductivity) * np.exp(-sorptive_number * height) - flux
    require(conductivity > 0, "height exceeds the limit of steady capillary rise")
    return np.log(conductivity / saturated_conductivity) / sorptive_number


@model
def green_ampt_rate(cumulative_infiltration, saturated_conductivity, suction_head, moisture_deficit):
    """Green-Ampt infiltration capacity of a ponded soil.

    A sharp wetting front moves down behind a ponded surface, pulled by capillary suction that matters less and less as the
    wetted zone deepens.

    Parameters
    ----------
    cumulative_infiltration : float or array_like
        m, depth infiltrated so far; constraint: positive.
    saturated_conductivity : float or array_like
        m/s; constraint: positive.
    suction_head : float or array_like
        m, capillary suction at the wetting front; constraint: positive.
    moisture_deficit : float or array_like
        1, saturated minus initial water content; constraint: between 0 (exclusive) and 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Green and Ampt (1911): a piston-like front in a uniform soil with a negligible ponding depth; the rate tends to the
    saturated conductivity.
    Equation: ``saturated_conductivity * (1 + suction_head * moisture_deficit / cumulative_infiltration)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    cumulative_infiltration = quantity("cumulative_infiltration", cumulative_infiltration, lower=0, strict_lower=True)
    saturated_conductivity = quantity("saturated_conductivity", saturated_conductivity, lower=0, strict_lower=True)
    suction_head = quantity("suction_head", suction_head, lower=0, strict_lower=True)
    moisture_deficit = quantity("moisture_deficit", moisture_deficit, lower=0, upper=1, strict_lower=True)
    return saturated_conductivity * (1 + suction_head * moisture_deficit / cumulative_infiltration)


@model
def green_ampt_time(cumulative_infiltration, saturated_conductivity, suction_head, moisture_deficit):
    """Time for a ponded Green-Ampt soil to absorb a given depth.

    Integrating the Green-Ampt rate gives the time in closed form, so the infiltrated depth follows by inverting it.

    Parameters
    ----------
    cumulative_infiltration : float or array_like
        m; constraint: nonnegative.
    saturated_conductivity : float or array_like
        m/s; constraint: positive.
    suction_head : float or array_like
        m; constraint: positive.
    moisture_deficit : float or array_like
        1; constraint: between 0 (exclusive) and 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s], since ponding at the surface began. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Green and Ampt (1911), ponded from time zero; early on the depth grows as the square root of time.
    Equation: ``(cumulative_infiltration - suction_head * moisture_deficit * np.log(1 + cumulative_infiltration /
    (suction_head * moisture_deficit))) / saturated_conductivity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    cumulative_infiltration = quantity("cumulative_infiltration", cumulative_infiltration, lower=0)
    saturated_conductivity = quantity("saturated_conductivity", saturated_conductivity, lower=0, strict_lower=True)
    suction_head = quantity("suction_head", suction_head, lower=0, strict_lower=True)
    moisture_deficit = quantity("moisture_deficit", moisture_deficit, lower=0, upper=1, strict_lower=True)
    storage = suction_head * moisture_deficit
    return (cumulative_infiltration - storage * np.log1p(cumulative_infiltration / storage)) / saturated_conductivity


@model
def horton_infiltration_rate(time, initial_rate, final_rate, decay_constant):
    """Horton's infiltration capacity during a storm.

    Infiltration capacity falls from a high initial value to a steady one as the soil wets, and rain beyond it runs off
    over the surface.

    Parameters
    ----------
    time : float or array_like
        s, since the soil began to infiltrate at capacity; constraint: nonnegative.
    initial_rate : float or array_like
        m/s; constraint: at least final_rate.
    final_rate : float or array_like
        m/s; constraint: nonnegative.
    decay_constant : float or array_like
        1/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Horton (1933, 1940), an empirical law for rain at or above capacity throughout.
    Equation: ``final_rate + (initial_rate - final_rate) * np.exp(-decay_constant * time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    time = quantity("time", time, lower=0)
    initial_rate = quantity("initial_rate", initial_rate, lower=0)
    final_rate = quantity("final_rate", final_rate, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    require(initial_rate >= final_rate, "initial rate must be at least the final rate")
    return final_rate + (initial_rate - final_rate) * np.exp(-decay_constant * time)


@model
def horton_cumulative_infiltration(time, initial_rate, final_rate, decay_constant):
    """Depth infiltrated under Horton's law.

    The time integral of Horton's capacity is the water the soil can take up before the steady rate dominates.

    Parameters
    ----------
    time : float or array_like
        s; constraint: nonnegative.
    initial_rate : float or array_like
        m/s; constraint: at least final_rate.
    final_rate : float or array_like
        m/s; constraint: nonnegative.
    decay_constant : float or array_like
        1/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Horton (1933, 1940), with rain at or above capacity from time zero.
    Equation: ``final_rate * time + (initial_rate - final_rate) * (1 - np.exp(-decay_constant * time)) / decay_constant``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    time = quantity("time", time, lower=0)
    initial_rate = quantity("initial_rate", initial_rate, lower=0)
    final_rate = quantity("final_rate", final_rate, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    require(initial_rate >= final_rate, "initial rate must be at least the final rate")
    return final_rate * time + (initial_rate - final_rate) * -np.expm1(-decay_constant * time) / decay_constant


@safe_arithmetic
def simulate_infiltration(saturated_conductivity, suction_head, moisture_deficit, rainfall_rate, time):
    """Infiltration and runoff during steady rain on a Green-Ampt soil.

    Parameters
    ----------
    saturated_conductivity : float
        Positive scalar [m/s].
    suction_head : float
        Positive scalar wetting-front suction [m].
    moisture_deficit : float
        Scalar saturated minus initial water content, in (0, 1] [1].
    rainfall_rate : float
        Nonnegative scalar rain intensity [m/s].
    time : array_like
        Nonnegative strictly increasing times; the rain starts at the first sample [s].

    Returns
    -------
    TimeSeries
        Cumulative infiltration [m], infiltration rate [m/s] and infiltration-excess runoff rate [m/s].

    Notes
    -----
    Mein and Larson (1973): the soil absorbs all the rain until the Green-Ampt capacity
    ``K (1 + psi dtheta / F)`` falls to the rain rate, at ``F = K psi dtheta / (r - K)``, and infiltrates at capacity
    afterwards, so dF/dt = min(r, capacity). Rain no faster than K never ponds. The excess runs off at once; surface
    storage, redistribution and a water table are neglected.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    conductivity = quantity("saturated_conductivity", saturated_conductivity, lower=0, strict_lower=True)
    suction = quantity("suction_head", suction_head, lower=0, strict_lower=True)
    deficit = quantity("moisture_deficit", moisture_deficit, lower=0, upper=1, strict_lower=True)
    rain = quantity("rainfall_rate", rainfall_rate, lower=0)
    require(all(v.ndim == 0 for v in (conductivity, suction, deficit, rain)), "infiltration parameters must be scalars")
    time = time_grid(time)
    storage = suction * deficit

    def rate(depth):
        depth = np.asarray(depth, dtype=float)
        ponded = depth * (rain - conductivity) > conductivity * storage
        capacity = conductivity * (1 + storage / np.where(ponded, depth, 1.0))
        return np.where(ponded, capacity, rain)

    depth = integrate_ode(lambda t, state: rate(state), [0.0], time - time[0], names=("F",), units=("m",)).values[:, 0]
    infiltration = rate(depth)
    values = np.column_stack([depth, infiltration, rain - infiltration])
    return TimeSeries(time, values, ("cumulative_infiltration", "infiltration_rate", "runoff_rate"), ("m", "m/s", "m/s"))
