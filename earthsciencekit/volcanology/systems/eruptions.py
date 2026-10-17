"""Eruptions: buoyant plume rise, eruption-column height, the Volcanic Explosivity Index and a rising column parcel.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.results import TimeSeries

__all__ = [
    "plume_rise_height",
    "column_height",
    "volcanic_explosivity_index",
    "simulate_eruption_column",
]


@model
def plume_rise_height(buoyancy_flux, buoyancy_frequency):
    """Maximum height of a buoyant plume in a stratified fluid.

    Turbulent entrainment dilutes a plume until the stratification removes its buoyancy, which caps its height at a value
    set by only the source buoyancy flux and the stratification.

    Parameters
    ----------
    buoyancy_flux : float or array_like
        m^4/s^3, source ``reduced_gravity * velocity * radius**2`` (the factor pi omitted); constraint: positive.
    buoyancy_frequency : float or array_like
        1/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Morton, Taylor and Turner (1956) for a steady point source of buoyancy in uniform stratification, with Briggs's
    constant 5.0 (entrainment about 0.1); source momentum, wind and moisture are neglected.
    Equation: ``5.0 * buoyancy_flux**0.25 * buoyancy_frequency**-0.75``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    buoyancy_flux = quantity("buoyancy_flux", buoyancy_flux, lower=0, strict_lower=True)
    buoyancy_frequency = quantity("buoyancy_frequency", buoyancy_frequency, lower=0, strict_lower=True)
    return 5.0 * buoyancy_flux**0.25 * buoyancy_frequency**-0.75


@model
def column_height(heat_flux):
    """Height of a Plinian eruption column from its heat flux.

    The thermal energy released at the vent, not the eruption velocity, controls how high a volcanic column rises.

    Parameters
    ----------
    heat_flux : float or array_like
        W, rate of heat release ``mass_eruption_rate * specific_heat * temperature_excess``; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Wilson, Sparks, Huang and Watkins (1978), JGR 83: 1829-1836, for a standard atmosphere and a strong steady column;
    wind bends weak plumes lower.
    Equation: ``8.2 * heat_flux**0.25``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    heat_flux = quantity("heat_flux", heat_flux, lower=0, strict_lower=True)
    return 8.2 * heat_flux**0.25


@model
def volcanic_explosivity_index(erupted_volume):
    """Volcanic Explosivity Index from the erupted tephra volume.

    A single logarithmic number lets eruptions of very different style and size be compared and catalogued.

    Parameters
    ----------
    erupted_volume : float or array_like
        m^3, bulk tephra volume; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], an integer from 0 to 8. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Newhall and Self (1982), volume criterion only: VEI 0 below 1e4 m^3, 1 below 1e6 m^3, then one step per factor of
    ten up to VEI 8 at 1e12 m^3 and above. Column height and other criteria are not used.
    Equation: ``np.where(erupted_volume < 1e4, 0, np.where(erupted_volume < 1e6, 1,
    np.clip(np.floor(np.log10(erupted_volume)) - 4, 2, 8)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    erupted_volume = quantity("erupted_volume", erupted_volume, lower=0)
    steps = np.clip(np.floor(np.log10(np.maximum(erupted_volume, 1e6))) - 4, 2, 8)
    return np.where(erupted_volume < 1e4, 0.0, np.where(erupted_volume < 1e6, 1.0, steps))[()]


@safe_arithmetic
def simulate_eruption_column(vent_radius, vent_velocity, reduced_gravity, buoyancy_frequency, time, *, entrainment=0.12):
    """Follow a parcel up the axis of a steady buoyant eruption column.

    Parameters
    ----------
    vent_radius : float
        Positive scalar source radius [m].
    vent_velocity : float
        Positive scalar source velocity [m/s].
    reduced_gravity : float
        Positive scalar source buoyancy ``g (rho_air - rho_plume) / rho_air`` [m/s^2].
    buoyancy_frequency : float
        Nonnegative scalar Brunt-Vaisala frequency of the atmosphere [1/s].
    time : array_like
        Nonnegative strictly increasing travel times since the parcel left the vent [s].
    entrainment : float
        Positive scalar top-hat entrainment coefficient [1]; 0.12 reproduces the constant 5.0 of plume_rise_height.

    Returns
    -------
    TimeSeries
        Parcel height [m], vertical velocity [m/s] and reduced gravity [m/s^2].

    Notes
    -----
    The Morton, Taylor and Turner (1956) top-hat plume equations for the fluxes Q = b**2 w, M = b**2 w**2 and
    F = b**2 w g' (pi omitted), carried by a parcel moving at w = M / Q: dz/dt = M / Q, dQ/dt = 2 alpha abs(M)**1.5 / Q,
    dM/dt = F and dF/dt = -N**2 M. In a uniform atmosphere the parcel stops near plume_rise_height. Past that height the
    steady theory ends; the parcel is continued as a fountain that falls back with the same entrainment law and oscillates
    about the level of neutral buoyancy, near three quarters of the maximum height, where the umbrella cloud spreads.
    This illustrative closure omits wind, the dense gas-thrust region, particle fallout and moisture.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    radius = quantity("vent_radius", vent_radius, lower=0, strict_lower=True)
    velocity = quantity("vent_velocity", vent_velocity, lower=0, strict_lower=True)
    buoyancy = quantity("reduced_gravity", reduced_gravity, lower=0, strict_lower=True)
    frequency = quantity("buoyancy_frequency", buoyancy_frequency, lower=0)
    alpha = quantity("entrainment", entrainment, lower=0, strict_lower=True)
    require(all(v.ndim == 0 for v in (radius, velocity, buoyancy, frequency, alpha)), "eruption column parameters must be scalars")
    time = time_grid(time)

    def rhs(t, state):
        _, volume, momentum, flux = state
        return [momentum / volume, 2 * alpha * np.abs(momentum) ** 1.5 / volume, flux, -(frequency**2) * momentum]

    initial = [0.0, radius**2 * velocity, radius**2 * velocity**2, radius**2 * velocity * buoyancy]
    fluxes = integrate_ode(rhs, initial, time - time[0], names=("height", "Q", "M", "F"), units=("m", "m^3/s", "m^4/s^2", "m^4/s^3"))
    _, volume, momentum, flux = fluxes.values.T
    values = np.column_stack([fluxes.values[:, 0], momentum / volume, flux / volume])
    return TimeSeries(time, values, ("height", "vertical_velocity", "reduced_gravity"), ("m", "m/s", "m/s^2"))
