"""Glaciers: mass balance, the equilibrium line, response time and a shallow-ice flowline glacier.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.constants import GRAVITY, ICE_DENSITY
from earthsciencekit.results import TimeSeries

__all__ = [
    "ela_shift",
    "linear_mass_balance",
    "accumulation_area_ratio",
    "glacier_response_time",
    "simulate_flowline_glacier",
]


@model
def ela_shift(temperature_change, lapse_rate):
    """Shift of the equilibrium-line altitude with a change of climate.

    A colder climate brings the snow line down the mountains, and the glaciers grow to meet it.

    Parameters
    ----------
    temperature_change : float or array_like
        K, change of summer air temperature at fixed altitude; constraint: real.
    lapse_rate : float or array_like
        K/m, decrease of temperature with height; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], positive upward, so cooling lowers the line. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Temperature alone controls ablation and the lapse rate is uniform; changes of snowfall also move the line.
    Equation: ``temperature_change / lapse_rate``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature_change = quantity("temperature_change", temperature_change)
    lapse_rate = quantity("lapse_rate", lapse_rate, lower=0, strict_lower=True)
    return temperature_change / lapse_rate


@model
def linear_mass_balance(elevation, equilibrium_line_altitude, balance_gradient):
    """Surface mass balance rising linearly with altitude.

    A glacier gains mass above its equilibrium line and loses it below, at a rate that grows steadily with the distance
    from the line.

    Parameters
    ----------
    elevation : float or array_like
        m; constraint: real.
    equilibrium_line_altitude : float or array_like
        m; constraint: real.
    balance_gradient : float or array_like
        1/s, metres of ice per second per metre of altitude; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s] of ice, positive for net accumulation. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Meier (1962): the balance gradient is about 0.003-0.01 m of water per year per metre on mid-latitude glaciers; real
    profiles flatten at high altitude.
    Equation: ``balance_gradient * (elevation - equilibrium_line_altitude)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    elevation = quantity("elevation", elevation)
    equilibrium_line_altitude = quantity("equilibrium_line_altitude", equilibrium_line_altitude)
    balance_gradient = quantity("balance_gradient", balance_gradient, lower=0)
    return balance_gradient * (elevation - equilibrium_line_altitude)


@model
def accumulation_area_ratio(accumulation_area, glacier_area):
    """Fraction of a glacier's area above its equilibrium line.

    A glacier in balance with its climate keeps a characteristic share of its area in the accumulation zone, about 0.6.

    Parameters
    ----------
    accumulation_area : float or array_like
        m^2; constraint: between 0 and glacier_area.
    glacier_area : float or array_like
        m^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Meier (1962); the steady-state ratio depends on the glacier's hypsometry and is not universal.
    Equation: ``accumulation_area / glacier_area``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    accumulation_area = quantity("accumulation_area", accumulation_area, lower=0)
    glacier_area = quantity("glacier_area", glacier_area, lower=0, strict_lower=True)
    require(accumulation_area <= glacier_area, "accumulation area exceeds the glacier area")
    return accumulation_area / glacier_area


@model
def glacier_response_time(thickness, terminus_ablation):
    """Volume response time of a glacier to a change of climate.

    A glacier adjusts to a new climate in about the time its terminus ablation takes to remove its own thickness.

    Parameters
    ----------
    thickness : float or array_like
        m, a characteristic thickness; constraint: positive.
    terminus_ablation : float or array_like
        m/s of ice, ablation rate at the terminus as a positive magnitude; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s], an e-folding time. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Johannesson, Raymond and Waddington (1989): a linear response of volume; decades for valley glaciers, millennia for
    ice sheets.
    Equation: ``thickness / terminus_ablation``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    terminus_ablation = quantity("terminus_ablation", terminus_ablation, lower=0, strict_lower=True)
    return thickness / terminus_ablation


@safe_arithmetic
def simulate_flowline_glacier(
    bed_elevation,
    spacing,
    equilibrium_line_altitude,
    balance_gradient,
    time,
    *,
    initial_thickness=None,
    rate_factor=2.4e-24,
    exponent=3.0,
    density=ICE_DENSITY,
    gravity=GRAVITY,
):
    """Shallow-ice flowline glacier under a linear mass balance.

    Parameters
    ----------
    bed_elevation : array_like, shape (m,), m >= 3
        Bed elevation at equally spaced nodes along the flowline, head first [m].
    spacing : float
        Positive scalar node spacing [m].
    equilibrium_line_altitude : float
        Scalar equilibrium-line altitude held from the first sample on [m].
    balance_gradient : float
        Positive scalar mass-balance gradient [1/s].
    time : array_like
        Nonnegative strictly increasing output times [s]; the initial state belongs to the first sample.
    initial_thickness : array_like, shape (m,), optional
        Nonnegative ice thickness at the first sample [m]; ice-free by default. Pass a steady state for one climate
        with the equilibrium line of another to see the response to a step.
    rate_factor : float, optional
        Positive scalar Glen rate factor [Pa^-n/s], 2.4e-24 for ice at 0 degC (Cuffey and Paterson 2010).
    exponent : float, optional
        Scalar Glen exponent of at least 1 [1].
    density, gravity : float, optional
        Positive scalar ice density [kg/m^3] and gravity [m/s^2].

    Returns
    -------
    TimeSeries
        Glacier length [m] (the spacing times the number of nodes with more than 1 m of ice), cross-sectional area per
        unit width [m^2], and the thickness at every node [m], labeled thickness_<node>.

    Notes
    -----
    Solves ``dH/dt = b(b_bed + H) - dq/dx`` with the shallow-ice flux ``q = -D ds/dx`` of
    :func:`shallow_ice_diffusivity` and the balance of :func:`linear_mass_balance`, by the method of lines on a
    staggered grid in Mahaffy's (1976) form, with ``integrators.integrate_ode`` (BDF) in time. The ends carry no
    flux, the thickness is kept nonnegative (ice-free nodes cannot ablate), and there is no sliding, calving,
    width variation or isostasy.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    bed = quantity("bed_elevation", bed_elevation)
    require(bed.ndim == 1 and bed.size >= 3, "bed_elevation must be a vector of at least three nodes")
    spacing = quantity("spacing", spacing, lower=0, strict_lower=True)
    ela = quantity("equilibrium_line_altitude", equilibrium_line_altitude)
    gradient = quantity("balance_gradient", balance_gradient, lower=0, strict_lower=True)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    scalars = (spacing, ela, gradient, rate_factor, exponent, density, gravity)
    require(all(v.ndim == 0 for v in scalars), "flowline parameters must be scalars")
    time = time_grid(time)
    initial = np.zeros_like(bed) if initial_thickness is None else quantity("initial_thickness", initial_thickness, lower=0)
    require(initial.shape == bed.shape, "initial_thickness must match bed_elevation")
    stiffness = 2 * rate_factor / (exponent + 2) * (density * gravity) ** exponent

    def rhs(t, state):
        thickness = np.maximum(state, 0.0)
        surface = bed + thickness
        slope = np.diff(surface) / spacing
        mean = 0.5 * (thickness[1:] + thickness[:-1])
        flux = -stiffness * mean ** (exponent + 2) * np.abs(slope) ** (exponent - 1) * slope
        divergence = np.diff(np.concatenate([[0.0], flux, [0.0]])) / spacing
        rate = gradient * (surface - ela) - divergence
        return np.where((state <= 0) & (rate < 0), 0.0, rate)

    result = integrate_ode(rhs, initial, time, method="BDF", rtol=1e-6, atol=1e-6)
    thickness = np.maximum(result.values, 0.0)
    length = spacing * np.count_nonzero(thickness > 1.0, axis=1)
    area = spacing * thickness.sum(axis=1)
    names = ("length", "area") + tuple(f"thickness_{i}" for i in range(bed.size))
    units = ("m", "m^2") + ("m",) * bed.size
    return TimeSeries(time, np.column_stack([length, area, thickness]), names, units)
