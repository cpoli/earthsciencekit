"""Bubbles: diffusive growth of a water bubble in a supersaturated melt, after Scriven (1959).

SI units unless stated otherwise. Water contents are mass fractions; the bubble holds water vapour as an ideal gas.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "diffusive_supersaturation",
    "scriven_growth_constant",
    "scriven_bubble_radius",
    "simulate_bubble_growth",
]


@model
def diffusive_supersaturation(initial_water, surface_water, melt_density, gas_density):
    """Dimensionless supersaturation that drives diffusive bubble growth.

    After decompression the melt holds more water than it can keep in solution, and the excess, measured in volumes of
    gas, sets how fast bubbles grow.

    Parameters
    ----------
    initial_water : float or array_like
        1, mass fraction dissolved far from the bubble; constraint: between surface_water and 1.
    surface_water : float or array_like
        1, mass fraction in equilibrium with the gas at the bubble wall; constraint: between 0 and 1.
    melt_density : float or array_like
        kg/m^3; constraint: positive.
    gas_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The mass-transfer analogue of the Jakob number of boiling, the input of :func:`scriven_growth_constant`.
    Equation: ``melt_density * (initial_water - surface_water) / gas_density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import diffusive_supersaturation
    >>> round(float(diffusive_supersaturation(0.05, 0.03, 2300, 100)), 6)
    0.46
    """
    initial_water = quantity("initial_water", initial_water, lower=0, upper=1)
    surface_water = quantity("surface_water", surface_water, lower=0, upper=1)
    melt_density = quantity("melt_density", melt_density, lower=0, strict_lower=True)
    gas_density = quantity("gas_density", gas_density, lower=0, strict_lower=True)
    require(initial_water >= surface_water, "initial_water must be at least surface_water")
    return melt_density * (initial_water - surface_water) / gas_density


def _scriven_supersaturation(beta, advection):
    """Scriven's supersaturation for a growth constant, integrated in units of the boundary-layer width."""
    from scipy.integrate import quad

    width = 1 / (1 + 2 * beta * (1 - advection))
    rate = 2 * beta * (1 - advection) * width

    def integrand(u):
        t = width * u
        x = beta + t
        return np.exp(-t * (2 * beta + t) + 2 * advection * beta**2 * t / x) * (beta / x) ** 2

    # The exponent stays below -(rate * u + (width * u)**2), so the integrand is under exp(-50) beyond this point.
    end = (np.sqrt(rate**2 + 200 * width**2) - rate) / (2 * width**2)
    edges = [0, 1, 10, end] if end > 10 else [0, 1, end] if end > 1 else [0, end]
    pieces = (quad(integrand, a, b, epsabs=0, epsrel=1e-12, limit=200)[0] for a, b in zip(edges[:-1], edges[1:]))
    return 2 * beta * width * sum(pieces)


@model
def scriven_growth_constant(supersaturation, density_ratio):
    """Growth constant of a bubble growing by diffusion from a supersaturated liquid.

    A bubble nucleated in a supersaturated melt grows as the square root of time, at a rate set by how much excess gas
    the melt holds and by how fast it diffuses.

    Parameters
    ----------
    supersaturation : float or array_like
        1, from :func:`diffusive_supersaturation`; constraint: positive and below ``1 / density_ratio``.
    density_ratio : float or array_like
        1, gas density over melt density; constraint: between 0 and 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], ``beta`` in ``R = 2 beta sqrt(D t)``. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Scriven (1959): the self-similar growth of a sphere from zero radius at constant supersaturation, including the
    radial flow the expanding bubble drives in the liquid. Solved by root finding on the integral equation, whose
    right-hand side rises from 0 to ``1 / density_ratio`` as ``beta`` grows; for small
    supersaturation ``beta = np.sqrt(supersaturation / 2)``, and for a light gas and large supersaturation
    ``beta = np.sqrt(3 / np.pi) * supersaturation``. Surface tension, viscosity and inertia are neglected.
    Equation: ``supersaturation = 2 * beta**3 * np.exp(beta**2 + 2 * e * beta**2) * integral(x**-2 * np.exp(-x**2 -
    2 * e * beta**3 / x), x from beta to infinity)`` with ``e = 1 - density_ratio``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import scriven_growth_constant
    >>> float(scriven_growth_constant(0.4842556877173758, 1))
    1.0
    """
    from scipy.optimize import brentq

    supersaturation = quantity("supersaturation", supersaturation, lower=0, strict_lower=True)
    density_ratio = quantity("density_ratio", density_ratio, lower=0, upper=1)
    require(supersaturation * density_ratio < 1, "supersaturation must be below 1 / density_ratio")
    supersaturation, density_ratio = np.broadcast_arrays(supersaturation, density_ratio)
    result = np.empty(supersaturation.shape)
    for index, (target, ratio) in enumerate(zip(supersaturation.flat, density_ratio.flat)):
        advection = 1 - ratio
        upper = 2 * max(np.sqrt(target), target)
        while _scriven_supersaturation(upper, advection) < target:
            require(upper < 1e4, "supersaturation is too large for this density_ratio")
            upper *= 2
        result.flat[index] = brentq(lambda b: _scriven_supersaturation(b, advection) - target, 0, upper, xtol=1e-15, rtol=1e-14)
    return result[()] if result.ndim == 0 else result


@model
def scriven_bubble_radius(time, diffusivity, growth_constant):
    """Radius of a bubble growing by diffusion at constant supersaturation.

    Diffusion-limited growth slows as the depleted shell around the bubble thickens, so the radius grows only as the
    square root of time.

    Parameters
    ----------
    time : float or array_like
        s, since nucleation; constraint: nonnegative.
    diffusivity : float or array_like
        m^2/s, of water in the melt (about 1e-11 for rhyolite near 850 degC); constraint: positive.
    growth_constant : float or array_like
        1, from :func:`scriven_growth_constant`; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Scriven (1959), growth from zero radius in an unbounded melt with constant diffusivity; a melt that runs out of
    water, as in :func:`simulate_bubble_growth`, slows the growth further.
    Equation: ``2 * growth_constant * np.sqrt(diffusivity * time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import scriven_bubble_radius
    >>> float(scriven_bubble_radius(4, 1, 0.5))
    2.0
    """
    time = quantity("time", time, lower=0)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    growth_constant = quantity("growth_constant", growth_constant, lower=0)
    return 2 * growth_constant * np.sqrt(diffusivity * time)


@safe_arithmetic
def simulate_bubble_growth(
    initial_radius,
    shell_radius,
    initial_water,
    pressure,
    temperature,
    time,
    *,
    solubility_constant=4.1e-6,
    diffusivity=1e-11,
    melt_density=2300.0,
):
    """Growth of one water bubble in its shell of melt after a sudden decompression.

    Parameters
    ----------
    initial_radius : float
        Positive scalar bubble radius at the first sample [m].
    shell_radius : float
        Scalar outer radius of the bubble's share of melt, above initial_radius [m]; for ``N`` bubbles per cubic
        metre it is about ``(3 / (4 pi N))**(1/3)``.
    initial_water : float
        Scalar dissolved water mass fraction before growth, in [0, 1) [1].
    pressure : float
        Positive scalar ambient pressure after decompression [Pa].
    temperature : float
        Positive scalar melt and gas temperature [K].
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    solubility_constant : float, optional
        Positive scalar of :func:`water_solubility` [1/sqrt(Pa)], 4.1e-6 for rhyolite.
    diffusivity : float, optional
        Positive scalar water diffusivity in the melt [m^2/s].
    melt_density : float, optional
        Positive scalar [kg/m^3].

    Returns
    -------
    TimeSeries
        Bubble radius [m], mean dissolved water in the shell [1] and vesicularity, the gas volume fraction [1].

    Notes
    -----
    The bubble wall holds the water solubility at the ambient pressure, and the gas is ideal water vapour of density
    ``pressure / (R_VAPOR * temperature)``. Growth is quasi-static diffusion, ``d(gas mass)/dt = 4 pi R melt_density D
    (c - c_s)``, which is Scriven's (1959) small-supersaturation limit, driven by the mean dissolved water ``c`` of a
    melt shell of fixed mass that the bubble depletes (a cell model, after Proussevitch et al. 1993). The bubble
    approaches the radius at which the melt is just saturated. Surface tension, viscous resistance, gas expansion
    during ascent and neighbouring bubbles' overlap are neglected.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.constants import R_VAPOR
    from earthsciencekit.integrators import integrate_ode
    from earthsciencekit.results import TimeSeries
    from earthsciencekit.volcanology.systems.magma import water_solubility

    radius = quantity("initial_radius", initial_radius, lower=0, strict_lower=True)
    shell = quantity("shell_radius", shell_radius, lower=0, strict_lower=True)
    water = quantity("initial_water", initial_water, lower=0, upper=1, strict_upper=True)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    diffusivity = quantity("diffusivity", diffusivity, lower=0, strict_lower=True)
    melt_density = quantity("melt_density", melt_density, lower=0, strict_lower=True)
    scalars = (radius, shell, water, pressure, temperature, diffusivity, melt_density)
    require(all(v.ndim == 0 for v in scalars), "bubble parameters must be scalars")
    require(shell > radius, "shell_radius must exceed initial_radius")
    time = time_grid(time)
    surface = water_solubility(pressure, solubility_constant)
    require(surface <= water, "the melt must be saturated or supersaturated at the new pressure")
    gas_density = pressure / (R_VAPOR * temperature)
    melt_volume = shell**3 - radius**3  # times 4 pi / 3

    def dissolved(r):
        return water - gas_density * (r**3 - radius**3) / (melt_density * melt_volume)

    def rhs(t, state):
        return [melt_density * diffusivity * (dissolved(state[0]) - surface) / (gas_density * state[0])]

    result = integrate_ode(rhs, [radius], time, names=("radius",), units=("m",), rtol=1e-10, atol=1e-16)
    r = result.values[:, 0]
    vesicularity = r**3 / (r**3 + melt_volume)
    values = np.column_stack([r, dissolved(r), vesicularity])
    return TimeSeries(time, values, ("radius", "dissolved_water", "vesicularity"), ("m", "1", "1"))
