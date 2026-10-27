"""Energy-balance models: the two-box (surface and deep ocean) response and North's latitudinal diffusive model.

SI units unless stated otherwise. Latitude enters the diffusive model through its sine ``x``, so equal steps in ``x``
cover equal areas of the sphere.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "two_box_timescales",
    "two_box_fast_warming",
    "insolation_distribution",
    "north_temperature_profile",
    "simulate_two_box",
    "simulate_latitudinal_energy_balance",
]


@model
def two_box_timescales(feedback, exchange, surface_heat_capacity, deep_heat_capacity):
    """Fast and slow e-folding times of the two-box climate.

    A shallow mixed layer coupled to a deep ocean responds on two timescales: years for the surface to adjust, and
    centuries for the deep ocean to fill with heat.

    Parameters
    ----------
    feedback : float or array_like
        W/(m^2 K), the climate feedback parameter lambda; constraint: positive.
    exchange : float or array_like
        W/(m^2 K), the heat-exchange coefficient gamma between the boxes; constraint: positive.
    surface_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.
    deep_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [s]: the fast and the slow timescale on the last axis.

    Notes
    -----
    Held et al. (2010): ``C dT/dt = F - lambda T - gamma (T - T0)`` and ``C0 dT0/dt = gamma (T - T0)``; the timescales
    are minus the reciprocals of the eigenvalues of that linear system. For ``C0 >> C`` they tend to
    ``C / (lambda + gamma)`` and ``C0 (lambda + gamma) / (lambda gamma)``.
    Equation: ``2 / (b -+ np.sqrt(b**2 - 4 * lambda * gamma / (C * C0)))`` with
    ``b = (lambda + gamma) / C + gamma / C0``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.climate import two_box_timescales
    >>> np.round(two_box_timescales(1, 1, 1, 1), 5).tolist()
    [0.38197, 2.61803]
    """
    feedback = quantity("feedback", feedback, lower=0, strict_lower=True)
    exchange = quantity("exchange", exchange, lower=0, strict_lower=True)
    surface = quantity("surface_heat_capacity", surface_heat_capacity, lower=0, strict_lower=True)
    deep = quantity("deep_heat_capacity", deep_heat_capacity, lower=0, strict_lower=True)
    trace = (feedback + exchange) / surface + exchange / deep
    determinant = feedback * exchange / (surface * deep)
    root = np.sqrt(trace**2 - 4 * determinant)
    # The slow rate is the small root; determinant / fast avoids cancellation in trace - root.
    fast = (trace + root) / 2
    return np.stack(np.broadcast_arrays(1 / fast, fast / determinant), axis=-1)


@model
def two_box_fast_warming(forcing, feedback, exchange):
    """Warming reached once the fast response is over, while the deep ocean is still cold.

    After a few years the surface sits on a plateau, held down by the heat the deep ocean still absorbs; the rest of
    the equilibrium warming arrives over centuries.

    Parameters
    ----------
    forcing : float or array_like
        W/m^2; constraint: real.
    feedback : float or array_like
        W/(m^2 K); constraint: positive.
    exchange : float or array_like
        W/(m^2 K); constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Held et al. (2010), the two-box surface temperature with the deep ocean still at its initial temperature; it
    approximates the transient climate response. Without exchange it is the equilibrium warming.
    Equation: ``forcing / (feedback + exchange)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import two_box_fast_warming
    >>> round(float(two_box_fast_warming(3.7, 1.2, 0.7)), 5)
    1.94737
    """
    forcing = quantity("forcing", forcing)
    feedback = quantity("feedback", feedback, lower=0, strict_lower=True)
    exchange = quantity("exchange", exchange, lower=0)
    return forcing / (feedback + exchange)


@model
def insolation_distribution(sine_latitude, second_legendre_coefficient):
    """Annual-mean insolation relative to its global mean.

    The tropics receive more than twice the annual sunshine of the poles, the contrast that drives poleward heat
    transport.

    Parameters
    ----------
    sine_latitude : float or array_like
        1; constraint: between -1 and 1.
    second_legendre_coefficient : float or array_like
        1, ``s2``, about -0.482 for the present Earth; constraint: between -1 and 2.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], with global mean 1. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    North (1975): the two-term Legendre fit to the annual mean for the present obliquity; it is within a few percent of
    the exact annual mean.
    Equation: ``1 + second_legendre_coefficient * (3 * sine_latitude**2 - 1) / 2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import insolation_distribution
    >>> float(insolation_distribution(1, -0.482))
    0.518
    """
    sine_latitude = quantity("sine_latitude", sine_latitude, lower=-1, upper=1)
    second_legendre_coefficient = quantity("second_legendre_coefficient", second_legendre_coefficient, lower=-1, upper=2)
    return 1 + second_legendre_coefficient * (3 * sine_latitude**2 - 1) / 2


@model
def north_temperature_profile(
    sine_latitude,
    solar_constant,
    albedo,
    outgoing_intercept,
    outgoing_slope,
    diffusion,
    second_legendre_coefficient,
):
    """Steady temperature of North's diffusive energy-balance model with a fixed albedo.

    Diffusive heat transport flattens the equator-to-pole contrast of absorbed sunshine, the more so the larger the
    diffusivity.

    Parameters
    ----------
    sine_latitude : float or array_like
        1; constraint: between -1 and 1.
    solar_constant : float or array_like
        W/m^2; constraint: positive.
    albedo : float or array_like
        1; constraint: fraction.
    outgoing_intercept : float or array_like
        W/m^2, ``A`` in the outgoing longwave ``A + B (T - 273.15)``; constraint: real.
    outgoing_slope : float or array_like
        W/(m^2 K), ``B``; constraint: positive.
    diffusion : float or array_like
        W/(m^2 K), ``D``; constraint: nonnegative.
    second_legendre_coefficient : float or array_like
        1, ``s2`` of :func:`insolation_distribution`; constraint: between -1 and 2.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    North (1975): ``Q s(x) (1 - albedo) = A + B (T - 273.15) - D d/dx[(1 - x**2) dT/dx]`` with ``Q = S / 4``. The
    second Legendre polynomial is an eigenfunction of the diffusion operator (eigenvalue ``-6 D``), so the solution has
    two modes; with a temperature-dependent albedo it no longer does.
    Equation: ``273.15 + (Q * (1 - albedo) - A) / B + Q * (1 - albedo) * s2 / (B + 6 * D) * (3 * x**2 - 1) / 2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import north_temperature_profile
    >>> float(north_temperature_profile(1, 1360, 0.3, 210, 2, 0.5, -0.482))
    264.2068
    """
    x = quantity("sine_latitude", sine_latitude, lower=-1, upper=1)
    absorbed = quantity("solar_constant", solar_constant, lower=0, strict_lower=True) / 4
    absorbed = absorbed * (1 - quantity("albedo", albedo, lower=0, upper=1))
    intercept = quantity("outgoing_intercept", outgoing_intercept)
    slope = quantity("outgoing_slope", outgoing_slope, lower=0, strict_lower=True)
    diffusion = quantity("diffusion", diffusion, lower=0)
    s2 = quantity("second_legendre_coefficient", second_legendre_coefficient, lower=-1, upper=2)
    return 273.15 + (absorbed - intercept) / slope + absorbed * s2 / (slope + 6 * diffusion) * (3 * x**2 - 1) / 2


@safe_arithmetic
def simulate_two_box(forcing, feedback, exchange, surface_heat_capacity, deep_heat_capacity, time, *, initial_anomaly=(0.0, 0.0)):
    """Surface and deep-ocean temperature anomalies of the two-box climate under a constant forcing.

    Parameters
    ----------
    forcing : float
        Constant signed radiative forcing from the first sample on [W/m^2].
    feedback : float
        Positive scalar climate feedback parameter lambda [W/(m^2 K)].
    exchange : float
        Positive scalar heat-exchange coefficient gamma [W/(m^2 K)].
    surface_heat_capacity, deep_heat_capacity : float
        Positive scalar heat capacities C and C0 [J/(m^2 K)].
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    initial_anomaly : sequence of two floats, optional
        Surface and deep anomalies at the first sample [K]. Start from a warmed state with zero forcing to follow
        Held's abrupt return to preindustrial forcing.

    Returns
    -------
    TimeSeries
        Surface and deep-ocean temperature anomalies [K].

    Notes
    -----
    Held et al. (2010): ``C dT/dt = F - lambda T - gamma (T - T0)`` and ``C0 dT0/dt = gamma (T - T0)``, integrated with
    ``integrators.integrate_ode``. The decay rates are those of :func:`two_box_timescales`. There is no efficacy factor,
    state-dependent feedback or ocean circulation change.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    forcing = quantity("forcing", forcing)
    feedback = quantity("feedback", feedback, lower=0, strict_lower=True)
    exchange = quantity("exchange", exchange, lower=0, strict_lower=True)
    surface = quantity("surface_heat_capacity", surface_heat_capacity, lower=0, strict_lower=True)
    deep = quantity("deep_heat_capacity", deep_heat_capacity, lower=0, strict_lower=True)
    initial = quantity("initial_anomaly", initial_anomaly)
    require(all(v.ndim == 0 for v in (forcing, feedback, exchange, surface, deep)), "two-box parameters must be scalars")
    require(initial.shape == (2,), "initial_anomaly must hold the surface and deep anomalies")
    time = time_grid(time)

    def rhs(t, state):
        flux = exchange * (state[0] - state[1])
        return [(forcing - feedback * state[0] - flux) / surface, flux / deep]

    return integrate_ode(rhs, initial, time, names=("surface_anomaly", "deep_anomaly"), units=("K", "K"))


@safe_arithmetic
def simulate_latitudinal_energy_balance(
    solar_constant,
    time,
    *,
    nodes=60,
    initial_temperature=288.0,
    heat_capacity=4.2e7,
    outgoing_intercept=210.0,
    outgoing_slope=2.0,
    diffusion=0.555,
    cold_albedo=0.62,
    warm_albedo=0.3,
    warm_albedo_gradient=0.078,
    freezing_temperature=263.15,
    transition_width=1.0,
    second_legendre_coefficient=-0.482,
):
    """North's one-dimensional diffusive energy-balance model with ice-albedo feedback.

    Parameters
    ----------
    solar_constant : float
        Positive scalar solar constant S [W/m^2]; the global-mean insolation is S / 4.
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    nodes : int, optional
        Number of equal-area cells between the equator and the pole, at least 3.
    initial_temperature : float or array_like, optional
        Positive temperature at the first sample, a scalar or one value per cell from the equator [K].
    heat_capacity : float, optional
        Positive scalar areal heat capacity [J/(m^2 K)], about 10 m of water by default.
    outgoing_intercept, outgoing_slope : float, optional
        Outgoing longwave ``A + B (T - 273.15)``: A [W/m^2] and positive B [W/(m^2 K)].
    diffusion : float, optional
        Nonnegative scalar diffusion coefficient D [W/(m^2 K)].
    cold_albedo, warm_albedo : float, optional
        Albedo of ice-covered cells and the global-mean albedo of ice-free ones [1].
    warm_albedo_gradient : float, optional
        Coefficient of ``(3 x**2 - 1) / 2`` in the ice-free albedo, which rises towards the poles as the Sun gets
        lower [1]. The ice-free albedo must stay within 0 and cold_albedo.
    freezing_temperature : float, optional
        Positive scalar temperature of the ice edge [K], -10 degC by default.
    transition_width : float, optional
        Positive scalar width of the smooth albedo step [K], as in :func:`ice_albedo`.
    second_legendre_coefficient : float, optional
        ``s2`` of :func:`insolation_distribution`, between -1 and 2 [1].

    Returns
    -------
    TimeSeries
        Global-mean temperature [K], ice-edge latitude [rad] (``pi / 2`` without ice, 0 for a snowball) and the
        temperature of every cell [K], labeled temperature_<cell> from the equator.

    Notes
    -----
    ``C dT/dt = (S / 4) s(x) (1 - albedo(T)) - (A + B (T - 273.15)) + D d/dx[(1 - x**2) dT/dx]`` for one hemisphere,
    symmetric about the equator, after Budyko (1969), Sellers (1969) and North (1975). Cells are equal steps in the sine
    of latitude, the diffusive flux vanishes at the equator and the pole, the albedo is :func:`ice_albedo` and the time
    integration is ``integrators.integrate_ode`` with BDF. The ice-free albedo has the latitude dependence of North
    et al. (1981), and the defaults put the ice edge near 71 degrees and the global mean near 287 K for
    ``S = 1361``. With
    ``cold_albedo == warm_albedo`` and no gradient the steady state is :func:`north_temperature_profile`. There are no
    seasons, land, ocean dynamics or cloud feedbacks.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.climate.systems.response import ice_albedo
    from earthsciencekit.integrators import integrate_ode
    from earthsciencekit.results import TimeSeries

    solar = quantity("solar_constant", solar_constant, lower=0, strict_lower=True)
    require(isinstance(nodes, (int, np.integer)) and not isinstance(nodes, bool) and nodes >= 3, "nodes must be an integer >= 3")
    capacity = quantity("heat_capacity", heat_capacity, lower=0, strict_lower=True)
    intercept = quantity("outgoing_intercept", outgoing_intercept)
    slope = quantity("outgoing_slope", outgoing_slope, lower=0, strict_lower=True)
    diffusion = quantity("diffusion", diffusion, lower=0)
    freezing = quantity("freezing_temperature", freezing_temperature, lower=0, strict_lower=True)
    scalars = (solar, capacity, intercept, slope, diffusion, freezing)
    require(all(v.ndim == 0 for v in scalars), "energy-balance parameters must be scalars")
    time = time_grid(time)
    initial = np.broadcast_to(quantity("initial_temperature", initial_temperature, lower=0, strict_lower=True), (nodes,)).copy()
    spacing = 1 / nodes
    x = (np.arange(nodes) + 0.5) * spacing
    faces = np.arange(1, nodes) * spacing
    absorbed_max = solar / 4 * insolation_distribution(x, second_legendre_coefficient)
    gradient = quantity("warm_albedo_gradient", warm_albedo_gradient)
    require(gradient.ndim == 0, "warm_albedo_gradient must be a scalar")
    ice_free = quantity("warm_albedo", warm_albedo, lower=0, upper=1) + gradient * (3 * x**2 - 1) / 2

    def albedo(temperature):
        return ice_albedo(temperature, cold_albedo, ice_free, freezing, transition_width)

    albedo(initial)  # validates the albedo parameters before integrating

    def rhs(t, temperature):
        flux = diffusion * (1 - faces**2) * np.diff(temperature) / spacing
        transport = np.diff(np.concatenate([[0.0], flux, [0.0]])) / spacing
        return (absorbed_max * (1 - albedo(temperature)) - intercept - slope * (temperature - 273.15) + transport) / capacity

    result = integrate_ode(rhs, initial, time, method="BDF", rtol=1e-8, atol=1e-6)
    temperature = result.values
    warm = temperature >= freezing
    # Interpolate the edge in x between the last warm cell and the first cold one, counting cells from the equator.
    first_cold = np.where(warm.all(axis=1), nodes, np.argmin(warm, axis=1))
    edge = np.empty(time.size)
    for row, index in enumerate(first_cold):
        if index == nodes:
            edge[row] = 1.0
        elif index == 0:
            edge[row] = 0.0
        else:
            above, below = temperature[row, index - 1], temperature[row, index]
            edge[row] = x[index - 1] + spacing * (above - freezing) / (above - below)
    values = np.column_stack([temperature.mean(axis=1), np.arcsin(edge), temperature])
    names = ("global_mean_temperature", "ice_edge") + tuple(f"temperature_{i}" for i in range(nodes))
    return TimeSeries(time, values, names, ("K", "rad") + ("K",) * nodes)
