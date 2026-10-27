"""Climate response: feedbacks, heat capacity, transient and stochastic one-box response, ice-albedo equilibria and the Gregory regression.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.constants import SIGMA_SB

__all__ = [
    "planck_feedback",
    "equilibrium_warming",
    "heat_capacity_column",
    "response_timescale",
    "step_temperature",
    "heat_uptake",
    "ice_albedo",
    "equilibrium_solar_constant",
    "red_noise_spectrum",
    "gregory_regression",
    "simulate_one_box",
]


@model
def planck_feedback(temperature, emissivity):
    """Blackbody restoring derivative.

    The temperature derivative of emission supplies a radiative restoring coefficient.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    emissivity : float or array_like
        1; constraint: open fraction.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Positive restoring magnitude in W per square metre per kelvin; other feedbacks are separate.
    Equation: ``4 * emissivity * SIGMA_SB * temperature**3``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import planck_feedback
    >>> round(float(planck_feedback(300, 1)), 5)
    6.124
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    emissivity = quantity("emissivity", emissivity, lower=0, upper=1, strict_lower=True)
    return 4 * emissivity * SIGMA_SB * temperature**3


@model
def equilibrium_warming(forcing, restoring):
    """Linear climate response.

    A restoring feedback maps steady forcing to equilibrium temperature change.

    Parameters
    ----------
    forcing : float or array_like
        W/m^2; constraint: real.
    restoring : float or array_like
        W/(m^2 K); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Stable linear response with a positive net restoring coefficient.
    Equation: ``forcing / restoring``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import equilibrium_warming
    >>> float(equilibrium_warming(4, 2))
    2.0
    """
    forcing = quantity("forcing", forcing)
    restoring = quantity("restoring", restoring, lower=0, strict_lower=True)
    return forcing / restoring


@model
def heat_capacity_column(density, heat_capacity, depth):
    """Mixed-layer areal heat capacity.

    Mixed-layer depth determines how much energy is needed for warming.

    Parameters
    ----------
    density : float or array_like
        kg/m^3; constraint: positive.
    heat_capacity : float or array_like
        J/(kg K); constraint: positive.
    depth : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uniform slab with constant material properties.
    Equation: ``density * heat_capacity * depth``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import heat_capacity_column
    >>> float(heat_capacity_column(1000, 4000, 50))
    200000000.0
    """
    density = quantity("density", density, lower=0, strict_lower=True)
    heat_capacity = quantity("heat_capacity", heat_capacity, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    return density * heat_capacity * depth


@model
def response_timescale(areal_heat_capacity, restoring):
    """Climate relaxation time.

    Heat storage delays the equilibrium response to forcing.

    Parameters
    ----------
    areal_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.
    restoring : float or array_like
        W/(m^2 K); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    One-box linear model; output is seconds.
    Equation: ``areal_heat_capacity / restoring``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import response_timescale
    >>> float(response_timescale(200000000.0, 2))
    100000000.0
    """
    areal_heat_capacity = quantity("areal_heat_capacity", areal_heat_capacity, lower=0, strict_lower=True)
    restoring = quantity("restoring", restoring, lower=0, strict_lower=True)
    return areal_heat_capacity / restoring


@model
def step_temperature(forcing, restoring, areal_heat_capacity, time):
    """Transient one-box warming.

    A step forcing exposes the lag between forcing and temperature adjustment.

    Parameters
    ----------
    forcing : float or array_like
        W/m^2; constraint: real.
    restoring : float or array_like
        W/(m^2 K); constraint: positive.
    areal_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Zero initial anomaly and constant parameters; no deep-ocean exchange.
    Equation: ``forcing / restoring * (-np.expm1(-restoring * time / areal_heat_capacity))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import step_temperature
    >>> float(step_temperature(4, 2, 2, 0))
    0.0
    """
    forcing = quantity("forcing", forcing)
    restoring = quantity("restoring", restoring, lower=0, strict_lower=True)
    areal_heat_capacity = quantity("areal_heat_capacity", areal_heat_capacity, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0)
    return forcing / restoring * (-np.expm1(-restoring * time / areal_heat_capacity))


@model
def heat_uptake(areal_heat_capacity, temperature_change):
    """Areal heat-content change.

    Energy accounting connects warming to ocean heat storage.

    Parameters
    ----------
    areal_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.
    temperature_change : float or array_like
        K; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    One uniform layer; temperature_change is an anomaly.
    Equation: ``areal_heat_capacity * temperature_change``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import heat_uptake
    >>> float(heat_uptake(200000000.0, 2))
    400000000.0
    """
    areal_heat_capacity = quantity("areal_heat_capacity", areal_heat_capacity, lower=0, strict_lower=True)
    temperature_change = quantity("temperature_change", temperature_change)
    return areal_heat_capacity * temperature_change


@model
def ice_albedo(temperature, cold_albedo, warm_albedo, transition_temperature, width):
    """Smooth illustrative ice-albedo transition.

    Temperature-dependent reflection illustrates ice-albedo feedback.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    cold_albedo : float or array_like
        1; constraint: fraction.
    warm_albedo : float or array_like
        1; constraint: fraction.
    transition_temperature : float or array_like
        K; constraint: positive.
    width : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Illustrative smooth closure, not a calibrated satellite albedo retrieval.
    Equation: ``warm_albedo + (cold_albedo - warm_albedo) * (1 - np.tanh((temperature - transition_temperature) / width)) / 2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import ice_albedo
    >>> float(ice_albedo(273, 0.6, 0.2, 273, 5))
    0.4
    """
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    cold_albedo = quantity("cold_albedo", cold_albedo, lower=0, upper=1)
    warm_albedo = quantity("warm_albedo", warm_albedo, lower=0, upper=1)
    transition_temperature = quantity("transition_temperature", transition_temperature, lower=0, strict_lower=True)
    width = quantity("width", width, lower=0, strict_lower=True)
    require(cold_albedo >= warm_albedo, "cold albedo must be at least warm albedo")
    return warm_albedo + (cold_albedo - warm_albedo) * (1 - np.tanh((temperature - transition_temperature) / width)) / 2


@model
def equilibrium_solar_constant(temperature, emissivity, cold_albedo, warm_albedo, transition_temperature, width):
    """Solar constant that holds a planet with ice-albedo feedback at a given temperature.

    Inverting the energy balance traces every equilibrium, including the unstable branch between a warm and an ice-covered state.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: positive.
    emissivity : float or array_like
        1; constraint: open fraction.
    cold_albedo : float or array_like
        1; constraint: fraction below one.
    warm_albedo : float or array_like
        1; constraint: fraction.
    transition_temperature : float or array_like
        K; constraint: positive.
    width : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Zero-dimensional greybody balance with the ice_albedo closure; equilibria with a negative slope dS/dT are unstable.
    Equation: ``4 * emissivity * SIGMA_SB * temperature**4 / (1 - albedo)``, with ``albedo`` from ice_albedo.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import equilibrium_solar_constant
    >>> round(float(equilibrium_solar_constant(288, 1, 0.6, 0.2, 288, 5)), 2)
    2600.7
    """
    emissivity = quantity("emissivity", emissivity, lower=0, upper=1, strict_lower=True)
    require(quantity("cold_albedo", cold_albedo) < 1, "cold albedo must be below one")
    albedo = ice_albedo(temperature, cold_albedo, warm_albedo, transition_temperature, width)
    return 4 * emissivity * SIGMA_SB * quantity("temperature", temperature) ** 4 / (1 - albedo)


@model
def red_noise_spectrum(frequency, forcing_spectrum, restoring, areal_heat_capacity):
    """Temperature spectrum of a one-box climate driven by white noise.

    Heat capacity integrates random weather forcing into slow red-noise climate variability.

    Parameters
    ----------
    frequency : float or array_like
        Hz; constraint: real.
    forcing_spectrum : float or array_like
        W^2/(m^4 Hz), two-sided; constraint: nonnegative.
    restoring : float or array_like
        W/(m^2 K); constraint: positive.
    areal_heat_capacity : float or array_like
        J/(m^2 K); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K^2/Hz], two-sided. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Linear C dT/dt = F - lambda T with white-noise F; integrating over all frequencies gives the
    variance forcing_spectrum / (2 * restoring * areal_heat_capacity).
    Equation: ``forcing_spectrum / (restoring**2 + (2 * np.pi * frequency * areal_heat_capacity)**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.climate import red_noise_spectrum
    >>> float(red_noise_spectrum(0, 4, 2, 1))
    1.0
    """
    frequency = quantity("frequency", frequency)
    forcing_spectrum = quantity("forcing_spectrum", forcing_spectrum, lower=0)
    restoring = quantity("restoring", restoring, lower=0, strict_lower=True)
    areal_heat_capacity = quantity("areal_heat_capacity", areal_heat_capacity, lower=0, strict_lower=True)
    return forcing_spectrum / (restoring**2 + (2 * np.pi * frequency * areal_heat_capacity) ** 2)


@model
def gregory_regression(temperature_change, imbalance):
    """Forcing and feedback from a regression of imbalance on warming.

    In a run with constant forcing, the top-of-atmosphere imbalance falls linearly with warming, so a straight-line fit gives both.

    Parameters
    ----------
    temperature_change : array_like, shape (n,), n >= 2
        K; constraint: real, not all equal.
    imbalance : array_like, shape (n,)
        W/m^2, net downward top-of-atmosphere flux; constraint: real.

    Returns
    -------
    numpy.ndarray, shape (2,)
        ``[forcing, restoring]``: the intercept [W/m^2] and minus the slope [W/(m^2 K)] of the least-squares line.

    Notes
    -----
    Assumes N = F - lambda T with constant F and lambda; equilibrium warming is forcing / restoring.
    Equation: ``imbalance = forcing - restoring * temperature_change`` (ordinary least squares).
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.climate import gregory_regression
    >>> np.round(gregory_regression([0, 1, 2], [4, 2.8, 1.6]), 5).tolist()
    [4.0, 1.2]
    """
    temperature_change = quantity("temperature_change", temperature_change)
    imbalance = quantity("imbalance", imbalance)
    require(temperature_change.ndim == 1 and temperature_change.size >= 2, "temperature_change must be a vector of at least two samples")
    require(imbalance.shape == temperature_change.shape, "imbalance must match temperature_change")
    anomaly = temperature_change - temperature_change.mean()
    require(np.any(anomaly != 0), "temperature_change must vary")
    slope = np.dot(anomaly, imbalance - imbalance.mean()) / np.dot(anomaly, anomaly)
    return np.array([imbalance.mean() - slope * temperature_change.mean(), -slope])


@safe_arithmetic
def simulate_one_box(forcing, restoring, areal_heat_capacity, time, *, initial_anomaly=0.0):
    """Integrate a linear one-box climate energy balance.

    Parameters
    ----------
    forcing : float
        Constant signed radiative perturbation [W/m^2].
    restoring : float
        Positive scalar net restoring coefficient [W/(m^2 K)].
    areal_heat_capacity : float
        Positive scalar effective heat capacity [J/(m^2 K)].
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    initial_anomaly : float
        Finite scalar temperature anomaly [K] at the first supplied time.

    Returns
    -------
    TimeSeries
        Temperature anomaly [K].

    Notes
    -----
    C dT/dt = F - lambda T. Compare with step_temperature for zero initial
    anomaly, using elapsed time since the first sample. This educational model
    omits deep-ocean exchange, nonlinear feedbacks and geographic structure.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    forcing = quantity("forcing", forcing)
    restoring = quantity("restoring", restoring, lower=0, strict_lower=True)
    capacity = quantity("areal_heat_capacity", areal_heat_capacity, lower=0, strict_lower=True)
    initial = quantity("initial_anomaly", initial_anomaly)
    require(all(v.ndim == 0 for v in (forcing, restoring, capacity, initial)), "climate parameters must be scalars")
    time = time_grid(time)
    return integrate_ode(lambda t, y: (forcing - restoring * y) / capacity, [initial], time, names=("temperature_anomaly",), units=("K",))
