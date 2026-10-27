"""Channels: uniform open-channel flow by Chezy and Manning, and Muskingum and kinematic-wave flood routing.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "hydraulic_radius",
    "chezy_velocity",
    "manning_velocity",
    "chezy_from_manning",
    "muskingum_coefficients",
    "muskingum_routing",
    "kinematic_wave_celerity",
    "simulate_kinematic_wave",
]


@model
def hydraulic_radius(width, depth):
    """Hydraulic radius of a rectangular channel.

    Friction acts on the wetted perimeter while gravity acts on the cross-section, so their ratio sets the resistance of
    a channel.

    Parameters
    ----------
    width : float or array_like
        m, bed width; constraint: positive.
    depth : float or array_like
        m, flow depth; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m], flow area over wetted perimeter. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rectangular section; tends to the depth for a wide channel.
    Equation: ``width * depth / (width + 2 * depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import hydraulic_radius
    >>> round(float(hydraulic_radius(4, 1)), 6)
    0.666667
    """
    width = quantity("width", width, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    return width * depth / (width + 2 * depth)


@model
def chezy_velocity(chezy_coefficient, hydraulic_radius, slope):
    """Chezy's mean velocity of uniform flow in an open channel.

    In uniform flow the weight of the water down the slope balances bed friction that grows as the square of the
    velocity.

    Parameters
    ----------
    chezy_coefficient : float or array_like
        m^0.5/s, typically 30-90; constraint: positive.
    hydraulic_radius : float or array_like
        m; constraint: positive.
    slope : float or array_like
        1, bed (energy) slope; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Chezy (1775): steady uniform turbulent flow; the coefficient is not a constant but depends on roughness and depth.
    Equation: ``chezy_coefficient * np.sqrt(hydraulic_radius * slope)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import chezy_velocity
    >>> float(chezy_velocity(50, 2, 0.02))
    10.0
    """
    chezy_coefficient = quantity("chezy_coefficient", chezy_coefficient, lower=0, strict_lower=True)
    hydraulic_radius = quantity("hydraulic_radius", hydraulic_radius, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0)
    return chezy_coefficient * np.sqrt(hydraulic_radius * slope)


@model
def manning_velocity(roughness, hydraulic_radius, slope):
    """Manning's mean velocity of uniform flow in an open channel.

    One roughness number per kind of channel bed predicts the velocity at any depth, which made Manning's formula the
    standard of channel design.

    Parameters
    ----------
    roughness : float or array_like
        s/m^(1/3), Manning's n (about 0.012 for concrete, 0.03-0.05 for natural rivers); constraint: positive.
    hydraulic_radius : float or array_like
        m; constraint: positive.
    slope : float or array_like
        1; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Manning (1891), SI form; fully rough turbulent flow, steady and uniform.
    Equation: ``hydraulic_radius**(2 / 3) * np.sqrt(slope) / roughness``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import manning_velocity
    >>> round(float(manning_velocity(0.04, 8, 0.0004)), 5)
    2.0
    """
    roughness = quantity("roughness", roughness, lower=0, strict_lower=True)
    hydraulic_radius = quantity("hydraulic_radius", hydraulic_radius, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0)
    return hydraulic_radius ** (2 / 3) * np.sqrt(slope) / roughness


@model
def chezy_from_manning(roughness, hydraulic_radius):
    """Chezy coefficient equivalent to a Manning roughness.

    Manning's formula is Chezy's with a coefficient that grows as the sixth root of the hydraulic radius.

    Parameters
    ----------
    roughness : float or array_like
        s/m^(1/3); constraint: positive.
    hydraulic_radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^0.5/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Manning (1891), SI units.
    Equation: ``hydraulic_radius**(1 / 6) / roughness``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import chezy_from_manning
    >>> float(chezy_from_manning(0.05, 64))
    40.0
    """
    roughness = quantity("roughness", roughness, lower=0, strict_lower=True)
    hydraulic_radius = quantity("hydraulic_radius", hydraulic_radius, lower=0, strict_lower=True)
    return hydraulic_radius ** (1 / 6) / roughness


@model
def muskingum_coefficients(storage_constant, weighting, time_step):
    """Muskingum routing coefficients.

    Storage in a river reach is a weighted mix of inflow and outflow, which makes the outflow a fixed linear combination of
    the latest flows.

    Parameters
    ----------
    storage_constant : float or array_like
        s, Muskingum K, about the travel time through the reach; constraint: positive.
    weighting : float or array_like
        1, Muskingum X; constraint: between 0 and 0.5.
    time_step : float or array_like
        s; constraint: positive.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [1]: ``[C0, C1, C2]``, which sum to one.

    Notes
    -----
    McCarthy (1938), storage ``K (X I + (1 - X) O)``; all coefficients are nonnegative when
    ``2 K X <= time_step <= 2 K (1 - X)``.
    Equation: ``[time_step - 2 K X, time_step + 2 K X, 2 K (1 - X) - time_step] / (2 K (1 - X) + time_step)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.hydrology import muskingum_coefficients
    >>> np.round(muskingum_coefficients(1, 0.25, 1), 6).tolist()
    [0.2, 0.6, 0.2]
    """
    storage_constant = quantity("storage_constant", storage_constant, lower=0, strict_lower=True)
    weighting = quantity("weighting", weighting, lower=0, upper=0.5)
    time_step = quantity("time_step", time_step, lower=0, strict_lower=True)
    denominator = 2 * storage_constant * (1 - weighting) + time_step
    return np.stack(
        np.broadcast_arrays(
            (time_step - 2 * storage_constant * weighting) / denominator,
            (time_step + 2 * storage_constant * weighting) / denominator,
            (2 * storage_constant * (1 - weighting) - time_step) / denominator,
        ),
        axis=-1,
    )


@model
def muskingum_routing(inflow, storage_constant, weighting, time_step):
    """Outflow hydrograph of a river reach by the Muskingum method.

    Routing a flood through a reach delays and flattens its peak, the effect forecasters need downstream.

    Parameters
    ----------
    inflow : array_like, shape (n,), n >= 2
        m^3/s, sampled every time step; constraint: nonnegative vector.
    storage_constant : float
        s; constraint: positive scalar.
    weighting : float
        1; constraint: scalar between 0 and 0.5.
    time_step : float
        s; constraint: positive scalar.

    Returns
    -------
    ndarray, shape (n,)
        Model value [m^3/s], starting from a steady state ``outflow[0] = inflow[0]``.

    Notes
    -----
    McCarthy (1938): ``O[j+1] = C0 I[j+1] + C1 I[j] + C2 O[j]``; choose the time step so that no coefficient is negative.
    Equation: ``O[j + 1] = C0 * I[j + 1] + C1 * I[j] + C2 * O[j]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.hydrology import muskingum_routing
    >>> np.round(muskingum_routing([1, 3, 2], 1, 0.25, 1), 5).tolist()
    [1.0, 1.4, 2.48]
    """
    inflow = quantity("inflow", inflow, lower=0)
    require(inflow.ndim == 1 and inflow.size >= 2, "inflow must be a vector with at least two samples")
    coefficients = muskingum_coefficients(storage_constant, weighting, time_step)
    require(coefficients.ndim == 1, "routing parameters must be scalars")
    c0, c1, c2 = coefficients
    outflow = np.empty_like(inflow)
    outflow[0] = inflow[0]
    for j in range(inflow.size - 1):
        outflow[j + 1] = c0 * inflow[j + 1] + c1 * inflow[j] + c2 * outflow[j]
    return outflow


@model
def kinematic_wave_celerity(velocity, exponent):
    """Speed at which a flood wave travels down a channel.

    A flood crest moves faster than the water in it, because deeper water flows faster: by 5/3 for Manning's law in a
    wide channel.

    Parameters
    ----------
    velocity : float or array_like
        m/s, the mean flow velocity; constraint: nonnegative.
    exponent : float or array_like
        1, ``m`` in the rating ``Q = a A**m`` (5/3 for Manning and 3/2 for Chezy in a wide channel); constraint:
        positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Lighthill and Whitham (1955): ``dQ/dA`` of a single-valued rating, valid where inertia and the pressure gradient
    are negligible against friction and the bed slope; kinematic waves translate without attenuating.
    Equation: ``exponent * velocity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import kinematic_wave_celerity
    >>> float(kinematic_wave_celerity(1.5, 1.6666666666666667))
    2.5
    """
    velocity = quantity("velocity", velocity, lower=0)
    exponent = quantity("exponent", exponent, lower=0, strict_lower=True)
    return exponent * velocity


@safe_arithmetic
def simulate_kinematic_wave(inflow, length, width, roughness, slope, time, *, nodes=100):
    """Route an inflow hydrograph down a wide channel as a kinematic wave.

    Parameters
    ----------
    inflow : array_like, shape (n,)
        Positive upstream discharge at each sample of time, interpolated linearly between them [m^3/s].
    length : float
        Positive scalar reach length [m].
    width : float
        Positive scalar channel width [m].
    roughness : float
        Positive scalar Manning n [s/m^(1/3)].
    slope : float
        Positive scalar bed slope [1].
    time : array_like, shape (n,)
        Nonnegative strictly increasing sampled times [s].
    nodes : int, optional
        Number of equal cells along the reach, at least 2.

    Returns
    -------
    TimeSeries
        Outflow at the downstream end [m^3/s], water stored in the reach [m^3] and the discharge of every cell
        [m^3/s], labeled discharge_<cell> from upstream.

    Notes
    -----
    Lighthill and Whitham (1955): ``dA/dt + dQ/dx = 0`` with Manning's uniform-flow rating for a wide rectangular
    channel, ``Q = np.sqrt(slope) / (roughness * width**(2 / 3)) * A**(5 / 3)`` (hydraulic radius equal to the depth),
    by the method of lines with first-order upwind differences and ``integrators.integrate_ode``. The reach starts in
    uniform flow at the first inflow. The upwind scheme adds numerical diffusion that attenuates the peak slightly;
    refine the cells to reduce it. There are no backwater effects, lateral inflow or floodplain storage.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode
    from earthsciencekit.results import TimeSeries

    time = time_grid(time)
    inflow = quantity("inflow", inflow, lower=0, strict_lower=True)
    require(inflow.shape == time.shape, "inflow must have one value per sample of time")
    length = quantity("length", length, lower=0, strict_lower=True)
    width = quantity("width", width, lower=0, strict_lower=True)
    roughness = quantity("roughness", roughness, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0, strict_lower=True)
    require(all(v.ndim == 0 for v in (length, width, roughness, slope)), "channel parameters must be scalars")
    require(isinstance(nodes, (int, np.integer)) and not isinstance(nodes, bool) and nodes >= 2, "nodes must be an integer >= 2")
    spacing = length / nodes
    rating = np.sqrt(slope) / (roughness * width ** (2 / 3))

    def discharge(area):
        return rating * np.maximum(area, 0.0) ** (5 / 3)

    def rhs(t, area):
        upstream = np.interp(t, time, inflow)
        return -np.diff(np.concatenate([[upstream], discharge(area)])) / spacing

    initial = np.full(nodes, (inflow[0] / rating) ** (3 / 5))
    area = integrate_ode(rhs, initial, time, rtol=1e-8, atol=1e-8 * initial[0]).values
    flow = discharge(area)
    names = ("outflow", "storage") + tuple(f"discharge_{i}" for i in range(nodes))
    units = ("m^3/s", "m^3") + ("m^3/s",) * nodes
    return TimeSeries(time, np.column_stack([flow[:, -1], area.sum(axis=1) * spacing, flow]), names, units)
