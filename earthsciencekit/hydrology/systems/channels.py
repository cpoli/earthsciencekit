"""Channels: uniform open-channel flow by Chezy and Manning, and Muskingum flood routing.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "hydraulic_radius",
    "chezy_velocity",
    "manning_velocity",
    "chezy_from_manning",
    "muskingum_coefficients",
    "muskingum_routing",
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
