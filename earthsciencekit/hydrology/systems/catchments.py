"""Catchments: the water balance, unit hydrographs, Horton's stream laws and TOPMODEL's topographic index.

SI units unless stated otherwise.
"""

import numpy as np
from scipy.special import gammaln

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "water_balance",
    "runoff_coefficient",
    "nash_unit_hydrograph",
    "unit_hydrograph_convolution",
    "stream_number",
    "stream_length",
    "topographic_index",
    "topmodel_deficit",
    "topmodel_baseflow",
]


@model
def water_balance(precipitation, evapotranspiration, runoff):
    """Rate of change of the water stored in a catchment.

    Whatever falls on a catchment and neither evaporates nor leaves by its river must be stored in its soil, aquifers and
    lakes.

    Parameters
    ----------
    precipitation : float or array_like
        m/s, depth per unit time over the catchment; constraint: nonnegative.
    evapotranspiration : float or array_like
        m/s; constraint: nonnegative.
    runoff : float or array_like
        m/s, river discharge divided by the catchment area; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], positive when storage grows. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Perrault (1674): a closed catchment with no groundwater exchange across its divide; any common depth rate (such as
    mm/day) may be used throughout.
    Equation: ``precipitation - evapotranspiration - runoff``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import water_balance
    >>> float(water_balance(3, 1, 0.5))
    1.5
    """
    precipitation = quantity("precipitation", precipitation, lower=0)
    evapotranspiration = quantity("evapotranspiration", evapotranspiration, lower=0)
    runoff = quantity("runoff", runoff, lower=0)
    return precipitation - evapotranspiration - runoff


@model
def runoff_coefficient(runoff, precipitation):
    """Fraction of the precipitation that leaves a catchment as river flow.

    Perrault's measurement that the Seine carries only about a sixth of the rain on its basin showed that rainfall is
    enough to feed rivers.

    Parameters
    ----------
    runoff : float or array_like
        m, runoff depth over a period; constraint: nonnegative.
    precipitation : float or array_like
        m, precipitation depth over the same period; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Both depths over the same area and period; long periods make the storage change negligible.
    Equation: ``runoff / precipitation``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import runoff_coefficient
    >>> round(float(runoff_coefficient(1, 6)), 6)
    0.166667
    """
    runoff = quantity("runoff", runoff, lower=0)
    precipitation = quantity("precipitation", precipitation, lower=0, strict_lower=True)
    return runoff / precipitation


@model
def nash_unit_hydrograph(time, reservoirs, storage_constant):
    """Instantaneous unit hydrograph of a cascade of linear reservoirs.

    A catchment that drains like a chain of equal linear reservoirs turns an instant of rain into a skewed, gamma-shaped
    pulse of river flow.

    Parameters
    ----------
    time : float or array_like
        s, since the instant of excess rainfall; constraint: nonnegative.
    reservoirs : float or array_like
        1, Nash's n; constraint: at least 1.
    storage_constant : float or array_like
        s, Nash's k of each reservoir; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1/s], which integrates to one; multiply by the catchment area for m^3/s per metre of excess
        rainfall. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nash (1957), a gamma density with mean lag ``reservoirs * storage_constant``; a smooth synthetic shape for Sherman's
    unit hydrograph, linear and time-invariant.
    Equation: ``time**(reservoirs - 1) * np.exp(-time / storage_constant) / (storage_constant**reservoirs *
    scipy.special.gamma(reservoirs))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import nash_unit_hydrograph
    >>> round(float(nash_unit_hydrograph(2, 2, 1)), 6)
    0.270671
    """
    time = quantity("time", time, lower=0)
    reservoirs = quantity("reservoirs", reservoirs, lower=1)
    storage_constant = quantity("storage_constant", storage_constant, lower=0, strict_lower=True)
    positive = np.where(time > 0, time, 1.0)
    log_density = (reservoirs - 1) * np.log(positive) - positive / storage_constant - reservoirs * np.log(storage_constant)
    density = np.exp(log_density - gammaln(reservoirs))
    return np.where(time > 0, density, np.where(reservoirs == 1, 1 / storage_constant, 0.0))[()]


@model
def unit_hydrograph_convolution(excess_rainfall, unit_hydrograph):
    """Storm hydrograph from excess rainfall and a unit hydrograph.

    If a catchment responds linearly, the flood from any storm is the sum of unit-hydrograph responses to its successive
    bursts of excess rainfall, each shifted and scaled.

    Parameters
    ----------
    excess_rainfall : array_like, shape (n,)
        m, excess (effective) rainfall depth in each interval; constraint: nonnegative vector.
    unit_hydrograph : array_like, shape (m,)
        m^3/s per metre of excess rainfall, ordinates at the same interval; constraint: nonnegative vector.

    Returns
    -------
    ndarray, shape (n + m - 1,)
        Model value [m^3/s], the direct runoff at the same interval.

    Notes
    -----
    Sherman (1932): linear superposition and time invariance; baseflow and losses are removed beforehand.
    Equation: ``np.convolve(excess_rainfall, unit_hydrograph)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.hydrology import unit_hydrograph_convolution
    >>> np.round(unit_hydrograph_convolution([1, 2], [3, 4, 5]), 4).tolist()
    [3.0, 10.0, 13.0, 10.0]
    """
    excess_rainfall = quantity("excess_rainfall", excess_rainfall, lower=0)
    unit_hydrograph = quantity("unit_hydrograph", unit_hydrograph, lower=0)
    require(excess_rainfall.ndim == 1 and unit_hydrograph.ndim == 1, "excess rainfall and unit hydrograph must be vectors")
    return np.convolve(excess_rainfall, unit_hydrograph)


@model
def stream_number(order, highest_order, bifurcation_ratio):
    """Number of streams of a given order in a drainage network.

    Horton found that stream counts fall geometrically with order, so a whole river network is summarized by one
    bifurcation ratio.

    Parameters
    ----------
    order : float or array_like
        1, stream order; constraint: at least 1.
    highest_order : float or array_like
        1, order of the main stream; constraint: at least order.
    bifurcation_ratio : float or array_like
        1, typically 3-5; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], one stream of the highest order. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Horton's law of stream numbers (1945), an empirical average; Strahler's (1952) ordering is the usual convention.
    Equation: ``bifurcation_ratio**(highest_order - order)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import stream_number
    >>> float(stream_number(1, 3, 4))
    16.0
    """
    order = quantity("order", order, lower=1)
    highest_order = quantity("highest_order", highest_order, lower=1)
    bifurcation_ratio = quantity("bifurcation_ratio", bifurcation_ratio, lower=0, strict_lower=True)
    require(highest_order >= order, "highest order must be at least the stream order")
    return bifurcation_ratio ** (highest_order - order)


@model
def stream_length(order, first_order_length, length_ratio):
    """Mean length of the streams of a given order.

    Streams grow longer geometrically with order, Horton's companion to his law of stream numbers.

    Parameters
    ----------
    order : float or array_like
        1, stream order; constraint: at least 1.
    first_order_length : float or array_like
        m, mean length of first-order streams; constraint: positive.
    length_ratio : float or array_like
        1, typically 1.5-3.5; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Horton's law of stream lengths (1945), an empirical average over a network.
    Equation: ``first_order_length * length_ratio**(order - 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import stream_length
    >>> float(stream_length(3, 100, 2))
    400.0
    """
    order = quantity("order", order, lower=1)
    first_order_length = quantity("first_order_length", first_order_length, lower=0, strict_lower=True)
    length_ratio = quantity("length_ratio", length_ratio, lower=0, strict_lower=True)
    return first_order_length * length_ratio ** (order - 1)


@model
def topographic_index(specific_area, slope):
    """TOPMODEL topographic wetness index.

    Ground that drains a large upslope area through a gentle slope is the first to saturate, so one number from the
    topography maps where runoff starts.

    Parameters
    ----------
    specific_area : float or array_like
        m, upslope area drained per unit contour width; constraint: positive.
    slope : float or array_like
        rad, local surface slope; constraint: between 0 and pi/2, exclusive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], with the specific area in metres. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Beven and Kirkby (1979): steady recharge, a water table parallel to the surface and transmissivity falling
    exponentially with deficit.
    Equation: ``np.log(specific_area / np.tan(slope))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import topographic_index
    >>> float(topographic_index(2.718281828459045, 0.7853981633974483))
    1.0
    """
    specific_area = quantity("specific_area", specific_area, lower=0, strict_lower=True)
    slope = quantity("slope", slope, lower=0, upper=np.pi / 2, strict_lower=True, strict_upper=True)
    return np.log(specific_area / np.tan(slope))


@model
def topmodel_deficit(mean_deficit, topographic_index, mean_index, scaling_depth):
    """Local storage deficit of TOPMODEL from the catchment mean.

    Each point's deficit departs from the catchment mean in proportion to how much wetter its topographic index is than
    average, so one state variable maps the saturated area.

    Parameters
    ----------
    mean_deficit : float or array_like
        m, catchment-mean storage deficit; constraint: nonnegative.
    topographic_index : float or array_like
        1, local ``ln(a / tan(beta))``; constraint: real.
    mean_index : float or array_like
        1, catchment mean of the index; constraint: real.
    scaling_depth : float or array_like
        m, TOPMODEL's m, the deficit that reduces transmissivity by a factor e; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]; zero or negative where the ground is saturated (negative values are return flow). Scalar input
        returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Beven and Kirkby (1979) with uniform transmissivity at saturation and quasi-steady drainage.
    Equation: ``mean_deficit + scaling_depth * (mean_index - topographic_index)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import topmodel_deficit
    >>> float(topmodel_deficit(0.1, 7, 6, 0.05))
    0.05
    """
    mean_deficit = quantity("mean_deficit", mean_deficit, lower=0)
    topographic_index = quantity("topographic_index", topographic_index)
    mean_index = quantity("mean_index", mean_index)
    scaling_depth = quantity("scaling_depth", scaling_depth, lower=0, strict_lower=True)
    return mean_deficit + scaling_depth * (mean_index - topographic_index)


@model
def topmodel_baseflow(saturated_transmissivity, mean_index, mean_deficit, scaling_depth):
    """TOPMODEL subsurface drainage from the catchment-mean deficit.

    Baseflow falls exponentially as the catchment dries, the recession behaviour TOPMODEL is built on.

    Parameters
    ----------
    saturated_transmissivity : float or array_like
        m^2/s, transmissivity when the soil is saturated to the surface; constraint: positive.
    mean_index : float or array_like
        1, catchment mean of ``ln(a / tan(beta))`` with a in metres; constraint: real.
    mean_deficit : float or array_like
        m; constraint: nonnegative.
    scaling_depth : float or array_like
        m, TOPMODEL's m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], discharge per unit catchment area. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Beven and Kirkby (1979): exponential transmissivity profile, water table parallel to the surface.
    Equation: ``saturated_transmissivity * np.exp(-mean_index) * np.exp(-mean_deficit / scaling_depth)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.hydrology import topmodel_baseflow
    >>> round(float(topmodel_baseflow(2, 0, 1, 1)), 6)
    0.735759
    """
    saturated_transmissivity = quantity("saturated_transmissivity", saturated_transmissivity, lower=0, strict_lower=True)
    mean_index = quantity("mean_index", mean_index)
    mean_deficit = quantity("mean_deficit", mean_deficit, lower=0)
    scaling_depth = quantity("scaling_depth", scaling_depth, lower=0, strict_lower=True)
    return saturated_transmissivity * np.exp(-mean_index) * np.exp(-mean_deficit / scaling_depth)
