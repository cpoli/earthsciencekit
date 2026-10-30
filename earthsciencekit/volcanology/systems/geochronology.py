"""Geochronology: radioactive decay, the age equation, radiocarbon, lead-lead ages, the concordia and the isochron.

SI units unless stated otherwise; ages and half-lives are seconds and decay constants 1/s.
"""

import numpy as np
from scipy.optimize import brentq

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import LAMBDA_U235, LAMBDA_U238, LIBBY_MEAN_LIFE, SECONDS_PER_YEAR, U238_U235

__all__ = [
    "decay_constant",
    "parent_fraction",
    "radiometric_age",
    "radiogenic_daughter",
    "radiocarbon_age",
    "radiogenic_lead_ratio",
    "lead_lead_age",
    "concordia_ratios",
    "isochron_fit",
]

_YOUNGEST, _OLDEST = SECONDS_PER_YEAR, 1e10 * SECONDS_PER_YEAR  # s, the bracket of lead_lead_age


@model
def decay_constant(half_life):
    """Decay constant from a half-life.

    Every radioactive nuclide has a fixed probability per unit time of decaying, set by its half-life.

    Parameters
    ----------
    half_life : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rutherford and Soddy (1902): first-order decay, independent of temperature, pressure and chemistry.
    Equation: ``np.log(2) / half_life``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import decay_constant
    >>> round(float(decay_constant(2)), 6)
    0.346574
    """
    half_life = quantity("half_life", half_life, lower=0, strict_lower=True)
    return np.log(2) / half_life


@model
def parent_fraction(decay_constant, time):
    """Fraction of radioactive parent atoms remaining.

    Exponential decay is the clock behind every radiometric age.

    Parameters
    ----------
    decay_constant : float or array_like
        1/s; constraint: positive.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rutherford and Soddy (1902): a closed system with no parent added after time zero.
    Equation: ``np.exp(-decay_constant * time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import parent_fraction
    >>> float(parent_fraction(0.6931471805599453, 1))
    0.5
    """
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0)
    return np.exp(-decay_constant * time)


@model
def radiometric_age(daughter_parent_ratio, decay_constant):
    """Age from a radiogenic daughter-to-parent ratio.

    The ratio of accumulated daughter atoms to surviving parent atoms grows with time, so measuring it dates the rock.

    Parameters
    ----------
    daughter_parent_ratio : float or array_like
        1, radiogenic daughter atoms per parent atom; constraint: nonnegative.
    decay_constant : float or array_like
        1/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Closed system with no initial daughter; Boltwood (1907) applied it to lead in uranium minerals.
    Equation: ``np.log1p(daughter_parent_ratio) / decay_constant``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import radiometric_age
    >>> float(radiometric_age(1.718281828459045, 2))
    0.5
    """
    daughter_parent_ratio = quantity("daughter_parent_ratio", daughter_parent_ratio, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    return np.log1p(daughter_parent_ratio) / decay_constant


@model
def radiogenic_daughter(parent, decay_constant, time):
    """Radiogenic daughter accumulated beside the surviving parent.

    Forward accumulation shows how the daughter-to-parent ratio of a rock records its age.

    Parameters
    ----------
    parent : float or array_like
        Surviving parent amount today (atoms, moles or any consistent unit); constraint: nonnegative.
    decay_constant : float or array_like
        1/s; constraint: positive.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [unit of parent]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Closed system with no initial daughter; the inverse of radiometric_age. Holmes (1913) used it to build the first
    radiometric time scale.
    Equation: ``parent * np.expm1(decay_constant * time)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import radiogenic_daughter
    >>> float(radiogenic_daughter(3, 0.6931471805599453, 1))
    3.0
    """
    parent = quantity("parent", parent, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    time = quantity("time", time, lower=0)
    return parent * np.expm1(decay_constant * time)


@model
def radiocarbon_age(fraction_modern):
    """Conventional radiocarbon age.

    Living matter exchanges carbon-14 with the atmosphere; after death the carbon-14 decays, dating organic remains.

    Parameters
    ----------
    fraction_modern : float or array_like
        1, sample 14C/12C relative to the modern standard; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]; negative for bomb carbon above the standard. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Libby et al. (1949) with the conventional 5568-year (Libby) half-life, mean life LIBBY_MEAN_LIFE; it is a radiocarbon
    age, not a calendar age, until calibrated.
    Equation: ``-LIBBY_MEAN_LIFE * np.log(fraction_modern)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import radiocarbon_age
    >>> round(float(radiocarbon_age(0.36787944117144233)), 0)
    253502200800.0
    """
    fraction_modern = quantity("fraction_modern", fraction_modern, lower=0, strict_lower=True)
    return -LIBBY_MEAN_LIFE * np.log(fraction_modern)


@model
def radiogenic_lead_ratio(time):
    """Radiogenic 207Pb/206Pb ratio produced over a given age.

    Two uranium isotopes decaying at different rates make the ratio of their lead daughters a clock that needs no uranium
    measurement.

    Parameters
    ----------
    time : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Patterson (1956): closed system with today's 238U/235U = U238_U235 and the Jaffey et al. (1971) decay constants. The
    ratio tends to LAMBDA_U235 / (U238_U235 * LAMBDA_U238) at zero age.
    Equation: ``np.expm1(LAMBDA_U235 * time) / (U238_U235 * np.expm1(LAMBDA_U238 * time))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import radiogenic_lead_ratio
    >>> round(float(radiogenic_lead_ratio(1.4100926004988576e+17)), 6)
    0.58385
    """
    time = quantity("time", time, lower=0, strict_lower=True)
    return np.expm1(LAMBDA_U235 * time) / (U238_U235 * np.expm1(LAMBDA_U238 * time))


@model
def lead_lead_age(lead_ratio):
    """Age from a radiogenic 207Pb/206Pb ratio.

    Patterson dated the Earth from the slope of a meteorite lead isochron, which equals this ratio.

    Parameters
    ----------
    lead_ratio : float or array_like
        1, radiogenic 207Pb*/206Pb* or the slope of a Pb-Pb isochron; constraint: between the ratios at one year
        (about 0.0461) and at 10 Ga.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Inverts radiogenic_lead_ratio by Brent's method between one year and 10 Ga; same closed-system assumptions.
    Equation: ``radiogenic_lead_ratio(age) = lead_ratio``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import lead_lead_age
    >>> round(float(lead_lead_age(0.5838503844760619)), -10)
    1.4100926e+17
    """
    lead_ratio = quantity("lead_ratio", lead_ratio)
    require(lead_ratio >= radiogenic_lead_ratio(_YOUNGEST), "lead_ratio is below the ratio at one year")
    require(lead_ratio <= radiogenic_lead_ratio(_OLDEST), "lead_ratio exceeds the ratio at 10 Ga")

    def solve(ratio):
        return brentq(lambda age: radiogenic_lead_ratio(age) - ratio, _YOUNGEST, _OLDEST, xtol=1e-6 * SECONDS_PER_YEAR, rtol=1e-14)

    return np.vectorize(solve, otypes=[float])(lead_ratio)[()]


@model
def concordia_ratios(time):
    """Point on the Wetherill U-Pb concordia.

    Both uranium-lead clocks agree only on this curve, so a zircon off it has lost lead and its chord dates the loss.

    Parameters
    ----------
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [1]: ``[207Pb*/235U, 206Pb*/238U]``, the horizontal and vertical concordia coordinates.

    Notes
    -----
    Wetherill (1956): closed system with the Jaffey et al. (1971) decay constants. Episodic lead loss puts samples on a
    straight chord (discordia) between the concordia points of formation and loss.
    Equation: ``[np.expm1(LAMBDA_U235 * time), np.expm1(LAMBDA_U238 * time)]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.volcanology import concordia_ratios
    >>> np.round(concordia_ratios(1.4100926004988576e+17), 4).tolist()
    [80.5013, 1.0]
    """
    time = quantity("time", time, lower=0)
    return np.stack([np.expm1(LAMBDA_U235 * time), np.expm1(LAMBDA_U238 * time)], axis=-1)


@model
def isochron_fit(parent_ratio, daughter_ratio, decay_constant):
    """Age and initial ratio from a whole-rock or mineral isochron.

    Cogenetic minerals with different parent contents line up on a straight line whose slope gives the age without knowing
    the initial daughter.

    Parameters
    ----------
    parent_ratio : array_like, shape (n,), n >= 2
        1, parent over a stable reference isotope (e.g. 87Rb/86Sr); constraint: nonnegative, not all equal.
    daughter_ratio : array_like, shape (n,)
        1, daughter over the same reference isotope (e.g. 87Sr/86Sr); constraint: nonnegative.
    decay_constant : float
        1/s; constraint: positive scalar.

    Returns
    -------
    numpy.ndarray, shape (2,)
        ``[age, initial_ratio]``: the age [s] from the least-squares slope and the intercept [1].

    Notes
    -----
    Nicolaysen (1961): every sample formed at the same time with the same initial ratio and stayed closed; ordinary least
    squares without measurement-error weighting (a York fit would weight them). The slope must be positive.
    Equation: ``daughter_ratio = initial_ratio + parent_ratio * np.expm1(decay_constant * age)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.volcanology import isochron_fit
    >>> np.round(isochron_fit([0, 1, 2], [0.7, 2.418281828459045, 4.13656365691809], 1), 5).tolist()
    [1.0, 0.7]
    """
    parent_ratio = quantity("parent_ratio", parent_ratio, lower=0)
    daughter_ratio = quantity("daughter_ratio", daughter_ratio, lower=0)
    decay_constant = quantity("decay_constant", decay_constant, lower=0, strict_lower=True)
    require(parent_ratio.ndim == 1 and parent_ratio.size >= 2, "parent_ratio must be a vector of at least two samples")
    require(daughter_ratio.shape == parent_ratio.shape, "daughter_ratio must match parent_ratio")
    require(decay_constant.ndim == 0, "decay_constant must be a scalar")
    anomaly = parent_ratio - parent_ratio.mean()
    require(np.any(anomaly != 0), "parent_ratio must vary")
    slope = np.dot(anomaly, daughter_ratio - daughter_ratio.mean()) / np.dot(anomaly, anomaly)
    require(slope > 0, "isochron slope must be positive")
    return np.array([np.log1p(slope) / decay_constant, daughter_ratio.mean() - slope * parent_ratio.mean()])
