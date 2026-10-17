"""Radiative forcing of the well-mixed greenhouse gases.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity

__all__ = [
    "co2_forcing",
    "methane_forcing",
    "nitrous_oxide_forcing",
]


@model
def co2_forcing(concentration, reference_concentration, coefficient):
    """Logarithmic CO2 forcing approximation.

    A logarithmic approximation connects concentration changes to a forcing perturbation.

    Parameters
    ----------
    concentration : float or array_like
        ppm; constraint: positive.
    reference_concentration : float or array_like
        ppm; constraint: positive.
    coefficient : float or array_like
        W/m^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    User-supplied coefficient; an educational parameterization, not a current assessed forcing product.
    Equation: ``coefficient * np.log(concentration / reference_concentration)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    concentration = quantity("concentration", concentration, lower=0, strict_lower=True)
    reference_concentration = quantity("reference_concentration", reference_concentration, lower=0, strict_lower=True)
    coefficient = quantity("coefficient", coefficient, lower=0, strict_lower=True)
    return coefficient * np.log(concentration / reference_concentration)


def _methane_nitrous_overlap(methane, nitrous_oxide):
    """Myhre et al. (1998) CH4-N2O band overlap in W/m^2; concentrations in ppb."""
    product = methane * nitrous_oxide
    return 0.47 * np.log(1 + 2.01e-5 * product**0.75 + 5.31e-15 * methane * product**1.52)


@model
def methane_forcing(concentration, reference_concentration, nitrous_oxide):
    """Methane radiative forcing (Myhre et al. 1998).

    Methane forcing grows as a square root of concentration, less the band it shares with nitrous oxide.

    Parameters
    ----------
    concentration : float or array_like
        ppb; constraint: positive.
    reference_concentration : float or array_like
        ppb; constraint: positive.
    nitrous_oxide : float or array_like
        ppb, reference N2O concentration; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Simplified fit to line-by-line calculations, valid near present-day concentrations; with
    ``f(M, N) = 0.47 * log(1 + 2.01e-5 * (M * N)**0.75 + 5.31e-15 * M * (M * N)**1.52)``.
    Equation: ``0.036 * (sqrt(M) - sqrt(M0)) - (f(M, N0) - f(M0, N0))``, with ``M = concentration``,
    ``M0 = reference_concentration`` and ``N0 = nitrous_oxide``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    concentration = quantity("concentration", concentration, lower=0, strict_lower=True)
    reference_concentration = quantity("reference_concentration", reference_concentration, lower=0, strict_lower=True)
    nitrous_oxide = quantity("nitrous_oxide", nitrous_oxide, lower=0, strict_lower=True)
    overlap = _methane_nitrous_overlap(concentration, nitrous_oxide) - _methane_nitrous_overlap(reference_concentration, nitrous_oxide)
    return 0.036 * (np.sqrt(concentration) - np.sqrt(reference_concentration)) - overlap


@model
def nitrous_oxide_forcing(concentration, reference_concentration, methane):
    """Nitrous oxide radiative forcing (Myhre et al. 1998).

    Nitrous oxide forcing grows as a square root of concentration, less the band it shares with methane.

    Parameters
    ----------
    concentration : float or array_like
        ppb; constraint: positive.
    reference_concentration : float or array_like
        ppb; constraint: positive.
    methane : float or array_like
        ppb, reference CH4 concentration; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Simplified fit to line-by-line calculations, valid near present-day concentrations; f is the overlap term of methane_forcing.
    Equation: ``0.12 * (sqrt(N) - sqrt(N0)) - (f(M0, N) - f(M0, N0))``, with ``N = concentration``,
    ``N0 = reference_concentration`` and ``M0 = methane``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    concentration = quantity("concentration", concentration, lower=0, strict_lower=True)
    reference_concentration = quantity("reference_concentration", reference_concentration, lower=0, strict_lower=True)
    methane = quantity("methane", methane, lower=0, strict_lower=True)
    overlap = _methane_nitrous_overlap(methane, concentration) - _methane_nitrous_overlap(methane, reference_concentration)
    return 0.12 * (np.sqrt(concentration) - np.sqrt(reference_concentration)) - overlap
