"""Geochemistry: Rayleigh fractionation, delta notation, melt crystallization, lattice strain and the oxygen-isotope thermometer.

SI units unless stated otherwise; isotope delta values are per mil.
"""

import numpy as np

from earthsciencekit._validation import model, quantity
from earthsciencekit.constants import AVOGADRO, R_GAS

__all__ = [
    "rayleigh_fractionation",
    "delta_value",
    "fractional_crystallization_melt",
    "equilibrium_crystallization_melt",
    "lattice_strain_partition",
    "isotope_temperature",
]


@model
def rayleigh_fractionation(initial_ratio, remaining_fraction, fractionation_factor):
    """Isotope ratio of a reservoir depleted by Rayleigh distillation.

    Removing a product that is always slightly enriched in one isotope drives the remaining reservoir further and further
    from its starting composition.

    Parameters
    ----------
    initial_ratio : float or array_like
        1, heavy-to-light isotope ratio of the reservoir at the start; constraint: positive.
    remaining_fraction : float or array_like
        1, fraction of the reservoir left; constraint: fraction above zero.
    fractionation_factor : float or array_like
        1, product ratio over reservoir ratio at each instant; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rayleigh (1896): the product leaves as soon as it forms, the reservoir stays well mixed and the factor is constant.
    Equation: ``initial_ratio * remaining_fraction**(fractionation_factor - 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import rayleigh_fractionation
    >>> float(rayleigh_fractionation(2, 0.25, 0.5))
    4.0
    """
    initial_ratio = quantity("initial_ratio", initial_ratio, lower=0, strict_lower=True)
    remaining_fraction = quantity("remaining_fraction", remaining_fraction, lower=0, upper=1, strict_lower=True)
    fractionation_factor = quantity("fractionation_factor", fractionation_factor, lower=0, strict_lower=True)
    return initial_ratio * remaining_fraction ** (fractionation_factor - 1)


@model
def delta_value(ratio, standard_ratio):
    """Isotope delta value relative to a standard.

    Isotope ratios vary in the fourth decimal place, so geochemists report deviations from a standard in per mil.

    Parameters
    ----------
    ratio : float or array_like
        1, heavy-to-light isotope ratio of the sample; constraint: positive.
    standard_ratio : float or array_like
        1, the same ratio in the standard (e.g. VSMOW, VPDB); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [per mil]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Linear in the ratio; delta values of different standards must be converted before they are compared.
    Equation: ``(ratio / standard_ratio - 1) * 1000``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import delta_value
    >>> round(float(delta_value(1.01, 1)), 4)
    10.0
    """
    ratio = quantity("ratio", ratio, lower=0, strict_lower=True)
    standard_ratio = quantity("standard_ratio", standard_ratio, lower=0, strict_lower=True)
    return (ratio / standard_ratio - 1) * 1000


@model
def fractional_crystallization_melt(initial_concentration, melt_fraction, partition_coefficient):
    """Trace-element concentration in a melt during fractional crystallization.

    Crystals removed as soon as they form leave the residual melt enriched in the elements they reject.

    Parameters
    ----------
    initial_concentration : float or array_like
        Concentration in the parent melt (any unit, e.g. ppm); constraint: nonnegative.
    melt_fraction : float or array_like
        1, fraction of melt remaining; constraint: fraction above zero.
    partition_coefficient : float or array_like
        1, bulk solid/melt concentration ratio; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [unit of initial_concentration]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rayleigh fractional crystallization (Bowen 1928; Gast 1968): crystals do not re-equilibrate and the bulk partition
    coefficient is constant.
    Equation: ``initial_concentration * melt_fraction**(partition_coefficient - 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import fractional_crystallization_melt
    >>> float(fractional_crystallization_melt(10, 0.25, 0.5))
    20.0
    """
    initial_concentration = quantity("initial_concentration", initial_concentration, lower=0)
    melt_fraction = quantity("melt_fraction", melt_fraction, lower=0, upper=1, strict_lower=True)
    partition_coefficient = quantity("partition_coefficient", partition_coefficient, lower=0)
    return initial_concentration * melt_fraction ** (partition_coefficient - 1)


@model
def equilibrium_crystallization_melt(initial_concentration, melt_fraction, partition_coefficient):
    """Trace-element concentration in a melt during equilibrium crystallization.

    Crystals that keep reacting with the melt, as in Bowen's continuous series, limit how far the melt can evolve.

    Parameters
    ----------
    initial_concentration : float or array_like
        Concentration in the parent melt (any unit, e.g. ppm); constraint: nonnegative.
    melt_fraction : float or array_like
        1, fraction of melt remaining; constraint: fraction.
    partition_coefficient : float or array_like
        1, bulk solid/melt concentration ratio; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [unit of initial_concentration]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Batch mass balance (Shaw 1970): solid and melt stay in equilibrium, so an incompatible element is enriched at most by
    1 / partition_coefficient.
    Equation: ``initial_concentration / (partition_coefficient + melt_fraction * (1 - partition_coefficient))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import equilibrium_crystallization_melt
    >>> round(float(equilibrium_crystallization_melt(10, 0.5, 0.5)), 4)
    13.3333
    """
    initial_concentration = quantity("initial_concentration", initial_concentration, lower=0)
    melt_fraction = quantity("melt_fraction", melt_fraction, lower=0, upper=1)
    partition_coefficient = quantity("partition_coefficient", partition_coefficient, lower=0, strict_lower=True)
    return initial_concentration / (partition_coefficient + melt_fraction * (1 - partition_coefficient))


@model
def lattice_strain_partition(ionic_radius, optimal_radius, peak_partition, youngs_modulus, temperature):
    """Mineral-melt partition coefficient from the strain of the crystal lattice.

    Ions close in size to the site they replace enter a crystal easily, which is Goldschmidt's first rule of substitution.

    Parameters
    ----------
    ionic_radius : float or array_like
        m; constraint: positive.
    optimal_radius : float or array_like
        m, radius that fits the site without strain; constraint: positive.
    peak_partition : float or array_like
        1, partition coefficient of an ion of the optimal radius; constraint: positive.
    youngs_modulus : float or array_like
        Pa, apparent Young's modulus of the site; constraint: positive.
    temperature : float or array_like
        K; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Blundy and Wood (1994) quantification of Goldschmidt (1937) for isovalent ions on one site; charge differences
    need a separate parabola.
    Equation: ``peak_partition * np.exp(-4 * np.pi * youngs_modulus * AVOGADRO * (optimal_radius / 2 * (ionic_radius - optimal_radius)**2
    + (ionic_radius - optimal_radius)**3 / 3) / (R_GAS * temperature))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import lattice_strain_partition
    >>> round(float(lattice_strain_partition(2e-10, 1e-10, 3, 1318422678.149293, 1000)), 5)
    1.10364
    """
    ionic_radius = quantity("ionic_radius", ionic_radius, lower=0, strict_lower=True)
    optimal_radius = quantity("optimal_radius", optimal_radius, lower=0, strict_lower=True)
    peak_partition = quantity("peak_partition", peak_partition, lower=0, strict_lower=True)
    youngs_modulus = quantity("youngs_modulus", youngs_modulus, lower=0, strict_lower=True)
    temperature = quantity("temperature", temperature, lower=0, strict_lower=True)
    misfit = ionic_radius - optimal_radius
    strain = 4 * np.pi * youngs_modulus * AVOGADRO * (optimal_radius / 2 * misfit**2 + misfit**3 / 3)
    return peak_partition * np.exp(-strain / (R_GAS * temperature))


@model
def isotope_temperature(delta_calcite, delta_water):
    """Calcite-water oxygen-isotope paleotemperature.

    The oxygen-isotope fractionation between a shell and the water it grew in falls with temperature, so fossil shells
    record ancient ocean temperatures.

    Parameters
    ----------
    delta_calcite : float or array_like
        per mil, delta18O of the carbonate (CO2 from the calcite, PDB scale); constraint: real.
    delta_water : float or array_like
        per mil, delta18O of the water on the same scale; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Urey (1947) predicted the effect; this is the empirical calibration of Epstein, Buchsbaum, Lowenstam and Urey (1953),
    valid for about 7-30 degC and equilibrium calcite.
    Equation: ``273.15 + 16.5 - 4.3 * d + 0.14 * d**2`` with ``d = delta_calcite - delta_water``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.volcanology import isotope_temperature
    >>> round(float(isotope_temperature(1, 0)), 3)
    285.49
    """
    delta_calcite = quantity("delta_calcite", delta_calcite)
    delta_water = quantity("delta_water", delta_water)
    difference = delta_calcite - delta_water
    return 273.15 + 16.5 - 4.3 * difference + 0.14 * difference**2
