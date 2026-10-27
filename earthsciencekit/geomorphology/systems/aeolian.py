"""Aeolian transport: the threshold of blown sand, saltation flux and dune migration.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = ["threshold_shear_velocity", "bagnold_sand_flux", "kawamura_sand_flux", "dune_migration_rate"]


@model
def threshold_shear_velocity(diameter, sediment_density, air_density, gravity, coefficient):
    """Fluid threshold shear velocity for wind to lift sand.

    Wind starts to move sand once the drag on a grain overcomes its weight, which takes a stronger wind for coarser
    grains.

    Parameters
    ----------
    diameter : float or array_like
        m; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than air_density.
    air_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    coefficient : float or array_like
        1, about 0.1 for the fluid threshold and 0.08 for the lower impact threshold; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Bagnold (1941): dry cohesionless sand coarser than about 0.1 mm; moisture, crusts and fine grains raise the
    threshold.
    Equation: ``coefficient * np.sqrt((sediment_density - air_density) / air_density * gravity * diameter)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import threshold_shear_velocity
    >>> round(float(threshold_shear_velocity(0.00025, 2650, 1.2, 9.81, 0.1)), 12)
    0.232668970858
    """
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    air_density = quantity("air_density", air_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    coefficient = quantity("coefficient", coefficient, lower=0, strict_lower=True)
    require(sediment_density > air_density, "sediment must be denser than the air")
    return coefficient * np.sqrt((sediment_density - air_density) / air_density * gravity * diameter)


@model
def bagnold_sand_flux(shear_velocity, diameter, air_density, gravity, coefficient):
    """Bagnold's saltation mass flux of blown sand.

    The mass of sand carried by the wind grows as the cube of the shear velocity, so a few storms move most of a
    desert's sand.

    Parameters
    ----------
    shear_velocity : float or array_like
        m/s; constraint: nonnegative.
    diameter : float or array_like
        m; constraint: positive.
    air_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    coefficient : float or array_like
        1, 1.5 for uniform, 1.8 for naturally graded and 2.8 for poorly sorted sand; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/(m s)] per unit width across the wind. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Bagnold (1941), referred to a 0.25 mm sand; well above the threshold, since it has no threshold term (see
    :func:`kawamura_sand_flux`).
    Equation: ``coefficient * np.sqrt(diameter / 2.5e-4) * air_density / gravity * shear_velocity**3``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import bagnold_sand_flux
    >>> round(float(bagnold_sand_flux(0.5, 0.00025, 1.2, 9.81, 1.8)), 13)
    0.0275229357798
    """
    shear_velocity = quantity("shear_velocity", shear_velocity, lower=0)
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    air_density = quantity("air_density", air_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    coefficient = quantity("coefficient", coefficient, lower=0, strict_lower=True)
    return coefficient * np.sqrt(diameter / 2.5e-4) * air_density / gravity * shear_velocity**3


@model
def kawamura_sand_flux(shear_velocity, threshold_shear_velocity, air_density, gravity, coefficient):
    """Kawamura's saltation mass flux with a threshold.

    Below the threshold the wind moves no sand at all, and just above it only a little, which Bagnold's cube law
    misses.

    Parameters
    ----------
    shear_velocity : float or array_like
        m/s; constraint: nonnegative.
    threshold_shear_velocity : float or array_like
        m/s, from :func:`threshold_shear_velocity`; constraint: nonnegative.
    air_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.
    coefficient : float or array_like
        1, about 2.78; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/(m s)], zero below the threshold. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Kawamura (1951): steady saltation over a flat dry sand bed; tends to the cube law far above the threshold.
    Equation: ``coefficient * air_density / gravity * np.maximum(u - ut, 0) * (u + ut)**2`` with
    ``u = shear_velocity`` and ``ut = threshold_shear_velocity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import kawamura_sand_flux
    >>> round(float(kawamura_sand_flux(0.5, 0.2, 1.2, 9.81, 2.78)), 13)
    0.0499889908257
    """
    shear_velocity = quantity("shear_velocity", shear_velocity, lower=0)
    threshold = quantity("threshold_shear_velocity", threshold_shear_velocity, lower=0)
    air_density = quantity("air_density", air_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    coefficient = quantity("coefficient", coefficient, lower=0, strict_lower=True)
    return coefficient * air_density / gravity * np.maximum(shear_velocity - threshold, 0) * (shear_velocity + threshold) ** 2


@model
def dune_migration_rate(sand_flux, dune_height, bulk_density):
    """Migration rate of a dune that keeps its shape.

    Sand blown over the crest settles on the lee face, so a dune advances by the sand flux divided by its height; small
    dunes outrun large ones.

    Parameters
    ----------
    sand_flux : float or array_like
        kg/(m s) across the crest; constraint: nonnegative.
    dune_height : float or array_like
        m; constraint: positive.
    bulk_density : float or array_like
        kg/m^3 of the deposited sand, about 1600; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Bagnold (1941): a two-dimensional dune of constant shape that traps all the sand crossing its crest.
    Equation: ``sand_flux / (bulk_density * dune_height)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import dune_migration_rate
    >>> float(dune_migration_rate(0.016, 10, 1600))
    1e-06
    """
    sand_flux = quantity("sand_flux", sand_flux, lower=0)
    dune_height = quantity("dune_height", dune_height, lower=0, strict_lower=True)
    bulk_density = quantity("bulk_density", bulk_density, lower=0, strict_lower=True)
    return sand_flux / (bulk_density * dune_height)
