"""Sea level: thermal expansion and the mass equivalent of added water.

SI units unless stated otherwise.
"""

from earthsciencekit._validation import model, quantity

__all__ = [
    "thermal_expansion_sea_level",
    "water_equivalent_sea_level",
]


@model
def thermal_expansion_sea_level(expansion_coefficient, depth, temperature_change):
    """Steric sea-level slab approximation.

    Ocean expansion connects heat uptake to a sea-level contribution.

    Parameters
    ----------
    expansion_coefficient : float or array_like
        1/K; constraint: positive.
    depth : float or array_like
        m; constraint: positive.
    temperature_change : float or array_like
        K; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant coefficient and uniform warming; salinity and pressure dependence are omitted.
    Equation: ``expansion_coefficient * depth * temperature_change``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    expansion_coefficient = quantity("expansion_coefficient", expansion_coefficient, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    temperature_change = quantity("temperature_change", temperature_change)
    return expansion_coefficient * depth * temperature_change


@model
def water_equivalent_sea_level(water_mass, ocean_area, water_density):
    """Mass-equivalent sea-level contribution.

    Added ocean water changes mean level through mass conservation.

    Parameters
    ----------
    water_mass : float or array_like
        kg; constraint: nonnegative.
    ocean_area : float or array_like
        m^2; constraint: positive.
    water_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [W/m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uniform ocean area and density; gravitational fingerprints and land motion are omitted.
    Equation: ``water_mass / (ocean_area * water_density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    water_mass = quantity("water_mass", water_mass, lower=0)
    ocean_area = quantity("ocean_area", ocean_area, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    return water_mass / (ocean_area * water_density)
