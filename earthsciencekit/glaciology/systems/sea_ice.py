"""Sea ice: Stefan's growth law and freeboard.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "stefan_ice_thickness",
    "stefan_growth_rate",
    "sea_ice_freeboard",
]


@model
def stefan_ice_thickness(freezing_degree_time, conductivity, density, latent_heat):
    """Stefan's thickness of ice grown by conduction through itself.

    The latent heat released at the base of the ice must be conducted up through ice already formed, so the thickness
    grows only as the square root of the accumulated cold.

    Parameters
    ----------
    freezing_degree_time : float or array_like
        K s, time integral of the freezing temperature minus the surface temperature; constraint: nonnegative.
    conductivity : float or array_like
        W/(m K), about 2.2 for ice; constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    latent_heat : float or array_like
        J/kg; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Stefan (1891): ice free of snow, with its surface at the air temperature, a linear temperature profile, no ocean
    heat flux and no initial ice; it overestimates growth under snow.
    Equation: ``np.sqrt(2 * conductivity * freezing_degree_time / (density * latent_heat))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import stefan_ice_thickness
    >>> round(float(stefan_ice_thickness(260000000.0, 2.2, 917, 334000.0)), 5)
    1.93266
    """
    freezing_degree_time = quantity("freezing_degree_time", freezing_degree_time, lower=0)
    conductivity = quantity("conductivity", conductivity, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    latent_heat = quantity("latent_heat", latent_heat, lower=0, strict_lower=True)
    return np.sqrt(2 * conductivity * freezing_degree_time / (density * latent_heat))


@model
def stefan_growth_rate(thickness, temperature_difference, conductivity, density, latent_heat):
    """Basal growth rate of ice in Stefan's problem.

    Thin ice grows fast and thick ice slowly, because the conductive heat flux falls with thickness.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: positive.
    temperature_difference : float or array_like
        K, freezing temperature minus surface temperature; constraint: nonnegative.
    conductivity : float or array_like
        W/(m K); constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    latent_heat : float or array_like
        J/kg; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Stefan (1891), quasi-steady conduction with no ocean heat flux.
    Equation: ``conductivity * temperature_difference / (density * latent_heat * thickness)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import stefan_growth_rate
    >>> float(stefan_growth_rate(1, 20, 2, 1000, 200000.0))
    2e-07
    """
    thickness = quantity("thickness", thickness, lower=0, strict_lower=True)
    temperature_difference = quantity("temperature_difference", temperature_difference, lower=0)
    conductivity = quantity("conductivity", conductivity, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    latent_heat = quantity("latent_heat", latent_heat, lower=0, strict_lower=True)
    return conductivity * temperature_difference / (density * latent_heat * thickness)


@model
def sea_ice_freeboard(thickness, ice_density, water_density):
    """Height of floating snow-free ice above the waterline.

    Satellite altimeters measure the freeboard, and Archimedes' principle turns it into a thickness about ten times
    larger.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: nonnegative.
    ice_density : float or array_like
        kg/m^3; constraint: less than water_density.
    water_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hydrostatic equilibrium of a uniform floe without snow load.
    Equation: ``thickness * (1 - ice_density / water_density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.glaciology import sea_ice_freeboard
    >>> round(float(sea_ice_freeboard(2, 900, 1000)), 6)
    0.2
    """
    thickness = quantity("thickness", thickness, lower=0)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    require(ice_density < water_density, "ice must be less dense than water to float")
    return thickness * (1 - ice_density / water_density)
