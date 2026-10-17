"""Ice sheets: driving stress, Nye's plastic and Vialov's shallow-ice profiles, Weertman sliding and the marine grounding line.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "driving_stress",
    "nye_plastic_profile",
    "weertman_sliding_velocity",
    "flotation_thickness",
    "ice_shelf_strain_rate",
    "shallow_ice_flux",
    "shallow_ice_diffusivity",
    "vialov_dome_height",
    "vialov_profile",
    "schoof_grounding_line_flux",
]


@model
def driving_stress(thickness, surface_slope, density, gravity):
    """Gravitational driving stress of a glacier or ice sheet.

    The weight of the ice resolved down the surface slope is what the bed and the ice must resist.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: nonnegative.
    surface_slope : float or array_like
        1, tangent of the surface slope; constraint: nonnegative.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa], typically 50-150 kPa under glaciers. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Shallow-ice balance: horizontal scales much larger than the thickness, small slopes, longitudinal stresses neglected.
    Equation: ``density * gravity * thickness * surface_slope``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return density * gravity * thickness * surface_slope


@model
def nye_plastic_profile(distance, yield_stress, density, gravity):
    """Surface height of a perfectly plastic ice sheet on a flat bed.

    If ice yields at a fixed stress, the ice sheet takes a parabolic profile whose height depends only on the distance
    from its margin.

    Parameters
    ----------
    distance : float or array_like
        m, inland from the margin; constraint: nonnegative.
    yield_stress : float or array_like
        Pa, about 100 kPa; constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nye (1951): a horizontal bed, basal shear stress everywhere equal to the yield stress and no isostatic depression.
    Equation: ``np.sqrt(2 * yield_stress * distance / (density * gravity))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    distance = quantity("distance", distance, lower=0)
    yield_stress = quantity("yield_stress", yield_stress, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return np.sqrt(2 * yield_stress * distance / (density * gravity))


@model
def weertman_sliding_velocity(basal_stress, sliding_coefficient, exponent):
    """Weertman's sliding velocity of temperate ice over a rough bed.

    Ice gets past bed obstacles by regelation around small ones and enhanced creep around large ones, which makes the
    sliding speed a power of the basal stress.

    Parameters
    ----------
    basal_stress : float or array_like
        Pa; constraint: nonnegative.
    sliding_coefficient : float or array_like
        m/(s Pa^m); constraint: positive.
    exponent : float or array_like
        1, m = (n + 1) / 2, about 2 for Glen's n = 3; constraint: at least 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Weertman (1957): a bed at the melting point, free of cavities and sediment; the coefficient lumps the bed roughness.
    Equation: ``sliding_coefficient * basal_stress**exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    basal_stress = quantity("basal_stress", basal_stress, lower=0)
    sliding_coefficient = quantity("sliding_coefficient", sliding_coefficient, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    return sliding_coefficient * basal_stress**exponent


@model
def flotation_thickness(bed_depth, ice_density, water_density):
    """Thickness at which ice resting on a submarine bed begins to float.

    The grounding line of a marine ice sheet sits where the ice is just thin enough to float, so its thickness is set by
    the water depth alone.

    Parameters
    ----------
    bed_depth : float or array_like
        m, below sea level; constraint: nonnegative.
    ice_density : float or array_like
        kg/m^3; constraint: positive.
    water_density : float or array_like
        kg/m^3; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Archimedes' principle for an ice column of uniform density.
    Equation: ``water_density * bed_depth / ice_density``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    bed_depth = quantity("bed_depth", bed_depth, lower=0)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    return water_density * bed_depth / ice_density


@model
def ice_shelf_strain_rate(thickness, rate_factor, exponent, ice_density, water_density, gravity):
    """Spreading rate of a freely floating ice shelf.

    A floating shelf spreads under the difference between its weight and the water pressure, faster the thicker it is.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    ice_density : float or array_like
        kg/m^3; constraint: less than water_density.
    water_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1/s], along-flow stretching. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Weertman (1957) for an unconfined shelf spreading in one direction, without buttressing or side drag.
    Equation: ``rate_factor * (ice_density * gravity * thickness * (1 - ice_density / water_density) / 4) ** exponent``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(ice_density < water_density, "ice must be less dense than water to float")
    return rate_factor * (ice_density * gravity * thickness * (1 - ice_density / water_density) / 4) ** exponent


@model
def shallow_ice_flux(thickness, surface_slope, rate_factor, exponent, density, gravity):
    """Ice flux per unit width in the shallow-ice approximation.

    In a thin ice sheet the flux at each point depends only on the local thickness and surface slope.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: nonnegative.
    surface_slope : float or array_like
        1, magnitude of the surface gradient; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s], down the surface gradient. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Hutter (1983) and Mahaffy (1976): Glen ice frozen to its bed, aspect ratio small, uniform rate factor.
    Equation: ``2 * rate_factor / (exponent + 2) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return 2 * rate_factor / (exponent + 2) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 2)


@model
def shallow_ice_diffusivity(thickness, surface_slope, rate_factor, exponent, density, gravity):
    """Nonlinear diffusivity of the shallow-ice equation.

    The shallow-ice approximation makes the ice thickness obey a diffusion equation whose diffusivity grows steeply with
    thickness and slope.

    Parameters
    ----------
    thickness : float or array_like
        m; constraint: nonnegative.
    surface_slope : float or array_like
        1, magnitude of the surface gradient; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s], so that the flux is minus the diffusivity times the surface gradient. Scalar input returns a
        NumPy scalar; arrays broadcast.

    Notes
    -----
    Mahaffy (1976): the shallow-ice approximation for Glen ice frozen to its bed. Below, A = rate_factor, n = exponent.
    Equation: ``2 * A / (n + 2) * (density * gravity) ** n * thickness ** (n + 2) * surface_slope ** (n - 1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    thickness = quantity("thickness", thickness, lower=0)
    surface_slope = quantity("surface_slope", surface_slope, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    stiffness = 2 * rate_factor / (exponent + 2) * (density * gravity) ** exponent
    return stiffness * thickness ** (exponent + 2) * surface_slope ** (exponent - 1)


@model
def vialov_dome_height(half_width, accumulation, rate_factor, exponent, density, gravity):
    """Central height of a steady shallow-ice sheet.

    An ice sheet in balance with a uniform snowfall has a dome height that grows only as the square root of its width.

    Parameters
    ----------
    half_width : float or array_like
        m, divide to margin; constraint: positive.
    accumulation : float or array_like
        m/s of ice; constraint: positive.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Vialov (1958): a flat bed, uniform accumulation up to a fixed margin, shallow-ice flow without sliding. Below,
    A = rate_factor, n = exponent.
    Equation: ``(2 * ((n + 2) * accumulation / (2 * A * (density * gravity) ** n)) ** (1 / n)) ** (n / (2 * n + 2)) * np.sqrt(half_width)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    half_width = quantity("half_width", half_width, lower=0, strict_lower=True)
    accumulation = quantity("accumulation", accumulation, lower=0, strict_lower=True)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    scale = ((exponent + 2) * accumulation / (2 * rate_factor * (density * gravity) ** exponent)) ** (1 / exponent)
    return (2 * scale) ** (exponent / (2 * exponent + 2)) * np.sqrt(half_width)


@model
def vialov_profile(distance, half_width, dome_height, exponent):
    """Surface profile of a steady shallow-ice sheet.

    Glen's law flattens the dome and steepens the margin of a steady ice sheet compared with Nye's parabola.

    Parameters
    ----------
    distance : float or array_like
        m, from the divide; constraint: between 0 and half_width.
    half_width : float or array_like
        m; constraint: positive.
    dome_height : float or array_like
        m; constraint: nonnegative.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Vialov (1958), see vialov_dome_height; the shallow-ice approximation fails within a few thicknesses of the divide
    and the margin.
    Equation: ``dome_height * (1 - (distance / half_width) ** ((exponent + 1) / exponent)) ** (exponent / (2 * exponent + 2))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    distance = quantity("distance", distance, lower=0)
    half_width = quantity("half_width", half_width, lower=0, strict_lower=True)
    dome_height = quantity("dome_height", dome_height, lower=0)
    exponent = quantity("exponent", exponent, lower=1)
    require(distance <= half_width, "distance exceeds the half-width")
    return dome_height * (1 - (distance / half_width) ** ((exponent + 1) / exponent)) ** (exponent / (2 * exponent + 2))


@model
def schoof_grounding_line_flux(
    grounding_thickness, rate_factor, exponent, friction_coefficient, friction_exponent, ice_density, water_density, gravity
):
    """Schoof's ice flux across the grounding line of a marine ice sheet.

    The flux out of a marine ice sheet is set by the ice thickness at its grounding line, and rises so steeply with it
    that a grounding line on a bed deepening inland cannot be steady.

    Parameters
    ----------
    grounding_thickness : float or array_like
        m, the flotation thickness at the grounding line; constraint: nonnegative.
    rate_factor : float or array_like
        Pa^-n/s; constraint: positive.
    exponent : float or array_like
        1, Glen's n; constraint: at least 1.
    friction_coefficient : float or array_like
        Pa (s/m)^m, the C of a basal stress C u^m; constraint: positive.
    friction_exponent : float or array_like
        1, the m of a basal stress C u^m, such as 1/3; constraint: positive.
    ice_density : float or array_like
        kg/m^3; constraint: less than water_density.
    water_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s], per unit width. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Schoof (2007): a flowline with a power-law sliding bed, an unbuttressed shelf and a boundary layer at the grounding
    line; steady states are stable only where the bed rises inland. Below, A = rate_factor, n = exponent, C =
    friction_coefficient, m = friction_exponent, rho_i = ice_density, rho_w = water_density, h = grounding_thickness.
    Equation: ``(A * (rho_i * gravity)**(n + 1) * (1 - rho_i / rho_w)**n / (4**n * C))**(1 / (m + 1)) * h**((m + n + 3) / (m + 1))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    grounding_thickness = quantity("grounding_thickness", grounding_thickness, lower=0)
    rate_factor = quantity("rate_factor", rate_factor, lower=0, strict_lower=True)
    exponent = quantity("exponent", exponent, lower=1)
    friction_coefficient = quantity("friction_coefficient", friction_coefficient, lower=0, strict_lower=True)
    friction_exponent = quantity("friction_exponent", friction_exponent, lower=0, strict_lower=True)
    ice_density = quantity("ice_density", ice_density, lower=0, strict_lower=True)
    water_density = quantity("water_density", water_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(ice_density < water_density, "ice must be less dense than water to float")
    m = friction_exponent
    factor = rate_factor * (ice_density * gravity) ** (exponent + 1) * (1 - ice_density / water_density) ** exponent
    return (factor / (4**exponent * friction_coefficient)) ** (1 / (m + 1)) * grounding_thickness ** ((m + exponent + 3) / (m + 1))
