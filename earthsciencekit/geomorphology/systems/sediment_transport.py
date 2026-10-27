"""Sediment transport: settling, bed shear, the threshold of motion, suspension and bedload.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import VON_KARMAN

__all__ = [
    "stokes_settling_velocity",
    "ferguson_church_settling_velocity",
    "bed_shear_stress",
    "shear_velocity",
    "shields_number",
    "dimensionless_grain_size",
    "critical_shields_number",
    "keulegan_mean_velocity",
    "rouse_number",
    "rouse_concentration",
    "meyer_peter_muller_transport",
    "bedload_flux",
]


@model
def stokes_settling_velocity(diameter, sediment_density, fluid_density, viscosity, gravity):
    """Stokes settling velocity of a small sphere.

    Fine grains fall so slowly that viscosity alone resists them, which is why silt and clay stay in suspension.

    Parameters
    ----------
    diameter : float or array_like
        m; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than fluid_density.
    fluid_density : float or array_like
        kg/m^3; constraint: positive.
    viscosity : float or array_like
        Pa s, dynamic viscosity of the fluid; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], positive downward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Stokes (1851): a smooth sphere at particle Reynolds number below about 1, roughly quartz finer than 0.1 mm in water.
    Equation: ``(sediment_density - fluid_density) * gravity * diameter**2 / (18 * viscosity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import stokes_settling_velocity
    >>> round(float(stokes_settling_velocity(0.0001, 2650, 1000, 0.001, 9.81)), 14)
    0.0089925
    """
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    fluid_density = quantity("fluid_density", fluid_density, lower=0, strict_lower=True)
    viscosity = quantity("viscosity", viscosity, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(sediment_density > fluid_density, "sediment must be denser than the fluid")
    return (sediment_density - fluid_density) * gravity * diameter**2 / (18 * viscosity)


@model
def ferguson_church_settling_velocity(diameter, sediment_density, fluid_density, kinematic_viscosity, gravity):
    """Settling velocity of natural sand grains at any size.

    One formula spans viscous settling of silt and drag-dominated settling of gravel, joining Stokes's law to the
    turbulent limit.

    Parameters
    ----------
    diameter : float or array_like
        m, sieve diameter; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than fluid_density.
    fluid_density : float or array_like
        kg/m^3; constraint: positive.
    kinematic_viscosity : float or array_like
        m^2/s; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s], positive downward. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Ferguson and Church (2004) with C1 = 18 and C2 = 1 for natural grains; Stokes's law for fine grains and a constant
    drag coefficient for coarse ones; spheres need C2 = 0.4.
    Equation: ``R * gravity * diameter**2 / (18 * kinematic_viscosity + np.sqrt(0.75 * R * gravity * diameter**3))`` with
    ``R = sediment_density / fluid_density - 1``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import ferguson_church_settling_velocity
    >>> round(float(ferguson_church_settling_velocity(0.0001, 2650, 1000, 1e-06, 9.81)), 14)
    0.00753413002789
    """
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    fluid_density = quantity("fluid_density", fluid_density, lower=0, strict_lower=True)
    kinematic_viscosity = quantity("kinematic_viscosity", kinematic_viscosity, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(sediment_density > fluid_density, "sediment must be denser than the fluid")
    submerged = sediment_density / fluid_density - 1
    return submerged * gravity * diameter**2 / (18 * kinematic_viscosity + np.sqrt(0.75 * submerged * gravity * diameter**3))


@model
def bed_shear_stress(hydraulic_radius, slope, density, gravity):
    """Mean bed shear stress of steady uniform flow.

    The weight of the water resolved down the slope is balanced by friction on the bed, and that friction moves the
    sediment.

    Parameters
    ----------
    hydraulic_radius : float or array_like
        m, the depth in a wide channel; constraint: nonnegative.
    slope : float or array_like
        1, tangent of the water-surface slope; constraint: nonnegative.
    density : float or array_like
        kg/m^3 of the fluid; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The depth-slope product of steady uniform flow; bank friction and bedforms take part of it in real channels.
    Equation: ``density * gravity * hydraulic_radius * slope``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import bed_shear_stress
    >>> float(bed_shear_stress(2, 0.001, 1000, 9.81))
    19.62
    """
    hydraulic_radius = quantity("hydraulic_radius", hydraulic_radius, lower=0)
    slope = quantity("slope", slope, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return density * gravity * hydraulic_radius * slope


@model
def shear_velocity(shear_stress, density):
    """Shear velocity of a flow over its bed.

    The bed shear stress expressed as a velocity sets the scale of near-bed turbulence.

    Parameters
    ----------
    shear_stress : float or array_like
        Pa; constraint: nonnegative.
    density : float or array_like
        kg/m^3 of the fluid; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    A definition; it applies to any turbulent boundary layer, in water or air.
    Equation: ``np.sqrt(shear_stress / density)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import shear_velocity
    >>> round(float(shear_velocity(19.62, 1000)), 12)
    0.140071410359
    """
    shear_stress = quantity("shear_stress", shear_stress, lower=0)
    density = quantity("density", density, lower=0, strict_lower=True)
    return np.sqrt(shear_stress / density)


@model
def shields_number(shear_stress, diameter, sediment_density, fluid_density, gravity):
    """Shields number, the dimensionless bed shear stress.

    The ratio of the force that drives a grain to the submerged weight that holds it decides whether the bed moves.

    Parameters
    ----------
    shear_stress : float or array_like
        Pa; constraint: nonnegative.
    diameter : float or array_like
        m; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than fluid_density.
    fluid_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Shields (1936): uniform cohesionless grains on a nearly flat bed.
    Equation: ``shear_stress / ((sediment_density - fluid_density) * gravity * diameter)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import shields_number
    >>> round(float(shields_number(19.62, 0.01, 2650, 1000, 9.81)), 12)
    0.121212121212
    """
    shear_stress = quantity("shear_stress", shear_stress, lower=0)
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    fluid_density = quantity("fluid_density", fluid_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(sediment_density > fluid_density, "sediment must be denser than the fluid")
    return shear_stress / ((sediment_density - fluid_density) * gravity * diameter)


@model
def dimensionless_grain_size(diameter, sediment_density, fluid_density, kinematic_viscosity, gravity):
    """Dimensionless grain size of a sediment in a fluid.

    One number combines grain size, density and viscosity, so that a single curve gives the threshold of motion for
    sand in water and in air.

    Parameters
    ----------
    diameter : float or array_like
        m; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than fluid_density.
    fluid_density : float or array_like
        kg/m^3; constraint: positive.
    kinematic_viscosity : float or array_like
        m^2/s; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Soulsby and Whitehouse (1997); about 25 per millimetre of quartz in water at 20 degC.
    Equation: ``diameter * ((sediment_density / fluid_density - 1) * gravity / kinematic_viscosity**2)**(1 / 3)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import dimensionless_grain_size
    >>> round(float(dimensionless_grain_size(0.001, 2650, 1000, 1e-06, 9.81)), 10)
    25.2959494002
    """
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    fluid_density = quantity("fluid_density", fluid_density, lower=0, strict_lower=True)
    kinematic_viscosity = quantity("kinematic_viscosity", kinematic_viscosity, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(sediment_density > fluid_density, "sediment must be denser than the fluid")
    return diameter * ((sediment_density / fluid_density - 1) * gravity / kinematic_viscosity**2) ** (1 / 3)


@model
def critical_shields_number(dimensionless_grain_size):
    """Critical Shields number for the threshold of motion.

    The Shields curve falls from fine to medium sand and levels off for gravel, where a grain moves once the bed
    stress reaches about 5% of its submerged weight per unit area.

    Parameters
    ----------
    dimensionless_grain_size : float or array_like
        1, from :func:`dimensionless_grain_size`; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Soulsby and Whitehouse (1997) fit to Shields's (1936) data, with the rise for fine grains; uniform cohesionless
    grains on a flat bed, with a scatter of about a factor of two.
    Equation: ``0.30 / (1 + 1.2 * D) + 0.055 * (1 - np.exp(-0.020 * D))`` with ``D = dimensionless_grain_size``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import critical_shields_number
    >>> round(float(critical_shields_number(10)), 13)
    0.0330467316576
    """
    size = quantity("dimensionless_grain_size", dimensionless_grain_size, lower=0, strict_lower=True)
    return 0.30 / (1 + 1.2 * size) + 0.055 * (1 - np.exp(-0.020 * size))


@model
def keulegan_mean_velocity(shear_velocity, depth, roughness_height):
    """Depth-averaged velocity of a rough turbulent channel flow.

    Integrating the logarithmic velocity profile over the depth turns a bed shear stress into the mean velocity a
    current meter would record.

    Parameters
    ----------
    shear_velocity : float or array_like
        m/s; constraint: nonnegative.
    depth : float or array_like
        m; constraint: greater than roughness_height / 11.
    roughness_height : float or array_like
        m, Nikuradse's equivalent sand roughness, about 2.5 grain diameters on a flat bed; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Keulegan (1938): fully rough flow in a wide channel with the von Karman constant 0.41; the depth should exceed
    the roughness several times.
    Equation: ``shear_velocity / VON_KARMAN * np.log(11 * depth / roughness_height)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import keulegan_mean_velocity
    >>> round(float(keulegan_mean_velocity(0.1, 2, 0.022)), 11)
    1.68481836073
    """
    shear_velocity = quantity("shear_velocity", shear_velocity, lower=0)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    roughness_height = quantity("roughness_height", roughness_height, lower=0, strict_lower=True)
    require(11 * depth > roughness_height, "depth must exceed roughness_height / 11")
    return shear_velocity / VON_KARMAN * np.log(11 * depth / roughness_height)


@model
def rouse_number(settling_velocity, shear_velocity):
    """Rouse number of suspended sediment.

    Settling pulls grains down while turbulence mixes them up; the ratio of the two sets how evenly they fill the
    flow.

    Parameters
    ----------
    settling_velocity : float or array_like
        m/s; constraint: nonnegative.
    shear_velocity : float or array_like
        m/s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]: below about 0.8 the load is nearly uniform, above about 2.5 it is bedload. Scalar input
        returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rouse (1937), with the eddy diffusivity of sediment equal to that of momentum and the von Karman constant 0.41.
    Equation: ``settling_velocity / (VON_KARMAN * shear_velocity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import rouse_number
    >>> float(rouse_number(0.041, 0.1))
    1.0
    """
    settling_velocity = quantity("settling_velocity", settling_velocity, lower=0)
    shear_velocity = quantity("shear_velocity", shear_velocity, lower=0, strict_lower=True)
    return settling_velocity / (VON_KARMAN * shear_velocity)


@model
def rouse_concentration(height, depth, reference_height, reference_concentration, rouse_number):
    """Rouse profile of suspended-sediment concentration.

    The balance of settling and turbulent mixing concentrates sediment near the bed, the more so the coarser the
    grains.

    Parameters
    ----------
    height : float or array_like
        m above the bed; constraint: between reference_height and depth.
    depth : float or array_like
        m; constraint: greater than reference_height.
    reference_height : float or array_like
        m, where the concentration is known; constraint: positive.
    reference_concentration : float or array_like
        Any concentration unit, at the reference height; constraint: nonnegative.
    rouse_number : float or array_like
        1, from :func:`rouse_number`; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value in the unit of reference_concentration. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Rouse (1937): steady uniform flow with a parabolic eddy diffusivity and dilute suspension; zero at the surface.
    Equation: ``reference_concentration * ((depth - height) / height * reference_height / (depth - reference_height))**rouse_number``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import rouse_concentration
    >>> round(float(rouse_concentration(0.5, 1, 0.05, 2, 1)), 12)
    0.105263157895
    """
    height = quantity("height", height, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    reference_height = quantity("reference_height", reference_height, lower=0, strict_lower=True)
    reference_concentration = quantity("reference_concentration", reference_concentration, lower=0)
    rouse_number = quantity("rouse_number", rouse_number, lower=0)
    require(reference_height < depth, "reference_height must lie below the surface")
    require((height >= reference_height) & (height <= depth), "height must lie between reference_height and depth")
    ratio = (depth - height) / height * reference_height / (depth - reference_height)
    return reference_concentration * ratio**rouse_number


@model
def meyer_peter_muller_transport(shields_number, critical_shields_number):
    """Meyer-Peter and Müller dimensionless bedload transport rate.

    Bedload grows as the 3/2 power of the excess of the bed stress over its threshold, so a modest flood moves far
    more gravel than a long ordinary flow.

    Parameters
    ----------
    shields_number : float or array_like
        1, from :func:`shields_number`; constraint: nonnegative.
    critical_shields_number : float or array_like
        1, 0.047 in the original fit; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1], the Einstein number, zero below the threshold. Scalar input returns a NumPy scalar; arrays
        broadcast.

    Notes
    -----
    Meyer-Peter and Müller (1948): flume data for well-sorted sand and gravel on a plane bed, with the skin friction
    part of the stress; Wong and Parker (2006) correct it to ``3.97 (theta - 0.0495)**1.5``.
    Equation: ``8 * np.maximum(shields_number - critical_shields_number, 0)**1.5``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import meyer_peter_muller_transport
    >>> round(float(meyer_peter_muller_transport(0.147, 0.047)), 12)
    0.252982212813
    """
    shields_number = quantity("shields_number", shields_number, lower=0)
    critical_shields_number = quantity("critical_shields_number", critical_shields_number, lower=0)
    return 8 * np.maximum(shields_number - critical_shields_number, 0) ** 1.5


@model
def bedload_flux(transport_rate, diameter, sediment_density, fluid_density, gravity):
    """Volumetric bedload flux per unit width from the Einstein number.

    The dimensionless transport rate scales back to a volume of grains crossing each metre of channel width per
    second.

    Parameters
    ----------
    transport_rate : float or array_like
        1, the Einstein number of :func:`meyer_peter_muller_transport`; constraint: nonnegative.
    diameter : float or array_like
        m; constraint: positive.
    sediment_density : float or array_like
        kg/m^3; constraint: greater than fluid_density.
    fluid_density : float or array_like
        kg/m^3; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s] of solid grains, without pore space. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Einstein (1950) scaling; uniform grains.
    Equation: ``transport_rate * np.sqrt((sediment_density / fluid_density - 1) * gravity * diameter) * diameter``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geomorphology import bedload_flux
    >>> round(float(bedload_flux(1, 0.01, 2650, 1000, 9.81)), 14)
    0.00402324495899
    """
    transport_rate = quantity("transport_rate", transport_rate, lower=0)
    diameter = quantity("diameter", diameter, lower=0, strict_lower=True)
    sediment_density = quantity("sediment_density", sediment_density, lower=0, strict_lower=True)
    fluid_density = quantity("fluid_density", fluid_density, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    require(sediment_density > fluid_density, "sediment must be denser than the fluid")
    return transport_rate * np.sqrt((sediment_density / fluid_density - 1) * gravity * diameter) * diameter
