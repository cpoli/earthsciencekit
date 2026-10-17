"""Figure of the Earth: rotational flattening, normal gravity, J2 from satellite orbits and the geoid.

SI units unless stated otherwise; latitude, inclination and spherical distance in radians.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "centrifugal_ratio",
    "homogeneous_flattening",
    "clairaut_gravity",
    "normal_gravity",
    "flattening_from_j2",
    "nodal_precession_rate",
    "stokes_function",
    "degree_geoid_height",
]

# GRS80 Somigliana constants: equatorial normal gravity, k = (b gamma_p) / (a gamma_e) - 1 and the first eccentricity squared.
_GAMMA_E = 9.7803267715
_SOMIGLIANA_K = 0.001931851353
_E2 = 0.00669438002290


@model
def centrifugal_ratio(angular_speed, radius, gravity):
    """Ratio of centrifugal to gravitational acceleration at the equator.

    Rotation flattens a planet in proportion to how strongly it competes with gravity.

    Parameters
    ----------
    angular_speed : float or array_like
        rad/s; constraint: nonnegative.
    radius : float or array_like
        m; constraint: positive. Equatorial radius.
    gravity : float or array_like
        m/s^2; constraint: positive. Equatorial gravity.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    The geodetic ``m`` uses ``GM / (a b)`` in place of gravity; the two differ at the 0.5% level.
    Equation: ``angular_speed**2 * radius / gravity``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    angular_speed = quantity("angular_speed", angular_speed, lower=0)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return angular_speed**2 * radius / gravity


@model
def homogeneous_flattening(centrifugal_ratio):
    """Flattening of a slowly rotating homogeneous fluid planet.

    A uniform self-gravitating fluid in rotation settles into an oblate spheroid.

    Parameters
    ----------
    centrifugal_ratio : float or array_like
        1; constraint: 0 <= value <= 0.1.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]: ``(a - c) / a``.

    Notes
    -----
    First order in rotation (Maclaurin spheroid). A centrally condensed Earth is less flattened than this.
    Equation: ``5 * centrifugal_ratio / 4``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    centrifugal_ratio = quantity("centrifugal_ratio", centrifugal_ratio, lower=0, upper=0.1)
    return 5 * centrifugal_ratio / 4


@model
def clairaut_gravity(latitude, equatorial_gravity, flattening, centrifugal_ratio):
    """Gravity on a rotating spheroid from Clairaut's theorem.

    The increase of gravity toward the poles measures the flattening without any surveying.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2. Geodetic latitude.
    equatorial_gravity : float or array_like
        m/s^2; constraint: positive.
    flattening : float or array_like
        1; constraint: 0 <= value < 1.
    centrifugal_ratio : float or array_like
        1; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    First order in flattening and rotation for a spheroid in hydrostatic equilibrium; the sin^4 term is omitted.
    Equation: ``equatorial_gravity * (1 + (5 * centrifugal_ratio / 2 - flattening) * np.sin(latitude) ** 2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    equatorial_gravity = quantity("equatorial_gravity", equatorial_gravity, lower=0, strict_lower=True)
    flattening = quantity("flattening", flattening, lower=0, upper=1, strict_upper=True)
    centrifugal_ratio = quantity("centrifugal_ratio", centrifugal_ratio, lower=0)
    return equatorial_gravity * (1 + (5 * centrifugal_ratio / 2 - flattening) * np.sin(latitude) ** 2)


@model
def normal_gravity(latitude):
    """Normal gravity on the GRS80 ellipsoid.

    The reference field that gravity anomalies are measured against.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2. Geodetic latitude.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Somigliana's closed form on the ellipsoid surface with the GRS80 constants
    ``gamma_e = 9.7803267715``, ``k = 0.001931851353`` and ``e2 = 0.00669438002290``; exact for the level ellipsoid.
    Equation: ``gamma_e * (1 + k * np.sin(latitude) ** 2) / np.sqrt(1 - e2 * np.sin(latitude) ** 2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    sine_squared = np.sin(latitude) ** 2
    return _GAMMA_E * (1 + _SOMIGLIANA_K * sine_squared) / np.sqrt(1 - _E2 * sine_squared)


@model
def flattening_from_j2(j2, centrifugal_ratio):
    """Flattening of a hydrostatic Earth from its J2.

    The oblateness term of the gravity field and the rotation fix the shape of the level surface.

    Parameters
    ----------
    j2 : float or array_like
        1; constraint: real. Unnormalized second zonal harmonic ``-C20``.
    centrifugal_ratio : float or array_like
        1; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    First order in flattening (Clairaut's relation for the external field); second-order terms change ``1/f`` by
    about 0.1%.
    Equation: ``3 * j2 / 2 + centrifugal_ratio / 2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    j2 = quantity("j2", j2)
    centrifugal_ratio = quantity("centrifugal_ratio", centrifugal_ratio, lower=0)
    return 3 * j2 / 2 + centrifugal_ratio / 2


@model
def nodal_precession_rate(semi_major_axis, eccentricity, inclination, j2, radius, gravitational_parameter):
    """Secular regression of a satellite orbit's node.

    The equatorial bulge torques an inclined orbit, so the node drifts at a rate that measures J2.

    Parameters
    ----------
    semi_major_axis : float or array_like
        m; constraint: positive.
    eccentricity : float or array_like
        1; constraint: 0 <= value < 1.
    inclination : float or array_like
        rad; constraint: 0 <= inclination <= pi.
    j2 : float or array_like
        1; constraint: real.
    radius : float or array_like
        m; constraint: positive. Equatorial reference radius of J2.
    gravitational_parameter : float or array_like
        m^3/s^2; constraint: positive. ``G M`` of the Earth.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]: negative (westward) for prograde orbits and positive J2.

    Notes
    -----
    First-order secular J2 perturbation, without drag or higher harmonics; perigee must clear the reference radius.
    Equation: ``-1.5 * n * j2 * (radius / semi_major_axis) ** 2 * np.cos(inclination) / (1 - eccentricity**2) ** 2``
    with mean motion ``n = np.sqrt(gravitational_parameter / semi_major_axis**3)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    semi_major_axis = quantity("semi_major_axis", semi_major_axis, lower=0, strict_lower=True)
    eccentricity = quantity("eccentricity", eccentricity, lower=0, upper=1, strict_upper=True)
    inclination = quantity("inclination", inclination, lower=0, upper=np.pi)
    j2 = quantity("j2", j2)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    gravitational_parameter = quantity("gravitational_parameter", gravitational_parameter, lower=0, strict_lower=True)
    require(semi_major_axis * (1 - eccentricity) > radius, "perigee must lie above the reference radius")
    mean_motion = np.sqrt(gravitational_parameter / semi_major_axis**3)
    return -1.5 * mean_motion * j2 * (radius / semi_major_axis) ** 2 * np.cos(inclination) / (1 - eccentricity**2) ** 2


@model
def stokes_function(spherical_distance):
    """Stokes's kernel for the geoid.

    Weighting gravity anomalies over the whole globe by spherical distance gives the geoid height.

    Parameters
    ----------
    spherical_distance : float or array_like
        rad; constraint: 0 < value <= pi.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Spherical approximation with no masses outside the geoid; singular as ``2 / psi`` at zero distance. Its Legendre
    expansion is ``sum((2 n + 1) / (n - 1) P_n(cos psi), n >= 2)``, so the geoid is
    ``radius / (4 pi gravity) * integral(anomaly * S dsigma)``.
    Equation: ``1 / s - 6 * s + 1 - 5 * np.cos(psi) - 3 * np.cos(psi) * np.log(s + s**2)`` with ``s = np.sin(psi / 2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    psi = quantity("spherical_distance", spherical_distance, lower=0, upper=np.pi, strict_lower=True)
    s = np.sin(psi / 2)
    return 1 / s - 6 * s + 1 - 5 * np.cos(psi) - 3 * np.cos(psi) * np.log(s + s**2)


@model
def degree_geoid_height(anomaly_coefficient, degree, radius, gravity):
    """Geoid height of one spherical-harmonic degree of the gravity anomaly.

    Stokes's integral divides each degree of the anomaly by ``n - 1``, so long wavelengths dominate the geoid.

    Parameters
    ----------
    anomaly_coefficient : float or array_like
        m/s^2; constraint: real. Coefficient of the degree-``n`` surface harmonic of the gravity anomaly.
    degree : float or array_like
        1; constraint: integer >= 2.
    radius : float or array_like
        m; constraint: positive.
    gravity : float or array_like
        m/s^2; constraint: positive. Mean normal gravity.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m]: coefficient of the same surface harmonic in the geoid height.

    Notes
    -----
    Spherical approximation of Stokes's formula; degrees 0 and 1 are excluded by construction.
    Equation: ``radius * anomaly_coefficient / ((degree - 1) * gravity)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    anomaly_coefficient = quantity("anomaly_coefficient", anomaly_coefficient)
    degree = quantity("degree", degree, lower=2)
    require(degree == np.round(degree), "degree must be an integer")
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return radius * anomaly_coefficient / ((degree - 1) * gravity)
