"""Geomagnetism and paleomagnetism: the dipole, Gauss coefficients, field directions, poles and crustal anomalies.

SI units unless stated otherwise; latitude, longitude, declination and inclination in radians.
"""

import numpy as np
from scipy.special import gammaln, lpmv

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import MU0

__all__ = [
    "dipole_field",
    "dipole_inclination",
    "paleolatitude",
    "zonal_field",
    "gauss_dipole_moment",
    "field_elements",
    "direction_angle",
    "virtual_geomagnetic_pole",
    "block_magnetic_anomaly",
    "schmidt_legendre",
    "geomagnetic_field",
]


def _latitude(name, value):
    return quantity(name, value, lower=-np.pi / 2, upper=np.pi / 2)


@model
def dipole_field(latitude, distance, moment):
    """Field of a geocentric axial dipole.

    A uniformly magnetized Earth gives a field that dips more steeply toward the poles.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2. Geocentric latitude.
    distance : float or array_like
        m; constraint: positive. Distance from the Earth's centre.
    moment : float or array_like
        A m^2; constraint: real. Positive for today's polarity, whose field points down in the northern hemisphere.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [T]: the northward and the downward component.

    Notes
    -----
    Pure axial dipole, outside its sources. The field at the poles is twice that at the equator.
    Equation: ``MU0 * moment / (4 * np.pi * distance**3) * [np.cos(latitude), 2 * np.sin(latitude)]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import dipole_field
    >>> np.round(dipole_field(0, 1, 10000000), 5).tolist()
    [1.0, 0.0]
    """
    latitude = _latitude("latitude", latitude)
    distance = quantity("distance", distance, lower=0, strict_lower=True)
    moment = quantity("moment", moment)
    scale = MU0 * moment / (4 * np.pi * distance**3)
    return np.stack(np.broadcast_arrays(scale * np.cos(latitude), 2 * scale * np.sin(latitude)), axis=-1)


@model
def dipole_inclination(latitude):
    """Inclination of an axial dipole field.

    The dip of a compass needle tells a navigator, or a rock, its magnetic latitude.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad]: positive downward, for today's polarity.

    Notes
    -----
    Geocentric axial dipole; ``tan(I) = 2 tan(latitude)``.
    Equation: ``np.arctan2(2 * np.sin(latitude), np.cos(latitude))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import dipole_inclination
    >>> round(float(dipole_inclination(0.7853981633974483)), 5)
    1.10715
    """
    latitude = _latitude("latitude", latitude)
    return np.arctan2(2 * np.sin(latitude), np.cos(latitude))


@model
def paleolatitude(inclination):
    """Paleolatitude from a remanent inclination.

    A rock's frozen-in dip records how far from the pole it formed.

    Parameters
    ----------
    inclination : float or array_like
        rad; constraint: -pi/2 <= inclination <= pi/2. Positive downward.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad]: the magnetic latitude of the site when the rock was magnetized.

    Notes
    -----
    Assumes a geocentric axial dipole averaged over secular variation; a reversed direction must be inverted first.
    Equation: ``np.arctan2(np.sin(inclination), 2 * np.cos(inclination))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import paleolatitude
    >>> round(float(paleolatitude(1.1071487177940906)), 6)
    0.785398
    """
    inclination = _latitude("inclination", inclination)
    return np.arctan2(np.sin(inclination), 2 * np.cos(inclination))


@model
def zonal_field(latitude, distance, degree, coefficient, reference_radius):
    """Field of one zonal Gauss coefficient.

    Expanding the potential in spherical harmonics separates the dipole from smaller, faster-decaying parts.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2. Geocentric latitude.
    distance : float or array_like
        m; constraint: positive. Distance from the Earth's centre.
    degree : float or array_like
        1; constraint: integer >= 1.
    coefficient : float or array_like
        T; constraint: real. Gauss coefficient ``g_n^0``; for zonal terms Schmidt and unnormalized Legendre functions agree.
    reference_radius : float or array_like
        m; constraint: positive. Usually 6371.2 km.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [T]: the northward and the downward component.

    Notes
    -----
    Internal sources only, potential ``V = reference_radius * sum((a / r) ** (n + 1) g_n^0 P_n(sin latitude))``.
    Sum several calls for a truncated zonal model; non-zonal terms (``m > 0``) are not included.
    Equation: ``(a / r) ** (n + 2) * coefficient * [lpmv(1, n, x), -(n + 1) * lpmv(0, n, x)]`` with
    ``x = np.sin(latitude)``, ``a = reference_radius``, ``r = distance`` and ``n = degree``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import zonal_field
    >>> np.round(zonal_field(0.5235987755982988, 1, 2, 1, 1), 5).tolist()
    [-1.29904, 0.375]
    """
    latitude = _latitude("latitude", latitude)
    distance = quantity("distance", distance, lower=0, strict_lower=True)
    degree = quantity("degree", degree, lower=1)
    require(degree == np.round(degree), "degree must be an integer")
    coefficient = quantity("coefficient", coefficient)
    reference_radius = quantity("reference_radius", reference_radius, lower=0, strict_lower=True)
    x = np.sin(latitude)
    scale = (reference_radius / distance) ** (degree + 2) * coefficient
    north = scale * lpmv(1, degree, x)
    down = -(degree + 1) * scale * lpmv(0, degree, x)
    return np.stack(np.broadcast_arrays(north, down), axis=-1)


@model
def gauss_dipole_moment(g10, g11, h11, reference_radius):
    """Dipole moment from the degree-1 Gauss coefficients.

    Gauss's analysis weighed the Earth's magnet from surface measurements alone.

    Parameters
    ----------
    g10, g11, h11 : float or array_like
        T; constraint: real. Schmidt semi-normalized degree-1 coefficients.
    reference_radius : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [A m^2]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Magnitude of the best-fitting geocentric (tilted) dipole; the tilt is ``arccos(abs(g10) / B0)``.
    Equation: ``4 * np.pi * reference_radius**3 / MU0 * np.sqrt(g10**2 + g11**2 + h11**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import gauss_dipole_moment
    >>> float(gauss_dipole_moment(3, 4, 12, 1))
    130000000.0
    """
    g10 = quantity("g10", g10)
    g11 = quantity("g11", g11)
    h11 = quantity("h11", h11)
    reference_radius = quantity("reference_radius", reference_radius, lower=0, strict_lower=True)
    return 4 * np.pi * reference_radius**3 / MU0 * np.sqrt(g10**2 + g11**2 + h11**2)


@model
def field_elements(north, east, down):
    """Declination, inclination and intensity of a field vector.

    Magnetic directions in rocks and observatories are reported as two angles and a strength.

    Parameters
    ----------
    north, east, down : float or array_like
        T; constraint: real, not all zero.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value: declination [rad] east of north in (-pi, pi], inclination [rad] positive down, intensity [T].

    Notes
    -----
    Local north-east-down frame. A reversed field has declination near pi and inclination of the opposite sign.
    Equation: ``[np.arctan2(east, north), np.arctan2(down, np.hypot(north, east)), np.sqrt(north**2 + east**2 + down**2)]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import field_elements
    >>> np.round(field_elements(1, 1, 1.4142135623730951), 5).tolist()
    [0.7854, 0.7854, 2.0]
    """
    north = quantity("north", north)
    east = quantity("east", east)
    down = quantity("down", down)
    intensity = np.sqrt(north**2 + east**2 + down**2)
    require(intensity > 0, "the field vector must be nonzero")
    declination = np.arctan2(east, north)
    inclination = np.arctan2(down, np.hypot(north, east))
    return np.stack(np.broadcast_arrays(declination, inclination, intensity), axis=-1)


@model
def direction_angle(declination_a, inclination_a, declination_b, inclination_b):
    """Angle between two magnetic directions.

    A remanence more than 90 degrees from the present field records a reversed polarity.

    Parameters
    ----------
    declination_a, declination_b : float or array_like
        rad; constraint: real.
    inclination_a, inclination_b : float or array_like
        rad; constraint: -pi/2 <= inclination <= pi/2.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad] in [0, pi]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Great-circle distance between unit vectors; intensities are ignored.
    Equation: ``np.arccos(np.sin(Ia) * np.sin(Ib) + np.cos(Ia) * np.cos(Ib) * np.cos(Da - Db))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import direction_angle
    >>> round(float(direction_angle(0, 0.5, 0, -0.5)), 6)
    1.0
    """
    declination_a = quantity("declination_a", declination_a)
    inclination_a = _latitude("inclination_a", inclination_a)
    declination_b = quantity("declination_b", declination_b)
    inclination_b = _latitude("inclination_b", inclination_b)
    cosine = np.sin(inclination_a) * np.sin(inclination_b) + np.cos(inclination_a) * np.cos(inclination_b) * np.cos(
        declination_a - declination_b
    )
    return np.arccos(np.clip(cosine, -1, 1))


@model
def virtual_geomagnetic_pole(site_latitude, site_longitude, declination, inclination):
    """Virtual geomagnetic pole of a site direction.

    Poles from rocks of one age on one continent should agree; when they wander, the continent moved.

    Parameters
    ----------
    site_latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2.
    site_longitude : float or array_like
        rad; constraint: real.
    declination : float or array_like
        rad; constraint: real.
    inclination : float or array_like
        rad; constraint: -pi/2 <= inclination <= pi/2.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [rad]: pole latitude and pole longitude in [0, 2 pi).

    Notes
    -----
    The pole of the geocentric dipole that gives the observed direction at the site; a reversed direction gives a
    southern pole. Its angular distance from the site is ``p = np.arctan2(2 * np.cos(I), np.sin(I))``, along the
    great circle leaving the site at azimuth ``D``.
    Equation: ``pole = site * np.cos(p) + (north * np.cos(D) + east * np.sin(D)) * np.sin(p)`` with the site's unit
    position vector and local north and east unit vectors; latitude ``np.arctan2(z, np.hypot(x, y))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import virtual_geomagnetic_pole
    >>> np.round(virtual_geomagnetic_pole(0, 0, 1.5707963267948966, 0), 5).tolist()
    [0.0, 1.5708]
    """
    site_latitude = _latitude("site_latitude", site_latitude)
    site_longitude = quantity("site_longitude", site_longitude)
    declination = quantity("declination", declination)
    inclination = _latitude("inclination", inclination)
    p = np.arctan2(2 * np.cos(inclination), np.sin(inclination))
    sin_lat, cos_lat, sin_lon, cos_lon = np.sin(site_latitude), np.cos(site_latitude), np.sin(site_longitude), np.cos(site_longitude)
    along = np.cos(p)
    toward_north, toward_east = np.sin(p) * np.cos(declination), np.sin(p) * np.sin(declination)
    x = cos_lat * cos_lon * along - sin_lat * cos_lon * toward_north - sin_lon * toward_east
    y = cos_lat * sin_lon * along - sin_lat * sin_lon * toward_north + cos_lon * toward_east
    z = sin_lat * along + cos_lat * toward_north
    pole_latitude = np.arctan2(z, np.hypot(x, y))
    pole_longitude = np.mod(np.arctan2(y, x), 2 * np.pi)
    return np.stack(np.broadcast_arrays(pole_latitude, pole_longitude), axis=-1)


@model
def block_magnetic_anomaly(x, left, right, top, bottom, magnetization):
    """Magnetic anomaly of a vertically magnetized two-dimensional block.

    Strips of normally and reversely magnetized crust give the alternating anomalies over mid-ocean ridges.

    Parameters
    ----------
    x : float or array_like
        m; constraint: real. Horizontal position of the observer, perpendicular to the strike.
    left, right : float or array_like
        m; constraint: left < right. Edges of the block.
    top, bottom : float or array_like
        m; constraint: 0 < top < bottom. Depths below the observation level.
    magnetization : float or array_like
        A/m; constraint: real. Positive downward (normal polarity at the north magnetic pole).

    Returns
    -------
    numpy scalar or ndarray
        Model value [T]: vertical field anomaly, equal to the total-field anomaly reduced to the pole.

    Notes
    -----
    Infinite strike length, uniform vertical magnetization and a vertical ambient field. A laterally infinite layer has
    no external field.
    Equation: ``MU0 * magnetization / (2 * np.pi) * (phi(top) - phi(bottom))`` with
    ``phi(z) = np.arctan((right - x) / z) - np.arctan((left - x) / z)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import block_magnetic_anomaly
    >>> round(float(block_magnetic_anomaly(0, -1, 1, 1, 1.7320508075688772, 9549296.58551372)), 6)
    1.0
    """
    x = quantity("x", x)
    left = quantity("left", left)
    right = quantity("right", right)
    top = quantity("top", top, lower=0, strict_lower=True)
    bottom = quantity("bottom", bottom, lower=0, strict_lower=True)
    magnetization = quantity("magnetization", magnetization)
    require(left < right, "left must lie west of right")
    require(top < bottom, "top must be shallower than bottom")

    def subtended(depth):
        return np.arctan((right - x) / depth) - np.arctan((left - x) / depth)

    return MU0 * magnetization / (2 * np.pi) * (subtended(top) - subtended(bottom))


def _ferrers(order, degree, x):
    """Associated Legendre function without the Condon-Shortley phase; zero for order > degree."""
    return (-1.0) ** order * lpmv(order, degree, x) if order <= degree else np.zeros_like(x)


def _schmidt_factor(degree, order):
    return np.sqrt(np.where(order == 0, 1.0, 2.0) * np.exp(gammaln(degree - order + 1) - gammaln(degree + order + 1)))


@model
def schmidt_legendre(degree, order, latitude):
    """Schmidt semi-normalized associated Legendre function of the sine of latitude.

    Geomagnetism normalizes its spherical harmonics so that every order of one degree carries the same mean square,
    which makes Gauss coefficients of different orders directly comparable.

    Parameters
    ----------
    degree : float or array_like
        1; constraint: integer >= 0.
    order : float or array_like
        1; constraint: integer, 0 <= order <= degree.
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2. Geocentric latitude.

    Returns
    -------
    numpy scalar or ndarray
        Model value [1]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Without the Condon-Shortley phase, as in the IGRF (Winch et al. 2005); over the sphere
    ``mean((S_n^m cos(m lon))**2) = 1 / (2 n + 1)``. Degrees above about 150 lose accuracy in ``lpmv``.
    Equation: ``np.sqrt((2 - (order == 0)) * factorial(n - m) / factorial(n + m)) * (-1)**m * lpmv(m, n, np.sin(latitude))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.geophysics import schmidt_legendre
    >>> round(float(schmidt_legendre(2, 1, 0.5235987755982988)), 6)
    0.75
    """
    degree = quantity("degree", degree, lower=0)
    order = quantity("order", order, lower=0)
    require((degree == np.round(degree)) & (order == np.round(order)), "degree and order must be integers")
    require(order <= degree, "order must not exceed degree")
    latitude = _latitude("latitude", latitude)
    return _schmidt_factor(degree, order) * (-1.0) ** order * lpmv(order, degree, np.sin(latitude))


@model
def geomagnetic_field(latitude, longitude, distance, g, h, reference_radius):
    """Magnetic field of an internal spherical-harmonic model from its Gauss coefficients.

    A table of Gauss coefficients, such as the International Geomagnetic Reference Field, fixes the main field and its
    declination and inclination everywhere above the core.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 < latitude < pi/2. Geocentric latitude.
    longitude : float or array_like
        rad, east; constraint: real.
    distance : float or array_like
        m, from the Earth's centre; constraint: positive.
    g, h : array_like, shape (N + 1, N + 1), N >= 1
        T, Gauss coefficients ``g[n, m]`` and ``h[n, m]``; constraint: real, zero where ``m > n``, in row ``n = 0``
        and in ``h[:, 0]``.
    reference_radius : float
        m; constraint: positive scalar, 6371.2 km for the IGRF.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [T]: the northward, eastward and downward components on the last axis.

    Notes
    -----
    Potential ``V = a * sum((a / r)**(n + 1) * (g[n, m] cos(m lon) + h[n, m] sin(m lon)) * S_n^m)`` with
    :func:`schmidt_legendre` ``S_n^m`` of the colatitude ``theta``, and ``B = -grad V``. Spherical Earth, internal
    sources only; the m = 0 terms equal :func:`zonal_field`. The east component divides by ``sin(theta)``, so the
    geographic poles are excluded.
    Equation: ``sum((a / r)**(n + 2) * [(g cos + h sin) dS/dtheta, m (g sin - h cos) S / sin(theta),
    -(n + 1) (g cos + h sin) S])``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.geophysics import geomagnetic_field
    >>> (np.round(geomagnetic_field(0, 1.5707963267948966, 1, [[0, 0], [0, 1]], [[0, 0], [0, 0]], 1), 5) + 0.0).tolist()
    [0.0, 1.0, 0.0]
    """
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2, strict_lower=True, strict_upper=True)
    longitude = quantity("longitude", longitude)
    distance = quantity("distance", distance, lower=0, strict_lower=True)
    g, h = quantity("g", g), quantity("h", h)
    reference_radius = quantity("reference_radius", reference_radius, lower=0, strict_lower=True)
    require(reference_radius.ndim == 0, "reference_radius must be a scalar")
    require(g.ndim == 2 and g.shape[0] == g.shape[1] >= 2 and h.shape == g.shape, "g and h must be square (N + 1, N + 1) tables")
    require((np.triu(g, 1) == 0).all() & (np.triu(h, 1) == 0).all(), "coefficients with m > n must be zero")
    require((g[0] == 0).all() & (h[0] == 0).all() & (h[:, 0] == 0).all(), "row n = 0 and h[:, 0] must be zero")
    latitude, longitude, distance = np.broadcast_arrays(latitude, longitude, distance)
    x, sine = np.sin(latitude), np.cos(latitude)
    north, east, down = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
    for n in range(1, g.shape[0]):
        scale = (reference_radius / distance) ** (n + 2)
        for m in range(n + 1):
            factor = _schmidt_factor(n, m)
            cosine, sinusoid = np.cos(m * longitude), np.sin(m * longitude)
            even = g[n, m] * cosine + h[n, m] * sinusoid
            odd = g[n, m] * sinusoid - h[n, m] * cosine
            value = factor * _ferrers(m, n, x)
            slope = -_ferrers(1, n, x) if m == 0 else (n + m) * (n - m + 1) * _ferrers(m - 1, n, x) / 2 - _ferrers(m + 1, n, x) / 2
            north += scale * even * factor * slope
            east += scale * m * odd * value / sine
            down -= scale * (n + 1) * even * value
    return np.stack([north, east, down], axis=-1)
