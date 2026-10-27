"""Finite rotations: rotation matrices, reconstructed positions, and composition and interpolation of Euler poles.

Angles and coordinates are radians. A finite rotation is a pole (latitude, longitude) and an angle, positive
counterclockwise seen from above the pole (the right-hand rule). Composite rotations are returned with the angle in
``[0, pi]`` and the pole flipped to match; a zero rotation is reported about the north pole.
"""

import numpy as np

from earthsciencekit._validation import model, quantity

__all__ = [
    "finite_rotation_matrix",
    "rotate_position",
    "compose_finite_rotations",
    "interpolate_finite_rotation",
]


def _unit_vector(latitude, longitude):
    components = np.cos(latitude) * np.cos(longitude), np.cos(latitude) * np.sin(longitude), np.sin(latitude)
    return np.stack(np.broadcast_arrays(*components), axis=-1)


def _rotation(prefix, latitude, longitude, angle):
    latitude = quantity(f"{prefix}latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    longitude = quantity(f"{prefix}longitude", longitude)
    angle = quantity("angle", angle)
    return latitude, longitude, angle


def _quaternion(latitude, longitude, angle):
    axis = _unit_vector(latitude, longitude)
    half = np.asarray(angle)[..., None] / 2
    scalar, vector = np.broadcast_arrays(np.cos(half), np.sin(half) * axis)
    return np.concatenate([scalar[..., :1], vector], axis=-1)


def _multiply(a, b):
    """Hamilton product: rotating by ``b`` and then by ``a``."""
    a, b = np.broadcast_arrays(a, b)
    w = a[..., 0] * b[..., 0] - np.sum(a[..., 1:] * b[..., 1:], axis=-1)
    v = a[..., :1] * b[..., 1:] + b[..., :1] * a[..., 1:] + np.cross(a[..., 1:], b[..., 1:])
    return np.concatenate([w[..., None], v], axis=-1)


def _pole(quaternion):
    """Pole latitude, longitude and angle in [0, pi] of a unit quaternion; the north pole for no rotation."""
    q = np.where(quaternion[..., :1] < 0, -quaternion, quaternion)
    sine = np.linalg.norm(q[..., 1:], axis=-1)
    angle = 2 * np.arctan2(sine, q[..., 0])
    axis = np.where(sine[..., None] > 0, q[..., 1:] / np.where(sine > 0, sine, 1.0)[..., None], [0.0, 0.0, 1.0])
    latitude = np.arcsin(np.clip(axis[..., 2], -1, 1))
    longitude = np.arctan2(axis[..., 1], axis[..., 0])
    return np.stack([latitude, longitude, angle], axis=-1)


@model
def finite_rotation_matrix(pole_latitude, pole_longitude, angle):
    """Rotation matrix of a finite rotation about an Euler pole.

    Any displacement of a rigid plate on a sphere is a rotation about one pole, so a matrix carries every point of a
    plate to its reconstructed position.

    Parameters
    ----------
    pole_latitude : float or array_like
        rad; constraint: -pi/2 <= pole_latitude <= pi/2.
    pole_longitude : float or array_like
        rad; constraint: real.
    angle : float or array_like
        rad, counterclockwise seen from above the pole; constraint: real.

    Returns
    -------
    ndarray, shape (..., 3, 3)
        Model value [1]: the orthogonal matrix acting on geocentric unit vectors (x to 0 degE on the equator, z north).

    Notes
    -----
    Rodrigues's formula, the finite-rotation form of Euler's theorem as used by Bullard et al. (1965) and Cox and Hart
    (1986).
    Equation: ``np.cos(angle) * I + np.sin(angle) * [n]x + (1 - np.cos(angle)) * n n^T`` with ``n`` the pole's unit
    vector and ``[n]x`` its cross-product matrix.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.plate_tectonics import finite_rotation_matrix
    >>> (np.round(finite_rotation_matrix(1.5707963267948966, 0, 1.5707963267948966), 5) + 0.0).tolist()
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    """
    latitude, longitude, angle = _rotation("pole_", pole_latitude, pole_longitude, angle)
    n = _unit_vector(latitude, longitude)
    angle = np.broadcast_to(angle, n.shape[:-1])[..., None, None]
    zero = np.zeros(n.shape[:-1])
    cross = np.stack(
        [
            np.stack([zero, -n[..., 2], n[..., 1]], axis=-1),
            np.stack([n[..., 2], zero, -n[..., 0]], axis=-1),
            np.stack([-n[..., 1], n[..., 0], zero], axis=-1),
        ],
        axis=-2,
    )
    return np.cos(angle) * np.eye(3) + np.sin(angle) * cross + (1 - np.cos(angle)) * n[..., :, None] * n[..., None, :]


@model
def rotate_position(latitude, longitude, pole_latitude, pole_longitude, angle):
    """Position of a point after a finite rotation.

    Rotating a plate's present coastline or magnetic anomalies back by its total reconstruction rotation shows where
    they lay in the past.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: -pi/2 <= latitude <= pi/2.
    longitude : float or array_like
        rad; constraint: real.
    pole_latitude : float or array_like
        rad; constraint: -pi/2 <= pole_latitude <= pi/2.
    pole_longitude : float or array_like
        rad; constraint: real.
    angle : float or array_like
        rad, counterclockwise seen from above the pole; constraint: real.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [rad]: the rotated latitude and longitude (in ``(-pi, pi]``) on the last axis.

    Notes
    -----
    Rigid rotation of a geocentric unit vector; positions on a sphere, not an ellipsoid.
    Equation: ``R(pole, angle) @ [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.plate_tectonics import rotate_position
    >>> np.round(rotate_position(0, 0, 1.5707963267948966, 0, 1.5707963267948966), 5).tolist()
    [0.0, 1.5708]
    """
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    longitude = quantity("longitude", longitude)
    matrix = finite_rotation_matrix(pole_latitude, pole_longitude, angle)
    point = _unit_vector(latitude, longitude)
    rotated = np.einsum("...ij,...j->...i", matrix, point)
    return np.stack([np.arcsin(np.clip(rotated[..., 2], -1, 1)), np.arctan2(rotated[..., 1], rotated[..., 0])], axis=-1)


@model
def compose_finite_rotations(pole_latitude_1, pole_longitude_1, angle_1, pole_latitude_2, pole_longitude_2, angle_2):
    """Single finite rotation equal to one rotation followed by another.

    Rotations add around a plate circuit only as finite rotations, not as angles: B relative to A followed by C relative
    to B gives C relative to A about a new pole.

    Parameters
    ----------
    pole_latitude_1, pole_longitude_1 : float or array_like
        rad, the first rotation's pole; constraint: latitude within [-pi/2, pi/2].
    angle_1 : float or array_like
        rad; constraint: real.
    pole_latitude_2, pole_longitude_2 : float or array_like
        rad, the second rotation's pole, applied after the first; constraint: latitude within [-pi/2, pi/2].
    angle_2 : float or array_like
        rad; constraint: real.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [rad]: pole latitude, pole longitude and angle (in ``[0, pi]``) on the last axis.

    Notes
    -----
    Finite rotations do not commute; the result is ``R2 @ R1``, computed with unit quaternions (Cox and Hart 1986).
    Compose a rotation with its own pole and minus its angle to get the identity, and ``R2 @ R1^-1`` (pass ``-angle_1``
    as the first rotation) to get the stage rotation between two reconstructions.
    Equation: ``pole_and_angle(q2 * q1)`` with ``q = [np.cos(angle / 2), np.sin(angle / 2) * n]``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.plate_tectonics import compose_finite_rotations
    >>> np.round(compose_finite_rotations(1.5707963267948966, 0, 1.5707963267948966, 1.5707963267948966, 0, 1.5707963267948966), 5).tolist()
    [1.5708, 0.0, 3.14159]
    """
    first = _quaternion(*_rotation("pole_", pole_latitude_1, pole_longitude_1, angle_1))
    second = _quaternion(*_rotation("pole_", pole_latitude_2, pole_longitude_2, angle_2))
    return _pole(_multiply(second, first))


@model
def interpolate_finite_rotation(pole_latitude_1, pole_longitude_1, angle_1, pole_latitude_2, pole_longitude_2, angle_2, fraction):
    """Finite rotation part of the way between two reconstructions.

    Between two dated reconstructions a plate is assumed to turn steadily about the stage pole, which gives its
    position at any intermediate age.

    Parameters
    ----------
    pole_latitude_1, pole_longitude_1 : float or array_like
        rad, the rotation at the start (for example the younger age); constraint: latitude within [-pi/2, pi/2].
    angle_1 : float or array_like
        rad; constraint: real.
    pole_latitude_2, pole_longitude_2 : float or array_like
        rad, the rotation at the end (the older age); constraint: latitude within [-pi/2, pi/2].
    angle_2 : float or array_like
        rad; constraint: real.
    fraction : float or array_like
        1, 0 at the first rotation and 1 at the second; constraint: between 0 and 1.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [rad]: pole latitude, pole longitude and angle (in ``[0, pi]``) on the last axis.

    Notes
    -----
    A constant-rate stage rotation about the fixed-frame stage pole of ``R2 @ R1^-1``, taking the shorter way round
    (quaternion slerp); for rotations about one pole it interpolates the angle linearly. Ages map to the fraction
    linearly, for example ``(age - age_1) / (age_2 - age_1)``.
    Equation: ``pole_and_angle((q2 * conj(q1))**fraction * q1)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> import numpy as np
    >>> from earthsciencekit.plate_tectonics import interpolate_finite_rotation
    >>> np.round(interpolate_finite_rotation(1.5707963267948966, 0, 0, 1.5707963267948966, 0, 1, 0.5), 5).tolist()
    [1.5708, 0.0, 0.5]
    """
    first = _quaternion(*_rotation("pole_", pole_latitude_1, pole_longitude_1, angle_1))
    second = _quaternion(*_rotation("pole_", pole_latitude_2, pole_longitude_2, angle_2))
    fraction = quantity("fraction", fraction, lower=0, upper=1)
    stage = _pole(_multiply(second, first * [1, -1, -1, -1]))
    partial = _quaternion(stage[..., 0], stage[..., 1], stage[..., 2] * fraction)
    return _pole(_multiply(partial, first))
