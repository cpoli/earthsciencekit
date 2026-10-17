"""Rigid-plate kinematics: Euler rotation, plate speed and displacement, relative Euler vectors and hotspot tracks.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require

__all__ = [
    "angular_speed",
    "plate_speed",
    "euler_velocity_x",
    "euler_velocity_y",
    "euler_velocity_z",
    "plate_displacement",
    "hotspot_age",
    "compose_euler_vectors",
]


@model
def angular_speed(rotation_angle, time):
    """Euler-pole angular speed.

    Rigid plate motion is rotation about an Euler pole.

    Parameters
    ----------
    rotation_angle : float or array_like
        rad; constraint: real.
    time : float or array_like
        s; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant angular rate; angular input is radians and output radians per second.
    Equation: ``rotation_angle / time``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    rotation_angle = quantity("rotation_angle", rotation_angle)
    time = quantity("time", time, lower=0, strict_lower=True)
    return rotation_angle / time


@model
def plate_speed(angular_speed, radius, angular_distance):
    """Rigid plate speed about an Euler pole.

    Distance from the rotation axis controls surface plate speed.

    Parameters
    ----------
    angular_speed : float or array_like
        rad/s; constraint: nonnegative.
    radius : float or array_like
        m; constraint: positive.
    angular_distance : float or array_like
        rad; constraint: width.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Spherical rigid plate, magnitude only; angular distance must be at most pi.
    Equation: ``angular_speed * radius * np.sin(angular_distance)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    angular_speed = quantity("angular_speed", angular_speed, lower=0)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    angular_distance = quantity("angular_distance", angular_distance, lower=0, upper=2 * np.pi)
    require(angular_distance <= np.pi, "angular distance exceeds pi")
    return angular_speed * radius * np.sin(angular_distance)


@model
def euler_velocity_x(omega_y, omega_z, position_y, position_z):
    """Euler rotation east-fixed x velocity.

    The cross product of rotation and position gives Cartesian plate velocity.

    Parameters
    ----------
    omega_y : float or array_like
        rad/s; constraint: real.
    omega_z : float or array_like
        rad/s; constraint: real.
    position_y : float or array_like
        m; constraint: real.
    position_z : float or array_like
        m; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Earth-fixed Cartesian axes; x is not the local east component.
    Equation: ``omega_y * position_z - omega_z * position_y``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    omega_y = quantity("omega_y", omega_y)
    omega_z = quantity("omega_z", omega_z)
    position_y = quantity("position_y", position_y)
    position_z = quantity("position_z", position_z)
    return omega_y * position_z - omega_z * position_y


@model
def euler_velocity_y(omega_z, omega_x, position_z, position_x):
    """Euler rotation y velocity.

    A rigid plate velocity field is tangent to its spherical path.

    Parameters
    ----------
    omega_z : float or array_like
        rad/s; constraint: real.
    omega_x : float or array_like
        rad/s; constraint: real.
    position_z : float or array_like
        m; constraint: real.
    position_x : float or array_like
        m; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Same Earth-fixed Cartesian convention as euler_velocity_x.
    Equation: ``omega_z * position_x - omega_x * position_z``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    omega_z = quantity("omega_z", omega_z)
    omega_x = quantity("omega_x", omega_x)
    position_z = quantity("position_z", position_z)
    position_x = quantity("position_x", position_x)
    return omega_z * position_x - omega_x * position_z


@model
def euler_velocity_z(omega_x, omega_y, position_x, position_y):
    """Euler rotation polar velocity.

    The polar velocity component completes the Euler-vector cross product.

    Parameters
    ----------
    omega_x : float or array_like
        rad/s; constraint: real.
    omega_y : float or array_like
        rad/s; constraint: real.
    position_x : float or array_like
        m; constraint: real.
    position_y : float or array_like
        m; constraint: real.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Earth-fixed components, not local vertical motion or deformation.
    Equation: ``omega_x * position_y - omega_y * position_x``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    omega_x = quantity("omega_x", omega_x)
    omega_y = quantity("omega_y", omega_y)
    position_x = quantity("position_x", position_x)
    position_y = quantity("position_y", position_y)
    return omega_x * position_y - omega_y * position_x


@model
def plate_displacement(speed, time):
    """Plate kinematic displacement.

    Velocity connects short instrumental records with geological displacement.

    Parameters
    ----------
    speed : float or array_like
        m/s; constraint: real.
    time : float or array_like
        s; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [rad/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Constant signed speed along a local coordinate; not a finite Euler rotation.
    Equation: ``speed * time``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    speed = quantity("speed", speed)
    time = quantity("time", time, lower=0)
    return speed * time


@model
def compose_euler_vectors(omega_ab, omega_bc):
    """Composition of relative Euler vectors around a plate circuit.

    Angular velocities of rigid plates add, so two relative motions predict the third.

    Parameters
    ----------
    omega_ab : array_like, shape (..., 3)
        rad/s; constraint: real. Earth-fixed Cartesian Euler vector of plate A relative to plate B.
    omega_bc : array_like, shape (..., 3)
        rad/s; constraint: real. Euler vector of plate B relative to plate C.

    Returns
    -------
    ndarray, shape (..., 3)
        Model value [rad/s]: the Euler vector of plate A relative to plate C. Leading axes broadcast.

    Notes
    -----
    Instantaneous (infinitesimal) rotations only; finite rotations do not commute and do not add.
    Closure of the circuit A-B-C means ``omega_ab + omega_bc + omega_ca = 0``.
    Equation: ``omega_ab + omega_bc``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    omega_ab = quantity("omega_ab", omega_ab)
    omega_bc = quantity("omega_bc", omega_bc)
    require(omega_ab.ndim >= 1 and omega_ab.shape[-1] == 3, "omega_ab must have three Cartesian components on its last axis")
    require(omega_bc.ndim >= 1 and omega_bc.shape[-1] == 3, "omega_bc must have three Cartesian components on its last axis")
    return omega_ab + omega_bc


@model
def hotspot_age(track_length, angular_speed, radius, pole_distance):
    """Age progression along a hotspot track.

    A plate moving over a fixed plume leaves volcanoes that age with distance along a small circle.

    Parameters
    ----------
    track_length : float or array_like
        m; constraint: nonnegative.
    angular_speed : float or array_like
        rad/s; constraint: positive.
    radius : float or array_like
        m; constraint: positive.
    pole_distance : float or array_like
        rad; constraint: open interval (0, pi).

    Returns
    -------
    numpy scalar or ndarray
        Model value [s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Fixed hotspot and constant Euler rotation of the plate over it; track length is measured along the small circle.
    Equation: ``track_length / (angular_speed * radius * np.sin(pole_distance))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    track_length = quantity("track_length", track_length, lower=0)
    angular_speed = quantity("angular_speed", angular_speed, lower=0, strict_lower=True)
    radius = quantity("radius", radius, lower=0, strict_lower=True)
    pole_distance = quantity("pole_distance", pole_distance, lower=0, upper=np.pi, strict_lower=True, strict_upper=True)
    return track_length / (angular_speed * radius * np.sin(pole_distance))
