"""Rotating dynamics: the Coriolis parameter, geostrophy, the Rossby number, Rossby waves and barotropic vorticity.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic
from earthsciencekit.constants import OMEGA_EARTH

__all__ = [
    "coriolis_parameter",
    "geostrophic_speed",
    "rossby_number",
    "rossby_wave_speed",
    "simulate_barotropic_vorticity",
]


@model
def coriolis_parameter(latitude):
    """Rotating-Earth Coriolis parameter.

    Latitude controls rotational deflection and geophysical balance.

    Parameters
    ----------
    latitude : float or array_like
        rad; constraint: latitude.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Spherical rotating-Earth approximation; input is radians.
    Equation: ``2 * OMEGA_EARTH * np.sin(latitude)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import coriolis_parameter
    >>> float(coriolis_parameter(0))
    0.0
    """
    latitude = quantity("latitude", latitude, lower=-np.pi / 2, upper=np.pi / 2)
    return 2 * OMEGA_EARTH * np.sin(latitude)


@model
def geostrophic_speed(pressure_gradient, density, coriolis):
    """Pressure-gradient geostrophic component.

    Pressure gradients and rotation can balance in large-scale flow.

    Parameters
    ----------
    pressure_gradient : float or array_like
        Pa/m; constraint: real.
    density : float or array_like
        kg/m^3; constraint: positive.
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Signed cross-gradient component; invalid at the equator where coriolis vanishes.
    Equation: ``pressure_gradient / (density * coriolis)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import geostrophic_speed
    >>> float(geostrophic_speed(0.001, 1, 0.0001))
    10.0
    """
    pressure_gradient = quantity("pressure_gradient", pressure_gradient)
    density = quantity("density", density, lower=0, strict_lower=True)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    return pressure_gradient / (density * coriolis)


@model
def rossby_number(speed, length, coriolis):
    """Rossby dynamical scaling.

    Rossby number measures inertia relative to rotational effects.

    Parameters
    ----------
    speed : float or array_like
        m/s; constraint: nonnegative.
    length : float or array_like
        m; constraint: positive.
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/m^3]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Magnitude scaling only; small values support but do not prove geostrophic balance.
    Equation: ``speed / (np.abs(coriolis) * length)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import rossby_number
    >>> float(rossby_number(10, 1000000.0, 0.0001))
    0.1
    """
    speed = quantity("speed", speed, lower=0)
    length = quantity("length", length, lower=0, strict_lower=True)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    return speed / (np.abs(coriolis) * length)


@model
def rossby_wave_speed(zonal_wind, beta, wavenumber):
    """Barotropic Rossby-wave zonal phase speed.

    The latitude gradient of the Coriolis parameter makes vorticity waves drift westward.

    Parameters
    ----------
    zonal_wind : float or array_like
        m/s; constraint: real. Uniform background eastward wind.
    beta : float or array_like
        1/(m s); constraint: nonnegative. Northward gradient of the Coriolis parameter.
    wavenumber : float or array_like
        rad/m; constraint: positive. Total horizontal wavenumber.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m/s]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Nondivergent barotropic flow on a beta plane; positive values are eastward.
    Equation: ``zonal_wind - beta / wavenumber**2``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.

    Examples
    --------
    >>> from earthsciencekit.atmosphere import rossby_wave_speed
    >>> float(rossby_wave_speed(10, 2e-11, 1e-06))
    -10.0
    """
    zonal_wind = quantity("zonal_wind", zonal_wind)
    beta = quantity("beta", beta, lower=0)
    wavenumber = quantity("wavenumber", wavenumber, lower=0, strict_lower=True)
    return zonal_wind - beta / wavenumber**2


@safe_arithmetic
def simulate_barotropic_vorticity(initial_vorticity, time, beta, spacing):
    """Integrate the barotropic vorticity equation on a periodic beta plane.

    Parameters
    ----------
    initial_vorticity : array_like, shape (ny, nx), ny, nx >= 4
        Relative vorticity [1/s] on an equally spaced doubly periodic grid,
        rows northward and columns eastward.
    time : array_like
        Strictly increasing nonnegative output times [s]. Initial state belongs
        to the first supplied time; values need not start at zero.
    beta : float
        Nonnegative scalar northward gradient of the Coriolis parameter [1/(m s)].
    spacing : float
        Positive scalar grid spacing [m], the same in both directions.

    Returns
    -------
    TimeSeries
        Relative vorticity at each grid point [1/s], flattened row by row.
        Columns are labeled cell_<row>_<column>.

    Notes
    -----
    Solves ``d(zeta)/dt + J(psi, zeta) + beta * d(psi)/dx = 0`` with
    ``laplacian(psi) = zeta``. Derivatives are pseudo-spectral (FFT) with the
    two-thirds dealiasing rule, and ``integrators.integrate_ode`` steps in time.
    The domain mean of vorticity is removed and conserved; energy and enstrophy
    are conserved up to the time-stepping error. A single Fourier mode is an
    exact Rossby wave. There is no forcing, friction, divergence or sphere.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode

    zeta = quantity("initial_vorticity", initial_vorticity)
    require(zeta.ndim == 2 and min(zeta.shape) >= 4, "initial_vorticity must be a grid of at least 4 x 4 cells")
    time = time_grid(time)
    beta = quantity("beta", beta, lower=0)
    spacing = quantity("spacing", spacing, lower=0, strict_lower=True)
    require(beta.ndim == spacing.ndim == 0, "beta and spacing must be scalars")
    shape = zeta.shape
    ky = 2 * np.pi * np.fft.fftfreq(shape[0], spacing)[:, None]
    kx = 2 * np.pi * np.fft.rfftfreq(shape[1], spacing)[None, :]
    k2 = kx**2 + ky**2
    inverse = np.divide(-1, k2, out=np.zeros_like(k2), where=k2 > 0)
    keep = (np.abs(ky) < 2 / 3 * np.abs(ky).max()) & (np.abs(kx) < 2 / 3 * kx.max())

    def rhs(t, state):
        spectrum = np.fft.rfft2(state.reshape(shape)) * keep
        psi = spectrum * inverse
        u = np.fft.irfft2(-1j * ky * psi, s=shape)
        v = np.fft.irfft2(1j * kx * psi, s=shape)
        zeta_x = np.fft.irfft2(1j * kx * spectrum, s=shape)
        zeta_y = np.fft.irfft2(1j * ky * spectrum, s=shape)
        advection = np.fft.rfft2(u * zeta_x + v * zeta_y) * keep
        return np.fft.irfft2(-advection - beta * 1j * kx * psi, s=shape).ravel()

    initial = np.fft.irfft2(np.fft.rfft2(zeta) * keep * (k2 > 0), s=shape).ravel()
    names = tuple(f"cell_{i}_{j}" for i in range(shape[0]) for j in range(shape[1]))
    return integrate_ode(rhs, initial, time, names=names, units=("1/s",) * initial.size)
