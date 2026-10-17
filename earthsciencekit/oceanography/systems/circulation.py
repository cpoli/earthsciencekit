"""Circulation: geostrophy, Ekman and Sverdrup transport, the Stommel gyre, the Rossby radius and thermohaline boxes.

SI units unless stated otherwise.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require, safe_arithmetic

__all__ = [
    "geostrophic_current",
    "ekman_transport",
    "rossby_radius",
    "volume_transport",
    "ekman_spiral_velocity",
    "sverdrup_transport",
    "stommel_streamfunction",
    "thermohaline_forcing",
    "simulate_thermohaline_boxes",
]


@model
def geostrophic_current(surface_slope, coriolis, gravity):
    """Sea-surface geostrophic component.

    Sea-surface tilt can reveal large-scale currents.

    Parameters
    ----------
    surface_slope : float or array_like
        1; constraint: real.
    coriolis : float or array_like
        1/s; constraint: nonzero.
    gravity : float or array_like
        m/s^2; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Signed cross-slope component away from the equator; no ageostrophic velocity.
    Equation: ``gravity * surface_slope / coriolis``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    surface_slope = quantity("surface_slope", surface_slope)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    gravity = quantity("gravity", gravity, lower=0, strict_lower=True)
    return gravity * surface_slope / coriolis


@model
def ekman_transport(wind_stress, density, coriolis):
    """Cross-wind Ekman transport magnitude.

    Wind stress and rotation drive vertically integrated transport.

    Parameters
    ----------
    wind_stress : float or array_like
        Pa; constraint: real.
    density : float or array_like
        kg/m^3; constraint: positive.
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Signed cross-wind component per unit width; direction requires a consistent vector convention.
    Equation: ``wind_stress / (density * coriolis)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wind_stress = quantity("wind_stress", wind_stress)
    density = quantity("density", density, lower=0, strict_lower=True)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    return wind_stress / (density * coriolis)


@model
def ekman_spiral_velocity(depth, wind_stress, density, eddy_viscosity, coriolis):
    """Velocity profile of the Ekman spiral.

    Wind stress diffused downward by friction turns and decays with depth under rotation.

    Parameters
    ----------
    depth : float or array_like
        m; constraint: nonnegative.
    wind_stress : float or array_like
        Pa; constraint: real. Signed stress along the wind axis.
    density : float or array_like
        kg/m^3; constraint: positive.
    eddy_viscosity : float or array_like
        m^2/s; constraint: positive. Constant vertical eddy viscosity.
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    ndarray, shape (..., 2)
        Model value [m/s]: the downwind component and the crosswind component, positive to the right of the wind.

    Notes
    -----
    Steady, infinitely deep ocean with constant eddy viscosity and no pressure gradient. The surface current is 45 degrees
    to the right of the wind in the northern hemisphere (left in the southern). The crosswind component integrated over
    depth equals :func:`ekman_transport`; the downwind component integrates to zero.
    Equation: ``v0 * np.exp(-s) * [np.sin(np.pi / 4 - s), np.sign(coriolis) * np.cos(np.pi / 4 - s)]``
    with ``s = depth / d``, ``d = np.sqrt(2 * eddy_viscosity / abs(coriolis))`` and
    ``v0 = wind_stress / (density * np.sqrt(eddy_viscosity * abs(coriolis)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    depth = quantity("depth", depth, lower=0)
    wind_stress = quantity("wind_stress", wind_stress)
    density = quantity("density", density, lower=0, strict_lower=True)
    eddy_viscosity = quantity("eddy_viscosity", eddy_viscosity, lower=0, strict_lower=True)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    scaled = depth / np.sqrt(2 * eddy_viscosity / np.abs(coriolis))
    surface = wind_stress / (density * np.sqrt(eddy_viscosity * np.abs(coriolis)))
    decay = surface * np.exp(-scaled)
    return np.stack([decay * np.sin(np.pi / 4 - scaled), np.sign(coriolis) * decay * np.cos(np.pi / 4 - scaled)], axis=-1)


@model
def volume_transport(speed, width, depth):
    """Section volume transport.

    Integrating velocity over a section produces volume transport.

    Parameters
    ----------
    speed : float or array_like
        m/s; constraint: real.
    width : float or array_like
        m; constraint: positive.
    depth : float or array_like
        m; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Uniform signed velocity and rectangular section; output m cubed per second, not Sverdrups.
    Equation: ``speed * width * depth``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    speed = quantity("speed", speed)
    width = quantity("width", width, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    return speed * width * depth


@model
def sverdrup_transport(wind_stress_curl, density, beta):
    """Sverdrup meridional transport.

    Wind-stress curl drives interior flow across planetary vorticity contours.

    Parameters
    ----------
    wind_stress_curl : float or array_like
        Pa/m; constraint: real.
    density : float or array_like
        kg/m^3; constraint: positive.
    beta : float or array_like
        1/(m s); constraint: positive. Northward gradient of the Coriolis parameter.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^2/s]: depth-integrated northward transport per unit zonal width.

    Notes
    -----
    Steady linear interior balance away from boundaries, friction and the equator's singular Ekman layer.
    Negative curl (subtropical gyres) drives equatorward flow.
    Equation: ``wind_stress_curl / (density * beta)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wind_stress_curl = quantity("wind_stress_curl", wind_stress_curl)
    density = quantity("density", density, lower=0, strict_lower=True)
    beta = quantity("beta", beta, lower=0, strict_lower=True)
    return wind_stress_curl / (density * beta)


@model
def stommel_streamfunction(x, y, wind_stress, basin_width, basin_length, density, beta, friction):
    """Transport streamfunction of Stommel's wind-driven gyre.

    Bottom friction and the beta effect confine the return flow to a western boundary current.

    Parameters
    ----------
    x : float or array_like
        m; constraint: 0 <= x <= basin_width. Eastward distance from the western wall.
    y : float or array_like
        m; constraint: 0 <= y <= basin_length. Northward distance from the southern wall.
    wind_stress : float or array_like
        Pa; constraint: real. Amplitude of the zonal stress ``-wind_stress * cos(pi y / basin_length)``.
    basin_width : float or array_like
        m; constraint: positive.
    basin_length : float or array_like
        m; constraint: positive.
    density : float or array_like
        kg/m^3; constraint: positive.
    beta : float or array_like
        1/(m s); constraint: nonnegative. Zero gives Stommel's symmetric uniform-rotation gyre.
    friction : float or array_like
        1/s; constraint: positive. Linear bottom-drag coefficient on the depth-integrated flow.

    Returns
    -------
    numpy scalar or ndarray
        Model value [m^3/s]. Northward transport is the eastward derivative.

    Notes
    -----
    Steady linear barotropic vorticity balance ``friction * laplacian(psi) + beta * dpsi/dx = curl(tau) / density`` in a
    rectangular basin with ``psi = 0`` on the walls. Positive wind stress gives a clockwise (subtropical) gyre whose
    western boundary layer has width ``friction / beta``.
    Equation: ``wind_stress * basin_length / (np.pi * density * friction) * np.sin(np.pi * y / basin_length) * X``
    with ``X = 1 + ((np.exp(m2 * L) - 1) * np.exp(m1 * (x - L)) + (np.exp(-m1 * L) - 1) * np.exp(m2 * x))
    / (1 - np.exp((m2 - m1) * L))``, ``L = basin_width``,
    ``m1, m2 = (-beta +/- np.sqrt(beta**2 + 4 * friction**2 * k**2)) / (2 * friction)`` and ``k = np.pi / basin_length``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    x = quantity("x", x, lower=0)
    y = quantity("y", y, lower=0)
    wind_stress = quantity("wind_stress", wind_stress)
    basin_width = quantity("basin_width", basin_width, lower=0, strict_lower=True)
    basin_length = quantity("basin_length", basin_length, lower=0, strict_lower=True)
    density = quantity("density", density, lower=0, strict_lower=True)
    beta = quantity("beta", beta, lower=0)
    friction = quantity("friction", friction, lower=0, strict_lower=True)
    require(x <= basin_width, "x must lie inside the basin")
    require(y <= basin_length, "y must lie inside the basin")
    k = np.pi / basin_length
    root = np.sqrt(beta**2 + 4 * friction**2 * k**2)
    m1 = 2 * friction * k**2 / (beta + root)
    m2 = -(beta + root) / (2 * friction)
    shape = 1 + ((np.exp(m2 * basin_width) - 1) * np.exp(m1 * (x - basin_width)) + (np.exp(-m1 * basin_width) - 1) * np.exp(m2 * x)) / (
        1 - np.exp((m2 - m1) * basin_width)
    )
    return wind_stress * basin_length / (np.pi * density * friction) * np.sin(np.pi * y / basin_length) * shape


@model
def rossby_radius(buoyancy_frequency, depth, coriolis):
    """Baroclinic deformation scale.

    Stratification and rotation set a horizontal adjustment scale.

    Parameters
    ----------
    buoyancy_frequency : float or array_like
        1/s; constraint: positive.
    depth : float or array_like
        m; constraint: positive.
    coriolis : float or array_like
        1/s; constraint: nonzero.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    First-mode constant-N rigid-boundary estimate N H over pi f; not a general layered solution.
    Equation: ``buoyancy_frequency * depth / (np.pi * np.abs(coriolis))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    buoyancy_frequency = quantity("buoyancy_frequency", buoyancy_frequency, lower=0, strict_lower=True)
    depth = quantity("depth", depth, lower=0, strict_lower=True)
    coriolis = quantity("coriolis", coriolis)
    require(coriolis != 0, "coriolis must be nonzero")
    return buoyancy_frequency * depth / (np.pi * np.abs(coriolis))


@model
def thermohaline_forcing(flow, thermal_flow, haline_coefficient):
    """Salt forcing that holds Stommel's two-box overturning steady.

    Salt advection by the overturning must balance the freshwater-driven salinity contrast.

    Parameters
    ----------
    flow : float or array_like
        m^3/s; constraint: at most thermal_flow. Overturning, positive when polar water sinks.
    thermal_flow : float or array_like
        m^3/s; constraint: positive. Overturning the temperature contrast alone would drive.
    haline_coefficient : float or array_like
        m^3/s per g/kg; constraint: positive. Overturning lost per unit salinity contrast.

    Returns
    -------
    numpy scalar or ndarray
        Model value [g/kg m^3/s]: the virtual salt flux into the equatorial box.

    Notes
    -----
    Temperatures held fixed by fast restoring; flow ``q = thermal_flow - haline_coefficient * dS`` and steady salt budget
    ``abs(q) * dS = forcing``. Read as forcing against flow, the curve folds at ``q = thermal_flow / 2``: for
    ``0 < forcing < thermal_flow**2 / (4 * haline_coefficient)`` a thermal and a haline state coexist.
    Equation: ``np.abs(flow) * (thermal_flow - flow) / haline_coefficient``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    flow = quantity("flow", flow)
    thermal_flow = quantity("thermal_flow", thermal_flow, lower=0, strict_lower=True)
    haline_coefficient = quantity("haline_coefficient", haline_coefficient, lower=0, strict_lower=True)
    require(flow <= thermal_flow, "flow above thermal_flow implies a negative salinity contrast")
    return np.abs(flow) * (thermal_flow - flow) / haline_coefficient


@safe_arithmetic
def simulate_thermohaline_boxes(thermal_flow, haline_coefficient, salt_forcing, volume, time, *, initial_salinity_difference=0.0):
    """Integrate Stommel's two-box thermohaline circulation.

    Parameters
    ----------
    thermal_flow : float
        Positive scalar overturning driven by the fixed temperature contrast [m^3/s].
    haline_coefficient : float
        Positive scalar overturning lost per unit salinity contrast [m^3/s per g/kg].
    salt_forcing : float
        Nonnegative scalar virtual salt flux into the equatorial box [g/kg m^3/s].
    volume : float
        Positive scalar volume of each box [m^3].
    time : array_like
        Nonnegative strictly increasing sampled times [s].
    initial_salinity_difference : float
        Equatorial minus polar salinity at the first time [g/kg].

    Returns
    -------
    TimeSeries
        Salinity difference [g/kg] and overturning [m^3/s], positive when polar water sinks.

    Notes
    -----
    ``volume * dS/dt = salt_forcing - abs(q) * dS`` with ``q = thermal_flow - haline_coefficient * dS``.
    Temperatures are held fixed by fast restoring. Steady states satisfy
    :func:`thermohaline_forcing`; below the fold forcing
    thermal_flow**2 / (4 * haline_coefficient) the initial state decides
    between the thermal (q > 0) and haline (q < 0) regimes.
    """
    from earthsciencekit._validation import time_grid
    from earthsciencekit.integrators import integrate_ode
    from earthsciencekit.oceanography.core.base import TimeSeries

    thermal_flow = quantity("thermal_flow", thermal_flow, lower=0, strict_lower=True)
    haline_coefficient = quantity("haline_coefficient", haline_coefficient, lower=0, strict_lower=True)
    salt_forcing = quantity("salt_forcing", salt_forcing, lower=0)
    volume = quantity("volume", volume, lower=0, strict_lower=True)
    initial_salinity_difference = quantity("initial_salinity_difference", initial_salinity_difference)
    scalars = (thermal_flow, haline_coefficient, salt_forcing, volume, initial_salinity_difference)
    require(all(v.ndim == 0 for v in scalars), "box parameters must be scalars")
    time = time_grid(time)

    def rhs(t, state):
        flow = thermal_flow - haline_coefficient * state[0]
        return [(salt_forcing - np.abs(flow) * state[0]) / volume]

    salinity = integrate_ode(rhs, [initial_salinity_difference], time).values[:, 0]
    values = np.column_stack([salinity, thermal_flow - haline_coefficient * salinity])
    return TimeSeries(time, values, ("salinity_difference", "overturning"), ("g/kg", "m^3/s"))
