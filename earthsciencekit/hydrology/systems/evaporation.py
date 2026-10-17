"""Evaporation: Dalton's law, Penman's combination equation and the Penman-Monteith canopy resistance.

SI units unless stated otherwise. Evaporation rates are kg/(m^2 s), which equals mm/s of liquid water.
"""

import numpy as np

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.atmosphere.systems.moisture import saturation_vapor_pressure
from earthsciencekit.constants import CP_DRY_AIR, EPSILON_WATER, LATENT_HEAT_VAPORIZATION, R_VAPOR, VON_KARMAN

__all__ = [
    "dalton_evaporation",
    "saturation_slope",
    "psychrometric_constant",
    "penman_evaporation",
    "aerodynamic_resistance",
    "bulk_surface_resistance",
    "penman_monteith_evaporation",
]


@model
def dalton_evaporation(saturation_vapor_pressure, vapor_pressure, transfer_coefficient):
    """Dalton's evaporation from a wet surface.

    Evaporation is proportional to how far the air is from saturation at the temperature of the surface, and stops when
    the air is saturated.

    Parameters
    ----------
    saturation_vapor_pressure : float or array_like
        Pa, at the surface temperature; constraint: nonnegative.
    vapor_pressure : float or array_like
        Pa, of the air; constraint: nonnegative.
    transfer_coefficient : float or array_like
        kg/(m^2 s Pa), a wind function that grows with wind speed; constraint: nonnegative.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/(m^2 s)], negative for condensation. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dalton (1802), a bulk mass-transfer law; the surface temperature it needs is rarely measured.
    Equation: ``transfer_coefficient * (saturation_vapor_pressure - vapor_pressure)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    saturation_vapor_pressure = quantity("saturation_vapor_pressure", saturation_vapor_pressure, lower=0)
    vapor_pressure = quantity("vapor_pressure", vapor_pressure, lower=0)
    transfer_coefficient = quantity("transfer_coefficient", transfer_coefficient, lower=0)
    return transfer_coefficient * (saturation_vapor_pressure - vapor_pressure)


@model
def saturation_slope(temperature):
    """Slope of the saturation vapour-pressure curve.

    The steepness of the Clausius-Clapeyron curve converts a surface-to-air temperature difference into a vapour-pressure
    difference.

    Parameters
    ----------
    temperature : float or array_like
        K; constraint: between 180 and 330.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa/K]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Derivative of :func:`earthsciencekit.atmosphere.saturation_vapor_pressure` (constant latent heat 2.5e6 J/kg).
    Equation: ``saturation_vapor_pressure(temperature) * LATENT_HEAT_VAPORIZATION / (R_VAPOR * temperature**2)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    temperature = quantity("temperature", temperature, lower=180, upper=330)
    return saturation_vapor_pressure(temperature) * LATENT_HEAT_VAPORIZATION / (R_VAPOR * temperature**2)


@model
def psychrometric_constant(pressure):
    """Psychrometric constant.

    It converts a temperature difference into the vapour-pressure difference that carries the same energy, as sensible
    rather than latent heat.

    Parameters
    ----------
    pressure : float or array_like
        Pa, air pressure; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [Pa/K], about 65 Pa/K at sea level. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Dry-air heat capacity and constant latent heat; equal transfer coefficients for heat and vapour.
    Equation: ``CP_DRY_AIR * pressure / (EPSILON_WATER * LATENT_HEAT_VAPORIZATION)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    return CP_DRY_AIR * pressure / (EPSILON_WATER * LATENT_HEAT_VAPORIZATION)


@model
def penman_evaporation(net_radiation, temperature, aerodynamic_evaporation, pressure):
    """Penman's combination equation for open-water evaporation.

    Penman eliminated the unknown surface temperature by combining the energy balance with Dalton's law, so evaporation
    follows from routine weather measurements.

    Parameters
    ----------
    net_radiation : float or array_like
        W/m^2, net radiation minus heat stored below the surface; constraint: real.
    temperature : float or array_like
        K, air temperature; constraint: between 180 and 330.
    aerodynamic_evaporation : float or array_like
        kg/(m^2 s), Dalton's evaporation evaluated with the air temperature, ``f(u) (e_s(T_air) - e_air)``;
        constraint: real.
    pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/(m^2 s)]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Penman (1948): the saturation curve is linearized between air and surface temperature, and heat and vapour share one
    transfer coefficient; a saturated surface (open water or wet grass).
    Equation: ``(saturation_slope(temperature) * net_radiation / LATENT_HEAT_VAPORIZATION + psychrometric_constant(pressure)
    * aerodynamic_evaporation) / (saturation_slope(temperature) + psychrometric_constant(pressure))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    net_radiation = quantity("net_radiation", net_radiation)
    temperature = quantity("temperature", temperature, lower=180, upper=330)
    aerodynamic_evaporation = quantity("aerodynamic_evaporation", aerodynamic_evaporation)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    slope, gamma = saturation_slope(temperature), psychrometric_constant(pressure)
    return (slope * net_radiation / LATENT_HEAT_VAPORIZATION + gamma * aerodynamic_evaporation) / (slope + gamma)


@model
def aerodynamic_resistance(wind_speed, measurement_height, displacement_height, momentum_roughness, vapor_roughness):
    """Aerodynamic resistance to vapour transfer above a canopy.

    Turbulence in the logarithmic wind profile carries vapour away from the canopy, faster in stronger wind over rougher
    surfaces.

    Parameters
    ----------
    wind_speed : float or array_like
        m/s, at the measurement height; constraint: positive.
    measurement_height : float or array_like
        m; constraint: positive.
    displacement_height : float or array_like
        m, zero-plane displacement (about two thirds of the canopy height); constraint: nonnegative.
    momentum_roughness : float or array_like
        m; constraint: positive.
    vapor_roughness : float or array_like
        m, often a tenth of the momentum roughness; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s/m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Neutral stability; the measurement height above the displacement height must exceed both roughness lengths. For
    FAO grass (0.12 m) at 2 m it is about ``208 / wind_speed``.
    Equation: ``np.log((measurement_height - displacement_height) / momentum_roughness) * np.log((measurement_height -
    displacement_height) / vapor_roughness) / (VON_KARMAN**2 * wind_speed)``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    wind_speed = quantity("wind_speed", wind_speed, lower=0, strict_lower=True)
    measurement_height = quantity("measurement_height", measurement_height, lower=0, strict_lower=True)
    displacement_height = quantity("displacement_height", displacement_height, lower=0)
    momentum_roughness = quantity("momentum_roughness", momentum_roughness, lower=0, strict_lower=True)
    vapor_roughness = quantity("vapor_roughness", vapor_roughness, lower=0, strict_lower=True)
    height = measurement_height - displacement_height
    require((height > momentum_roughness) & (height > vapor_roughness), "measurement height must lie above the roughness layer")
    return np.log(height / momentum_roughness) * np.log(height / vapor_roughness) / (VON_KARMAN**2 * wind_speed)


@model
def bulk_surface_resistance(stomatal_resistance, active_leaf_area_index):
    """Bulk surface (canopy) resistance from the leaf stomatal resistance.

    The stomata of all leaves act as resistances in parallel, so a dense canopy offers less resistance than a single leaf.

    Parameters
    ----------
    stomatal_resistance : float or array_like
        s/m, of a single well-illuminated leaf (about 100); constraint: nonnegative.
    active_leaf_area_index : float or array_like
        1, sunlit leaf area per ground area (about half the leaf area index); constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [s/m]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Monteith (1965) big-leaf canopy with parallel leaves; FAO grass has about 70 s/m.
    Equation: ``stomatal_resistance / active_leaf_area_index``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    stomatal_resistance = quantity("stomatal_resistance", stomatal_resistance, lower=0)
    active_leaf_area_index = quantity("active_leaf_area_index", active_leaf_area_index, lower=0, strict_lower=True)
    return stomatal_resistance / active_leaf_area_index


@model
def penman_monteith_evaporation(
    available_energy, temperature, vapor_pressure_deficit, air_density, aerodynamic_resistance, surface_resistance, pressure
):
    """Penman-Monteith evaporation from a vegetated surface.

    Monteith added the resistance of the stomata to Penman's equation, so the same formula covers open water, wet leaves
    and transpiring crops.

    Parameters
    ----------
    available_energy : float or array_like
        W/m^2, net radiation minus ground heat flux; constraint: real.
    temperature : float or array_like
        K, air temperature; constraint: between 180 and 330.
    vapor_pressure_deficit : float or array_like
        Pa, ``e_s(T_air) - e_air``; constraint: real.
    air_density : float or array_like
        kg/m^3; constraint: positive.
    aerodynamic_resistance : float or array_like
        s/m; constraint: positive.
    surface_resistance : float or array_like
        s/m; constraint: nonnegative.
    pressure : float or array_like
        Pa; constraint: positive.

    Returns
    -------
    numpy scalar or ndarray
        Model value [kg/(m^2 s)]. Scalar input returns a NumPy scalar; arrays broadcast.

    Notes
    -----
    Monteith (1965): a big-leaf canopy at one temperature, neutral stability; zero surface resistance recovers Penman
    with ``aerodynamic_evaporation = air_density * CP_DRY_AIR * vapor_pressure_deficit / (psychrometric_constant *
    LATENT_HEAT_VAPORIZATION * aerodynamic_resistance)``.
    Equation: ``(Delta * available_energy + air_density * CP_DRY_AIR * vapor_pressure_deficit / aerodynamic_resistance) /
    (LATENT_HEAT_VAPORIZATION * (Delta + gamma * (1 + surface_resistance / aerodynamic_resistance)))``.
    Invalid, empty, nonreal or nonfinite inputs raise ValueError.
    NumPy broadcasting applies; singular and overflowing outputs raise ValueError.
    """
    available_energy = quantity("available_energy", available_energy)
    temperature = quantity("temperature", temperature, lower=180, upper=330)
    vapor_pressure_deficit = quantity("vapor_pressure_deficit", vapor_pressure_deficit)
    air_density = quantity("air_density", air_density, lower=0, strict_lower=True)
    aerodynamic_resistance = quantity("aerodynamic_resistance", aerodynamic_resistance, lower=0, strict_lower=True)
    surface_resistance = quantity("surface_resistance", surface_resistance, lower=0)
    pressure = quantity("pressure", pressure, lower=0, strict_lower=True)
    slope, gamma = saturation_slope(temperature), psychrometric_constant(pressure)
    numerator = slope * available_energy + air_density * CP_DRY_AIR * vapor_pressure_deficit / aerodynamic_resistance
    return numerator / (LATENT_HEAT_VAPORIZATION * (slope + gamma * (1 + surface_resistance / aerodynamic_resistance)))
