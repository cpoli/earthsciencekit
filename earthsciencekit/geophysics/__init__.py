"""Inspectable geophysics models: gravity, the figure of the Earth and geomagnetism; SI units unless stated otherwise."""

from earthsciencekit.geophysics.systems.deformation import (
    equivalent_water_height,
    love_number_h2,
    love_number_k2,
    tidal_gravimetric_factor,
)
from earthsciencekit.geophysics.systems.dynamo import (
    magnetic_decay_time,
    magnetic_reynolds_number,
    rikitake_fixed_point,
    simulate_rikitake_dynamo,
)
from earthsciencekit.geophysics.systems.figure import (
    centrifugal_ratio,
    clairaut_gravity,
    degree_geoid_height,
    flattening_from_j2,
    homogeneous_flattening,
    nodal_precession_rate,
    normal_gravity,
    stokes_function,
)
from earthsciencekit.geophysics.systems.gravity import (
    bouguer_anomaly,
    bouguer_correction,
    earth_mean_density,
    free_air_correction,
    sphere_gravity,
    torsion_balance_gravitational_constant,
)
from earthsciencekit.geophysics.systems.magnetism import (
    block_magnetic_anomaly,
    dipole_field,
    dipole_inclination,
    direction_angle,
    field_elements,
    gauss_dipole_moment,
    paleolatitude,
    virtual_geomagnetic_pole,
    zonal_field,
)

__all__ = [
    "sphere_gravity",
    "bouguer_correction",
    "free_air_correction",
    "bouguer_anomaly",
    "torsion_balance_gravitational_constant",
    "earth_mean_density",
    "centrifugal_ratio",
    "homogeneous_flattening",
    "clairaut_gravity",
    "normal_gravity",
    "flattening_from_j2",
    "nodal_precession_rate",
    "stokes_function",
    "degree_geoid_height",
    "dipole_field",
    "dipole_inclination",
    "paleolatitude",
    "zonal_field",
    "gauss_dipole_moment",
    "field_elements",
    "direction_angle",
    "virtual_geomagnetic_pole",
    "block_magnetic_anomaly",
    "magnetic_reynolds_number",
    "magnetic_decay_time",
    "rikitake_fixed_point",
    "love_number_h2",
    "love_number_k2",
    "tidal_gravimetric_factor",
    "equivalent_water_height",
    "simulate_rikitake_dynamo",
]
