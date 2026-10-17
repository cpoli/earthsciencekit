"""Inspectable plate tectonics models; SI units unless stated otherwise."""

from earthsciencekit.plate_tectonics.systems.boundaries import (
    convergence_rate,
    fault_slip_rate,
    full_spreading_rate,
    slab_pull_force,
    spreading_age,
    subduction_travel_time,
)
from earthsciencekit.plate_tectonics.systems.cooling import (
    diffuse_temperature,
    halfspace_heat_flux,
    halfspace_temperature,
    mantle_rayleigh_number,
    ocean_depth,
)
from earthsciencekit.plate_tectonics.systems.isostasy import (
    airy_root,
    flexural_deflection,
    flexural_parameter,
    flexural_rigidity,
)
from earthsciencekit.plate_tectonics.systems.kinematics import (
    angular_speed,
    compose_euler_vectors,
    euler_velocity_x,
    euler_velocity_y,
    euler_velocity_z,
    hotspot_age,
    plate_displacement,
    plate_speed,
)

__all__ = [
    "angular_speed",
    "plate_speed",
    "euler_velocity_x",
    "euler_velocity_y",
    "euler_velocity_z",
    "plate_displacement",
    "spreading_age",
    "full_spreading_rate",
    "convergence_rate",
    "fault_slip_rate",
    "subduction_travel_time",
    "slab_pull_force",
    "airy_root",
    "flexural_rigidity",
    "flexural_parameter",
    "halfspace_temperature",
    "halfspace_heat_flux",
    "mantle_rayleigh_number",
    "flexural_deflection",
    "hotspot_age",
    "ocean_depth",
    "compose_euler_vectors",
    "diffuse_temperature",
]
