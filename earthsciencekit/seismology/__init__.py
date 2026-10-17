"""Inspectable seismology models; SI units unless stated otherwise."""

from earthsciencekit.seismology.systems.interfaces import (
    reflection_coefficient,
    seismic_impedance,
    transmission_coefficient,
)
from earthsciencekit.seismology.systems.rays import (
    critical_angle,
    crossover_distance,
    head_wave_travel_time,
    ray_parameter,
    travel_time,
    vertical_travel_time,
)
from earthsciencekit.seismology.systems.sources import (
    gutenberg_richter,
    local_magnitude,
    moment_from_magnitude,
    moment_magnitude,
    omori_rate,
    seismic_moment,
    simulate_oscillator,
)
from earthsciencekit.seismology.systems.waves import (
    attenuated_amplitude,
    love_wave_angular_frequency,
    p_wave_speed,
    poisson_ratio,
    rayleigh_wave_speed,
    ricker_wavelet,
    s_wave_speed,
)

__all__ = [
    "p_wave_speed",
    "s_wave_speed",
    "poisson_ratio",
    "travel_time",
    "seismic_impedance",
    "reflection_coefficient",
    "transmission_coefficient",
    "critical_angle",
    "seismic_moment",
    "moment_magnitude",
    "moment_from_magnitude",
    "gutenberg_richter",
    "omori_rate",
    "ricker_wavelet",
    "attenuated_amplitude",
    "rayleigh_wave_speed",
    "head_wave_travel_time",
    "crossover_distance",
    "love_wave_angular_frequency",
    "ray_parameter",
    "local_magnitude",
    "vertical_travel_time",
    "simulate_oscillator",
]
