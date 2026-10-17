"""Cross-parameter relations that simple scalar bounds cannot establish."""

import pytest

from earthsciencekit import atmosphere as a
from earthsciencekit import climate as c
from earthsciencekit import plate_tectonics as p
from earthsciencekit import seismology as s


@pytest.mark.parametrize(
    "function,args",
    [
        (p.plate_speed, [1, 1, 4]),
        (p.subduction_travel_time, [1, 1, 0]),
        (p.airy_root, [1, 3, 2]),
        (p.halfspace_temperature, [1, 1, 1, 300, 200]),
        (p.halfspace_heat_flux, [1, 1, 1, 300, 200]),
        (p.ocean_depth, [1, 2500, 1000, 1000, 1, 1, 300, 1300]),
        (p.ocean_depth, [1, 2500, 2, 1, 1, 1, 300, 200]),
        (s.poisson_ratio, [1, 1]),
        (s.critical_angle, [2, 1]),
        (s.rayleigh_wave_speed, [1, 1]),
        (s.head_wave_travel_time, [4, 1, 2, 1]),
        (s.head_wave_travel_time, [1, 1, 1, 2]),
        (s.crossover_distance, [1, 2, 2]),
        (s.love_wave_angular_frequency, [1, 1, 1, 2, 1, 1]),
        (s.love_wave_angular_frequency, [2, 1, 1, 2, 1, 1]),
        (a.mixing_ratio, [1000, 100]),
        (a.two_layer_temperature, [0, 288, 0.0065, 50000]),
        (c.ice_albedo, [273, 0.2, 0.6, 273, 5]),
    ],
)
def test_cross_parameter_rejections(function, args):
    with pytest.raises(ValueError):
        function(*args)
