"""Independent identities, conservation checks and failure contracts."""

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from earthsciencekit import __version__, atmosphere, oceanography, seismology
from earthsciencekit._validation import model, quantity, require, time_grid
from earthsciencekit.integrators import integrate_ode
from earthsciencekit.results import TimeSeries
from earthsciencekit.units import available_units, convert
from earthsciencekit.visualizers import animate_series, plot_series


def test_version():
    assert __version__ == "0.2.0"


def test_bounds_and_scalar_array_contract():
    np.testing.assert_array_equal(quantity("x", [0, 1], lower=0, upper=1), [0, 1])
    assert quantity("x", 1).ndim == 0
    for kwargs, value in [
        ({"lower": 0}, -1),
        ({"lower": 0, "strict_lower": True}, 0),
        ({"upper": 1}, 2),
        ({"upper": 1, "strict_upper": True}, 1),
    ]:
        with pytest.raises(ValueError):
            quantity("x", value, **kwargs)
    with pytest.raises(ValueError):
        require([True, False], "invalid")
    require([True, True], "valid")


def test_arithmetic_rejection():
    with pytest.raises(ValueError, match="arithmetic"):
        seismology.moment_from_magnitude(1000)

    @model
    def nonfinite():
        return float("inf")

    with pytest.raises(ValueError, match="nonfinite"):
        nonfinite()


@pytest.mark.parametrize("bad", [[], [1], [0, 0], [1, 0], [-1, 1], [[0, 1]], [0, np.nan]])
def test_time_grid_bad(bad):
    with pytest.raises(ValueError):
        time_grid(bad)


def test_units():
    assert convert(1, "km", "m") == 1000
    assert convert(0, "degC", "K") == 273.15
    assert convert(273.15, "K", "degC") == 0
    assert convert(1, "Sv", "m^3/s") == 1e6
    assert convert(1, "year", "day") == 365.25
    assert convert(1, "Ga", "Ma") == 1000 and convert(1, "kyr", "year") == 1000 and convert(1, "GPa", "MPa") == 1000
    np.testing.assert_allclose(convert(convert([1, 5], "cm/year", "m/s"), "m/s", "cm/year"), [1, 5])
    assert convert(180, "degree", "rad") == np.pi
    assert np.isclose(convert(1, "mGal", "m/s^2"), 1e-5) and np.isclose(convert(5e4, "nT", "T"), 5e-5)
    units = available_units()
    units.pop("K")
    assert "K" in available_units()
    for source, target, value in [("m", "K", 1), ("bad", "m", 1), ([], "m", 1), ("K", "K", -1)]:
        with pytest.raises(ValueError):
            convert(value, source, target)


def test_result_copy_summary_and_plots(tmp_path):
    time = np.array([0.0, 1.0, 2.0])
    values = np.array([[1.0, 2.0], [0.0, 3.0], [-1.0, 4.0]])
    result = TimeSeries(time, values, ("a", "b"), ("m", "K"))
    time[0] = -100
    values[0, 0] = 999
    assert result.time[0] == 0 and result.values[0, 0] == 1
    assert not result.time.flags.writeable and not result.values.flags.writeable
    assert result.summary()["a"] == {"minimum": -1.0, "maximum": 1.0, "final": -1.0, "unit": "m"}
    figure, axes = result.plot()
    assert len(axes) == 2 and axes[0].get_ylabel() == "a [m]"
    plt.close(figure)
    animation = result.animate(interval=15)
    animation.save(tmp_path / "history.gif", writer="pillow")
    assert (tmp_path / "history.gif").stat().st_size > 100
    plt.close(animation._fig)
    with pytest.raises(ValueError):
        plot_series([1])
    with pytest.raises(ValueError):
        animate_series(result, interval=0)
    with pytest.raises(ValueError):
        animate_series(result, interval=[1, 2])


@pytest.mark.parametrize(
    "time, values, names, units",
    [
        ([0], [[1]], ("x",), ("m",)),
        ([0, 0], [[1], [2]], ("x",), ("m",)),
        ([0, 1], [1, 2], ("x",), ("m",)),
        ([0, 1], [[1]], ("x",), ("m",)),
        ([0, 1], [[1], [2]], (), ()),
        ([0, 1], [[1], [2]], ("",), ("m",)),
        ([0, 1], [[1, 2], [3, 4]], ("x", "x"), ("m", "m")),
    ],
)
def test_result_rejects_bad_shape_and_labels(time, values, names, units):
    with pytest.raises(ValueError):
        TimeSeries(time, values, names, units)


def test_ode_reference_and_defaults():
    time = np.linspace(0, 3, 40)
    result = integrate_ode(lambda t, y: -y, [2], time)
    np.testing.assert_allclose(result.values[:, 0], 2 * np.exp(-time), rtol=3e-8)
    assert result.names == ("x0",) and result.units == ("1",)
    reference = integrate_ode(lambda t, y: [y[1], -y[0]], [1, 0], time, names=("position", "velocity"), units=("m", "m/s"), method="DOP853")
    np.testing.assert_allclose(reference.values[:, 0], np.cos(time), atol=2e-8)
    np.testing.assert_allclose(np.sum(reference.values**2, axis=1), 1, atol=4e-8)


@pytest.mark.parametrize(
    "kwargs", [{"rtol": 0}, {"atol": np.nan}, {"atol": [1, 2]}, {"method": "BAD"}, {"names": ["a", "b"]}, {"units": [""]}]
)
def test_ode_invalid_options(kwargs):
    with pytest.raises(ValueError):
        integrate_ode(lambda t, y: y, [1], [0, 1], **kwargs)


def test_ode_failures(monkeypatch):
    for rhs, initial, time in [
        (0, [1], [0, 1]),
        (lambda t, y: y, [[1]], [0, 1]),
        (lambda t, y: y, [1], [0]),
        (lambda t, y: y, [1], [1, 0]),
        (lambda t, y: [1, 2], [1], [0, 1]),
        (lambda t, y: [np.nan], [1], [0, 1]),
    ]:
        with pytest.raises(ValueError):
            integrate_ode(rhs, initial, time)
    monkeypatch.setattr("earthsciencekit.integrators.solve_ivp", lambda *a, **k: SimpleNamespace(success=False, message="test failure"))
    with pytest.raises(ValueError, match="test failure"):
        integrate_ode(lambda t, y: -y, [1], [0, 1])


def test_humidity_and_wave_limits():
    ratio = atmosphere.mixing_ratio(1000, 100000)
    humidity = atmosphere.specific_humidity(ratio)
    assert np.isclose(humidity / (1 - humidity), ratio)
    assert atmosphere.coriolis_parameter(-np.pi / 4) < 0 < atmosphere.coriolis_parameter(np.pi / 4)
    k, h, g = 1e-6, 10.0, 9.81
    omega = oceanography.wave_angular_frequency(k, h, g)
    assert np.isclose(omega / k, oceanography.shallow_wave_speed(h, g), rtol=1e-9)
    k, h = 1.0, 100.0
    omega = oceanography.wave_angular_frequency(k, h, g)
    assert np.isclose(omega / k, oceanography.deep_wave_speed(2 * np.pi / k, g))


def test_dynamic_overflow_is_normalized():
    with pytest.raises(ValueError, match="arithmetic"):
        seismology.simulate_oscillator(1e308, 10, [0, 1])


@pytest.mark.parametrize("name", ["plate_tectonics", "seismology", "atmosphere", "climate", "oceanography", "geophysics", "volcanology"])
def test_domain_family_layout(name):
    import importlib
    import inspect

    domain = importlib.import_module(f"earthsciencekit.{name}")
    core = importlib.import_module(f"earthsciencekit.{name}.core")
    plots = importlib.import_module(f"earthsciencekit.{name}.visualizers")
    assert core.TimeSeries is TimeSeries
    assert plots.plot_series is plot_series and plots.animate_series is animate_series
    for public in domain.__all__:
        function = getattr(domain, public)
        assert inspect.getmodule(function).__name__.startswith(f"earthsciencekit.{name}.systems."), public
