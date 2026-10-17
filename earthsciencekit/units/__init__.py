"""Small explicit unit conversions; raw model arguments carry no unit metadata."""

from earthsciencekit._validation import model, quantity, require
from earthsciencekit.constants import SECONDS_PER_YEAR

# scale and offset convert each unit to its SI base; temperature offsets are exact.
_UNITS = {
    "m": ("length", 1, 0),
    "km": ("length", 1000, 0),
    "cm": ("length", 0.01, 0),
    "s": ("time", 1, 0),
    "day": ("time", 86400, 0),
    "year": ("time", SECONDS_PER_YEAR, 0),
    "kyr": ("time", 1e3 * SECONDS_PER_YEAR, 0),
    "Ma": ("time", 1e6 * SECONDS_PER_YEAR, 0),
    "Ga": ("time", 1e9 * SECONDS_PER_YEAR, 0),
    "Pa": ("pressure", 1, 0),
    "hPa": ("pressure", 100, 0),
    "MPa": ("pressure", 1e6, 0),
    "GPa": ("pressure", 1e9, 0),
    "K": ("temperature", 1, 0),
    "degC": ("temperature", 1, 273.15),
    "m/s": ("speed", 1, 0),
    "cm/year": ("speed", 0.01 / SECONDS_PER_YEAR, 0),
    "m/s^2": ("acceleration", 1, 0),
    "mGal": ("acceleration", 1e-5, 0),
    "T": ("magnetic field", 1, 0),
    "nT": ("magnetic field", 1e-9, 0),
    "rad": ("angle", 1, 0),
    "degree": ("angle", 0.017453292519943295, 0),
    "m^3/s": ("volume flow", 1, 0),
    "Sv": ("volume flow", 1e6, 0),
}


@model
def convert(value, source, target):
    """Convert a finite scalar or array between listed compatible units.

    Parameters
    ----------
    value : float or array_like
        Quantity in source units. Temperature conversions require nonnegative kelvin.
    source, target : str
        Exact case-sensitive unit labels from available_units().

    Returns
    -------
    numpy scalar or ndarray
        Quantity in target units. Unknown units or incompatible dimensions fail.
    """
    value = quantity("value", value)
    require(isinstance(source, str) and isinstance(target, str), "unit labels must be strings")
    require(source in _UNITS and target in _UNITS, "unsupported unit label")
    dimension, scale, offset = _UNITS[source]
    target_dimension, target_scale, target_offset = _UNITS[target]
    require(dimension == target_dimension, "incompatible unit dimensions")
    base = value * scale + offset
    if dimension == "temperature":
        require(base >= 0, "temperature is below absolute zero")
    return (base - target_offset) / target_scale


def available_units():
    """Return a new mapping of supported unit labels to dimension names."""
    return {name: spec[0] for name, spec in _UNITS.items()}
