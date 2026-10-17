"""Immutable sampled histories shared by the numerical models."""

from dataclasses import dataclass

import numpy as np

from ._validation import quantity, require


@dataclass(frozen=True)
class TimeSeries:
    """A time grid in seconds and one or more state columns.

    Parameters
    ----------
    time : array_like, shape (n,)
        Finite, strictly increasing seconds; at least two samples.
    values : array_like, shape (n, m)
        Finite real state columns in the units given by ``units``.
    names, units : tuple of str, length m
        Unique state names and explicit output unit labels.

    Notes
    -----
    Input data are copied and made read-only. Samples do not imply a continuous
    interpolation or independently validate the governing model.
    """

    time: np.ndarray
    values: np.ndarray
    names: tuple[str, ...]
    units: tuple[str, ...]

    def __post_init__(self):
        time = quantity("time", self.time)
        values = quantity("values", self.values)
        require(time.ndim == 1 and time.size >= 2, "time must be a vector with at least two samples")
        require(np.all(np.diff(time) > 0), "time must strictly increase")
        require(values.ndim == 2 and values.shape[0] == time.size, "values must have shape (time, state)")
        names, units = tuple(self.names), tuple(self.units)
        require(len(names) == values.shape[1] == len(units), "one name and unit per state required")
        require(all(isinstance(s, str) and s for s in names + units), "names and units must be nonempty strings")
        require(len(set(names)) == len(names), "state names must be unique")
        time.setflags(write=False)
        values.setflags(write=False)
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "units", units)

    def summary(self):
        """Return each state's minimum, maximum, final value and unit label."""
        return {
            name: {
                "minimum": float(self.values[:, i].min()),
                "maximum": float(self.values[:, i].max()),
                "final": float(self.values[-1, i]),
                "unit": self.units[i],
            }
            for i, name in enumerate(self.names)
        }

    def plot(self):
        """Create one labeled Matplotlib axis per state; return figure and axes."""
        from .visualizers import plot_series

        return plot_series(self)

    def animate(self, *, interval=40):
        """Create a history-reveal animation; interval is milliseconds per frame."""
        from .visualizers import animate_series

        return animate_series(self, interval=interval)
