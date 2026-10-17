"""Labeled figures and lightweight animations for sampled histories."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from earthsciencekit._validation import quantity, require
from earthsciencekit.results import TimeSeries

__all__ = ["plot_series", "animate_series"]


def plot_series(series):
    """Plot a TimeSeries with one state and explicit unit label per axis."""
    require(isinstance(series, TimeSeries), "series must be a TimeSeries")
    figure, axes = plt.subplots(len(series.names), 1, squeeze=False, sharex=True)
    axes = axes[:, 0]
    for i, axis in enumerate(axes):
        axis.plot(series.time, series.values[:, i])
        axis.set_ylabel(f"{series.names[i]} [{series.units[i]}]")
        axis.grid(alpha=0.2)
    axes[-1].set_xlabel("time [s]")
    figure.tight_layout()
    return figure, axes


def animate_series(series, *, interval=40):
    """Reveal a sampled history progressively; return a Matplotlib animation.

    Keep the returned object alive when displaying it, or call ``save`` with an
    appropriate Matplotlib writer. Interval is positive scalar milliseconds.
    """
    interval = quantity("interval", interval, lower=0, strict_lower=True)
    require(interval.ndim == 0, "interval must be a scalar")
    figure, axes = plot_series(series)
    lines = [axis.lines[0] for axis in axes]
    for line in lines:
        line.set_data([], [])

    def update(frame):
        for i, line in enumerate(lines):
            line.set_data(series.time[: frame + 1], series.values[: frame + 1, i])
        return lines

    return FuncAnimation(figure, update, frames=np.arange(series.time.size), interval=float(interval), blit=True)
