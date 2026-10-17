"""Plots for earthsciencekit.geophysics.

The dynamic models return a TimeSeries, whose ``plot`` and ``animate`` methods
call these package-wide helpers.
"""

from earthsciencekit.visualizers import animate_series, plot_series

__all__ = ["plot_series", "animate_series"]
