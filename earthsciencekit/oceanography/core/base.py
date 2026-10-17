"""Result containers for earthsciencekit.oceanography.

Algebraic models return NumPy scalars or arrays. The sampled dynamic models
return the package-wide immutable :class:`~earthsciencekit.results.TimeSeries`,
which is re-exported here as this domain's result type.
"""

from earthsciencekit.results import TimeSeries

__all__ = ["TimeSeries"]
