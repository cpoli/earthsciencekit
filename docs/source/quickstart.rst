Quickstart
==========

Install from PyPI or an editable checkout::

    python -m pip install earthsciencekit
    # From a source checkout:
    python -m pip install -e '.[dev]'

The import name is ``earthsciencekit``. Version 0.3.0 contains nine domain
subpackages, and ``geomorphology`` is unreleased; shared utilities support them.

Plate motion with explicit units
--------------------------------

.. doctest::

    >>> from earthsciencekit.units import convert
    >>> from earthsciencekit.plate_tectonics import plate_displacement
    >>> speed = convert(5, "cm/year", "m/s")
    >>> elapsed = convert(1e6, "year", "s")
    >>> round(float(plate_displacement(speed, elapsed)))
    50000

Seismic magnitude convention
----------------------------

.. doctest::

    >>> from earthsciencekit.seismology import moment_magnitude, moment_from_magnitude
    >>> round(float(moment_magnitude(moment_from_magnitude(6))), 10)
    6.0

A sampled climate history
-------------------------

.. code-block:: python

    import numpy as np
    from earthsciencekit.climate import simulate_one_box
    from earthsciencekit.units import convert

    time = convert(np.linspace(0, 100, 201), "year", "s")
    result = simulate_one_box(4, 1.2, 2e8, time)
    result.summary()
    figure, axes = result.plot()
    animation = result.animate(interval=40)
    # Optional: animation.save("warming.gif", writer="pillow")

Input errors
------------

.. doctest::

    >>> from earthsciencekit.atmosphere import dry_air_density
    >>> dry_air_density(101325, 0)
    Traceback (most recent call last):
        ...
    ValueError: temperature violates its lower bound

There is no automatic unit inference, data retrieval or operational forecast.
Select a guide before applying a model to observations.

The same steps are in a `quickstart notebook <https://colab.research.google.com/github/cpoli/earthsciencekit/blob/main/notebooks/quickstart.ipynb>`__ to run in Colab.
Every gallery example also has a **JupyterLite** button that runs it in the
browser, with no installation.
