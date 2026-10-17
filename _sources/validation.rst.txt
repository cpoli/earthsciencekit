Validation policy
=================

Validation has three separate goals: reject invalid inputs, verify numerical
implementation and explain physical applicability. None substitutes for another.

Tests include fixed reference values for all 162 algebraic models, elementwise
physical bounds, finite/real contracts and broadcasting. Cross-model tests cover
Euler-vector tangency and circuit closure, the load balance of a flexed plate,
ocean subsidence against integrated heat loss, wave-interface energy conservation, magnitude inversion,
moisture-ratio inversion and shallow/deep wave limits. Climate tests check the
top-of-atmosphere balance of a greenhouse layer, the global mean of daily
insolation, closure of every ice-albedo equilibrium, the variance of the red-noise
spectrum and recovery of one-box parameters by the Gregory regression.
Oceanography tests check that the integrated Ekman spiral equals the Ekman
transport, that the Stommel gyre interior matches Sverdrup transport, that the
abyssal profile solves its advection-diffusion equation, and that the simplified
seawater density reduces to the linear law near its reference and cabbels.
Geophysics tests check that the Bouguer anomaly recovers a buried body under
topography, that Clairaut's formula matches Somigliana's normal gravity, that
the numerical Stokes integral of a zonal anomaly gives the ``1 / (n - 1)`` law,
that the degree-1 zonal field is the axial dipole, that virtual geomagnetic
poles of a tilted dipole return its pole and that an infinite magnetized layer
has no external field.
Volcanology tests check that decay conserves parent plus daughter atoms and the
age equation inverts accumulation, that the Pb-Pb age inverts the radiogenic lead
ratio, that both concordia clocks agree and lead loss falls on a chord, that the
isochron recovers its age, that fractional crystallization conserves mass, that
lattice strain peaks at the optimal radius, that Shaw's lines share their pivot
and that the Mogi uplift integrates to ``2 (1 - nu) dV``.

Dynamic tests compare the seismic oscillator with its damped analytic solution,
climate adjustment with its closed form, dry ascent with conserved potential
temperature, inertial motion with analytic rotation/decay, the thermohaline
boxes with both steady states of their equilibrium curve, the Rikitake dynamo
with its steady states and spontaneous reversals, the eruption column with the
pure-plume similarity solution and the Morton-Taylor-Turner height, and periodic diffusion
with a discrete Fourier eigenmode and conserved mean. Plot labels, result
immutability, summaries and GIF export are exercised.

Coverage includes every production module, including private validation and
plotting utilities. Only test modules are omitted. The test command enforces
100% production statement coverage with no production exclusion pragmas.
Branch coverage is not currently configured. Statement execution is neither a
proof of scientific correctness nor evidence that a model applies to a real site.

Reproduce from the source checkout::

    python -m pip install -e '.[dev]'
    MPLBACKEND=Agg python -m pytest --cov=earthsciencekit --cov-report=term-missing
    python -m ruff check .
    python -m ruff format --check .
    MPLBACKEND=Agg python -m sphinx -E -W --keep-going -b html docs/source docs/_build/html
    python -m sphinx -W -b doctest docs/source docs/_build/doctest
    python -m build

The release's exact test counts, installed-wheel checks and link checks are
recorded in VERIFICATION.md. Examples execute during the HTML build and abort
on an error. Scientific teaching references are listed in :doc:`references`.
