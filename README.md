# EarthScienceKit

Inspectable Earth-system computations for learning and teaching.

[![PyPI](https://img.shields.io/pypi/v/earthsciencekit)](https://pypi.org/project/earthsciencekit/)
[![CI](https://github.com/cpoli/earthsciencekit/actions/workflows/ci.yml/badge.svg)](https://github.com/cpoli/earthsciencekit/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-cpoli.github.io-blue)](https://cpoli.github.io/earthsciencekit/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cpoli/earthsciencekit/blob/main/notebooks/quickstart.ipynb)
[![JupyterLite](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://cpoli.github.io/earthsciencekit/lite/lab/)

EarthScienceKit follows the MathematicsKit family structure with nine
subpackages and 248 domain models: 237 algebraic equations and eleven sampled
dynamic models. Release **0.2.0** has the first seven; `hydrology` and `glaciology` are unreleased. Every model documents its units, input bounds and applicability.
NumPy and SciPy supply mature numerical foundations.

| Subpackage | Models | Focus |
| --- | ---: | --- |
| `plate_tectonics` | 23 | Euler rotation, spreading, subduction, isostasy, flexure, convection, hotspots, cooling |
| `seismology` | 23 | Elastic, surface and head waves, rays, interfaces, magnitudes, pulses, rebound |
| `atmosphere` | 22 | Hydrostatics, dry thermodynamics, humidity, stratification, rotation, Rossby waves, barotropic vorticity |
| `climate` | 24 | Radiation, greenhouse layers, insolation, forcing, feedback, ice-albedo hysteresis, red noise, heat storage, carbon and sea level |
| `oceanography` | 24 | Tides, trochoidal and dispersive waves, density and TEOS-10 cabbeling, dynamic method, Ekman spiral, gyres, two-box overturning, abyssal mixing |
| `geophysics` | 31 | Gravity anomalies, normal gravity and flattening, satellite J2, Stokes's geoid, dipole and Gauss coefficients, paleomagnetic poles, reversal anomalies, Love numbers, GRACE water, Rikitake dynamo |
| `volcanology` | 24 | Decay and radiometric ages, radiocarbon, Pb-Pb, concordia and isochron dating, Rayleigh and crystal fractionation, lattice strain, isotope thermometry, magma viscosity and ascent, Mogi inflation, plume and column heights, VEI |
| `hydrology` | 39 | Water balance, Darcy and Dupuit flow, Theis and Cooper-Jacob wells, Richards, Green-Ampt and Horton infiltration, Chézy and Manning channels, unit hydrographs, Muskingum routing, Penman-Monteith evaporation, Horton's stream laws, TOPMODEL |
| `glaciology` | 38 | Glen's flow law, Forbes's viscous slab, Nye's plastic ice sheet and crevasses, Weertman sliding, shallow-ice and Vialov profiles, flowline glacier, mass balance, Stefan sea ice, glacial isostasy and Haskell rebound, ice-core dating, gas-age offset, Schoof grounding line |

## Install and explore

```bash
python -m pip install earthsciencekit
# Or from a checkout:
python -m pip install -e '.[dev]'
```

```python
import numpy as np
from earthsciencekit.climate import simulate_one_box
from earthsciencekit.units import convert

time = convert(np.linspace(0, 100, 201), "year", "s")
result = simulate_one_box(4, 1.2, 2e8, time)
print(result.summary())
figure, axes = result.plot()
animation = result.animate(interval=40)
```

SI units apply unless an argument explicitly states otherwise. Latitude and
angles use radians. Algebraic functions accept scalars and broadcast arrays;
finite real quantities and physical bounds are validated. Unit conversions are
explicit. Sampled histories copy their arrays, retain labels and support plots
and animation export.

## Documentation and validation

The Sphinx site has Subpackages, History, Examples and API entry points,
MathematicsKit-style cards/hubs and cross-links. Nine history pages connect
at least fifteen dated breakthroughs each to the APIs; nine guides explain equations,
conventions and independent checks. One executed gallery script per breakthrough
provides figures and Python/notebook downloads.

```bash
MPLBACKEND=Agg python -m pytest --cov=earthsciencekit --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
MPLBACKEND=Agg python -m sphinx -E -W --keep-going -b html docs/source docs/_build/html
python -m build
```

Coverage enforces 100% of production statements, including validators and
visualization; test files alone are omitted. This metric does not establish
physical applicability. Tests include independent analytic solutions,
conservation balances, limiting cases and invalid-input contracts. Exact release
results are in VERIFICATION.md; model limits are in docs/source/limitations.rst.
The documentation is hosted at https://cpoli.github.io/earthsciencekit/, and
`notebooks/quickstart.ipynb` runs in Colab from the badge above.

## Contribute

See CONTRIBUTING.md. Educational models should have explicit equations, units,
assumptions, independent benchmarks and an executable example. The MIT license
is in LICENSE. Please cite this software using CITATION.cff.
