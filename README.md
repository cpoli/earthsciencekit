# EarthScienceKit

Inspectable Earth-system computations for learning and teaching.

[![PyPI](https://img.shields.io/pypi/v/earthsciencekit)](https://pypi.org/project/earthsciencekit/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23465493.svg)](https://doi.org/10.5281/zenodo.23465493)
[![CI](https://github.com/cpoli/earthsciencekit/actions/workflows/ci.yml/badge.svg)](https://github.com/cpoli/earthsciencekit/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-cpoli.github.io-blue)](https://cpoli.github.io/earthsciencekit/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cpoli/earthsciencekit/blob/main/notebooks/quickstart.ipynb)
[![JupyterLite](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://cpoli.github.io/earthsciencekit/lite/lab/)

EarthScienceKit follows the MathematicsKit family structure with ten
subpackages and 323 domain models: 305 algebraic equations and 18 sampled
dynamic models. Release **0.3.0** has the first nine; `geomorphology` is unreleased. Every model documents its units, input bounds and applicability.
NumPy and SciPy supply mature numerical foundations.

| Subpackage | Models | Focus |
| --- | ---: | --- |
| `plate_tectonics` | 27 | Euler rotation, finite rotations and their composition, spreading, subduction, isostasy, flexure, convection, hotspots, cooling |
| `seismology` | 25 | Elastic, surface and head waves, rays and travel times in a layered Earth, interfaces, magnitudes, pulses, rebound |
| `atmosphere` | 30 | Hydrostatics, dry and moist thermodynamics, condensation level, CAPE and CIN, humidity, stratification, rotation, Rossby waves, barotropic vorticity |
| `climate` | 30 | Radiation, greenhouse layers, insolation, forcing, feedback, two-box response, diffusive energy balance, ice-albedo hysteresis, red noise, heat storage, carbon and sea level |
| `oceanography` | 25 | Tides, trochoidal and dispersive waves, density, the 75-term TEOS-10 polynomial and cabbeling, dynamic method, Ekman spiral, gyres, two-box overturning, abyssal mixing |
| `geophysics` | 33 | Gravity anomalies, normal gravity and flattening, satellite J2, Stokes's geoid, dipole and Gauss coefficients, Schmidt harmonics and field synthesis, paleomagnetic poles, reversal anomalies, Love numbers, GRACE water, Rikitake dynamo |
| `volcanology` | 28 | Decay and radiometric ages, radiocarbon, Pb-Pb, concordia and isochron dating, Rayleigh and crystal fractionation, lattice strain, isotope thermometry, magma viscosity and ascent, diffusive bubble growth, Mogi inflation, plume and column heights, VEI |
| `hydrology` | 45 | Water balance, Darcy and Dupuit flow, Theis and Cooper-Jacob wells, steady and transient Richards, Green-Ampt and Horton infiltration, Chézy and Manning channels, unit hydrographs, Muskingum and kinematic-wave routing, degree-day snowmelt, Penman-Monteith evaporation, Horton's stream laws, TOPMODEL |
| `glaciology` | 43 | Glen's flow law, Forbes's viscous slab, Nye's plastic ice sheet and crevasses, Weertman sliding, shallow-ice and Vialov profiles, flowline glacier, mass balance, Stefan sea ice, glacial isostasy and Haskell rebound, ice-core dating, Herron-Langway firn and gas-age offset, positive-degree-day melt, Schoof grounding line |
| `geomorphology` | 37 | Stokes and Ferguson-Church settling, the Shields threshold, Rouse suspension, Meyer-Peter and Müller bedload, hillslope creep and fault scarps, soil production, stream power, graded and transient river profiles, chi, Hack, Flint and Strahler, hydraulic geometry, effective discharge, blown sand and dunes, cosmogenic erosion rates |

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
MathematicsKit-style cards/hubs and cross-links. Ten history pages connect
at least fifteen dated breakthroughs each to the APIs; ten guides explain equations,
conventions and independent checks. One executed gallery script per breakthrough
provides figures and Python/notebook downloads.

```bash
MPLBACKEND=Agg python -m pytest --cov=earthsciencekit --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
cd docs && make html
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
