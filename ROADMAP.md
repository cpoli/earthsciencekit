# Roadmap

What earthsciencekit does not do yet, in priority order.

earthsciencekit 0.1.0 shipped five domains. Unreleased work since then
has moved every domain to the family layout (item 2), rewritten the
history pages as dated breakthroughs with one example each (item 3), and
added `geophysics` (item 6) and `volcanology` (item 7). The package now
has seven domains, 162 algebraic and nine dynamic models, 100% statement
coverage and an executed gallery (see the table at the end). Those four
items are summarized under [Done](#done); item numbers are kept because
CHANGELOG.md cites them.

## Still to do

| # | Item | State |
|---|---|---|
| 1 | Version control, CI and releases | In progress: repository, workflows, JupyterLite and Colab done; the 0.2.0 tag, Zenodo DOI and conda-forge recipe remain |
| 4 | `hydrology` domain | Not started |
| 5 | `glaciology` domain | Not started |
| 8 | Deeper models in the existing domains | Not started, apart from a simplified TEOS-10 fit |
| 9 | Course material | Not started; tutorials depend on items 4 and 5 |
| 10 | Typing and coverage | Not started |

Each item lists why it matters, the current state, and a sketch of the
approach. The candidate breakthroughs are a draft: verify each date
against its primary source while writing the history page, and when the
match between a breakthrough and its example is unclear, change the
example rather than the breakthrough.

The project conventions apply to all of them:
- the family layout (item 2), with 100% statement coverage and no
  `# pragma: no cover`;
- the `@model` / `quantity` / `require` validation pattern, SI units,
  radians and explicit `units.convert`;
- NumPy-style docstrings with the governing equation, every unit and the
  applicability limit, citing their source;
- a golden reference value plus an independent identity or limit in the
  tests, and the entry in `model_catalog.json`;
- runtime dependencies stay NumPy, SciPy and Matplotlib;
- a history page with at least 15 dated breakthroughs, each with its own
  gallery example, titled after it and linked from no other entry.

A new domain also needs `SUBPACKAGES` in `docs/source/conf.py`, the docs
index toctrees, the root `__init__` import and `__all__`, the
`test_domain_family_layout` parametrize list, a `model_catalog.json` key,
and the shared counts in README, CLAUDE.md, AGENTS.md, `validation.rst`,
CHANGELOG and the table below.

## 1. Version control, CI and releases

**Why.** Without a release the package cannot be installed from PyPI or
cited with a DOI.

**Done.** The package is the public GitHub repository
`cpoli/earthsciencekit` (author `cpoli374 <cpoli374@gmail.com>`), with
the stray `.DS_Store`, `dist/`, `coverage.json` and
`installed-coverage.json` removed and ignored. `ci.yml` runs ruff,
pytest with the coverage gate on Python 3.10 to 3.13, the narrow mypy
check (item 10 widens it) and `build` plus `twine check`; `docs.yml`
builds the HTML and doctests. `release.yml`, adapted from
mathematicskit, checks a `vX.Y.Z` tag against `__version__`,
`CITATION.cff` and a `## [X.Y.Z] - date` CHANGELOG section, then
publishes to PyPI through trusted publishing, creates the GitHub Release
and deploys `gh-pages`. The first release is 0.2.0 with the seven
domains. Every gallery example has a JupyterLite button, and
`notebooks/quickstart.ipynb` has a Colab badge in the README.

**Still to do.**
- One-time setup: a pending trusted publisher on PyPI (owner `cpoli`,
  repository `earthsciencekit`, workflow `release.yml`, environment
  `pypi`), the `pypi` environment and Pages from `gh-pages` on GitHub,
  and the repository enabled on Zenodo. Then push the tag `v0.2.0`.
- Add that release's own **version DOI** (never the concept DOI) to
  `CITATION.cff` and as a README badge, and update both on every
  release.
- Add a conda-forge `recipe.yaml` after the first PyPI release, with the
  package name written inline rather than in `context`.

## 4. `hydrology`: water on and under land

**Why.** None of the domains models land-surface water. Flow through
soil and aquifers, runoff and evaporation make up a course of their own,
and the domain gives the package its first diffusion problem with a
closed-form well solution.

**Now.** Nothing. `atmosphere` has the vapour-pressure and humidity
helpers that evaporation needs.

**Approach.** Darcy flow and hydraulic head; steady unconfined flow;
infiltration (Green-Ampt, Horton); the Theis and Cooper-Jacob well
solutions (`scipy.special.exp1`); open-channel flow (Chézy, Manning);
unit-hydrograph convolution and Muskingum routing; Penman-Monteith
evaporation. Dynamic model: `simulate_infiltration` (Green-Ampt) or
reservoir routing of a storm hydrograph. Candidate breakthroughs:

1. 1674: Perrault's catchment water balance
2. 1775: Chézy's open-channel formula
3. 1802: Dalton's evaporation law
4. 1856: Darcy's law
5. 1863: Dupuit's unconfined flow
6. 1889: Manning's formula
7. 1911: Green and Ampt's infiltration
8. 1931: Richards's equation
9. 1932: Sherman's unit hydrograph
10. 1933: Horton's infiltration capacity
11. 1935: Theis's well solution
12. 1938: McCarthy's Muskingum routing
13. 1945: Horton's laws of stream numbers
14. 1946: Cooper and Jacob's straight-line method
15. 1948: Penman's evaporation
16. 1965: Monteith's canopy resistance
17. 1979: Beven and Kirkby's TOPMODEL

## 5. `glaciology`: ice sheets, glaciers and sea ice

**Why.** Ice links the `climate` sea-level helpers with the isostasy in
`plate_tectonics`, but the package has no ice dynamics. Glen's law and
the perfectly plastic ice sheet are short closed-form results that
students can check against real ice-sheet profiles.

**Now.** `climate.ice_albedo` and `climate.water_equivalent_sea_level`;
`plate_tectonics.airy_root` and `flexural_deflection` for ice loads;
`geophysics.equivalent_water_height` for GRACE ice-mass change. No ice
flow, mass balance or rebound in time.

**Approach.** Glen's flow law and deformation velocity; Nye's plastic
profile and the shallow-ice approximation; Weertman sliding; mass
balance and equilibrium-line altitude; crevasse depth; Stefan sea-ice
growth; postglacial rebound with a viscous half-space relaxation time;
the Dansgaard-Johnsen age-depth model. Dynamic model: a 1D flowline glacier
responding to a step in the equilibrium line, with
`integrators.integrate_ode` (method of lines). Candidate breakthroughs:

1. 1840: Agassiz's ice age
2. 1843: Forbes's viscous glacier flow
3. 1865: Jamieson's glacial isostasy
4. 1891: Stefan's sea-ice growth
5. 1935: Haskell's mantle viscosity from rebound
6. 1951: Nye's perfectly plastic ice sheet
7. 1955: Glen's flow law
8. 1955: Nye's crevasse depth
9. 1957: Weertman's sliding law
10. 1962: Meier's glacier mass balance
11. 1964: Dansgaard's isotope thermometer
12. 1969: Dansgaard and Johnsen's age-depth model
13. 1974: Weertman's marine ice-sheet instability
14. 1976: Mahaffy's shallow-ice model
15. 1987: Barnola's Vostok CO2 record
16. 2007: Schoof's grounding-line flux

## 8. Deeper models in the existing domains

**Why.** `docs/source/limitations.rst` lists what each domain leaves out.
The gaps that come up most in a course are moist convection and the real
seawater equation of state.

**Now.** `oceanography.simplified_seawater_density` (the Roquet et al.
2015 fit, added for item 3) approximates TEOS-10 but is not TEOS-10.
Nothing else on this list exists.

**Approach.**
- `atmosphere`: a moist pseudo-adiabat and lifting condensation level,
  CAPE from a sounding, and moist ascent in `simulate_dry_ascent`'s
  successor.
- `oceanography`: a TEOS-10 density (the 75-term polynomial, without a
  `gsw` dependency) alongside `linear_density` and
  `simplified_seawater_density`.
- `climate`: a two-box (surface and deep ocean) response, and a 1D
  latitudinal energy-balance model.
- `seismology`: ray tracing in a radially layered Earth, with travel-time
  curves.
- `plate_tectonics`: finite rotations and the composition of Euler
  poles.
- `geophysics`: non-zonal (m ≠ 0) spherical harmonics, since
  `zonal_field` covers m = 0 only.
- `volcanology`: bubble growth (only water solubility is modelled).

Where one of these extends a scientific concept, it gets a history entry
and its own example.

## 9. Course material

**Why.** The package is "for learning and teaching", but it has no
exercises, and the only cross-domain material is the README example.

**Approach.** As in mathematicskit item 1:
- an exercise page per domain, five to ten problems built on gallery
  examples, with solutions in `sphinx-design` dropdowns that are executed
  in the docs build;
- cross-domain tutorials, such as "the Earth's energy budget"
  (`climate`, `atmosphere`, `oceanography`), "isostasy everywhere"
  (`plate_tectonics`, `glaciology`, `geophysics`) and "dating the Earth"
  (`volcanology`, `plate_tectonics`, `geophysics`). The first and last
  can be written now; "isostasy everywhere" waits for item 5;
- a quickstart notebook in JupyterLite (item 1).

## 10. Typing and coverage

- **Type-check the whole package**: the mypy job in `ci.yml` checks
  only `results.py` and `integrators`. Extend it to the whole package
  with `check_untyped_defs` and fix the errors.
- **Branch coverage**: statement coverage is at 100%. Turn on
  `branch = true` and close the gap.
- **Locally verify every CI Python version**: `VERIFICATION.md` covers
  3.14 only. CI runs 3.10 to 3.13; add 3.14, as the other kits do.

## Done

Finished items move here with the release they shipped in, as in
CHANGELOG.md.

### Unreleased

- **2. The family layout.** Each domain has four or five `systems`
  modules (for example `plate_tectonics`: `kinematics`, `boundaries`,
  `isostasy`, `cooling`), with `core/`, `visualizers/` and `tests/`. No
  domain has its own result class or plot yet, so `core/base.py` and
  `visualizers/` re-export the shared `TimeSeries`, `plot_series` and
  `animate_series`, and there is no `utils/`. `_validation.py` and
  `results.py` stay at the root, as in engineeringkit and physicskit;
  `constants` and `units` became packages. No import path changed.
  Physics checks for one domain live in `<domain>/tests/test_physics.py`.
  It landed after item 3 rather than before it.
- **3. History pages as dated breakthroughs.** All five original domains
  have 16 dated entries with one example each (17 for `plate_tectonics`)
  and the three shared scripts per domain are gone. New models added
  along the way: five in `plate_tectonics`, seven in `seismology`, five
  plus `simulate_barotropic_vorticity` in `atmosphere`, eight in
  `climate`, and seven plus `simulate_thermohaline_boxes` in
  `oceanography`. Deviations from the draft list: Richardson's 1922
  forecast is shown by the surface pressure tendency, with the Richardson
  number as its own 1920 entry; Revelle-Suess and Keeling are computed
  from `carbon_airborne_fraction` and `carbon_decay`; the 1968 internal
  waves entry became Ekman's dead water (1904); and TEOS-10 is shown with
  the Roquet et al. (2015) fit, leaving the full polynomial to item 8.
- **6. `geophysics`.** 30 algebraic models in `gravity`, `figure`,
  `magnetism`, `dynamo` and `deformation`, the dynamic
  `simulate_rikitake_dynamo`, and 16 dated entries from Gilbert (1600)
  to GRACE (2002). The dipole uses `lpmv` through `zonal_field` (m = 0
  only); Matuyama is shown by `direction_angle`, Larmor by the free decay
  time and magnetic Reynolds number, and the polarity timescale by summed
  `block_magnetic_anomaly`. New units `m/s^2`, `mGal`, `T` and `nT`.
- **7. `volcanology`.** Kept as one domain, since neither volcanology nor
  geochronology reaches 15 breakthroughs alone: 23 algebraic models in
  `geochronology`, `geochemistry`, `magma` and `eruptions`, the dynamic
  `simulate_eruption_column`, and 17 dated entries from Rayleigh (1896)
  to the VEI (1982). The 1978 column-height entry is named after Wilson,
  Sparks, Huang and Watkins. Goldschmidt uses the Blundy-Wood
  lattice-strain model, Urey the Epstein et al. (1953) calibration. New
  units `kyr`, `Ma`, `Ga` and `GPa`; uranium and radiocarbon decay
  constants in `constants`.

### 0.1.0

The initial educational release: five domains, 77 algebraic and five
dynamic models, shared integration, immutable results, explicit units,
labelled plots and animations, MathematicsKit-style Sphinx hubs and 15
executed examples, with a 100% statement-coverage gate.

## Where earthsciencekit stands

| Domain | Models | History entries | Dated breakthroughs | Gallery examples |
|---|---:|---:|---:|---:|
| `plate_tectonics` | 23 | 16 | 16 | 17 |
| `seismology` | 23 | 16 | 16 | 16 |
| `atmosphere` | 22 | 16 | 16 | 16 |
| `climate` | 24 | 16 | 16 | 16 |
| `oceanography` | 24 | 16 | 16 | 16 |
| `geophysics` | 31 | 16 | 16 | 16 |
| `volcanology` | 24 | 17 | 17 | 17 |
| `hydrology` (item 4) | — | — | — | — |
| `glaciology` (item 5) | — | — | — | — |
