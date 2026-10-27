# Roadmap

What earthsciencekit does not do yet, in priority order.

earthsciencekit 0.1.0 shipped five domains. 0.2.0 (2026-10-17) moved them
to the family layout (item 2), rewrote the history pages as dated
breakthroughs with one example each (item 3), and added `geophysics`
(item 6) and `volcanology` (item 7). 0.3.0 (tagged 2026-10-18, on PyPI
since 2026-10-19) added `hydrology` (item 4) and `glaciology` (item 5).
Item 8 then deepened all nine domains, and item 11 added `geomorphology`.
The package now has ten domains, 305 algebraic and 18 dynamic models, 100%
statement coverage and an executed gallery of 178 examples (see the table at the end). Finished items are summarized under
[Done](#done); item numbers are kept because CHANGELOG.md cites them, and
new items take the next free number.

## Still to do

| # | Item | State |
|---|---|---|
| 1 | Version control, CI and releases | 0.3.0 released with its DOI recorded, and later DOIs recorded by `release.yml`; conda-forge recipe awaits review |
| 9 | Course material | Not started |
| 12 | On the way to 1.0 | Not started |

Each item lists why it matters, the current state, and a sketch of the
approach. [Tasks](#tasks) breaks them into pieces that can each be given
to one prompt. The candidate breakthroughs are a draft: verify each date
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

A new model also changes the shared counts in README, CLAUDE.md,
AGENTS.md, `validation.rst`, `VERIFICATION.md`, CHANGELOG `[Unreleased]`
and the table at the end, and usually a sentence in `limitations.rst`.
Several prompts may run at once, so re-read those lines just before
editing them, edit by exact replacement, and add your delta to the
current totals rather than rewriting from an older copy.

A new domain also needs `SUBPACKAGES` in `docs/source/conf.py`, the four
docs index toctrees (api, history, examples, subpackages), the root
`__init__` import and `__all__`, the `test_domain_family_layout`
parametrize list and a `model_catalog.json` key.

## 1. Version control, CI and releases

**Why.** Without a release the package cannot be installed or cited with
a DOI; without automation every release leaves hand-work behind.

**Done.** The public repository is `cpoli/earthsciencekit` (author
`cpoli374 <cpoli374@gmail.com>`). `ci.yml` runs ruff, pytest with the
coverage gate on Python 3.10 to 3.13, a narrow mypy check (item 10) and
`build` plus `twine check`; `docs.yml` builds the HTML and doctests.
`release.yml`, adapted from mathematicskit, checks a `vX.Y.Z` tag against
`__version__`, `CITATION.cff` and a `## [X.Y.Z] - date` CHANGELOG section,
publishes to PyPI by trusted publishing, creates the GitHub Release and
deploys `gh-pages` (https://cpoli.github.io/earthsciencekit/); its
`zenodo-doi` job then polls Zenodo for the version DOI and commits it to
`CITATION.cff` and the README badge on `main`. Dependabot updates the
pinned actions weekly, and `pyproject.toml` has PyPI `classifiers` and
`keywords`. Every gallery example has a JupyterLite button, and
`notebooks/quickstart.ipynb` has a Colab badge. Zenodo concept record
23441033 holds 0.2.0 as 10.5281/zenodo.23441034 and 0.3.0 as
10.5281/zenodo.23465493, which `CITATION.cff` and the README badge cite.
The conda-forge recipe is conda-forge/staged-recipes#35239 (branch
`earthsciencekit` of `cpoli/staged-recipes`), with the name written
inline; the maintainer confirmation was posted and its checklist box
ticked, and on 2026-10-20 the recipe moved to 0.3.0 (commit c771edd),
with nine domains in the recipe and PR descriptions. No reviewer has
responded yet.

**Still to do.**
- Merge of conda-forge/staged-recipes#35239: answer review comments once
  a reviewer responds, then
  add a conda-forge badge and `conda install` line to the README as in
  mathematicskit.

## 9. Course material

**Why.** The package is "for learning and teaching", but it has no
exercises, and the only cross-domain material is the README example.

**Approach.** As in mathematicskit item 1:
- an exercise page per domain (`docs/source/exercises/<domain>.rst`),
  five to ten problems built on gallery examples (change a parameter,
  predict, check), with solutions in `sphinx-design` dropdowns that are
  executed in the docs build so an API change cannot break them
  silently;
- cross-domain tutorials (`docs/source/tutorials/`): "the Earth's energy
  budget" (`climate`, `atmosphere`, `oceanography`), "isostasy
  everywhere" (`plate_tectonics`, `glaciology`, `geophysics`) and "dating
  the Earth" (`volcanology`, `plate_tectonics`, `geophysics`), and the
  water cycle (`hydrology`, `atmosphere`, `glaciology`). All four can be
  written now;
- the quickstart notebook in JupyterLite: it is not in
  `docs/source/jupyterlite_contents/`, so only Colab opens it.

## 12. On the way to 1.0

**Why.** 1.0 freezes the public API, so the review and the policy come
before it. `CONTRIBUTING.md` says nothing yet about stability or
deprecation.

- **API review**: one pass over every domain's `__all__` for consistent
  names and argument order, with deprecations for whatever changes. For
  example `time` is the last positional argument of
  `simulate_dry_ascent` and `simulate_infiltration` but sits in the
  middle of `simulate_one_box` and `simulate_oscillator`, and
  `simulate_rikitake_dynamo` takes a bare `k`.
- **Stability policy**: a section in `CONTRIBUTING.md` and the docs on
  what 1.0 freezes, and a `DeprecationWarning` helper in `_validation.py`
  with one minor release of notice.
- **Cross-domain consistency of constants**: one source for shared values
  in `constants`, with a test that no domain hard-codes them. For example
  `atmosphere.saturation_vapor_pressure` writes `2.5e6` inline although
  `LATENT_HEAT_VAPORIZATION` exists.

## Scope decisions

What earthsciencekit deliberately does not do, so prompts do not drift
into it:
- no bundled observational datasets, downloads or network access; data
  are synthetic or tabulated constants with a citation;
- no wrapper dependencies (`gsw`, `pyproj`, `obspy`, `xarray`); a
  standard such as TEOS-10 is reimplemented and checked against its
  published values;
- no gridded 2D or 3D Earth-system models; dynamic models stay 0D or 1D
  and sampled through the shared integrators;
- no unit-carrying arrays; floats in SI with explicit `units.convert`.

## Tasks

Each task fits one prompt. Start the prompt with "In earthsciencekit, do
ROADMAP task T*n*" and select this line. All tasks follow the project
conventions above, add a CHANGELOG `[Unreleased]` line, and finish by
moving their bullet to [Done](#done) and updating the item's state in the
table. Unless "Needs" says otherwise, a task does not depend on another.

**Releases (item 1)**
- **T3. conda-forge follow-up.** Check conda-forge/staged-recipes#35239
  and answer any review. Once merged, add the badge
  and `conda install` line to the README. Repeat until merged.

**Course material (item 9)**
- **T14. Exercise infrastructure and `plate_tectonics` exercises**: the
  `docs/source/exercises/` index and toctree, the dropdown pattern with
  executed solutions, and five to ten problems. Done when the docs build
  with `-W`.
- **T15 to T22. Exercises for one domain each** (`seismology`,
  `atmosphere`, `climate`, `oceanography`, `geophysics`, `volcanology`,
  `hydrology`, `glaciology`): five to ten problems in the T14 pattern.
  Needs: T14.
- **T23. Tutorial "The Earth's energy budget"**, creating
  `docs/source/tutorials/` and its toctree if absent.
- **T24. Tutorial "Isostasy everywhere"**, same.
- **T25. Tutorial "Dating the Earth"**, same.
- **T26. Tutorial "The water cycle"**, same.
- **T27. Quickstart in JupyterLite**: add the notebook to the JupyterLite
  contents and a button next to the Colab badge.

**1.0 (item 12)**, in order
- **T35. API review report**: list inconsistent names across `__all__`
  with a proposed rename for each, without changing code.
- **T36. Stability policy and deprecation helper**, then apply the
  renames T35 agreed, with deprecated aliases. Needs: T35.
- **T37. Shared constants**: move duplicated physical values into
  `constants` and add the no-redefinition test.

## Done

Finished items move here with the release they shipped in, as in
CHANGELOG.md.

### Unreleased

- **11. `geomorphology`, T34.** 36 algebraic models in `sediment_transport`,
  `hillslopes`, `rivers`, `aeolian` and `dating`, the dynamic
  `simulate_river_profile` (detachment-limited stream power by the method of
  lines with BDF), the constant `LAMBDA_BE10`, and 17 dated entries from
  Playfair (1802) to Lal (1991). Horton's laws (1945) were left out, since
  they are already a `hydrology` entry with their own example; every other
  drafted date checked out. Playfair is shown by the chi coordinate (Perron
  and Royden 2013), Stokes alongside Ferguson and Church (2004), Gilbert and
  Davis by the river-profile simulation (Davis with slope exponents 1 and 2
  and `plate_tectonics.airy_root`), Shields by the Soulsby-Whitehouse fit,
  Culling by Hanks et al.'s (1984) scarp and morphologic age. Hjulström's
  diagram is rebuilt without cohesion, so its rise for mud is missing.
  `hydrology.hydraulic_radius` and `manning_velocity` feed the bedload and
  hydraulic-geometry examples. Hack's coefficient is 0.32 in SI.
- **8. Deeper models, T5 to T13.** 32 algebraic and six dynamic models,
  each with a catalog entry, golden and identity tests, a doctest and a
  dated history entry with its own example (15 entries in all), and the
  matching `limitations.rst` sentences rewritten. `atmosphere` (T5):
  Bolton's (1980) condensation level rather than Romps's exact form, the
  AMS pseudo-adiabat, parcel profiles, CAPE, CIN and
  `simulate_moist_ascent`, as Espy (1841), Hertz (1884) and Moncrieff and
  Miller (1976). `oceanography` (T6): `teos10_density`, matching the GSW-C
  check values, compared with the two fits in the TEOS-10 example rather
  than a new entry. `climate` (T7): Held et al.'s (2010) two-box response
  and North's (1975) diffusive energy balance with its ice-albedo
  hysteresis. `seismology` (T8): layered-sphere rays as Jeffreys and
  Bullen (1940), whose example shows the core shadow zone.
  `plate_tectonics` (T9): finite rotations as Bullard's fit (1965).
  `geophysics` (T10): Schmidt harmonics and `geomagnetic_field` as the IGRF
  (1968), with `zonal_field` as its m = 0 case. `volcanology` (T11):
  Scriven (1959) coupled to the water solubility. `hydrology` (T12): Celia
  et al.'s (1990) mixed-form Richards column, checked against the steady
  Gardner profiles, Lighthill and Whitham's (1955) kinematic waves and
  Finsterwalder and Schunk's (1887) degree-day snowpack. `glaciology`
  (T13): Herron and Langway (1980), whose close-off depth and mean density
  now compute Barnola's Vostok gas-age difference, and Reeh's (1991, volume
  dated 1989) positive-degree-day melt with the Calov-Greve closed form.
- **T32 (item 10).** The 32 models above have catalog-generated doctests;
  Sphinx doctest covers every algebraic model (572 tests, zero failures).
  Item 10 is done.
- **T1, T2, T4 (item 1).** `CITATION.cff` has `doi:
  10.5281/zenodo.23465493` and the README badge points at it. `release.yml`
  has tbkit's `zenodo-doi` job (concept record 23441033, pushing to
  `main`); its two steps, run locally against v0.3.0, found and wrote
  23465493. CLAUDE.md and AGENTS.md no longer list the manual step.
  `.github/dependabot.yml` updates GitHub Actions weekly in one grouped PR,
  and `pyproject.toml` has 13 keywords and 17 classifiers, all valid Trove
  names; `twine check --strict` passes.
- **T28 to T31, T33 and most of T32 (item 10).** mypy with
  `check-untyped-defs` passes on the whole package; the 140 errors were
  all in tests (mutated `args` lists, one shadowed loop variable, one
  `Optional` module). `branch = true` needed no new tests: branch coverage
  was already 100%. CI tests Python 3.14 and runs a Python 3.10 job at
  NumPy 1.24.0, SciPy 1.10.0 and Matplotlib 3.6.0 (pyparsing < 3.1, since
  newer releases warn inside Matplotlib 3.6). No floor needed raising:
  only tests and gallery examples used `numpy.trapezoid` (NumPy 2.0), and
  they now use `scipy.integrate.trapezoid`. `tests/test_catalog.py`
  checks the catalog's names per domain, parameter names and `expected`
  values; the only mismatch was a zero component that needed the golden
  tests' `atol`. 227 algebraic models have a generated `Examples` doctest
  from their catalog reference; Sphinx doctest runs 489 tests with zero
  failures. `linkcheck.yml` runs weekly; it found two permanent redirects,
  now updated, and volcano.si.edu, which answers bots with 403 and is
  ignored.

### 0.3.0 (2026-10-18)

- **5. `glaciology`.** 37 algebraic models in `rheology`, `ice_sheets`,
  `glaciers`, `sea_ice`, `rebound` and `ice_cores`, the dynamic
  `simulate_flowline_glacier` (shallow ice under a linear mass balance, by
  the method of lines with BDF), and 16 dated entries from Agassiz (1840) to
  Schoof (2007), all as drafted. Agassiz is shown by the equilibrium-line
  shift; Glen's entry carries Nye's (1952) laminar-flow profile; Weertman
  (1974) by the flotation thickness and Weertman's (1957) ice-shelf
  spreading; Mahaffy by Vialov's (1958) steady profile and the flowline
  glacier; Dansgaard-Johnsen with Nye's (1963) age for comparison; Barnola by
  the gas-age difference with `climate.co2_forcing`. New constants
  `ICE_DENSITY` and `LATENT_HEAT_FUSION`.

- **4. `hydrology`.** 38 algebraic models in `catchments`, `groundwater`,
  `soil_water`, `channels` and `evaporation`, the dynamic
  `simulate_infiltration` (Green-Ampt under steady rain, after Mein and
  Larson 1973), and 17 dated entries from Perrault (1674) to TOPMODEL
  (1979), all as drafted. Richards is shown by steady profiles in a Gardner
  (1958) soil rather than a transient solver, Sherman by convolution with a
  Nash (1957) cascade, and Monteith with the FAO-56 grass resistances.
  Evaporation reuses `atmosphere.saturation_vapor_pressure` and its 2.5e6
  J/kg latent heat; new constants `LATENT_HEAT_VAPORIZATION` and
  `VON_KARMAN`. Manning's paper was read in 1889 and printed in 1891.

### 0.2.0 (2026-10-17)

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
| `plate_tectonics` | 27 | 17 | 17 | 18 |
| `seismology` | 25 | 17 | 17 | 17 |
| `atmosphere` | 30 | 19 | 19 | 19 |
| `climate` | 30 | 18 | 18 | 18 |
| `oceanography` | 25 | 16 | 16 | 16 |
| `geophysics` | 33 | 17 | 17 | 17 |
| `volcanology` | 28 | 18 | 18 | 18 |
| `hydrology` | 45 | 20 | 20 | 20 |
| `glaciology` | 43 | 18 | 18 | 18 |
| `geomorphology` | 37 | 17 | 17 | 17 |
