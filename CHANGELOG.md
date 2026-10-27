# Changelog

## [Unreleased]

- Releases (ROADMAP item 1): `CITATION.cff` and the README badge cite the
  0.3.0 version DOI 10.5281/zenodo.23465493; a `zenodo-doi` job in
  `release.yml` records each later release's version DOI automatically;
  Dependabot keeps the GitHub Actions current; `pyproject.toml` has PyPI
  classifiers and keywords.
- Checks (ROADMAP item 10): mypy covers the whole package; coverage is
  measured on branches as well as statements, still gated at 100%; CI tests
  Python 3.14 and, in a separate job, Python 3.10 at the exact dependency
  floors; `tests/test_catalog.py` ties `model_catalog.json` to each domain's
  algebraic `__all__` and golden values; every algebraic model docstring has
  a doctested `Examples` section, run by `docs.yml`; a weekly `linkcheck.yml`
  checks the docs' external links. Tests and gallery examples use
  `scipy.integrate.trapezoid` so that they run on NumPy 1.24.
- Deeper models (ROADMAP item 8): 32 algebraic and six dynamic models, each
  documented, tested, catalogued and with a doctested `Examples` section, and
  15 new dated history entries with their own examples.
  - `atmosphere`: Bolton's condensation level (`lcl_temperature`,
    `lcl_pressure`, `lcl_height`), the pseudo-adiabatic lapse rate, parcel
    profiles, CAPE and CIN, and `simulate_moist_ascent`.
  - `oceanography`: the 75-term TEOS-10 density (`teos10_density`), checked
    against the published check values.
  - `climate`: the two-box response (`two_box_timescales`,
    `two_box_fast_warming`, `simulate_two_box`) and North's diffusive
    energy-balance model (`insolation_distribution`,
    `north_temperature_profile`, `simulate_latitudinal_energy_balance`).
  - `seismology`: rays in a radially layered sphere (`layered_ray_distance`,
    `layered_ray_travel_time`), with the core shadow zone.
  - `plate_tectonics`: finite rotations, their composition and interpolation.
  - `geophysics`: Schmidt semi-normalized harmonics of any order and field
    synthesis from Gauss coefficients (`geomagnetic_field`).
  - `volcanology`: Scriven's diffusive bubble growth and
    `simulate_bubble_growth`.
  - `hydrology`: a transient Richards column (`simulate_richards_column`),
    kinematic-wave routing and degree-day snowmelt.
  - `glaciology`: Herron-Langway firn densification, which now computes the
    Vostok gas-age difference, and Reeh's positive-degree-day melt
    (`expected_positive_degree_days`, `pdd_surface_melt`).
- New `geomorphology` domain (ROADMAP item 11): Stokes and Ferguson-Church
  settling, bed shear stress and shear velocity, the Shields number,
  dimensionless grain size and Soulsby-Whitehouse threshold, Keulegan's mean
  velocity, the Rouse number and profile, Meyer-Peter and Müller bedload;
  linear and nonlinear hillslope creep, Hanks's fault scarp and its steepest
  gradient, soil production and hilltop curvature; stream-power erosion, the
  steady-state slope, knickpoint celerity, Flint's slope-area law, the chi
  coordinate, Hack's law, Strahler order, downstream hydraulic geometry and
  the effective discharge; Bagnold's threshold and sand flux, Kawamura's flux
  and dune migration; cosmogenic production, accumulation and steady-state
  erosion rates, morphologic and luminescence ages: 36 algebraic models and
  the dynamic `simulate_river_profile`. Its history is a dated chronology of
  17 breakthroughs, from Playfair's accordant junctions (1802) to Lal's
  cosmogenic erosion rates (1991), each with its own executed example.
- New constant `LAMBDA_BE10`.

## [0.3.0] - 2026-10-18

- New `hydrology` domain (ROADMAP item 4): the catchment water balance and
  runoff coefficient, Darcy flow and hydraulic head, Dupuit's unconfined
  aquifer, the Theis well function and drawdown, the Cooper-Jacob
  approximation and transmissivity, Gardner's conductivity and steady Richards
  profiles, Green-Ampt and Horton infiltration, Chézy and Manning channel flow,
  Nash and Sherman unit hydrographs, Muskingum routing, Dalton, Penman and
  Penman-Monteith evaporation with the aerodynamic and surface resistances,
  Horton's stream laws and TOPMODEL's index, deficit and baseflow: 38
  algebraic models and the dynamic `simulate_infiltration` (Green-Ampt under
  steady rain). Its history is a dated chronology of 17 breakthroughs, from
  Perrault's catchment water balance (1674) to TOPMODEL (1979), each with its
  own executed example.
- New constants `LATENT_HEAT_VAPORIZATION` and `VON_KARMAN`.
- New `glaciology` domain (ROADMAP item 5): Forbes's viscous slab, Glen's
  flow law with its rate factor, effective viscosity and laminar velocity
  profile, Nye's plastic ice sheet and dry and water-filled crevasse depths,
  Weertman sliding, the driving stress, the shallow-ice flux and diffusivity
  with Vialov's steady profile, linear mass balance, the accumulation-area
  ratio, equilibrium-line shift and response time, Stefan sea-ice growth and
  freeboard, glacial isostasy and Haskell's viscous rebound, Dansgaard's
  isotope thermometer, the Nye and Dansgaard-Johnsen time scales, the gas-age
  difference, flotation thickness, ice-shelf spreading and Schoof's
  grounding-line flux: 37 algebraic models and the dynamic
  `simulate_flowline_glacier` (a shallow-ice flowline under a linear mass
  balance, by the method of lines). Its history is a dated chronology of 16
  breakthroughs, from Agassiz's ice age (1840) to Schoof's grounding-line flux
  (2007), each with its own executed example.
- New constants `ICE_DENSITY` and `LATENT_HEAT_FUSION`.
- The docs landing figure shows one panel per domain, all nine, and is
  regenerated by `docs/make_readme_figure.py`.

## [0.2.0] - 2026-10-17

- First release on PyPI, with the git history on GitHub, CI, tag-driven
  releases (`release.yml`), docs hosted on GitHub Pages (`docs.yml`), a mypy
  step in CI, JupyterLite buttons on the gallery and a Colab quickstart
  notebook (ROADMAP item 1).
- New `geophysics` domain (ROADMAP item 6): gravity anomalies and reductions,
  normal gravity, Clairaut's theorem and the flattening, satellite J2, Stokes's
  geoid kernel, the axial dipole and zonal Gauss coefficients, field elements,
  virtual geomagnetic poles, magnetic anomalies of reversed blocks, the dynamo
  scales, Love numbers and GRACE equivalent water: 30 algebraic models and the
  dynamic `simulate_rikitake_dynamo`. Its history is a dated chronology of 16
  breakthroughs, from Gilbert's magnetic Earth (1600) to GRACE (2002), each with
  its own executed example.
- New units `m/s^2`, `mGal`, `T` and `nT`.
- New `volcanology` domain (ROADMAP item 7): radioactive decay and the age
  equation, radiocarbon, Pb-Pb ages, the U-Pb concordia and the isochron,
  Rayleigh fractionation and delta values, fractional and equilibrium
  crystallization, lattice-strain partitioning, the oxygen-isotope
  thermometer, Arrhenius and Shaw viscosity, water solubility, dike ascent,
  the Mogi source, plume and eruption-column heights and the VEI: 23
  algebraic models and the dynamic `simulate_eruption_column`. Its history is
  a dated chronology of 17 breakthroughs, from Rayleigh's fractional
  distillation (1896) to the Volcanic Explosivity Index (1982), each with its
  own executed example.
- New units `kyr`, `Ma`, `Ga` and `GPa`; new constants `AVOGADRO`,
  `LAMBDA_U238`, `LAMBDA_U235`, `U238_U235` and `LIBBY_MEAN_LIFE`.
- Family layout: each domain is a package of `systems/<topic>.py` modules with
  `core/`, `visualizers/` and `tests/`; `constants` and `units` are packages and
  the cross-cutting tests moved from `tests/` to `earthsciencekit/tests/`. Every
  public import path, such as `from earthsciencekit.climate import simulate_one_box`,
  is unchanged.
- `plate_tectonics` history is a dated chronology of 16 breakthroughs, from Euler's
  rotation theorem (1776) to NUVEL-1 (1990), each with its own executed example;
  the three shared plate_tectonics gallery scripts are removed.
- New `plate_tectonics` models: `mantle_rayleigh_number`, `flexural_deflection`,
  `hotspot_age`, `ocean_depth` and `compose_euler_vectors`.
- `seismology` history is a dated chronology of 16 breakthroughs, from Poisson's
  two elastic waves (1829) to PREM (1981), each with its own executed example;
  the three shared seismology gallery scripts are removed.
- New `seismology` models: `rayleigh_wave_speed`, `head_wave_travel_time`,
  `crossover_distance`, `love_wave_angular_frequency`, `ray_parameter`,
  `local_magnitude` and `vertical_travel_time`.
- `atmosphere` history is a dated chronology of 16 breakthroughs, from Torricelli's
  barometer (1643) to the ENIAC forecast (1950), each with its own executed example;
  the three shared atmosphere gallery scripts are removed.
- New `atmosphere` models: `barometer_column_height`, `two_layer_temperature`,
  `richardson_number`, `surface_pressure_tendency`, `rossby_wave_speed` and the
  dynamic `simulate_barotropic_vorticity`.
- `climate` history is a dated chronology of 16 breakthroughs, from Fourier's
  greenhouse (1824) to the sea-level budget (2005), each with its own executed
  example; the three shared climate gallery scripts are removed.
- New `climate` models: `greenhouse_surface_temperature`, `daily_insolation`,
  `radiative_convective_temperature`, `equilibrium_solar_constant`,
  `red_noise_spectrum`, `methane_forcing`, `nitrous_oxide_forcing` and
  `gregory_regression`.
- `oceanography` history is a dated chronology of 16 breakthroughs, from Newton's
  equilibrium tide (1687) to TEOS-10 (2010), each with its own executed example;
  the three shared oceanography gallery scripts are removed.
- New `oceanography` models: `gerstner_orbit_radius`, `ekman_spiral_velocity`,
  `sverdrup_transport`, `stommel_streamfunction`, `thermohaline_forcing`,
  `abyssal_temperature`, `simplified_seawater_density` and the dynamic
  `simulate_thermohaline_boxes`.

## [0.1.0] - 2026-10-04

Initial educational release.

- Five domains: plate_tectonics, seismology, atmosphere, climate, oceanography.
- 77 validated algebraic models and five sampled dynamical models.
- Shared integration, immutable results, explicit units, labeled plots and animations.
- MathematicsKit-style Sphinx hubs, histories, guides and 15 executed examples.
- Independent reference, conservation, limit and invalid-input tests; coverage gate.
