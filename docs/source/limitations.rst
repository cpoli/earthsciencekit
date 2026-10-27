Scope and limitations
=====================

Version 0.3.0 contains plate_tectonics, seismology, atmosphere, climate,
oceanography, geophysics, volcanology, hydrology and glaciology; ``geomorphology`` is unreleased. Its data are deterministic synthetic examples, not bundled
observational products. It provides educational baselines with documented
assumptions, not an integrated Earth-system model.

- Plate kinematics use rigid spherical geometry; finite rotations are composed
  and interpolated at a constant stage rate. There is no plate catalogue,
  rotation-model data, fitting of reconstructions, uncertainty propagation or
  full mantle-convection model.
- Seismology has homogeneous and single-layer elastic helpers, a ray parameter,
  a vertical travel-time integral, ray tracing for surface sources in a
  radially layered sphere, and a single oscillator. There are no source-depth
  corrections, amplitudes, diffracted or reflected phases, 3D wave
  propagation, rupture inversion or earthquake forecasting.
- Parcel ascent is dry or pseudo-adiabatic, with Bolton's condensation level,
  no entrainment, ice or condensate loading, and saturation over liquid water
  only. There is no weather assimilation. The barotropic vorticity model is a periodic
  beta-plane teaching model, not numerical weather prediction.
- Climate uses grey radiation, linear one-box and two-box responses, and a
  zonally averaged, annual-mean diffusive energy-balance model without
  seasons, land or ocean dynamics. The carbon and sea-level helpers are
  deliberately simplified budget examples.
- Ocean density is locally linear, the simplified Roquet et al. (2015) fit or
  the 75-term TEOS-10 polynomial, which holds only in the oceanographic range
  and has no other TEOS-10 properties. Currents omit full primitive-equation
  dynamics and site-specific boundary conditions; the gyre is Stommel's linear
  rectangular-basin solution and the thermohaline model has two boxes.
- Geophysics uses first-order, spherical or planar reductions without terrain
  corrections, internal spherical-harmonic fields on a spherical Earth from
  user-supplied Gauss coefficients, a geocentric axial dipole for
  paleomagnetic directions, two-dimensional vertically magnetized blocks and the
  Love numbers of a uniform elastic sphere. There is no gravity or magnetic
  inversion, bundled IGRF coefficient set, secular variation, external field,
  geodetic coordinate conversion or geoid computation from data, and the
  Rikitake dynamo is a chaotic analogue, not a model of the core.
- Volcanology dating assumes closed systems with an unweighted isochron fit and
  uncalibrated radiocarbon ages. Geochemistry uses constant fractionation and
  partition coefficients. Magma viscosity and dike flow are for a crystal- and
  bubble-free Newtonian melt in laminar flow. Bubbles grow by diffusion alone at
  constant pressure and diffusivity, without surface tension, viscous
  resistance, nucleation or coalescence. The Mogi source is a point in an
  elastic half-space, and the eruption column is a windless top-hat plume without a gas-thrust region,
  particle fallout or moisture. There is no thermodynamic melting model,
  conduit-flow solver, tephra dispersal or hazard assessment.
- Hydrology uses homogeneous, isotropic media and lumped or one-dimensional
  forms: Dupuit flow without recharge, Theis wells in an infinite confined
  aquifer, steady Richards profiles in a Gardner soil, Green-Ampt and Horton
  infiltration in a uniform soil, steady uniform channel flow, a linear
  time-invariant unit hydrograph, Muskingum routing and a big-leaf
  Penman-Monteith canopy in neutral stability. The transient Richards column is
  one-dimensional in a homogeneous exponential soil above a fixed water table,
  kinematic-wave routing uses a wide channel without backwater, and snowmelt
  is a single-layer degree-day pack. There is no groundwater model, dynamic
  wave routing, energy-balance snowmelt, catchment calibration or
  flood-frequency analysis.
- Glaciology uses isothermal Glen ice with one rate factor, shallow-ice and
  plastic profiles, a power-law sliding law, a flowline glacier of uniform
  width with no sliding or calving, local isostasy with a single viscous
  half-space relaxation, Stefan growth without snow or ocean heat, and steady
  one-dimensional ice-core flow. Firn densifies by the empirical, steady-state
  Herron-Langway model, and melt is a one-layer positive-degree-day scheme
  without rain, superimposed ice or an energy balance. There is no
  thermomechanical coupling, higher-order or full-Stokes flow, ice-shelf
  buttressing, transient firn model, gas diffusion in the firn or
  three-dimensional ice sheet.
- Geomorphology uses uniform cohesionless grains on a flat bed in steady
  uniform flow, linear and Roering-type nonlinear creep, scarps degraded by
  linear diffusion from a single rupture, steady spallation-only cosmogenic
  production, and a one-dimensional detachment-limited stream-power river
  with a fixed outlet. There is no cohesion of mud, mixed grain sizes,
  bedforms, sediment-flux-dependent incision, erosion threshold, hillslope
  and channel coupling, landsliding, flexural isostasy or two-dimensional
  landscape evolution model.

Physical constants and sample coefficients are documented references, not a
service that automatically tracks revised measurements. Raw floats do not carry
unit metadata. Algebraic broadcasting can allocate large arrays; callers manage
memory. The shared tolerance wrapper uses scalar tolerances. Python 3.14 with
current NumPy, SciPy and Matplotlib, and Python 3.10 with the minimum supported
versions (NumPy 1.24, SciPy 1.10, Matplotlib 3.6), are the locally verified
runtimes; CI runs Python 3.10-3.14 and the same minimum-version job.
