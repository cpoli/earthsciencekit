Scope and limitations
=====================

Version 0.3.0 contains plate_tectonics, seismology, atmosphere, climate,
oceanography, geophysics, volcanology, hydrology and glaciology. Its data are deterministic synthetic examples, not bundled
observational products. It provides educational baselines with documented
assumptions, not an integrated Earth-system model.

- Plate kinematics use rigid spherical geometry. There is no plate catalogue,
  reconstruction inversion or full mantle-convection model.
- Seismology has homogeneous and single-layer elastic helpers, a ray parameter
  and a vertical travel-time integral, and a single oscillator. There is no ray
  tracer, 3D wave propagation, rupture inversion or earthquake forecasting.
- Atmospheric ascent is dry and idealized. There is no moisture-adjusted parcel
  solver or weather assimilation. The barotropic vorticity model is a periodic
  beta-plane teaching model, not numerical weather prediction.
- Climate uses grey radiation and linear one-box response. The carbon and
  sea-level helpers are deliberately simplified budget examples.
- Ocean density is locally linear or the simplified Roquet et al. (2015) fit,
  not the full TEOS-10 polynomial. Currents omit full primitive-equation
  dynamics and site-specific boundary conditions; the gyre is Stommel's linear
  rectangular-basin solution and the thermohaline model has two boxes.
- Geophysics uses first-order, spherical or planar reductions without terrain
  corrections, zonal (axisymmetric) harmonics only, a geocentric axial dipole for
  paleomagnetic directions, two-dimensional vertically magnetized blocks and the
  Love numbers of a uniform elastic sphere. There is no gravity or magnetic
  inversion, full IGRF evaluation or geoid computation from data, and the
  Rikitake dynamo is a chaotic analogue, not a model of the core.
- Volcanology dating assumes closed systems with an unweighted isochron fit and
  uncalibrated radiocarbon ages. Geochemistry uses constant fractionation and
  partition coefficients. Magma is a crystal- and bubble-free Newtonian melt in
  laminar dike flow, the Mogi source is a point in an elastic half-space, and the
  eruption column is a windless top-hat plume without a gas-thrust region,
  particle fallout or moisture. There is no thermodynamic melting model,
  conduit-flow solver, tephra dispersal or hazard assessment.
- Hydrology uses homogeneous, isotropic media and lumped or one-dimensional
  forms: Dupuit flow without recharge, Theis wells in an infinite confined
  aquifer, steady Richards profiles in a Gardner soil, Green-Ampt and Horton
  infiltration in a uniform soil, steady uniform channel flow, a linear
  time-invariant unit hydrograph, Muskingum routing and a big-leaf
  Penman-Monteith canopy in neutral stability. There is no transient Richards
  solver, groundwater model, kinematic or dynamic wave routing, snowmelt,
  catchment calibration or flood-frequency analysis.
- Glaciology uses isothermal Glen ice with one rate factor, shallow-ice and
  plastic profiles, a power-law sliding law, a flowline glacier of uniform
  width with no sliding or calving, local isostasy with a single viscous
  half-space relaxation, Stefan growth without snow or ocean heat, and steady
  one-dimensional ice-core flow. There is no thermomechanical coupling,
  higher-order or full-Stokes flow, ice-shelf buttressing, firn
  densification model, positive-degree-day melt or three-dimensional ice
  sheet.

Physical constants and sample coefficients are documented references, not a
service that automatically tracks revised measurements. Raw floats do not carry
unit metadata. Algebraic broadcasting can allocate large arrays; callers manage
memory. The shared tolerance wrapper uses scalar tolerances. Python 3.12 is the
locally verified runtime; the CI matrix is configured for 3.10-3.13.
