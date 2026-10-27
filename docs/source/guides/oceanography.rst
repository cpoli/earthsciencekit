Oceanography: units, equations and validation
=============================================

.. include:: /_generated/nav/oceanography.rst

Conventions and model limits
----------------------------

Depth is positive downward, but the density gradient for N squared is defined
upward. Hydrostatic pressure is gauge pressure and omits atmospheric pressure.
The linear density law is a local reference-state approximation, not TEOS-10.
Choose anomalies and coefficients that keep density positive and the
linearization small.

Geostrophic and Ekman helpers return signed scalar cross-gradient/cross-wind
components; constructing a vector requires a consistent orientation. Coriolis
parameter is signed, and zero is rejected. Volume transport is m^3/s; divide by
1e6 or convert to Sv for Sverdrups.

Surface-wave formulas assume small amplitude, no current and no capillarity.
Finite-depth dispersion approaches shallow-water speed sqrt(g h) and deep-water
speed sqrt(g/k). Deep-water group speed is half phase speed. The Rossby radius
uses the first-mode constant-N estimate ``N H/(pi abs(f))``, which differs from a
reduced-gravity estimate. A harmonic tide with user parameters is not a local
operational prediction.

Gerstner orbits use the rest depth of the orbit centre, positive downward, and
require ``amplitude * wavenumber <= 1``. The Ekman spiral returns the downwind
and crosswind components on the last axis, with crosswind positive to the right
of the wind; its depth integral equals ``ekman_transport`` with the same sign
convention. The Stommel streamfunction uses ``x`` east from the western wall and
``y`` north from the southern wall, with zonal stress
``-wind_stress * cos(pi y / basin_length)``; positive stress gives a clockwise
gyre and northward transport ``dpsi/dx``. Sverdrup transport is per unit zonal
width (m^2/s). In the two-box thermohaline model, positive overturning means
polar sinking and the salt forcing is a virtual salt flux in g/kg m^3/s. The
abyssal profile measures height upward from the bottom of the column. The
simplified seawater density is the Roquet et al. (2015) S-EOS with NEMO's
default coefficients; it reproduces TEOS-10 gradients, including cabbeling and
thermobaricity, but not TEOS-10 absolute density. ``teos10_density`` is the
75-term polynomial of the GSW toolbox, which matches TEOS-10 itself; it takes
the sea pressure (absolute pressure minus 101325 Pa) in pascals rather than a
depth.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. No unit
metadata are carried in raw numbers. Bounds protect mathematical domains and
obvious physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

hydrostatic_pressure
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.hydrostatic_pressure`

Output: **Pa**. Gauge pressure with constant density; atmospheric pressure is excluded.

.. code-block:: python

    density * gravity * depth

linear_density
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.linear_density`

Output: **kg/m^3**. Local linearization around a supplied reference; not TEOS-10 and no pressure dependence.

.. code-block:: python

    reference_density * (1 - thermal_expansion * temperature_anomaly + haline_contraction * salinity_anomaly)

buoyancy_frequency_squared
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.buoyancy_frequency_squared`

Output: **1/s^2**. Coordinate increases upward; negative output denotes unstable stratification.

.. code-block:: python

    -gravity * upward_density_gradient / density

geostrophic_current
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.geostrophic_current`

Output: **m/s**. Signed cross-slope component away from the equator; no ageostrophic velocity.

.. code-block:: python

    gravity * surface_slope / coriolis

ekman_transport
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.ekman_transport`

Output: **m^2/s**. Signed cross-wind component per unit width; direction requires a consistent vector convention.

.. code-block:: python

    wind_stress / (density * coriolis)

inertial_period
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.inertial_period`

Output: **s**. No pressure gradient or friction; invalid at zero Coriolis parameter.

.. code-block:: python

    2 * np.pi / np.abs(coriolis)

shallow_wave_speed
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.shallow_wave_speed`

Output: **m/s**. Wavelength much greater than depth, small amplitude, no dispersion.

.. code-block:: python

    np.sqrt(gravity * depth)

deep_wave_speed
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.deep_wave_speed`

Output: **m/s**. Deep-water small-amplitude gravity waves without capillarity.

.. code-block:: python

    np.sqrt(gravity * wavelength / (2 * np.pi))

wave_angular_frequency
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.wave_angular_frequency`

Output: **rad/s**. Linear inviscid surface gravity wave, no current or surface tension.

.. code-block:: python

    np.sqrt(gravity * wavenumber * np.tanh(wavenumber * depth))

deep_group_speed
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.deep_group_speed`

Output: **m/s**. Deep-water gravity waves; group speed is half phase speed.

.. code-block:: python

    0.5 * np.sqrt(gravity * wavelength / (2 * np.pi))

rossby_radius
~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.rossby_radius`

Output: **m**. First-mode constant-N rigid-boundary estimate N H over pi f; not a general layered solution.

.. code-block:: python

    buoyancy_frequency * depth / (np.pi * np.abs(coriolis))

reduced_gravity
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.reduced_gravity`

Output: **m/s^2**. Boussinesq small density contrast; output is acceleration.

.. code-block:: python

    gravity * density_difference / reference_density

internal_wave_speed
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.internal_wave_speed`

Output: **m/s**. Hydrostatic two-layer fluid under a rigid lid.

.. code-block:: python

    np.sqrt(reduced_gravity * upper_depth * lower_depth / (upper_depth + lower_depth))

volume_transport
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.volume_transport`

Output: **m^3/s**. Uniform signed velocity and rectangular section; output m cubed per second, not Sverdrups.

.. code-block:: python

    speed * width * depth

tidal_elevation
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.tidal_elevation`

Output: **m**. One supplied constituent; not a site-specific tide prediction or navigational product.

.. code-block:: python

    amplitude * np.cos(2 * np.pi * time / period + phase)

gerstner_orbit_radius
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.gerstner_orbit_radius`

Output: **m**. Exact deep-water trochoidal wave; requires ``amplitude * wavenumber <= 1``.

.. code-block:: python

    amplitude * np.exp(-wavenumber * rest_depth)

ekman_spiral_velocity
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.ekman_spiral_velocity`

Output: **m/s**, shape (..., 2). Steady, infinitely deep, constant eddy viscosity; crosswind positive to the right of the wind.

.. code-block:: python

    d = np.sqrt(2 * eddy_viscosity / abs(coriolis))
    v0 = wind_stress / (density * np.sqrt(eddy_viscosity * abs(coriolis)))
    s = depth / d
    v0 * np.exp(-s) * [np.sin(np.pi / 4 - s), np.sign(coriolis) * np.cos(np.pi / 4 - s)]

sverdrup_transport
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.sverdrup_transport`

Output: **m^2/s**. Steady linear interior balance away from boundaries and the equator.

.. code-block:: python

    wind_stress_curl / (density * beta)

stommel_streamfunction
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.stommel_streamfunction`

Output: **m^3/s**. Steady linear barotropic gyre with bottom friction in a rectangular basin.

.. code-block:: python

    k = np.pi / basin_length
    m1, m2 = (-beta + np.sqrt(beta**2 + 4 * friction**2 * k**2)) / (2 * friction), (-beta - ...) / (2 * friction)
    L = basin_width
    X = 1 + ((np.exp(m2 * L) - 1) * np.exp(m1 * (x - L)) + (np.exp(-m1 * L) - 1) * np.exp(m2 * x)) / (1 - np.exp((m2 - m1) * L))
    wind_stress * basin_length / (np.pi * density * friction) * np.sin(np.pi * y / basin_length) * X

thermohaline_forcing
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.thermohaline_forcing`

Output: **g/kg m^3/s**. Stommel two-box steady state with fixed temperatures; folds at ``flow = thermal_flow / 2``.

.. code-block:: python

    np.abs(flow) * (thermal_flow - flow) / haline_coefficient

abyssal_temperature
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.abyssal_temperature`

Output: **K**. Steady one-dimensional upwelling-diffusion balance with fixed end temperatures.

.. code-block:: python

    h = diffusivity / upwelling
    bottom_temperature + (top_temperature - bottom_temperature) * (np.exp((height - column_height) / h) - np.exp(-column_height / h)) / (1 - np.exp(-column_height / h))

simplified_seawater_density
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.simplified_seawater_density`

Output: **kg/m^3**. Roquet et al. (2015) fit to TEOS-10 gradients; not TEOS-10 absolute density.

.. code-block:: python

    t, s = conservative_temperature - 283.15, absolute_salinity - 35
    1026 - a0 * (1 + lambda1 * t / 2 + mu1 * depth) * t + b0 * (1 - lambda2 * s / 2 - mu2 * depth) * s - nu * t * s

teos10_density
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.teos10_density`

Output: **kg/m^3**. Roquet et al. (2015) 75-term polynomial for the specific volume, as in the GSW toolbox's ``gsw_rho``; it matches the TEOS-10 Gibbs function within the oceanographic range (about 0-42 g/kg, -2-40 degC and 0-80 MPa).

.. code-block:: python

    ys = (conservative_temperature - 273.15) / 40
    xs = np.sqrt(0.0248826675584615 * absolute_salinity + 0.5971840214030754)
    z = sea_pressure / 1e8
    1 / sum(v_ijk * ys**i * xs**j * z**k)  # the 75 coefficients v_ijk of GSW-C

Thermohaline boxes
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.oceanography.simulate_thermohaline_boxes` integrates
``volume * dS/dt = salt_forcing - abs(q) * dS`` with
``q = thermal_flow - haline_coefficient * dS`` and returns the salinity
difference and the overturning. Its steady states satisfy
``thermohaline_forcing``; below the fold the initial salinity difference
selects the thermal or the haline state.

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The shared suite adds cross-model identities, energy/conservation
checks and independent analytic solutions for dynamics. ``teos10_density``
reproduces the TEOS-10 check values of the GSW-C toolbox to 1e-13, stays within
3 kg/m^3 of the simplified fit, and expands more with temperature at depth
(thermobaricity). See :doc:`/validation`
for the exact evidence and :doc:`/numerical_methods` for numerical conventions.
