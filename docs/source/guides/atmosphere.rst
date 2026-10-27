Atmosphere: units, equations and validation
===========================================

.. include:: /_generated/nav/atmosphere.rst

Conventions and model limits
----------------------------

Pressure is Pa, temperature K, height m and latitude radians. The isothermal
barometric model assumes constant temperature; the dry-ascent model has a
changing temperature and should not be compared with that profile directly.
Potential temperature uses a caller-supplied reference pressure.

Mixing ratio is water mass per dry-air mass; specific humidity is water mass
per total moist-air mass. Relative humidity is a fraction and can exceed one.
Saturation vapor pressure is a constant-latent-heat Clausius-Clapeyron
approximation over liquid water, restricted to 180-330 K. It is not saturation
over ice. Dry parcel ascent excludes condensation.

Moist convection uses the same saturation formula. The lifting condensation
level is Bolton's (1980) fit, within about 0.1 K, and its height assumes dry
adiabatic ascent. The moist lapse rate is pseudo-adiabatic: condensate falls
out at once and its heat capacity is ignored. ``parcel_temperature_profile``
takes a strictly decreasing pressure vector that starts at the parcel's level.
CAPE and CIN take heights and virtual temperatures on the same levels and
integrate buoyancy exactly for a profile linear between levels. CAPE is
nonnegative. CIN is nonpositive and counts only the layers below the level of
free convection, where the parcel first becomes buoyant.

The two-layer temperature profile has a constant lapse rate up to the
tropopause and is isothermal above. Richardson number and the surface pressure
tendency are signed. Rossby-wave speed takes the total horizontal wavenumber
in rad/m and returns the eastward phase speed. The barotropic vorticity model
uses a doubly periodic beta plane, so its zonal and meridional boundaries are
artificial.

N squared is signed: negative values indicate instability. Geostrophic speed
and Rossby scaling reject exactly zero Coriolis parameter, but the caller must
also recognize their poor applicability near the equator.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. No unit
metadata are carried in raw numbers. Bounds protect mathematical domains and
obvious physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

dry_air_density
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.dry_air_density`

Output: **kg/m^3**. Dry ideal air with fixed specific gas constant.

.. code-block:: python

    pressure / (R_DRY_AIR * temperature)

scale_height
~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.scale_height`

Output: **m**. Dry isothermal atmosphere with constant gravity.

.. code-block:: python

    R_DRY_AIR * temperature / gravity

barometric_pressure
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.barometric_pressure`

Output: **Pa**. Isothermal column; signed height is relative to the reference level.

.. code-block:: python

    surface_pressure * np.exp(-gravity * height / (R_DRY_AIR * temperature))

potential_temperature
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.potential_temperature`

Output: **K**. Dry ideal gas with constant heat capacity; not equivalent potential temperature.

.. code-block:: python

    temperature * (reference_pressure / pressure)**(R_DRY_AIR / CP_DRY_AIR)

temperature_from_potential
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.temperature_from_potential`

Output: **K**. Uses the same dry-air heat capacity and reference pressure convention.

.. code-block:: python

    potential_temperature * (pressure / reference_pressure)**(R_DRY_AIR / CP_DRY_AIR)

dry_adiabatic_lapse_rate
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.dry_adiabatic_lapse_rate`

Output: **K/m**. Returned positive lapse magnitude in K per metre; temperature decreases upward.

.. code-block:: python

    gravity / heat_capacity

saturation_vapor_pressure
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.saturation_vapor_pressure`

Output: **Pa**. Liquid-water reference with constant latent heat; restricted to 180-330 K, not an ice saturation formula.

.. code-block:: python

    611.2 * np.exp((2.5e6 / R_VAPOR) * (1 / 273.15 - 1 / temperature))

mixing_ratio
~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.mixing_ratio`

Output: **kg/kg dry air**. Ideal vapor mixture; vapor pressure must be below total pressure.

.. code-block:: python

    EPSILON_WATER * vapor_pressure / (total_pressure - vapor_pressure)

specific_humidity
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.specific_humidity`

Output: **kg/kg moist air**. Output is water mass per total moist-air mass.

.. code-block:: python

    mixing_ratio / (1 + mixing_ratio)

virtual_temperature
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.virtual_temperature`

Output: **K**. Ideal vapor mixture without suspended liquid or ice.

.. code-block:: python

    temperature * (1 + (1 / EPSILON_WATER - 1) * specific_humidity)

relative_humidity
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.relative_humidity`

Output: **1**. Returned fraction may exceed one during supersaturation; no clipping.

.. code-block:: python

    vapor_pressure / saturation_pressure

coriolis_parameter
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.coriolis_parameter`

Output: **1/s**. Spherical rotating-Earth approximation; input is radians.

.. code-block:: python

    2 * OMEGA_EARTH * np.sin(latitude)

brunt_vaisala_squared
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.brunt_vaisala_squared`

Output: **1/s^2**. Returns signed N squared; negative values represent instability rather than an imaginary frequency.

.. code-block:: python

    gravity * vertical_gradient / potential_temperature

geostrophic_speed
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.geostrophic_speed`

Output: **m/s**. Signed cross-gradient component; invalid at the equator where coriolis vanishes.

.. code-block:: python

    pressure_gradient / (density * coriolis)

rossby_number
~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.rossby_number`

Output: **1**. Magnitude scaling only; small values support but do not prove geostrophic balance.

.. code-block:: python

    speed / (np.abs(coriolis) * length)

barometer_column_height
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.barometer_column_height`

Output: **m**. Incompressible liquid with a vacuum above the column; vapour pressure and capillarity are neglected.

.. code-block:: python

    pressure / (liquid_density * gravity)

two_layer_temperature
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.two_layer_temperature`

Output: **K**. Constant-lapse troposphere under an isothermal lower stratosphere; the stratospheric warming above about 20 km is omitted.

.. code-block:: python

    surface_temperature - lapse_rate * np.minimum(height, tropopause_height)

richardson_number
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.richardson_number`

Output: **1**. Signed; stratified shear flow can become unstable only where it falls below 1/4 (Miles-Howard).

.. code-block:: python

    brunt_vaisala_squared / shear**2

surface_pressure_tendency
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.surface_pressure_tendency`

Output: **Pa/s**. Hydrostatic column with no surface mass flux; divergence at one level is far larger than the column mean.

.. code-block:: python

    -surface_pressure * mean_divergence

rossby_wave_speed
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.rossby_wave_speed`

Output: **m/s**. Nondivergent barotropic flow on a beta plane; positive values are eastward.

.. code-block:: python

    zonal_wind - beta / wavenumber**2

simulate_barotropic_vorticity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.simulate_barotropic_vorticity`

Output: **1/s**, relative vorticity on a doubly periodic grid, flattened row by
row. Pseudo-spectral derivatives with two-thirds dealiasing; the domain mean is
removed. No forcing, friction, divergence or sphere.

.. code-block:: python

    d(zeta)/dt = -J(psi, zeta) - beta * d(psi)/dx,  laplacian(psi) = zeta

lcl_temperature
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.lcl_temperature`

Output: **K**. Bolton (1980, equation 15), accurate to about 0.1 K; it does not depend on pressure.

.. code-block:: python

    1 / (1 / (dewpoint - 56) + np.log(temperature / dewpoint) / 800) + 56

lcl_pressure
~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.lcl_pressure`

Output: **Pa**. Poisson's dry adiabat with constant heat capacity.

.. code-block:: python

    pressure * (lcl_temperature / temperature)**(CP_DRY_AIR / R_DRY_AIR)

lcl_height
~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.lcl_height`

Output: **m**. Dry adiabatic ascent without mixing, about 125 m per kelvin of dew-point depression.

.. code-block:: python

    CP_DRY_AIR * (temperature - lcl_temperature) / GRAVITY

moist_adiabatic_lapse_rate
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.moist_adiabatic_lapse_rate`

Output: **K/m**, a positive magnitude. Pseudo-adiabatic, 180-330 K; ``r_s`` is the saturation mixing ratio and ``L`` the latent heat of vaporization.

.. code-block:: python

    GRAVITY * (1 + L * r_s / (R_DRY_AIR * T)) / (CP_DRY_AIR + L**2 * r_s * EPSILON_WATER / (R_DRY_AIR * T**2))

parcel_temperature_profile
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.parcel_temperature_profile`

Output: **K** at each pressure. Dry adiabat to Bolton's lifting condensation level, then the pseudo-adiabat integrated in ``ln p`` with ``integrators.integrate_ode``.

.. code-block:: python

    T = temperature * (pressure / pressure[0])**(R_DRY_AIR / CP_DRY_AIR)   # below cloud base
    dT/d(ln p) = moist_adiabatic_lapse_rate(T, p) * R_DRY_AIR * T / GRAVITY  # above

convective_available_potential_energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.convective_available_potential_energy`

Output: **J/kg**. Moncrieff and Miller (1976): every positively buoyant layer, exact for buoyancy linear between levels.

.. code-block:: python

    sum(integral of max(GRAVITY * (Tv_parcel - Tv_env) / Tv_env, 0) dz)

convective_inhibition
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.convective_inhibition`

Output: **J/kg**, nonpositive. Negative buoyancy below the level of free convection, or over the whole profile if the parcel is never buoyant.

.. code-block:: python

    -sum(integral of max(-GRAVITY * (Tv_parcel - Tv_env) / Tv_env, 0) dz)  # below the first buoyant level

simulate_moist_ascent
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.atmosphere.simulate_moist_ascent`

Output: temperature **K**, pressure **Pa**, vapour mixing ratio and condensed
water **kg/kg**. The parcel rises at constant speed. The dry phase is closed
form, the moment of saturation is found by root finding, and the moist phase
follows the pseudo-adiabat with hydrostatic pressure. There is no entrainment
or ice.

.. code-block:: python

    dT/dt = -w * (GRAVITY / CP_DRY_AIR if unsaturated else moist_adiabatic_lapse_rate(T, p))
    dp/dt = -GRAVITY * p * w / (R_DRY_AIR * T)

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The shared suite adds cross-model identities, energy/conservation
checks and independent analytic solutions for dynamics. The moist-convection
tests check that a parcel lifted to Bolton's condensation level is saturated
there, that the moist lapse rate tends to the dry rate in cold air, that
``simulate_moist_ascent`` keeps its potential temperature below cloud base,
conserves total water and follows ``parcel_temperature_profile`` above it, and
that CAPE and CIN split the buoyancy integral of a sounding. See :doc:`/validation`
for the exact evidence and :doc:`/numerical_methods` for numerical conventions.
