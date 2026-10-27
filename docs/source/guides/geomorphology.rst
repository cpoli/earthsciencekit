Geomorphology: units, equations and validation
==============================================

.. include:: /_generated/nav/geomorphology.rst

Conventions and model limits
----------------------------

Lengths, depths and grain diameters are metres, areas square metres, times
seconds and stresses pascals. Erosion, uplift and soil-production rates are
metres of rock per second; divide a rate in m/year by ``SECONDS_PER_YEAR``.
Densities are arguments (2650 kg/m^3 for quartz, 1000 for water, about 1.2
for air); viscosities are dynamic (Pa s) in Stokes's law and kinematic
(m^2/s) elsewhere. Gradients are signed ``dz/dx`` and fluxes positive toward
increasing ``x``; slopes are nonnegative tangents. The von Karman constant is
``constants.VON_KARMAN = 0.41``.

Sediment transport assumes uniform, cohesionless grains on a flat bed in
steady uniform flow: the bed shear stress is the depth-slope product, the
threshold is the Soulsby-Whitehouse fit to Shields's data, suspension
follows Rouse's parabolic diffusivity and bedload Meyer-Peter and Müller's
law with the skin-friction stress. Sand fluxes are mass per unit width
(kg/(m s)) and bedload fluxes volume of grains per unit width (m^2/s).

Stream power is ``E = K A**m S**n`` with drainage area as a proxy for
discharge; the units of ``K`` are m^(1 - 2m)/s when ``n = 1``. Hack's
coefficient is about 0.32 in SI for an exponent of 0.6 (Hack's 1.4 in
miles). Strahler orders take a vector of downstream link indices with -1 at
outlets. Cosmogenic concentrations are atoms/kg, production rates
atoms/(kg s) and attenuation lengths kg/m^2 (1600 for spallation);
``constants.LAMBDA_BE10`` is the beryllium-10 decay constant for a 1.387 Myr
half-life. Morphologic ages are diffusivity times age in m^2.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments.
Sediment must be denser than its fluid, a Rouse height must lie between the
reference height and the surface, a nonlinear-creep gradient must be gentler
than the critical one, and a cosmogenic concentration cannot exceed its
saturation. ``chi_coordinate`` needs strictly increasing distances and
``strahler_order`` a network without cycles. No unit metadata are carried in
raw numbers. Bounds protect mathematical domains and obvious physical
contradictions, but do not establish model applicability.

Equation and output map
-----------------------

stokes_settling_velocity
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.stokes_settling_velocity`

Output: **m/s**. Stokes (1851): a smooth sphere at particle Reynolds number below about 1, roughly quartz finer than 0.1 mm in water.

.. code-block:: python

    (sediment_density - fluid_density) * gravity * diameter**2 / (18 * viscosity)

ferguson_church_settling_velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.ferguson_church_settling_velocity`

Output: **m/s**. Ferguson and Church (2004) with C1 = 18 and C2 = 1 for natural grains; Stokes's law for fine grains and a constant drag coefficient for coarse ones; spheres need C2 = 0.4.

.. code-block:: python

    R * gravity * diameter**2 / (18 * kinematic_viscosity + np.sqrt(0.75 * R * gravity * diameter**3))
    R = sediment_density / fluid_density - 1

bed_shear_stress
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.bed_shear_stress`

Output: **Pa**. The depth-slope product of steady uniform flow; bank friction and bedforms take part of it in real channels.

.. code-block:: python

    density * gravity * hydraulic_radius * slope

shear_velocity
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.shear_velocity`

Output: **m/s**. A definition; it applies to any turbulent boundary layer, in water or air.

.. code-block:: python

    np.sqrt(shear_stress / density)

shields_number
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.shields_number`

Output: **1**. Shields (1936): uniform cohesionless grains on a nearly flat bed.

.. code-block:: python

    shear_stress / ((sediment_density - fluid_density) * gravity * diameter)

dimensionless_grain_size
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.dimensionless_grain_size`

Output: **1**. Soulsby and Whitehouse (1997); about 25 per millimetre of quartz in water at 20 degC.

.. code-block:: python

    diameter * ((sediment_density / fluid_density - 1) * gravity / kinematic_viscosity**2)**(1 / 3)

critical_shields_number
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.critical_shields_number`

Output: **1**. Soulsby and Whitehouse (1997) fit to Shields's (1936) data, with the rise for fine grains; uniform cohesionless grains on a flat bed, with a scatter of about a factor of two.

.. code-block:: python

    0.30 / (1 + 1.2 * D) + 0.055 * (1 - np.exp(-0.020 * D))
    D = dimensionless_grain_size

keulegan_mean_velocity
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.keulegan_mean_velocity`

Output: **m/s**. Keulegan (1938): fully rough flow in a wide channel with the von Karman constant 0.41; the depth should exceed the roughness several times.

.. code-block:: python

    shear_velocity / VON_KARMAN * np.log(11 * depth / roughness_height)

rouse_number
~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.rouse_number`

Output: **1**. Rouse (1937), with the eddy diffusivity of sediment equal to that of momentum and the von Karman constant 0.41.

.. code-block:: python

    settling_velocity / (VON_KARMAN * shear_velocity)

rouse_concentration
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.rouse_concentration`

Output: **same as reference_concentration**. Rouse (1937): steady uniform flow with a parabolic eddy diffusivity and dilute suspension; zero at the surface.

.. code-block:: python

    reference_concentration * ((depth - height) / height * reference_height / (depth - reference_height))**rouse_number

meyer_peter_muller_transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.meyer_peter_muller_transport`

Output: **1**. Meyer-Peter and Müller (1948): flume data for well-sorted sand and gravel on a plane bed, with the skin friction part of the stress; Wong and Parker (2006) correct it to ``3.97 (theta - 0.0495)**1.5``.

.. code-block:: python

    8 * np.maximum(shields_number - critical_shields_number, 0)**1.5

bedload_flux
~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.bedload_flux`

Output: **m^2/s**. Einstein (1950) scaling; uniform grains.

.. code-block:: python

    transport_rate * np.sqrt((sediment_density / fluid_density - 1) * gravity * diameter) * diameter

diffusive_hillslope_flux
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.diffusive_hillslope_flux`

Output: **m^2/s**. Gilbert (1909) and Culling (1960): gentle soil-mantled slopes; the flux grows faster than linearly near the angle of repose (see :func:`nonlinear_hillslope_flux`).

.. code-block:: python

    -diffusivity * gradient

nonlinear_hillslope_flux
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.nonlinear_hillslope_flux`

Output: **m^2/s**. Roering, Kirchner and Dietrich (1999): reduces to :func:`diffusive_hillslope_flux` on gentle slopes; landslides take over at the critical gradient.

.. code-block:: python

    -diffusivity * gradient / (1 - (gradient / critical_gradient)**2)

scarp_elevation
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.scarp_elevation`

Output: **m**. Culling (1960) and Hanks et al. (1984): a vertical step on a uniform slope degraded by linear diffusion without base-level change; a fresh step (``diffusivity * time = 0``) is returned exactly.

.. code-block:: python

    half_offset * erf(distance / (2 * np.sqrt(diffusivity * time))) + far_field_gradient * distance

scarp_maximum_gradient
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.scarp_maximum_gradient`

Output: **1**. Hanks et al. (1984), the derivative of :func:`scarp_elevation` at its midpoint; same assumptions.

.. code-block:: python

    half_offset / np.sqrt(np.pi * diffusivity * time) + far_field_gradient

soil_production_rate
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.soil_production_rate`

Output: **m/s**. Heimsath et al. (1997), from cosmogenic nuclides under soils of a few metres at most; some sites show a humped function with a maximum under a thin soil.

.. code-block:: python

    maximum_rate * np.exp(-soil_depth / decay_depth)

hilltop_curvature
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.hilltop_curvature`

Output: **1/m**. Roering, Perron and Kirchner (2007): steady state with linear creep at the divide, where the gradient is small.

.. code-block:: python

    -rock_density * erosion_rate / (soil_density * diffusivity)

stream_power_erosion
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.stream_power_erosion`

Output: **m/s**. Howard and Kerby (1983): detachment-limited bedrock channels, with discharge and width folded into powers of the area; no threshold, sediment cover or tools effect.

.. code-block:: python

    erodibility * drainage_area**area_exponent * slope**slope_exponent

steady_state_slope
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.steady_state_slope`

Output: **1**. Steady state of :func:`stream_power_erosion` under uniform uplift; the slope follows Flint's law with concavity ``area_exponent / slope_exponent``.

.. code-block:: python

    (uplift_rate / (erodibility * drainage_area**area_exponent))**(1 / slope_exponent)

knickpoint_celerity
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.knickpoint_celerity`

Output: **m/s**. Whipple and Tucker (1999), the characteristic speed of :func:`stream_power_erosion`; independent of the slope only when ``slope_exponent = 1``, where knickpoints keep their shape.

.. code-block:: python

    erodibility * drainage_area**area_exponent * slope**(slope_exponent - 1)

flint_slope
~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.flint_slope`

Output: **1**. Flint (1974): below the colluvial-fluvial transition at about 1 km^2, in rivers without knickpoints or changes of lithology.

.. code-block:: python

    steepness_index * drainage_area**(-concavity)

chi_coordinate
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.chi_coordinate`

Output: **m**. Perron and Royden (2013), integrated by the trapezoidal rule between nodes; in steady state with :func:`steady_state_slope` the elevation is linear in chi.

.. code-block:: python

    chi = integral of (reference_area / drainage_area)**concavity d distance

hack_length
~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.hack_length`

Output: **m**. Hack (1957), from basins in Virginia and Maryland; the exponent drifts toward 0.5 in the largest basins.

.. code-block:: python

    coefficient * drainage_area**exponent

strahler_order
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.strahler_order`

Output: **1**. Strahler (1952, 1957): a link takes the highest order among those draining into it, plus one when two or more share that highest order; several outlets are allowed.

.. code-block:: python

    order = max(upstream orders) + (count of that maximum >= 2)

downstream_hydraulic_geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.downstream_hydraulic_geometry`

Output: **width m; depth m; velocity m/s**. Leopold and Maddock (1953): power laws fitted to gauging stations, with the velocity from continuity so that its coefficient and exponent are ``1 / (a c)`` and ``1 - b - f``.

.. code-block:: python

    width = a * Q**b
    depth = c * Q**f
    velocity = Q / (width * depth)

effective_discharge
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.effective_discharge`

Output: **m^3/s**. Wolman and Miller (1960), with a lognormal frequency of discharge and a power-law sediment rating; the maximum of ``Q**b`` times the lognormal density.

.. code-block:: python

    median_discharge * np.exp((rating_exponent - 1) * log_standard_deviation**2)

threshold_shear_velocity
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.threshold_shear_velocity`

Output: **m/s**. Bagnold (1941): dry cohesionless sand coarser than about 0.1 mm; moisture, crusts and fine grains raise the threshold.

.. code-block:: python

    coefficient * np.sqrt((sediment_density - air_density) / air_density * gravity * diameter)

bagnold_sand_flux
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.bagnold_sand_flux`

Output: **kg/(m s)**. Bagnold (1941), referred to a 0.25 mm sand; well above the threshold, since it has no threshold term (see :func:`kawamura_sand_flux`).

.. code-block:: python

    coefficient * np.sqrt(diameter / 2.5e-4) * air_density / gravity * shear_velocity**3

kawamura_sand_flux
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.kawamura_sand_flux`

Output: **kg/(m s)**. Kawamura (1951): steady saltation over a flat dry sand bed; tends to the cube law far above the threshold.

.. code-block:: python

    coefficient * air_density / gravity * np.maximum(u - ut, 0) * (u + ut)**2
    u = shear_velocity
    ut = threshold_shear_velocity

dune_migration_rate
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.dune_migration_rate`

Output: **m/s**. Bagnold (1941): a two-dimensional dune of constant shape that traps all the sand crossing its crest.

.. code-block:: python

    sand_flux / (bulk_density * dune_height)

cosmogenic_production_at_depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.cosmogenic_production_at_depth`

Output: **atoms/(kg s)**. Lal (1991): neutron spallation only, on a flat surface; muons dominate below a few metres.

.. code-block:: python

    surface_production_rate * np.exp(-density * depth / attenuation_length)

cosmogenic_concentration
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.cosmogenic_concentration`

Output: **atoms/kg**. Lal (1991): steady erosion and production, spallation only, no inheritance; ``production_rate * time`` when neither decay nor erosion removes nuclides.

.. code-block:: python

    production_rate / r * (1 - np.exp(-r * time))
    r = decay_constant + density * erosion_rate / attenuation_length

steady_state_erosion_rate
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.steady_state_erosion_rate`

Output: **m/s**. Lal (1991), the inverse of :func:`cosmogenic_concentration` at long times; valid once the surface has been exposed for several times ``attenuation_length / (density * erosion_rate)``.

.. code-block:: python

    (production_rate / concentration - decay_constant) * attenuation_length / density

morphologic_age
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.morphologic_age`

Output: **m^2**. Hanks et al. (1984), inverting :func:`scarp_maximum_gradient`: one rupture, linear diffusion and no base-level change at the foot.

.. code-block:: python

    half_offset**2 / (np.pi * (maximum_gradient - far_field_gradient)**2)

luminescence_age
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.luminescence_age`

Output: **s**. Aitken (1985): complete bleaching at deposition, a constant dose rate and no fading or saturation of the signal.

.. code-block:: python

    equivalent_dose / dose_rate

River profile
~~~~~~~~~~~~~

:py:func:`earthsciencekit.geomorphology.simulate_river_profile` integrates
the elevation ``z`` of a river from a fixed outlet upstream,
``dz/dt = U - K A**m max(dz/dx, 0)**n``, by the method of lines with
first-order upwind slopes (toward each node's downstream neighbour) and
SciPy's BDF integrator. Its discrete steady state is the cumulative sum of
:func:`~earthsciencekit.geomorphology.steady_state_slope` with the area of
the upstream node, which ``chi_coordinate`` matches to discretization error.
Pass a steady profile as ``initial_elevation`` with a new uplift rate to send
a knickpoint upstream, or with zero uplift to watch relief decay. It returns
the relief, the mean elevation and the elevation of every node; there is no
sediment flux, threshold, hillslope coupling or isostasy.

Independent validation
----------------------

Every algebraic API has a fixed-value check, bound checks and invalid-input
contracts. The geomorphology tests also check that the Ferguson-Church
velocity tends to Stokes's law for fine grains and to a square-root law for
coarse ones, that the Shields curve has its minimum in the sand sizes and
levels at 0.055, that Keulegan's velocity is the depth average of the log
profile, that the Rouse profile balances settling against turbulent mixing,
that nonlinear creep is linear on gentle slopes, that the degraded scarp
solves the diffusion equation and its morphologic age inverts the steepest
gradient, that the hilltop curvature balances steady lowering, that the
steady stream-power slope is Flint's law and erodes at the uplift rate, that
chi reduces to distance at the reference area, that Strahler orders balanced
and unbalanced trees, that the hydraulic geometry conserves discharge, that
the effective discharge is the numerical maximum of rate times frequency,
that Kawamura's flux is Bagnold's without a threshold, and that cosmogenic
concentrations accumulate, saturate and invert to the erosion rate.
``simulate_river_profile`` reaches its discrete steady state, sends a
knickpoint upstream at the stream-power celerity and loses relief ever more
slowly with a slope exponent of 2. See :doc:`/validation` for the exact
evidence and :doc:`/numerical_methods` for numerical conventions.
