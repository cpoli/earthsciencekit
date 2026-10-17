Glaciology: units, equations and validation
===========================================

.. include:: /_generated/nav/glaciology.rst

Conventions and model limits
----------------------------

Lengths and thicknesses are metres, times seconds and stresses pascals.
Accumulation, ablation and mass-balance rates are metres of ice per second;
divide a rate in m/year by ``SECONDS_PER_YEAR`` and convert ages with
``units.convert(value, "s", "kyr")``. The balance gradient is in 1/s, metres
of ice per second per metre of altitude. Temperatures are kelvin; Dansgaard's
relation takes kelvin and returns delta-18O in permil, with the slope in
permil/K and the offset at 0 degC. :mod:`earthsciencekit.constants` holds
``ICE_DENSITY = 917`` kg/m^3 and ``LATENT_HEAT_FUSION = 3.34e5`` J/kg; the
models take densities as arguments.

Glen's rate factor ``A`` is in Pa^-n/s and its exponent ``n`` is at least 1;
with ``n = 1`` Glen's law is a Newtonian fluid of viscosity ``1 / (2 A)``,
Forbes's slab. Slopes are tangents (``surface_slope`` is a nonnegative
magnitude) and velocities and fluxes are positive down the surface slope.
The shallow-ice models (``driving_stress``, ``glen_velocity_profile``,
``shallow_ice_flux``, ``shallow_ice_diffusivity``, Vialov's profile and the
flowline glacier) assume ice much thinner than it is wide, frozen to its bed,
with a uniform rate factor. Weertman's sliding coefficient is in
m/(s Pa^m). Schoof's friction law is ``tau_b = C u**m`` with ``C`` in
Pa (s/m)^m; the grounding-line flux takes the flotation thickness, which
``flotation_thickness`` returns from the water depth.

Heights in ice are measured up from the bed and depths in ice cores down from
the surface, both in metres of ice equivalent; ages diverge at a frozen bed.
Isostatic depressions are positive downward and assume local isostasy;
``viscous_relaxation_time`` is for a uniform Newtonian half-space without a
lithosphere. Stefan's law takes the freezing degree-time in K s, the time
integral of the freezing point minus the surface temperature.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. Ice
must be less dense than water where it floats, heights cannot exceed the ice
thickness, and core depths must lie above the bed. No unit metadata are
carried in raw numbers. Bounds protect mathematical domains and obvious
physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

ela_shift
~~~~~~~~~

:py:func:`earthsciencekit.glaciology.ela_shift`

Output: **m**. Temperature alone controls ablation and the lapse rate is uniform; changes of snowfall also move the line.

.. code-block:: python

    temperature_change / lapse_rate

driving_stress
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.driving_stress`

Output: **Pa**. Shallow-ice balance: horizontal scales much larger than the thickness, small slopes, longitudinal stresses neglected.

.. code-block:: python

    density * gravity * thickness * surface_slope

viscous_slab_velocity
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.viscous_slab_velocity`

Output: **m/s**. Forbes (1843): a wide parallel-sided slab with a linear viscosity, no sliding and small slope.

.. code-block:: python

    density * gravity * surface_slope * (2 * thickness * height - height**2) / (2 * viscosity)

isostatic_depression
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.isostatic_depression`

Output: **m**. Jamieson (1865): local isostasy for a load much wider than the flexural wavelength of the lithosphere.

.. code-block:: python

    ice_thickness * ice_density / mantle_density

stefan_ice_thickness
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.stefan_ice_thickness`

Output: **m**. Stefan (1891): ice free of snow, with its surface at the air temperature, a linear temperature profile, no ocean heat flux and no initial ice; it overestimates growth under snow.

.. code-block:: python

    np.sqrt(2 * conductivity * freezing_degree_time / (density * latent_heat))

stefan_growth_rate
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.stefan_growth_rate`

Output: **m/s**. Stefan (1891), quasi-steady conduction with no ocean heat flux.

.. code-block:: python

    conductivity * temperature_difference / (density * latent_heat * thickness)

sea_ice_freeboard
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.sea_ice_freeboard`

Output: **m**. Hydrostatic equilibrium of a uniform floe without snow load.

.. code-block:: python

    thickness * (1 - ice_density / water_density)

viscous_relaxation_time
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.viscous_relaxation_time`

Output: **s**. Haskell (1935): a uniform Newtonian half-space without an elastic lithosphere, which slows the shortest wavelengths.

.. code-block:: python

    4 * np.pi * viscosity / (density * gravity * wavelength)

haskell_viscosity
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.haskell_viscosity`

Output: **Pa s**. Haskell (1935), the inverse of viscous_relaxation_time for a uniform half-space.

.. code-block:: python

    density * gravity * wavelength * relaxation_time / (4 * np.pi)

remaining_depression
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.remaining_depression`

Output: **m**. Haskell (1935): one wavelength on a viscous half-space, with the load removed at once.

.. code-block:: python

    initial_depression * np.exp(-time / relaxation_time)

rebound_uplift_rate
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.rebound_uplift_rate`

Output: **m/s**. Haskell (1935), the time derivative of remaining_depression.

.. code-block:: python

    initial_depression / relaxation_time * np.exp(-time / relaxation_time)

nye_plastic_profile
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.nye_plastic_profile`

Output: **m**. Nye (1951): a horizontal bed, basal shear stress everywhere equal to the yield stress and no isostatic depression.

.. code-block:: python

    np.sqrt(2 * yield_stress * distance / (density * gravity))

glen_strain_rate
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glen_strain_rate`

Output: **1/s**. Glen (1955): secondary creep of isotropic ice at stresses of about 0.05-1 MPa.

.. code-block:: python

    rate_factor * stress**exponent

glen_stress
~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glen_stress`

Output: **Pa**. Glen (1955), inverted; secondary creep of isotropic ice.

.. code-block:: python

    (strain_rate / rate_factor) ** (1 / exponent)

glen_rate_factor
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glen_rate_factor`

Output: **Pa^-n/s**. An Arrhenius law with one activation energy, about 60 kJ/mol below -10 degC and 115 kJ/mol above (Cuffey and Paterson 2010), for temperate-ice water contents neglected.

.. code-block:: python

    reference_factor * np.exp(-activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))

effective_viscosity
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.effective_viscosity`

Output: **Pa s**. Glen (1955), written as a stress-dependent Newtonian viscosity; with n = 1 it is the viscosity 1 / (2 A).

.. code-block:: python

    1 / (2 * rate_factor * stress ** (exponent - 1))

glen_velocity_profile
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glen_velocity_profile`

Output: **m/s**. Nye (1952) laminar flow with Glen (1955) rheology: a wide slab, uniform rate factor and no sliding.

.. code-block:: python

    2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * (thickness ** (exponent + 1) - (thickness - height) ** (exponent + 1))

glen_surface_velocity
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glen_surface_velocity`

Output: **m/s**. Nye (1952) laminar flow with Glen (1955) rheology, no sliding.

.. code-block:: python

    2 * rate_factor / (exponent + 1) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 1)

nye_crevasse_depth
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.nye_crevasse_depth`

Output: **m**. Nye (1955): a field of closely spaced crevasses in uniform ice of density independent of depth; fracture toughness is neglected.

.. code-block:: python

    2 * (strain_rate / rate_factor) ** (1 / exponent) / (density * gravity)

water_filled_crevasse_depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.water_filled_crevasse_depth`

Output: **m**. Nye (1955) with the water pressure added (Benn et al. 2007); the result is meaningful only while it is less than the ice thickness and at least the water depth.

.. code-block:: python

    2 * (strain_rate / rate_factor) ** (1 / exponent) / (ice_density * gravity) + water_density / ice_density * water_depth

weertman_sliding_velocity
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.weertman_sliding_velocity`

Output: **m/s**. Weertman (1957): a bed at the melting point, free of cavities and sediment; the coefficient lumps the bed roughness.

.. code-block:: python

    sliding_coefficient * basal_stress**exponent

linear_mass_balance
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.linear_mass_balance`

Output: **m/s**. Meier (1962): the balance gradient is about 0.003-0.01 m of water per year per metre on mid-latitude glaciers; real profiles flatten at high altitude.

.. code-block:: python

    balance_gradient * (elevation - equilibrium_line_altitude)

accumulation_area_ratio
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.accumulation_area_ratio`

Output: **1**. Meier (1962); the steady-state ratio depends on the glacier's hypsometry and is not universal.

.. code-block:: python

    accumulation_area / glacier_area

dansgaard_delta
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.dansgaard_delta`

Output: **permil**. Dansgaard (1964): a spatial relation for mid- and high-latitude coastal and polar stations; the temporal slope at one site can differ.

.. code-block:: python

    slope * (temperature - 273.15) + offset

dansgaard_temperature
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.dansgaard_temperature`

Output: **K**. Dansgaard (1964), inverted; the spatial slope is an approximation of the temporal one.

.. code-block:: python

    273.15 + (delta - offset) / slope

nye_age
~~~~~~~

:py:func:`earthsciencekit.glaciology.nye_age`

Output: **s**. Nye (1963): a steady divide with uniform vertical strain rate, no basal melting or sliding; the age diverges at the bed.

.. code-block:: python

    thickness / accumulation * np.log(thickness / (thickness - depth))

dansgaard_johnsen_downward_velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.dansgaard_johnsen_downward_velocity`

Output: **m/s**. Dansgaard and Johnsen (1969): steady state, constant horizontal strain rate above the kink and a linear decrease below it, no basal melting. Below, z = thickness - depth and h = kink_height.

.. code-block:: python

    np.where(z >= h, accumulation * (2 * z - h) / (2 * thickness - h), accumulation * z**2 / (h * (2 * thickness - h)))

dansgaard_johnsen_age
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.dansgaard_johnsen_age`

Output: **s**. Dansgaard and Johnsen (1969), the integral of the reciprocal of dansgaard_johnsen_downward_velocity; it tends to nye_age as the kink height goes to zero. Below, z = thickness - depth, h = kink_height, a = accumulation and s = 2 * thickness - h.

.. code-block:: python

    s / (2 * a) * np.log(s / (2 * np.maximum(z, h) - h)) + h * s / a * (1 / np.minimum(z, h) - 1 / h)

flotation_thickness
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.flotation_thickness`

Output: **m**. Archimedes' principle for an ice column of uniform density.

.. code-block:: python

    water_density * bed_depth / ice_density

ice_shelf_strain_rate
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.ice_shelf_strain_rate`

Output: **1/s**. Weertman (1957) for an unconfined shelf spreading in one direction, without buttressing or side drag.

.. code-block:: python

    rate_factor * (ice_density * gravity * thickness * (1 - ice_density / water_density) / 4) ** exponent

shallow_ice_flux
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.shallow_ice_flux`

Output: **m^2/s**. Hutter (1983) and Mahaffy (1976): Glen ice frozen to its bed, aspect ratio small, uniform rate factor.

.. code-block:: python

    2 * rate_factor / (exponent + 2) * (density * gravity * surface_slope) ** exponent * thickness ** (exponent + 2)

shallow_ice_diffusivity
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.shallow_ice_diffusivity`

Output: **m^2/s**. Mahaffy (1976): the shallow-ice approximation for Glen ice frozen to its bed. Below, A = rate_factor, n = exponent.

.. code-block:: python

    2 * A / (n + 2) * (density * gravity) ** n * thickness ** (n + 2) * surface_slope ** (n - 1)

vialov_dome_height
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.vialov_dome_height`

Output: **m**. Vialov (1958): a flat bed, uniform accumulation up to a fixed margin, shallow-ice flow without sliding. Below, A = rate_factor, n = exponent.

.. code-block:: python

    (2 * ((n + 2) * accumulation / (2 * A * (density * gravity) ** n)) ** (1 / n)) ** (n / (2 * n + 2)) * np.sqrt(half_width)

vialov_profile
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.vialov_profile`

Output: **m**. Vialov (1958), see vialov_dome_height; the shallow-ice approximation fails within a few thicknesses of the divide and the margin.

.. code-block:: python

    dome_height * (1 - (distance / half_width) ** ((exponent + 1) / exponent)) ** (exponent / (2 * exponent + 2))

glacier_response_time
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.glacier_response_time`

Output: **s**. Johannesson, Raymond and Waddington (1989): a linear response of volume; decades for valley glaciers, millennia for ice sheets.

.. code-block:: python

    thickness / terminus_ablation

gas_age_difference
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.gas_age_difference`

Output: **s**. Schwander and Stauffer (1984), as used for Vostok by Barnola et al. (1987): steady firn densification with no thinning above close-off; diffusion through the open pores, which makes the gas slightly older, is neglected.

.. code-block:: python

    close_off_depth * mean_firn_density / (ice_density * accumulation)

schoof_grounding_line_flux
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.schoof_grounding_line_flux`

Output: **m^2/s**. Schoof (2007): a flowline with a power-law sliding bed, an unbuttressed shelf and a boundary layer at the grounding line; steady states are stable only where the bed rises inland. Below, A = rate_factor, n = exponent, C = friction_coefficient, m = friction_exponent, rho_i = ice_density, rho_w = water_density, h = grounding_thickness.

.. code-block:: python

    (A * (rho_i * gravity)**(n + 1) * (1 - rho_i / rho_w)**n / (4**n * C))**(1 / (m + 1)) * h**((m + n + 3) / (m + 1))

Flowline glacier
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.glaciology.simulate_flowline_glacier` integrates
the ice thickness ``H`` on equally spaced nodes over a given bed ``b``,
``dH/dt = beta (b + H - ELA) - dq/dx`` with the shallow-ice flux
``q = -D ds/dx`` evaluated between nodes, as in Mahaffy (1976), by the method
of lines with SciPy's BDF integrator. The ends carry no flux and ice-free
nodes cannot ablate. Pass the steady thickness for one equilibrium line as
``initial_thickness`` with another to follow the response to a step. It
returns the length (nodes with more than 1 m of ice), the area per unit width
and the thickness at every node.

Independent validation
----------------------

Every algebraic API has a fixed-value check, bound checks and invalid-input
contracts. The glaciology tests also check that Glen's profile with
``n = 1`` is Forbes's viscous slab, that Glen's stress, strain rate and
effective viscosity are consistent, that the shallow-ice flux is the depth
integral of Glen's velocity profile and equals the diffusivity times the
slope, that Nye's plastic profile carries the yield stress everywhere, that
Vialov's profile carries the snowfall accumulated upstream, that water deepens
a crevasse, that Stefan's thickness grows at Stefan's rate, that Haskell's
viscosity inverts the relaxation time and the uplift rate integrates to the
recovered depression, that Dansgaard's thermometer inverts, that the
Dansgaard-Johnsen age integrates the reciprocal of its velocity and tends to
Nye's age, that Schoof's flux rises as the 4.75 power of the grounding
thickness, and that ``simulate_flowline_glacier`` reaches a steady state whose
flux equals the mass balance accumulated upstream, shrinks when the
equilibrium line rises and stays ice-free above the summit. See
:doc:`/validation` for the exact evidence and :doc:`/numerical_methods` for
numerical conventions.
