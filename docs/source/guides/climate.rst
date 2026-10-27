Climate: units, equations and validation
========================================

.. include:: /_generated/nav/climate.rst

Conventions and model limits
----------------------------

Solar and radiative fluxes are W/m^2. Global absorbed sunlight divides the
solar constant by four through sphere geometry. Effective radiative temperature
is not greenhouse surface temperature. A positive restoring coefficient means
a stable negative net feedback in the convention C dT/dt = F - lambda T.

Heat capacity is J/(m^2 K); relaxation time is seconds. Temperature anomalies
may be negative, unlike absolute temperatures used in Stefan-Boltzmann emission.
CO2 forcing uses concentration ratios with both values in the same ppm units
and a caller-supplied coefficient. The 5.35 example is an educational older
logarithmic parameterization, not a claim about a current assessed best estimate.

The single-exponential carbon reservoir omits the multiple removal timescales
of atmospheric carbon dioxide. Steric sea level is a constant-coefficient slab
approximation; added water mass is distributed over a supplied fixed area.

The two-box model of Held et al. (2010) couples a surface layer of heat
capacity ``C`` to a deep ocean ``C0`` through an exchange coefficient
``gamma``, both feedback and exchange in W/(m^2 K). ``two_box_timescales``
returns the fast and slow timescales on the last axis. The latitudinal
energy-balance model is North's (1975) diffusive model for one hemisphere, on
equal steps of the sine of latitude ``x``. Its outgoing longwave is linear,
``A + B (T - 273.15)``, its ice edge is where the temperature falls to the
freezing temperature (-10 degC by default), and ``ice_edge`` is a latitude in
radians: ``pi / 2`` when there is no ice and 0 for a snowball.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. No unit
metadata are carried in raw numbers. Bounds protect mathematical domains and
obvious physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

absorbed_solar
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.absorbed_solar`

Output: **W/m^2**. Global mean rather than local insolation; fixed bolometric albedo.

.. code-block:: python

    solar_constant * (1 - albedo) / 4

effective_temperature
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.effective_temperature`

Output: **K**. Single emitting level in equilibrium; not the surface temperature of a greenhouse atmosphere.

.. code-block:: python

    (absorbed_flux / (emissivity * SIGMA_SB))**0.25

outgoing_longwave
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.outgoing_longwave`

Output: **W/m^2**. Constant emissivity; not a spectral atmospheric radiative-transfer solver.

.. code-block:: python

    emissivity * SIGMA_SB * temperature**4

planck_feedback
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.planck_feedback`

Output: **W/(m^2 K)**. Positive restoring magnitude in W per square metre per kelvin; other feedbacks are separate.

.. code-block:: python

    4 * emissivity * SIGMA_SB * temperature**3

co2_forcing
~~~~~~~~~~~

:py:func:`earthsciencekit.climate.co2_forcing`

Output: **W/m^2**. User-supplied coefficient; an educational parameterization, not a current assessed forcing product.

.. code-block:: python

    coefficient * np.log(concentration / reference_concentration)

equilibrium_warming
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.equilibrium_warming`

Output: **K**. Stable linear response with a positive net restoring coefficient.

.. code-block:: python

    forcing / restoring

heat_capacity_column
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.heat_capacity_column`

Output: **J/(m^2 K)**. Uniform slab with constant material properties.

.. code-block:: python

    density * heat_capacity * depth

response_timescale
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.response_timescale`

Output: **s**. One-box linear model; output is seconds.

.. code-block:: python

    areal_heat_capacity / restoring

step_temperature
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.step_temperature`

Output: **K**. Zero initial anomaly and constant parameters; no deep-ocean exchange.

.. code-block:: python

    forcing / restoring * (-np.expm1(-restoring * time / areal_heat_capacity))

heat_uptake
~~~~~~~~~~~

:py:func:`earthsciencekit.climate.heat_uptake`

Output: **J/m^2**. One uniform layer; temperature_change is an anomaly.

.. code-block:: python

    areal_heat_capacity * temperature_change

ice_albedo
~~~~~~~~~~

:py:func:`earthsciencekit.climate.ice_albedo`

Output: **1**. Illustrative smooth closure, not a calibrated satellite albedo retrieval.

.. code-block:: python

    warm_albedo + (cold_albedo - warm_albedo) * (1 - np.tanh((temperature - transition_temperature) / width)) / 2

carbon_airborne_fraction
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.carbon_airborne_fraction`

Output: **1**. Same time interval and carbon mass units; feedback-driven values above one are not clipped.

.. code-block:: python

    atmospheric_increase / emissions

carbon_decay
~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.carbon_decay`

Output: **kg C**. Not a realistic complete atmospheric CO2 impulse response, which has multiple timescales.

.. code-block:: python

    initial_excess * np.exp(-time / timescale)

thermal_expansion_sea_level
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.thermal_expansion_sea_level`

Output: **m**. Constant coefficient and uniform warming; salinity and pressure dependence are omitted.

.. code-block:: python

    expansion_coefficient * depth * temperature_change

water_equivalent_sea_level
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.water_equivalent_sea_level`

Output: **m**. Uniform ocean area and density; gravitational fingerprints and land motion are omitted.

.. code-block:: python

    water_mass / (ocean_area * water_density)

greenhouse_surface_temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.greenhouse_surface_temperature`

Output: **K**. Blackbody surface, a grey isothermal layer transparent to sunlight; the layer is colder by 2**0.25.

.. code-block:: python

    (2 * absorbed_flux / ((2 - layer_emissivity) * SIGMA_SB))**0.25

daily_insolation
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.daily_insolation`

Output: **W/m^2**. Daily mean at the top of the atmosphere; declination fixed over the day; latitude and declination in radians.

.. code-block:: python

    h0 = arccos(clip(-tan(latitude) * tan(declination), -1, 1))
    solar_constant * distance_ratio**2 / pi * (h0 * sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * sin(h0))

radiative_convective_temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.radiative_convective_temperature`

Output: **K**. Fixed tropospheric lapse rate through the emission level under an isothermal stratosphere at the skin temperature.

.. code-block:: python

    np.maximum(effective_temperature + lapse_rate * (emission_height - height), effective_temperature / 2**0.25)

equilibrium_solar_constant
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.equilibrium_solar_constant`

Output: **W/m^2**. Zero-dimensional greybody balance with the ice_albedo closure; branches with a negative slope are unstable.

.. code-block:: python

    4 * emissivity * SIGMA_SB * temperature**4 / (1 - ice_albedo(temperature, cold_albedo, warm_albedo, transition_temperature, width))

red_noise_spectrum
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.red_noise_spectrum`

Output: **K^2/Hz**. Two-sided spectrum of a linear one-box model under white-noise forcing; frequency in cycles per second.

.. code-block:: python

    forcing_spectrum / (restoring**2 + (2 * np.pi * frequency * areal_heat_capacity)**2)

methane_forcing
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.methane_forcing`

Output: **W/m^2**. Myhre et al. (1998) fit, concentrations in ppb, valid near present-day values.

.. code-block:: python

    f(M, N) = 0.47 * log(1 + 2.01e-5 * (M * N)**0.75 + 5.31e-15 * M * (M * N)**1.52)
    0.036 * (sqrt(M) - sqrt(M0)) - (f(M, N0) - f(M0, N0))

nitrous_oxide_forcing
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.nitrous_oxide_forcing`

Output: **W/m^2**. Myhre et al. (1998) fit, concentrations in ppb, valid near present-day values.

.. code-block:: python

    0.12 * (sqrt(N) - sqrt(N0)) - (f(M0, N) - f(M0, N0))

gregory_regression
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.gregory_regression`

Output: **[W/m^2, W/(m^2 K)]**. Least-squares line through (temperature_change, imbalance) vectors; assumes constant forcing and feedback.

.. code-block:: python

    imbalance = forcing - restoring * temperature_change

two_box_timescales
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.two_box_timescales`

Output: **s**, the fast and the slow timescale on the last axis. Held et al. (2010) linear two-box model; minus the reciprocal eigenvalues of the system.

.. code-block:: python

    b = (feedback + exchange) / C + exchange / C0
    determinant = feedback * exchange / (C * C0)
    fast = (b + np.sqrt(b**2 - 4 * determinant)) / 2
    [1 / fast, fast / determinant]

two_box_fast_warming
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.two_box_fast_warming`

Output: **K**. The surface warming after the fast adjustment with the deep ocean still cold; approximates the transient climate response.

.. code-block:: python

    forcing / (feedback + exchange)

insolation_distribution
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.insolation_distribution`

Output: **1**, with global mean 1. North's (1975) two-term Legendre fit to the annual mean.

.. code-block:: python

    1 + second_legendre_coefficient * (3 * sine_latitude**2 - 1) / 2

north_temperature_profile
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.north_temperature_profile`

Output: **K**. North (1975) with a fixed albedo; the second Legendre polynomial is an eigenfunction of the diffusion operator.

.. code-block:: python

    Q = solar_constant / 4
    273.15 + (Q * (1 - albedo) - A) / B + Q * (1 - albedo) * s2 / (B + 6 * D) * (3 * x**2 - 1) / 2

Two-box response
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.simulate_two_box` integrates the surface and
deep anomalies under a constant forcing. Start from a warmed state with zero
forcing to repeat Held's abrupt return to preindustrial forcing.

.. code-block:: python

    C dT/dt = F - feedback * T - exchange * (T - T0)
    C0 dT0/dt = exchange * (T - T0)

Latitudinal energy balance
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.climate.simulate_latitudinal_energy_balance`
integrates the temperature of equal-area cells from the equator to the pole
with SciPy's BDF integrator. The diffusive flux vanishes at both ends, the
albedo is the smooth ``ice_albedo`` step between ice and an ice-free albedo
``warm_albedo + warm_albedo_gradient * P2(x)``, and the defaults give an ice
edge near 71 degrees for today's Sun. Pass the final state of one run as
``initial_temperature`` of the next to trace the hysteresis loop.

.. code-block:: python

    C dT/dt = (S / 4) s(x) (1 - albedo(T)) - (A + B (T - 273.15)) + D d/dx[(1 - x**2) dT/dx]

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The climate tests also check that one greenhouse layer balances the
top of the atmosphere, that daily insolation averages to a quarter of the solar
constant, that the radiative-convective profile passes through the effective
temperature, that every ice-albedo equilibrium closes the energy balance, that
the red-noise spectrum integrates to its variance, that the Myhre formulas
reproduce the IPCC values for 1750-1998, and that the Gregory regression
recovers the parameters of a one-box run. The two-box run matches the matrix
exponential and its timescales the eigenvalues, the fast plateau is
``F / (lambda + gamma)``, the latitudinal model without ice feedback converges
to North's two-mode solution, closes its global energy balance, and has a warm
and a snowball state for the same Sun. The shared suite adds cross-model identities, energy/conservation
checks and independent analytic solutions for dynamics. See :doc:`/validation`
for the exact evidence and :doc:`/numerical_methods` for numerical conventions.
