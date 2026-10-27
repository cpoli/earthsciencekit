Volcanology: units, equations and validation
============================================

.. include:: /_generated/nav/volcanology.rst

Conventions and model limits
----------------------------

Ages, half-lives and travel times are seconds and decay constants 1/s; convert
with ``units.convert(value, "Ma", "s")`` (``kyr``, ``Ma`` and ``Ga`` are
available). The uranium decay constants ``LAMBDA_U238`` and ``LAMBDA_U235``
(Jaffey et al. 1971), the present ``U238_U235 = 137.88`` and Libby's mean life
``LIBBY_MEAN_LIFE`` (5568-year half-life) live in :mod:`earthsciencekit.constants`.
Every dating model assumes a closed system: no parent or daughter gained or
lost after time zero, except where the concordia chord models episodic lead
loss. ``radiometric_age`` and ``radiogenic_daughter`` assume no initial
daughter; ``isochron_fit`` removes that assumption with an unweighted
least-squares line and requires a positive slope. ``radiocarbon_age`` is a
conventional radiocarbon age, not a calibrated calendar age. ``lead_lead_age``
solves for ages between one year and 10 Ga. ``concordia_ratios`` returns
``[207Pb*/235U, 206Pb*/238U]`` on the last axis, the horizontal and vertical
coordinates of Wetherill's diagram.

Isotope ratios are heavy over light; delta values are per mil relative to a
caller-supplied standard. Fractionation factors are product over reservoir.
``isotope_temperature`` is the Epstein et al. (1953) calibration, both delta
values on the same scale, valid near 7-30 degC. Crystallization concentrations
carry the unit of the initial concentration, and partition coefficients are
bulk solid over melt. ``lattice_strain_partition`` takes radii in metres and
the apparent Young's modulus of the site in pascals.

Viscosities are Pa s for a crystal- and bubble-free Newtonian melt. Shaw's
slope ``s`` is computed from the molar composition (Shaw 1972, Table 1) and is
supplied by the caller; his law is Arrhenius with activation energy
``s * 1e4 * R_GAS``. Water solubility is a mass fraction with a square-root
pressure law and no CO2. Dike flow is laminar plane Poiseuille flow with a
positive density contrast meaning buoyant magma. The Mogi source returns
``[radial, vertical]`` displacements on the last axis, outward and upward
positive, for a source small compared with its depth.

Bubble growth treats water vapour as the only volatile. The supersaturation
is ``melt_density * (initial_water - surface_water) / gas_density``, a ratio of
volumes. Scriven's growth constant ``beta`` gives ``R = 2 beta sqrt(D t)`` for
growth from zero radius in an unbounded melt. It needs the supersaturation
below ``1 / density_ratio``, where ``density_ratio`` is the gas density over
the melt density. ``simulate_bubble_growth`` follows one bubble in its own
shell of melt at a fixed pressure after a sudden decompression. The wall stays
at ``water_solubility`` and the gas density is ``pressure / (R_VAPOR *
temperature)``.

Plume buoyancy flux is ``reduced_gravity * velocity * radius**2`` without the
factor pi, as in Morton, Taylor and Turner (1956) and Briggs (1969).
``column_height`` takes the heat flux in watts. The Volcanic Explosivity Index
uses the bulk tephra volume only.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments, except
for the vector inputs of ``isochron_fit``. No unit metadata are carried in raw
numbers. Bounds protect mathematical domains and obvious physical
contradictions, but do not establish model applicability.

Equation and output map
-----------------------

rayleigh_fractionation
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.rayleigh_fractionation`

Output: **1**. Rayleigh (1896): the product leaves as soon as it forms, the reservoir stays well mixed and the factor is constant.

.. code-block:: python

    initial_ratio * remaining_fraction**(fractionation_factor - 1)

delta_value
~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.delta_value`

Output: **per mil**. Linear in the ratio; delta values of different standards must be converted before they are compared.

.. code-block:: python

    (ratio / standard_ratio - 1) * 1000

decay_constant
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.decay_constant`

Output: **1/s**. Rutherford and Soddy (1902): first-order decay, independent of temperature, pressure and chemistry.

.. code-block:: python

    np.log(2) / half_life

parent_fraction
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.parent_fraction`

Output: **1**. Rutherford and Soddy (1902): a closed system with no parent added after time zero.

.. code-block:: python

    np.exp(-decay_constant * time)

radiometric_age
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.radiometric_age`

Output: **s**. Closed system with no initial daughter; Boltwood (1907) applied it to lead in uranium minerals.

.. code-block:: python

    np.log1p(daughter_parent_ratio) / decay_constant

radiogenic_daughter
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.radiogenic_daughter`

Output: **unit of parent**. Closed system with no initial daughter; the inverse of radiometric_age. Holmes (1913) used it to build the first radiometric time scale.

.. code-block:: python

    parent * np.expm1(decay_constant * time)

fractional_crystallization_melt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.fractional_crystallization_melt`

Output: **unit of initial_concentration**. Rayleigh fractional crystallization (Bowen 1928; Gast 1968): crystals do not re-equilibrate and the bulk partition coefficient is constant.

.. code-block:: python

    initial_concentration * melt_fraction**(partition_coefficient - 1)

equilibrium_crystallization_melt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.equilibrium_crystallization_melt`

Output: **unit of initial_concentration**. Batch mass balance (Shaw 1970): solid and melt stay in equilibrium, so an incompatible element is enriched at most by 1 / partition_coefficient.

.. code-block:: python

    initial_concentration / (partition_coefficient + melt_fraction * (1 - partition_coefficient))

lattice_strain_partition
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.lattice_strain_partition`

Output: **1**. Blundy and Wood (1994) quantification of Goldschmidt (1937) for isovalent ions on one site; charge differences need a separate parabola.

.. code-block:: python

    peak_partition * np.exp(-4 * np.pi * youngs_modulus * AVOGADRO * (optimal_radius / 2 * (ionic_radius - optimal_radius)**2 + (ionic_radius - optimal_radius)**3 / 3) / (R_GAS * temperature))

isotope_temperature
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.isotope_temperature`

Output: **K**. Urey (1947) predicted the effect; this is the empirical calibration of Epstein, Buchsbaum, Lowenstam and Urey (1953), valid for about 7-30 degC and equilibrium calcite.

.. code-block:: python

    273.15 + 16.5 - 4.3 * d + 0.14 * d**2

radiocarbon_age
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.radiocarbon_age`

Output: **s**. Libby et al. (1949) with the conventional 5568-year (Libby) half-life, mean life LIBBY_MEAN_LIFE; it is a radiocarbon age, not a calendar age, until calibrated.

.. code-block:: python

    -LIBBY_MEAN_LIFE * np.log(fraction_modern)

radiogenic_lead_ratio
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.radiogenic_lead_ratio`

Output: **1**. Patterson (1956): closed system with today's 238U/235U = U238_U235 and the Jaffey et al. (1971) decay constants. The ratio tends to LAMBDA_U235 / (U238_U235 * LAMBDA_U238) at zero age.

.. code-block:: python

    np.expm1(LAMBDA_U235 * time) / (U238_U235 * np.expm1(LAMBDA_U238 * time))

lead_lead_age
~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.lead_lead_age`

Output: **s**. Inverts radiogenic_lead_ratio by Brent's method between one year and 10 Ga; same closed-system assumptions.

.. code-block:: python

    radiogenic_lead_ratio(age) = lead_ratio

concordia_ratios
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.concordia_ratios`

Output: **[1, 1]**. Wetherill (1956): closed system with the Jaffey et al. (1971) decay constants. Episodic lead loss puts samples on a straight chord (discordia) between the concordia points of formation and loss.

.. code-block:: python

    [np.expm1(LAMBDA_U235 * time), np.expm1(LAMBDA_U238 * time)]

plume_rise_height
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.plume_rise_height`

Output: **m**. Morton, Taylor and Turner (1956) for a steady point source of buoyancy in uniform stratification, with Briggs's constant 5.0 (entrainment about 0.1); source momentum, wind and moisture are neglected.

.. code-block:: python

    5.0 * buoyancy_flux**0.25 * buoyancy_frequency**-0.75

mogi_displacement
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.mogi_displacement`

Output: **[m, m]**. Mogi (1958) point source in an elastic half-space; the chamber radius must be small compared with its depth. The uplift halves at ``radial_distance = 0.766 * depth``.

.. code-block:: python

    (1 - poisson_ratio) * volume_change / np.pi * [radial_distance, depth] / (radial_distance**2 + depth**2)**1.5

isochron_fit
~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.isochron_fit`

Output: **[s, 1]**. Nicolaysen (1961): every sample formed at the same time with the same initial ratio and stayed closed; ordinary least squares without measurement-error weighting (a York fit would weight them). The slope must be positive.

.. code-block:: python

    daughter_ratio = initial_ratio + parent_ratio * np.expm1(decay_constant * age)

arrhenius_viscosity
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.arrhenius_viscosity`

Output: **Pa s**. Constant activation energy over the range of interest; crystal- and bubble-free Newtonian melt.

.. code-block:: python

    reference_viscosity * np.exp(activation_energy / R_GAS * (1 / temperature - 1 / reference_temperature))

shaw_viscosity
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.shaw_viscosity`

Output: **Pa s**. Shaw (1972), Am. J. Sci. 272: 870-893; all lines pass near ln(eta [poise]) = -6.4 at 1e4/T = 1.5. Crystal- and bubble-free Newtonian melt; about 0.2 log units of accuracy.

.. code-block:: python

    0.1 * np.exp(slope * (1e4 / temperature) - 1.5 * slope - 6.4)

column_height
~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.column_height`

Output: **m**. Wilson, Sparks, Huang and Watkins (1978), JGR 83: 1829-1836, for a standard atmosphere and a strong steady column; wind bends weak plumes lower.

.. code-block:: python

    8.2 * heat_flux**0.25

dike_ascent_velocity
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.dike_ascent_velocity`

Output: **m/s**. Wilson and Head (1981): plane Poiseuille flow between parallel walls driven by buoyancy alone; check that the Reynolds number ``magma_density * velocity * width / viscosity`` stays below about 2000.

.. code-block:: python

    width**2 * density_contrast * gravity / (12 * viscosity)

water_solubility
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.water_solubility`

Output: **1 (mass fraction)**. Square-root (Henry-type) law for pure water vapour, valid up to a few hundred MPa; ignores CO2 and temperature.

.. code-block:: python

    solubility_constant * np.sqrt(pressure)

diffusive_supersaturation
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.diffusive_supersaturation`

Output: **1**. The mass-transfer analogue of the Jakob number.

.. code-block:: python

    melt_density * (initial_water - surface_water) / gas_density

scriven_growth_constant
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.scriven_growth_constant`

Output: **1**. Scriven (1959), with the radial flow driven by the expanding bubble (``e = 1 - density_ratio``), solved by root finding; ``np.sqrt(phi / 2)`` for small supersaturation ``phi`` and ``np.sqrt(3 / np.pi) * phi`` for a light gas at large ``phi``.

.. code-block:: python

    phi = 2 * beta**3 * np.exp(beta**2 + 2 * e * beta**2) * integral(x**-2 * np.exp(-x**2 - 2 * e * beta**3 / x), beta, inf)

scriven_bubble_radius
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.scriven_bubble_radius`

Output: **m**. Growth from zero radius in an unbounded melt at constant supersaturation.

.. code-block:: python

    2 * growth_constant * np.sqrt(diffusivity * time)

volcanic_explosivity_index
~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.volcanic_explosivity_index`

Output: **1**. Newhall and Self (1982), volume criterion only: VEI 0 below 1e4 m^3, 1 below 1e6 m^3, then one step per factor of ten up to VEI 8 at 1e12 m^3 and above. Column height and other criteria are not used.

.. code-block:: python

    np.where(erupted_volume < 1e4, 0, np.where(erupted_volume < 1e6, 1, np.clip(np.floor(np.log10(erupted_volume)) - 4, 2, 8)))

Eruption column
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.simulate_eruption_column` integrates the
top-hat plume fluxes ``Q = b**2 w``, ``M = b**2 w**2`` and ``F = b**2 w g'``
along the travel time of a parcel on the axis:
``dz/dt = M / Q``, ``dQ/dt = 2 alpha abs(M)**1.5 / Q``, ``dM/dt = F`` and
``dF/dt = -N**2 M``. It returns the parcel height, vertical velocity and
reduced gravity. With ``N = 0`` it follows the pure-plume similarity solution;
in a uniform atmosphere its top lies near ``plume_rise_height`` for the default
entrainment 0.12. Past the top it continues as a fountain that falls back and
oscillates about the neutral buoyancy level, an illustrative closure.

Bubble growth
~~~~~~~~~~~~~

:py:func:`earthsciencekit.volcanology.simulate_bubble_growth` integrates the
radius of one bubble in a melt shell of fixed mass. The growth is quasi-static
diffusion driven by the shell's mean dissolved water ``c``, which falls as the
bubble fills. It returns the radius, ``c`` and the vesicularity. Without
depletion it gives ``R**2 = R0**2 + 2 phi D t``, Scriven's small-supersaturation
limit, which underestimates the full solution by about 15 % at ``phi = 0.05``.
The bubble stops growing when ``c`` reaches the solubility.

.. code-block:: python

    d(gas_density * R**3)/dt = 3 * R * melt_density * D * (c - c_s)

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting or bound checks and
invalid-input contracts. The volcanology tests also check that parent and
daughter atoms are conserved and that the age equation inverts accumulation,
that Libby's mean life gives a 5568-year half-life, that the lead ratio tends
to its zero-age limit and ``lead_lead_age`` inverts it, that both concordia
clocks return the same age and lead loss falls on a straight chord, that the
isochron recovers its age and initial ratio, that fractional crystallization
conserves mass and shares its initial slope with equilibrium crystallization,
that lattice strain peaks at the optimal radius, that Shaw's lines share their
pivot, that the Mogi uplift integrates to ``2 (1 - nu) dV``, and that the
eruption column follows the pure-plume similarity solution and reaches the
Morton-Taylor-Turner height. Scriven's growth constant matches its closed form
without advection and the small- and large-supersaturation limits, and the
bubble cell follows the quasi-static law, conserves water and saturates the
melt. See :doc:`/validation` for the exact evidence and
:doc:`/numerical_methods` for numerical conventions.
