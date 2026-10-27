Seismology: units, equations and validation
===========================================

.. include:: /_generated/nav/seismology.rst

Conventions and model limits
----------------------------

Moduli are Pa, density kg/m^3, distances m and wave speeds m/s. Reflection and
transmission coefficients refer to acoustic pressure amplitudes at normal
incidence. Transmitted energy includes an impedance ratio; pressure amplitudes
alone need not sum in quadrature to one. Oblique elastic mode conversion is
outside this model.

Moment is N m and magnitude uses (2/3)(log10(M0)-9.1). Switching to dyne cm
without changing the intercept is an error. Ricker time is centred at zero and
frequency is Hz. Constant-Q attenuation excludes geometrical spreading.
Oscillator Q gives damping ratio 1/(2Q); its natural frequency is in Hz.

Layered models are flat (head waves, Love waves) or spherically symmetric
(ray parameter, in s/rad with incidence measured from the radius). Local
magnitude takes Wood-Anderson amplitude in m and hypocentral distance in m.
Vertical travel time integrates a tabulated profile exactly for speeds that
vary linearly between its depths. The layered ray tracer takes the same kind
of table, starting at the surface, but lets the speed follow the Mohorovicic
law ``v = a r**b`` between depths, which makes each layer's ray integrals exact.
``radius / speed`` must decrease with depth within every layer; it may jump up
at a discontinuity, as at the core-mantle boundary. Distances are surface to
surface, in radians, and are not folded past the antipode.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. No unit
metadata are carried in raw numbers. Bounds protect mathematical domains and
obvious physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

p_wave_speed
~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.p_wave_speed`

Output: **m/s**. Homogeneous isotropic small-strain medium; shear modulus may vanish in a fluid.

.. code-block:: python

    np.sqrt((bulk_modulus + 4 * shear_modulus / 3) / density)

s_wave_speed
~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.s_wave_speed`

Output: **m/s**. Homogeneous isotropic medium; zero shear modulus returns zero speed.

.. code-block:: python

    np.sqrt(shear_modulus / density)

poisson_ratio
~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.poisson_ratio`

Output: **1**. Requires positive bulk modulus: vp squared must exceed 4 vs squared over 3.

.. code-block:: python

    (p_speed**2 - 2 * s_speed**2) / (2 * (p_speed**2 - s_speed**2))

travel_time
~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.travel_time`

Output: **s**. Straight path in a homogeneous medium; not a ray tracer.

.. code-block:: python

    distance / speed

seismic_impedance
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.seismic_impedance`

Output: **kg/(m^2 s)**. Scalar acoustic impedance; elastic mode conversion is excluded.

.. code-block:: python

    density * speed

reflection_coefficient
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.reflection_coefficient`

Output: **1**. Pressure amplitude convention at normal incidence; coefficient can be negative.

.. code-block:: python

    (impedance2 - impedance1) / (impedance2 + impedance1)

transmission_coefficient
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.transmission_coefficient`

Output: **1**. Normal-incidence acoustic pressure amplitudes, not displacement amplitudes.

.. code-block:: python

    2 * impedance2 / (impedance2 + impedance1)

critical_angle
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.critical_angle`

Output: **rad**. Incidence from the slower medium; second speed must exceed first.

.. code-block:: python

    np.arcsin(speed1 / speed2)

seismic_moment
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.seismic_moment`

Output: **N m**. Scalar mean-slip estimate; no spatial rupture distribution.

.. code-block:: python

    rigidity * area * slip

moment_magnitude
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.moment_magnitude`

Output: **1**. Uses the 9.1 intercept with moment in N m, not dyne cm.

.. code-block:: python

    (2 / 3) * (np.log10(moment) - 9.1)

moment_from_magnitude
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.moment_from_magnitude`

Output: **N m**. Same SI convention as moment_magnitude; very large inputs can overflow.

.. code-block:: python

    10 ** (1.5 * magnitude + 9.1)

gutenberg_richter
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.gutenberg_richter`

Output: **count**. Expected cumulative count for the fitted catalogue interval; completeness matters.

.. code-block:: python

    10 ** (a - b * magnitude)

omori_rate
~~~~~~~~~~

:py:func:`earthsciencekit.seismology.omori_rate`

Output: **count/s**. Empirical fitted rate; no earthquake timing or location prediction.

.. code-block:: python

    productivity / (offset + time)**exponent

ricker_wavelet
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.ricker_wavelet`

Output: **1**. Unit peak at zero time; frequency is the peak-spectrum convention.

.. code-block:: python

    (1 - 2 * (np.pi * frequency * time)**2) * np.exp(-(np.pi * frequency * time)**2)

attenuated_amplitude
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.attenuated_amplitude`

Output: **same amplitude units as input**. Constant-Q narrowband attenuation; geometric spreading is separate.

.. code-block:: python

    amplitude * np.exp(-np.pi * frequency * time / quality_factor)

rayleigh_wave_speed
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.rayleigh_wave_speed`

Output: **m/s**. Homogeneous half-space with a free surface; requires positive bulk modulus.

.. code-block:: python

    s_speed * np.sqrt(xi)  # xi in (0, 1): xi**3 - 8 xi**2 + (24 - 16 gamma) xi - 16 (1 - gamma) = 0, gamma = (s_speed / p_speed)**2

head_wave_travel_time
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.head_wave_travel_time`

Output: **s**. Flat layer over a faster half-space; distance must reach the critical distance.

.. code-block:: python

    distance / speed2 + 2 * thickness * np.sqrt(speed2**2 - speed1**2) / (speed1 * speed2)

crossover_distance
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.crossover_distance`

Output: **m**. Flat layer over a faster half-space; direct and head-wave times are equal there.

.. code-block:: python

    2 * thickness * np.sqrt((speed2 + speed1) / (speed2 - speed1))

love_wave_angular_frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.love_wave_angular_frequency`

Output: **rad/s**. Fundamental mode of one layer over a half-space; phase speed between the two shear speeds.

.. code-block:: python

    phase_speed * np.arctan(halfspace_rigidity * s2 / (layer_rigidity * s1)) / (thickness * s1)
    # s1 = np.sqrt(phase_speed**2 / layer_s_speed**2 - 1), s2 = np.sqrt(1 - phase_speed**2 / halfspace_s_speed**2)

ray_parameter
~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.ray_parameter`

Output: **s/rad**. Spherically symmetric Earth; incidence angle from the radial direction, conserved along a ray.

.. code-block:: python

    radius * np.sin(incidence_angle) / speed

local_magnitude
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.local_magnitude`

Output: **1**. Wood-Anderson amplitude in m, hypocentral distance in m; Hutton-Boore correction for about 10-700 km.

.. code-block:: python

    np.log10(amplitude / 1e-3) + 1.11 * np.log10(distance / 1e5) + 1.89e-6 * (distance - 1e5) + 3

vertical_travel_time
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.vertical_travel_time`

Output: **s**. Tabulated profile, speed linear between points; a repeated depth marks a discontinuity.

.. code-block:: python

    cumsum(dz * np.log(v2 / v1) / (v2 - v1))  # dz / v1 where v2 == v1

layered_ray_distance
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.layered_ray_distance`

Output: **rad**. Mohorovicic-law layers; a ray turns where ``eta = r / v`` falls to the ray parameter, or at the top of a layer it cannot enter.

.. code-block:: python

    2 * sum(np.arccos(p / eta_top) - np.arccos(p / eta_bottom)) / (1 - b)  # eta_bottom = p in the turning layer

layered_ray_travel_time
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.seismology.layered_ray_travel_time`

Output: **s**. Same rays; the slope of the travel-time curve is the ray parameter.

.. code-block:: python

    2 * sum(np.sqrt(eta_top**2 - p**2) - np.sqrt(eta_bottom**2 - p**2)) / (1 - b)

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The layered rays are chords of a uniform sphere, agree with direct
quadrature of the ray integrals, obey ``dT/d(distance) = p`` and reproduce the
two-shell core shadow. The shared suite adds cross-model identities, energy/conservation
checks and independent analytic solutions for dynamics. See :doc:`/validation`
for the exact evidence and :doc:`/numerical_methods` for numerical conventions.
