Plate Tectonics: units, equations and validation
================================================

.. include:: /_generated/nav/plate_tectonics.rst

Conventions and model limits
----------------------------

Rates use m/s and times use seconds. Convert cm/year explicitly with
``units.convert``. Euler-vector components and positions are Earth-fixed
Cartesian coordinates; their cross product gives Cartesian velocity, not local
north/east components. Scalar convergence is positive toward the supplied
boundary normal. Half-spreading rate is one-sided.

Cooling profiles use depth positive downward, absolute temperatures and
positive diffusivity. Heat flux is positive upward. The half-space solution
holds the surface fixed while periodic diffusion conserves the mean; those
boundary conditions describe different experiments. Flexural rigidity uses an
ideal stable isotropic elastic plate with -1 < nu < 0.5. Flexural deflection
and ocean depth are positive downward. Hotspot track length is measured along
the small circle about the Euler pole, and ``pole_distance`` is the angular
distance of the hotspot from that pole. Euler vectors passed to
``compose_euler_vectors`` carry their three Cartesian components on the last
axis; ``omega_ab`` is the motion of plate A relative to plate B.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments. No unit
metadata are carried in raw numbers. Bounds protect mathematical domains and
obvious physical contradictions, but do not establish model applicability.

Equation and output map
-----------------------

angular_speed
~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.angular_speed`

Output: **rad/s**. Constant angular rate; angular input is radians and output radians per second.

.. code-block:: python

    rotation_angle / time

plate_speed
~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.plate_speed`

Output: **m/s**. Spherical rigid plate, magnitude only; angular distance must be at most pi.

.. code-block:: python

    angular_speed * radius * np.sin(angular_distance)

euler_velocity_x
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.euler_velocity_x`

Output: **m/s**. Earth-fixed Cartesian axes; x is not the local east component.

.. code-block:: python

    omega_y * position_z - omega_z * position_y

euler_velocity_y
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.euler_velocity_y`

Output: **m/s**. Same Earth-fixed Cartesian convention as euler_velocity_x.

.. code-block:: python

    omega_z * position_x - omega_x * position_z

euler_velocity_z
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.euler_velocity_z`

Output: **m/s**. Earth-fixed components, not local vertical motion or deformation.

.. code-block:: python

    omega_x * position_y - omega_y * position_x

plate_displacement
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.plate_displacement`

Output: **m**. Constant signed speed along a local coordinate; not a finite Euler rotation.

.. code-block:: python

    speed * time

spreading_age
~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.spreading_age`

Output: **s**. Distance on one side and half-spreading rate; changes in rate are omitted.

.. code-block:: python

    distance / half_spreading_rate

full_spreading_rate
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.full_spreading_rate`

Output: **m/s**. Symmetric spreading; asymmetric velocities require their own relative-vector calculation.

.. code-block:: python

    2 * half_spreading_rate

convergence_rate
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.convergence_rate`

Output: **m/s**. Positive normal component denotes convergence; negative denotes divergence.

.. code-block:: python

    relative_speed * np.cos(obliquity)

fault_slip_rate
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.fault_slip_rate`

Output: **m/s**. Supplied obliquity is oriented from the positive convergence normal.

.. code-block:: python

    relative_speed * np.sin(obliquity)

subduction_travel_time
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.subduction_travel_time`

Output: **s**. Straight slab at constant speed; zero dip cannot reach depth and is rejected.

.. code-block:: python

    vertical_depth / (slab_speed * np.sin(dip))

slab_pull_force
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.slab_pull_force`

Output: **N**. Simple buoyancy force magnitude; drag, phase changes and mantle flow are excluded.

.. code-block:: python

    density_excess * slab_volume * gravity

airy_root
~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.airy_root`

Output: **m**. Local compensation without plate flexure; mantle must be denser than crust.

.. code-block:: python

    elevation * crust_density / (mantle_density - crust_density)

flexural_rigidity
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.flexural_rigidity`

Output: **N m**. Thin isotropic elastic plate with stable Poisson ratio, not a yielding viscoelastic lithosphere.

.. code-block:: python

    young_modulus * elastic_thickness**3 / (12 * (1 - poisson_ratio**2))

flexural_parameter
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.flexural_parameter`

Output: **m**. Infinite beam-like plate convention alpha to the fourth equals 4D over density contrast times g.

.. code-block:: python

    (4 * rigidity / (density_contrast * gravity))**0.25

halfspace_temperature
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.halfspace_temperature`

Output: **K**. Semi-infinite conduction with fixed surface temperature; no hydrothermal or advective heat transfer.

.. code-block:: python

    surface_temperature + (mantle_temperature - surface_temperature) * erf(depth / (2 * np.sqrt(diffusivity * age)))

halfspace_heat_flux
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.halfspace_heat_flux`

Output: **W/m^2**. Positive upward flux; semi-infinite conductive model excludes the zero-age singularity.

.. code-block:: python

    conductivity * (mantle_temperature - surface_temperature) / np.sqrt(np.pi * diffusivity * age)

mantle_rayleigh_number
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.mantle_rayleigh_number`

Output: **1**. Layer heated from below with uniform properties; convection starts above about 657.5 (free slip) or 1708 (rigid walls).

.. code-block:: python

    thermal_expansivity * density * gravity * temperature_difference * thickness**3 / (diffusivity * viscosity)

flexural_deflection
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.flexural_deflection`

Output: **m**. Unbroken infinite thin elastic plate, deflection positive downward; density contrast is mantle minus infill.

.. code-block:: python

    alpha = (4 * rigidity / (density_contrast * gravity))**0.25
    x = abs(distance) / alpha
    line_load * alpha**3 / (8 * rigidity) * np.exp(-x) * (np.cos(x) + np.sin(x))

hotspot_age
~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.hotspot_age`

Output: **s**. Fixed hotspot and constant Euler rotation of the plate over it; track length is measured along the small circle.

.. code-block:: python

    track_length / (angular_speed * radius * np.sin(pole_distance))

ocean_depth
~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.ocean_depth`

Output: **m**. Isostatic half-space cooling, depth positive downward; it overestimates depth on sea floor older than about 70 million years.

.. code-block:: python

    ridge_depth + 2 * mantle_density * thermal_expansivity * (mantle_temperature - surface_temperature)
    / (mantle_density - water_density) * np.sqrt(diffusivity * age / np.pi)

compose_euler_vectors
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.plate_tectonics.compose_euler_vectors`

Output: **rad/s**. Instantaneous rotations only; inputs and output have three Cartesian components on the last axis.

.. code-block:: python

    omega_ab + omega_bc

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The shared suite adds cross-model identities, energy/conservation
checks and independent analytic solutions for dynamics. See :doc:`/validation`
for the exact evidence and :doc:`/numerical_methods` for numerical conventions.
