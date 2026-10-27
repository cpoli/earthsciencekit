Geophysics: units, equations and validation
===========================================

.. include:: /_generated/nav/geophysics.rst

Conventions and model limits
----------------------------

Latitude, longitude, inclination, declination and spherical distance are in
radians. Gravity and its anomalies are in m/s^2 (``units.convert`` turns them
into mGal, 1 mGal = 1e-5 m/s^2), magnetic fields in tesla (convert to nT), and
magnetization in A/m. Depth is positive downward. Field vectors use the local
north-east-down frame: ``dipole_field`` and ``zonal_field`` return the
northward and downward components on the last axis, and ``field_elements``
returns declination east of north, inclination positive downward and intensity.
A positive dipole moment is today's polarity, whose field points down in the
northern hemisphere; Gauss's ``g10`` is then negative. ``geomagnetic_field``
returns the northward, eastward and downward components from square tables
``g[n, m]`` and ``h[n, m]`` of Schmidt semi-normalized Gauss coefficients in
tesla. Entries with ``m > n``, the row ``n = 0`` and ``h[:, 0]`` must be zero.
Its latitude excludes the geographic poles, where the east component is
undefined.

The gravity reductions are first order and planar: ``free_air_correction`` is
``2 g h / R`` and ``bouguer_correction`` is an infinite slab, with no terrain
correction. ``sphere_gravity`` requires the observer outside the sphere.
``clairaut_gravity`` and ``flattening_from_j2`` are first order in flattening;
``normal_gravity`` is Somigliana's exact closed form with the GRS80 constants.
``homogeneous_flattening`` is the slow-rotation Maclaurin limit and accepts a
centrifugal ratio of at most 0.1. ``nodal_precession_rate`` is the secular J2
term only. ``stokes_function`` is the spherical kernel, singular at zero
distance, and ``degree_geoid_height`` is its eigenvalue form for degree 2 and
up.

The paleomagnetic helpers assume a geocentric axial dipole. ``paleolatitude``
expects a normal-polarity inclination; a reversed direction gives a southern
virtual geomagnetic pole. ``block_magnetic_anomaly`` is a two-dimensional
block with vertical magnetization in a vertical field, which is the anomaly
reduced to the pole; sum blocks to build a profile. The Rikitake dynamo is
nondimensional: ``time_scale`` converts its time unit to seconds, and its
currents are scaled. The Love numbers are those of Kelvin's uniform
incompressible elastic sphere. ``equivalent_water_height`` uses fully
normalized Stokes coefficients and the load Love number of the same degree.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments.
Harmonic degrees must be integers. No unit metadata are carried in raw numbers.
Bounds protect mathematical domains and obvious physical contradictions, but do
not establish model applicability.

Equation and output map
-----------------------

.. rubric:: Gravity

sphere_gravity
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.sphere_gravity`

Output: **m/s^2**. Observer outside the sphere, which attracts as a point mass; last axis holds the horizontal (toward +x) and downward components.

.. code-block:: python

    M = 4 / 3 * np.pi * radius**3 * density_contrast
    G * M * [-x, depth] / (x**2 + depth**2) ** 1.5

bouguer_correction
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.bouguer_correction`

Output: **m/s^2**. Flat, uniform, laterally infinite slab; no terrain or curvature correction.

.. code-block:: python

    2 * np.pi * G * density * thickness

free_air_correction
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.free_air_correction`

Output: **m/s^2**. First order in height / radius for a spherical point-mass Earth; no rock between station and reference.

.. code-block:: python

    2 * gravity * height / radius

bouguer_anomaly
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.bouguer_anomaly`

Output: **m/s^2**. Free-air and infinite-slab reductions only, without terrain correction.

.. code-block:: python

    observed_gravity - normal_gravity + 2 * normal_gravity * height / radius - 2 * np.pi * G * density * height

torsion_balance_gravitational_constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.torsion_balance_gravitational_constant`

Output: **m^3/(kg s^2)**. Point masses, a massless rod and attraction by the near large ball only.

.. code-block:: python

    2 * np.pi**2 * arm_length * separation**2 * deflection / (attracting_mass * period**2)

earth_mean_density
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.earth_mean_density`

Output: **kg/m^3**. Spherically symmetric, nonrotating Earth; gravitational_constant defaults to constants.G.

.. code-block:: python

    3 * gravity / (4 * np.pi * gravitational_constant * radius)

.. rubric:: Figure of the Earth

centrifugal_ratio
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.centrifugal_ratio`

Output: **1**. Equatorial ratio with surface gravity; the geodetic m uses GM / (a b) instead.

.. code-block:: python

    angular_speed**2 * radius / gravity

homogeneous_flattening
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.homogeneous_flattening`

Output: **1**. First order in rotation (Maclaurin spheroid); centrifugal_ratio at most 0.1.

.. code-block:: python

    5 * centrifugal_ratio / 4

clairaut_gravity
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.clairaut_gravity`

Output: **m/s^2**. First order in flattening and rotation for a hydrostatic spheroid; no sin^4 term.

.. code-block:: python

    equatorial_gravity * (1 + (5 * centrifugal_ratio / 2 - flattening) * np.sin(latitude) ** 2)

normal_gravity
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.normal_gravity`

Output: **m/s^2**. Somigliana's closed form on the GRS80 ellipsoid surface.

.. code-block:: python

    gamma_e, k, e2 = 9.7803267715, 0.001931851353, 0.00669438002290
    gamma_e * (1 + k * np.sin(latitude) ** 2) / np.sqrt(1 - e2 * np.sin(latitude) ** 2)

flattening_from_j2
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.flattening_from_j2`

Output: **1**. First order in flattening; second-order terms change 1/f by about 0.1%.

.. code-block:: python

    3 * j2 / 2 + centrifugal_ratio / 2

nodal_precession_rate
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.nodal_precession_rate`

Output: **rad/s**. First-order secular J2 perturbation without drag or higher harmonics; perigee above the reference radius.

.. code-block:: python

    n = np.sqrt(gravitational_parameter / semi_major_axis**3)
    -1.5 * n * j2 * (radius / semi_major_axis) ** 2 * np.cos(inclination) / (1 - eccentricity**2) ** 2

stokes_function
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.stokes_function`

Output: **1**. Spherical approximation with no masses outside the geoid; singular at zero distance.

.. code-block:: python

    s = np.sin(psi / 2)
    1 / s - 6 * s + 1 - 5 * np.cos(psi) - 3 * np.cos(psi) * np.log(s + s**2)

degree_geoid_height
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.degree_geoid_height`

Output: **m**. Spherical approximation of Stokes's formula for degrees of at least 2.

.. code-block:: python

    radius * anomaly_coefficient / ((degree - 1) * gravity)

.. rubric:: Geomagnetism and paleomagnetism

dipole_field
~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.dipole_field`

Output: **T**. Pure axial dipole outside its sources; last axis holds the northward and downward components.

.. code-block:: python

    MU0 * moment / (4 * np.pi * distance**3) * [np.cos(latitude), 2 * np.sin(latitude)]

dipole_inclination
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.dipole_inclination`

Output: **rad**. Geocentric axial dipole; positive downward for today's polarity.

.. code-block:: python

    np.arctan2(2 * np.sin(latitude), np.cos(latitude))

paleolatitude
~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.paleolatitude`

Output: **rad**. Geocentric axial dipole averaged over secular variation.

.. code-block:: python

    np.arctan2(np.sin(inclination), 2 * np.cos(inclination))

zonal_field
~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.zonal_field`

Output: **T**. Internal zonal sources only; last axis holds the northward and downward components.

.. code-block:: python

    x, n, a, r = np.sin(latitude), degree, reference_radius, distance
    (a / r) ** (n + 2) * coefficient * [lpmv(1, n, x), -(n + 1) * lpmv(0, n, x)]

schmidt_legendre
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.schmidt_legendre`

Output: **1**. Schmidt semi-normalized, without the Condon-Shortley phase, as in the IGRF.

.. code-block:: python

    n, m = degree, order
    np.sqrt((2 - (m == 0)) * factorial(n - m) / factorial(n + m)) * (-1)**m * lpmv(m, n, np.sin(latitude))

geomagnetic_field
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.geomagnetic_field`

Output: **T**, north, east and down on the last axis. Internal sources on a spherical Earth; ``S`` is ``schmidt_legendre`` of the colatitude ``theta``, and the ``m = 0`` terms equal ``zonal_field``.

.. code-block:: python

    V = a * sum((a / r)**(n + 1) * (g[n, m] * cos(m lon) + h[n, m] * sin(m lon)) * S)
    north, east, down = sum((a / r)**(n + 2) * [(g cos + h sin) dS/dtheta,
                                                m (g sin - h cos) S / sin(theta),
                                                -(n + 1) (g cos + h sin) S])

gauss_dipole_moment
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.gauss_dipole_moment`

Output: **A m^2**. Magnitude of the best-fitting geocentric tilted dipole.

.. code-block:: python

    4 * np.pi * reference_radius**3 / MU0 * np.sqrt(g10**2 + g11**2 + h11**2)

field_elements
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.field_elements`

Output: **rad, rad, T**. Local north-east-down frame; the field vector must be nonzero.

.. code-block:: python

    [np.arctan2(east, north), np.arctan2(down, np.hypot(north, east)), np.sqrt(north**2 + east**2 + down**2)]

direction_angle
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.direction_angle`

Output: **rad**. Great-circle distance between unit vectors; intensities are ignored.

.. code-block:: python

    np.arccos(np.sin(Ia) * np.sin(Ib) + np.cos(Ia) * np.cos(Ib) * np.cos(Da - Db))

virtual_geomagnetic_pole
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.virtual_geomagnetic_pole`

Output: **rad**. Geocentric dipole pole; last axis holds pole latitude and longitude in [0, 2 pi).

.. code-block:: python

    p = np.arctan2(2 * np.cos(I), np.sin(I))
    pole = site * np.cos(p) + (north * np.cos(D) + east * np.sin(D)) * np.sin(p)
    [np.arctan2(pole_z, np.hypot(pole_x, pole_y)), np.arctan2(pole_y, pole_x) % (2 * np.pi)]

block_magnetic_anomaly
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.block_magnetic_anomaly`

Output: **T**. Infinite strike, vertical magnetization and ambient field (reduced to the pole); left < right and top < bottom.

.. code-block:: python

    phi = lambda z: np.arctan((right - x) / z) - np.arctan((left - x) / z)
    MU0 * magnetization / (2 * np.pi) * (phi(top) - phi(bottom))

.. rubric:: Geodynamo

magnetic_reynolds_number
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.magnetic_reynolds_number`

Output: **1**. Order-of-magnitude scale; dynamos need values above a threshold of order 10 to 100.

.. code-block:: python

    MU0 * conductivity * speed * length

magnetic_decay_time
~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.magnetic_decay_time`

Output: **s**. Uniform stationary conductor surrounded by an insulator.

.. code-block:: python

    MU0 * conductivity * radius**2 / np.pi**2

rikitake_fixed_point
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.rikitake_fixed_point`

Output: **1**. Positive state of the nondimensional Rikitake system; the other is (-x, -y, z).

.. code-block:: python

    [k, 1 / k, mu * k**2]

.. rubric:: Deformation and mass change

love_number_h2
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.love_number_h2`

Output: **1**. Uniform incompressible self-gravitating elastic sphere (Kelvin's solution).

.. code-block:: python

    2.5 / (1 + 19 * rigidity / (2 * density * gravity * radius))

love_number_k2
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.love_number_k2`

Output: **1**. Uniform incompressible self-gravitating elastic sphere; the fluid limit is 3/2.

.. code-block:: python

    1.5 / (1 + 19 * rigidity / (2 * density * gravity * radius))

tidal_gravimetric_factor
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.tidal_gravimetric_factor`

Output: **1**. Degree-2 tide on a spherically symmetric elastic Earth without ocean loading.

.. code-block:: python

    1 + h2 - 1.5 * k2

equivalent_water_height
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.equivalent_water_height`

Output: **m**. Thin surface layer on an elastic Earth (Wahr et al., 1998), fully normalized coefficients.

.. code-block:: python

    radius * mean_density * (2 * degree + 1) * coefficient_change / (3 * water_density * (1 + load_love_number))

Rikitake dynamo
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.geophysics.simulate_rikitake_dynamo` integrates
``dx/dt = -mu x + y z``, ``dy/dt = -mu y + (z - A) x`` and ``dz/dt = 1 - x y``
with ``A = mu (k**2 - 1 / k**2)``, in units of ``time_scale`` seconds, and
returns the scaled currents ``x`` and ``y`` and the angular speed ``z``. The
steady states are ``rikitake_fixed_point`` and its mirror ``(-x, -y, z)``;
the sign of ``x`` is the polarity.

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting check and invalid-input
contracts. The domain suite adds physical checks: the Bouguer anomaly recovers a
buried body from gravity on topography, Clairaut's formula matches Somigliana's
normal gravity, the numerical Stokes integral of zonal anomalies gives the
``1 / (n - 1)`` law, the degree-1 zonal field equals the axial dipole, the
Schmidt functions have the IGRF mean square ``1 / (2 n + 1)``, the full field
is minus the numerical gradient of its potential and reduces to the zonal
field and the axial dipole, the
virtual geomagnetic pole of a tilted dipole returns its pole, an infinite
magnetized layer has no external field, and the Rikitake fixed points are
steady while the free flow reverses. See :doc:`/validation` for the exact
evidence and :doc:`/numerical_methods` for numerical conventions.
