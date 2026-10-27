Hydrology: units, equations and validation
==========================================

.. include:: /_generated/nav/hydrology.rst

Conventions and model limits
----------------------------

Lengths and depths are metres, times seconds and rates m/s; divide a rain
rate in mm/h by ``3.6e6`` and convert days with
``units.convert(value, "day", "s")``. The water balance accepts any
common depth rate, such as mm/day, used consistently. Evaporation rates are
kg/(m^2 s), which equals mm/s of liquid water; multiply by 86400 for mm/day.
Vapour pressures are pascals, and ``saturation_slope`` is the derivative of
:func:`earthsciencekit.atmosphere.saturation_vapor_pressure` with the same
constant latent heat ``LATENT_HEAT_VAPORIZATION = 2.5e6`` J/kg; the von
Karman constant ``VON_KARMAN = 0.41`` sets the neutral aerodynamic resistance.
Both live in :mod:`earthsciencekit.constants`.

Darcy fluxes are positive down the head gradient, and ``hydraulic_head``
adds the elevation to the gauge pressure head. Aquifers are homogeneous and
isotropic: Dupuit's water table has no recharge and a horizontal base, and
the Theis and Cooper-Jacob drawdowns assume a fully penetrating well of
negligible radius in an infinite confined aquifer. Cooper-Jacob is the
late-time limit, within 1% of Theis for ``u < 0.01``.

In the unsaturated zone heights and fluxes are positive upward and pressure
heads negative under suction. ``richards_steady_head`` solves the steady
Richards equation for Gardner's exponential conductivity above a water table
at height zero; it rejects downward fluxes at or beyond the saturated
conductivity and heights above the limit of capillary supply. Green-Ampt and
Horton rates are infiltration capacities of a uniform soil; Green-Ampt times
start when the surface ponds, and the suction head is a positive magnitude.
``exponential_water_content`` is the retention curve that shares Gardner's
exponent, and ``simulate_richards_column`` solves the transient Richards
equation in that soil between the water table and a surface with a prescribed
flux.

Channel formulas are for steady uniform flow in SI units: Chézy's
coefficient is in m^0.5/s and Manning's ``n`` in s/m^(1/3). The unit
hydrograph is linear and time-invariant: ``unit_hydrograph_convolution``
takes excess rainfall depths per interval and ordinates in m^3/s per metre of
excess at the same interval, and the Nash shape integrates to one per second.
``muskingum_coefficients`` returns ``[C0, C1, C2]`` on the last axis;
``muskingum_routing`` starts from a steady state and needs scalar parameters.
``simulate_kinematic_wave`` routes an inflow hydrograph, sampled on its time
grid, down a wide rectangular channel with Manning's rating. The kinematic
celerity takes the exponent ``m`` of the rating ``Q = a A**m``. Degree-day
factors are metres of water equivalent per kelvin per second, so
``3e-3 / 86400`` is 3 mm per degC per day. ``snowpack_water_equivalent``
returns the snow water equivalent and the water released in each step on the
last axis.
Stream orders follow Strahler (1952). TOPMODEL indices use the upslope area
per unit contour width in metres, and nonpositive local deficits mark
saturated ground.

Input contracts
---------------

Algebraic models accept finite real scalars or nonempty arrays. Booleans,
complex inputs, nonnumeric values and NaN/Inf are rejected. Physical bounds
apply element by element; NumPy broadcasting applies between arguments, except
for the vector inputs of ``unit_hydrograph_convolution`` and
``muskingum_routing``. No unit metadata are carried in raw numbers. Bounds
protect mathematical domains and obvious physical contradictions, but do not
establish model applicability.

Equation and output map
-----------------------

water_balance
~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.water_balance`

Output: **m/s**. Perrault (1674): a closed catchment with no groundwater exchange across its divide; any common depth rate (such as mm/day) may be used throughout.

.. code-block:: python

    precipitation - evapotranspiration - runoff

runoff_coefficient
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.runoff_coefficient`

Output: **1**. Both depths over the same area and period; long periods make the storage change negligible.

.. code-block:: python

    runoff / precipitation

hydraulic_radius
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.hydraulic_radius`

Output: **m**. Rectangular section; tends to the depth for a wide channel.

.. code-block:: python

    width * depth / (width + 2 * depth)

chezy_velocity
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.chezy_velocity`

Output: **m/s**. Chezy (1775): steady uniform turbulent flow; the coefficient is not a constant but depends on roughness and depth.

.. code-block:: python

    chezy_coefficient * np.sqrt(hydraulic_radius * slope)

dalton_evaporation
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.dalton_evaporation`

Output: **kg/(m^2 s)**. Dalton (1802), a bulk mass-transfer law; the surface temperature it needs is rarely measured.

.. code-block:: python

    transfer_coefficient * (saturation_vapor_pressure - vapor_pressure)

hydraulic_head
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.hydraulic_head`

Output: **m**. Slow flow with negligible velocity head and uniform fluid density.

.. code-block:: python

    elevation + pressure / (density * gravity)

darcy_flux
~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.darcy_flux`

Output: **m/s**. Darcy (1856): laminar flow with a pore Reynolds number below about 1-10 in a saturated isotropic medium.

.. code-block:: python

    -hydraulic_conductivity * head_gradient

seepage_velocity
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.seepage_velocity`

Output: **m/s**. Average linear velocity; dispersion spreads tracers about it.

.. code-block:: python

    darcy_flux / porosity

dupuit_head
~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.dupuit_head`

Output: **m**. Dupuit (1863): no recharge, a horizontal base, and slopes small enough that the head is uniform on each vertical.

.. code-block:: python

    np.sqrt(upstream_head**2 - (upstream_head**2 - downstream_head**2) * distance / length)

dupuit_discharge
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.dupuit_discharge`

Output: **m^2/s**. Dupuit (1863), the Dupuit-Forchheimer discharge formula, exact for the discharge even where the head profile is not.

.. code-block:: python

    hydraulic_conductivity * (upstream_head**2 - downstream_head**2) / (2 * length)

manning_velocity
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.manning_velocity`

Output: **m/s**. Manning (1891), SI form; fully rough turbulent flow, steady and uniform.

.. code-block:: python

    hydraulic_radius**(2 / 3) * np.sqrt(slope) / roughness

chezy_from_manning
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.chezy_from_manning`

Output: **m^0.5/s**. Manning (1891), SI units.

.. code-block:: python

    hydraulic_radius**(1 / 6) / roughness

green_ampt_rate
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.green_ampt_rate`

Output: **m/s**. Green and Ampt (1911): a piston-like front in a uniform soil with a negligible ponding depth; the rate tends to the saturated conductivity.

.. code-block:: python

    saturated_conductivity * (1 + suction_head * moisture_deficit / cumulative_infiltration)

green_ampt_time
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.green_ampt_time`

Output: **s**. Green and Ampt (1911), ponded from time zero; early on the depth grows as the square root of time.

.. code-block:: python

    (cumulative_infiltration - suction_head * moisture_deficit * np.log(1 + cumulative_infiltration / (suction_head * moisture_deficit))) / saturated_conductivity

gardner_conductivity
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.gardner_conductivity`

Output: **m/s**. Gardner (1958), a one-parameter fit without hysteresis; it makes the steady Richards equation solvable in closed form.

.. code-block:: python

    saturated_conductivity * np.exp(sorptive_number * pressure_head)

unsaturated_darcy_flux
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.unsaturated_darcy_flux`

Output: **m/s**. Richards (1931): the total head is ``pressure_head + z`` and air stays at atmospheric pressure.

.. code-block:: python

    -conductivity * (pressure_head_gradient + 1)

richards_steady_head
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.richards_steady_head`

Output: **m**. Richards (1931) steady equation ``flux = -K(psi) (dpsi/dz + 1)`` with Gardner's (1958) conductivity; an upward flux can only be sustained below ``np.log(1 + saturated_conductivity / flux) / sorptive_number``.

.. code-block:: python

    np.log(((flux + saturated_conductivity) * np.exp(-sorptive_number * height) - flux) / saturated_conductivity) / sorptive_number

nash_unit_hydrograph
~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.nash_unit_hydrograph`

Output: **1/s**. Nash (1957), a gamma density with mean lag ``reservoirs * storage_constant``; a smooth synthetic shape for Sherman's unit hydrograph, linear and time-invariant.

.. code-block:: python

    time**(reservoirs - 1) * np.exp(-time / storage_constant) / (storage_constant**reservoirs * scipy.special.gamma(reservoirs))

unit_hydrograph_convolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.unit_hydrograph_convolution`

Output: **m^3/s**. Sherman (1932): linear superposition and time invariance; baseflow and losses are removed beforehand.

.. code-block:: python

    np.convolve(excess_rainfall, unit_hydrograph)

horton_infiltration_rate
~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.horton_infiltration_rate`

Output: **m/s**. Horton (1933, 1940), an empirical law for rain at or above capacity throughout.

.. code-block:: python

    final_rate + (initial_rate - final_rate) * np.exp(-decay_constant * time)

horton_cumulative_infiltration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.horton_cumulative_infiltration`

Output: **m**. Horton (1933, 1940), with rain at or above capacity from time zero.

.. code-block:: python

    final_rate * time + (initial_rate - final_rate) * (1 - np.exp(-decay_constant * time)) / decay_constant

well_function
~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.well_function`

Output: **1**. The exponential integral E1, ``-0.5772 - ln(u) + u - ...`` for small u.

.. code-block:: python

    scipy.special.exp1(u)

theis_drawdown
~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.theis_drawdown`

Output: **m**. Theis (1935): constant rate, a fully penetrating well of negligible radius in an infinite, homogeneous, isotropic confined aquifer with water released instantly from storage.

.. code-block:: python

    pumping_rate / (4 * np.pi * transmissivity) * scipy.special.exp1(radius**2 * storativity / (4 * transmissivity * time))

muskingum_coefficients
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.muskingum_coefficients`

Output: **1**. McCarthy (1938), storage ``K (X I + (1 - X) O)``; all coefficients are nonnegative when ``2 K X <= time_step <= 2 K (1 - X)``.

.. code-block:: python

    [time_step - 2 K X, time_step + 2 K X, 2 K (1 - X) - time_step] / (2 K (1 - X) + time_step)

muskingum_routing
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.muskingum_routing`

Output: **m^3/s**. McCarthy (1938): ``O[j+1] = C0 I[j+1] + C1 I[j] + C2 O[j]``; choose the time step so that no coefficient is negative.

.. code-block:: python

    O[j + 1] = C0 * I[j + 1] + C1 * I[j] + C2 * O[j]

kinematic_wave_celerity
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.kinematic_wave_celerity`

Output: **m/s**. Lighthill and Whitham (1955): ``dQ/dA`` of a rating ``Q = a A**m``; 5/3 for Manning and 3/2 for Chézy in a wide channel.

.. code-block:: python

    exponent * velocity

degree_day_melt
~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.degree_day_melt`

Output: **m/s** of water equivalent. Finsterwalder and Schunk (1887): periods of a day or longer with snow present.

.. code-block:: python

    degree_day_factor * np.maximum(temperature - threshold_temperature, 0)

snowpack_water_equivalent
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.snowpack_water_equivalent`

Output: **m**, snow water equivalent and released water on the last axis. One threshold separates snowfall from rain and starts melt; no refreezing or liquid storage.

.. code-block:: python

    melt = min(swe + snowfall, degree_day_factor * max(T - T0, 0) * time_step)
    swe = swe + snowfall - melt

exponential_water_content
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.exponential_water_content`

Output: **1**. Russo's (1988) retention for Gardner's soil; no air entry or hysteresis.

.. code-block:: python

    residual_water_content + (saturated_water_content - residual_water_content) * np.exp(sorptive_number * pressure_head)

stream_number
~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.stream_number`

Output: **1**. Horton's law of stream numbers (1945), an empirical average; Strahler's (1952) ordering is the usual convention.

.. code-block:: python

    bifurcation_ratio**(highest_order - order)

stream_length
~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.stream_length`

Output: **m**. Horton's law of stream lengths (1945), an empirical average over a network.

.. code-block:: python

    first_order_length * length_ratio**(order - 1)

cooper_jacob_drawdown
~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.cooper_jacob_drawdown`

Output: **m**. Cooper and Jacob (1946): the first two terms of the Theis series, within 1% for ``u < 0.01`` and meaningless (it goes negative) for u above about 0.56.

.. code-block:: python

    pumping_rate / (4 * np.pi * transmissivity) * (-np.euler_gamma - np.log(radius**2 * storativity / (4 * transmissivity * time)))

cooper_jacob_transmissivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.cooper_jacob_transmissivity`

Output: **m^2/s**. Cooper and Jacob (1946), from late-time data where ``u < 0.01``; the storativity follows from the zero-drawdown intercept ``t0`` as ``2.25 * transmissivity * t0 / radius**2``.

.. code-block:: python

    np.log(10) * pumping_rate / (4 * np.pi * drawdown_per_log_cycle)

saturation_slope
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.saturation_slope`

Output: **Pa/K**. Derivative of :func:`earthsciencekit.atmosphere.saturation_vapor_pressure` (constant latent heat 2.5e6 J/kg).

.. code-block:: python

    saturation_vapor_pressure(temperature) * LATENT_HEAT_VAPORIZATION / (R_VAPOR * temperature**2)

psychrometric_constant
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.psychrometric_constant`

Output: **Pa/K**. Dry-air heat capacity and constant latent heat; equal transfer coefficients for heat and vapour.

.. code-block:: python

    CP_DRY_AIR * pressure / (EPSILON_WATER * LATENT_HEAT_VAPORIZATION)

penman_evaporation
~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.penman_evaporation`

Output: **kg/(m^2 s)**. Penman (1948): the saturation curve is linearized between air and surface temperature, and heat and vapour share one transfer coefficient; a saturated surface (open water or wet grass).

.. code-block:: python

    (saturation_slope(temperature) * net_radiation / LATENT_HEAT_VAPORIZATION + psychrometric_constant(pressure) * aerodynamic_evaporation) / (saturation_slope(temperature) + psychrometric_constant(pressure))

aerodynamic_resistance
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.aerodynamic_resistance`

Output: **s/m**. Neutral stability; the measurement height above the displacement height must exceed both roughness lengths. For FAO grass (0.12 m) at 2 m it is about ``208 / wind_speed``.

.. code-block:: python

    np.log((measurement_height - displacement_height) / momentum_roughness) * np.log((measurement_height - displacement_height) / vapor_roughness) / (VON_KARMAN**2 * wind_speed)

bulk_surface_resistance
~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.bulk_surface_resistance`

Output: **s/m**. Monteith (1965) big-leaf canopy with parallel leaves; FAO grass has about 70 s/m.

.. code-block:: python

    stomatal_resistance / active_leaf_area_index

penman_monteith_evaporation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.penman_monteith_evaporation`

Output: **kg/(m^2 s)**. Monteith (1965): a big-leaf canopy at one temperature, neutral stability; zero surface resistance recovers Penman with ``aerodynamic_evaporation = air_density * CP_DRY_AIR * vapor_pressure_deficit / (psychrometric_constant * LATENT_HEAT_VAPORIZATION * aerodynamic_resistance)``.

.. code-block:: python

    (Delta * available_energy + air_density * CP_DRY_AIR * vapor_pressure_deficit / aerodynamic_resistance) / (LATENT_HEAT_VAPORIZATION * (Delta + gamma * (1 + surface_resistance / aerodynamic_resistance)))

topographic_index
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.topographic_index`

Output: **1**. Beven and Kirkby (1979): steady recharge, a water table parallel to the surface and transmissivity falling exponentially with deficit.

.. code-block:: python

    np.log(specific_area / np.tan(slope))

topmodel_deficit
~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.topmodel_deficit`

Output: **m**. Beven and Kirkby (1979) with uniform transmissivity at saturation and quasi-steady drainage.

.. code-block:: python

    mean_deficit + scaling_depth * (mean_index - topographic_index)

topmodel_baseflow
~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.topmodel_baseflow`

Output: **m/s**. Beven and Kirkby (1979): exponential transmissivity profile, water table parallel to the surface.

.. code-block:: python

    saturated_transmissivity * np.exp(-mean_index) * np.exp(-mean_deficit / scaling_depth)

Infiltration under steady rain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.simulate_infiltration` integrates the
cumulative infiltration ``F`` of a Green-Ampt soil under rain of rate ``r``,
``dF/dt = min(r, K (1 + psi dtheta / F))``, from ``F = 0`` at the first
sample. The surface ponds when ``F = K psi dtheta / (r - K)``, at
``t_p = K psi dtheta / (r (r - K))``; rain no faster than ``K`` never ponds.
It returns the cumulative infiltration, the infiltration rate and the
infiltration-excess runoff rate ``r - f``. Surface storage, redistribution and
a water table are neglected.

Transient Richards column
~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.simulate_richards_column` integrates the
water content of equal cells between a water table (zero head) and a surface
with a prescribed flux, in the mixed form of Celia et al. (1990), by the
method of lines with SciPy's BDF integrator. In the exponential soil the flux
is linear in ``u = exp(alpha h)``. The face fluxes are exponentially fitted and
exact for steady flow, so the column relaxes to ``richards_steady_head``. It
starts from hydrostatic equilibrium unless ``initial_head`` is given, and
returns the water stored above the residual content, the flux across the water
table and the head of every cell.

.. code-block:: python

    d(theta)/dt = -dq/dz,  q = -(K / alpha) du/dz - K u,  theta = theta_r + (theta_s - theta_r) u

Kinematic-wave routing
~~~~~~~~~~~~~~~~~~~~~~

:py:func:`earthsciencekit.hydrology.simulate_kinematic_wave` integrates the
flow area of equal cells along a reach with first-order upwind differences,
starting in uniform flow at the first inflow. It returns the outflow, the
stored volume and the discharge of every cell. Upwinding adds a little
numerical diffusion, which attenuates the peak; more cells reduce it.

.. code-block:: python

    dA/dt + dQ/dx = 0,  Q = np.sqrt(slope) / (roughness * width**(2 / 3)) * A**(5 / 3)

Independent validation
----------------------

Every algebraic API has a fixed-value check, broadcasting or bound checks and
invalid-input contracts. The hydrology tests also check that the water balance
closes over a cycle, that Manning's formula is Chézy's with ``R**(1/6) / n``
and the hydraulic radius tends to the depth, that every section of Dupuit's
water table carries the Dupuit discharge, that the steady Richards profile
carries a constant flux, is hydrostatic without flow and tends to unit
gradient, that the Green-Ampt time inverts the rate and grows as the square
of the depth early on, that the Nash hydrograph integrates to one and the
storm hydrograph conserves volume and superposes, that Horton's cumulative
infiltration integrates the rate, that Theis tends to Cooper-Jacob and its
log-cycle slope returns the transmissivity, that Muskingum coefficients sum
to one and routing conserves volume while delaying and lowering the peak,
that Penman-Monteith with zero surface resistance is Penman and
``saturation_slope`` is the derivative of the saturation curve, that TOPMODEL
deficits average to the mean and baseflow falls by ``e`` per scaling depth,
that ``simulate_infiltration`` matches the Mein-Larson closed form after
ponding, that the transient Richards column relaxes to the steady Gardner
profile and conserves water, that the kinematic wave is steady under steady
inflow, conserves water and carries a flood at its celerity, and that the
snowpack conserves water and melts at the degree-day rate. See :doc:`/validation` for the exact evidence and
:doc:`/numerical_methods` for numerical conventions.
