earthsciencekit: Earth-System Models in Python
==============================================

.. include:: /_generated/vars.rst

**earthsciencekit** makes the equations behind |num_subpackages| Earth-science
domains computable, inspectable and visualizable. Choose a domain to find its
foundational ideas, runnable examples and full API in one place. NumPy and SciPy
provide mature numerical solvers; EarthScienceKit exposes Earth-system equations,
SI conventions and intermediate results for learning and teaching.

.. image:: _static/images/readme_hero.png
   :alt: Streamlines of Stommel's wind-driven gyre with its western boundary current, two Rikitake dynamo trajectories started 0.001 apart, and the IGRF-13 field intensity showing the South Atlantic Anomaly and the magnetic equator, computed with EarthScienceKit
   :width: 100%

Each domain at a glance
-----------------------

- :mod:`earthsciencekit.atmosphere` -- Hydrostatic balance, dry thermodynamics, moisture, rotation and parcel ascent.

  .. image:: _static/images/readme_atmosphere.png
     :alt: Saturation vapour pressure against temperature, lifted parcels following dry and moist adiabats, and Rossby wave phase speeds
     :width: 100%

- :mod:`earthsciencekit.climate` -- Radiative balance, forcing, heat storage, carbon baselines and climate adjustment.

  .. image:: _static/images/readme_climate.png
     :alt: Ice-albedo bistability of equilibrium temperature against solar constant, two-box surface and deep-ocean warming, and daily insolation by latitude
     :width: 100%

- :mod:`earthsciencekit.geomorphology` -- Settling, the Shields threshold, suspension and bedload, hillslope creep and scarps, stream power and river profiles, network scaling, blown sand and cosmogenic erosion rates.

  .. image:: _static/images/readme_geomorphology.png
     :alt: The Shields curve, a diffusing fault scarp, and a stream-power river profile rising toward steady state
     :width: 100%

- :mod:`earthsciencekit.geophysics` -- Gravity anomalies, the figure of the Earth, the geoid, the geomagnetic field, paleomagnetism and the dynamo.

  .. image:: _static/images/readme_geophysics.png
     :alt: Gravity anomalies of a buried sphere, the magnetic anomaly of a magnetized block, and the Rikitake dynamo attractor
     :width: 100%

- :mod:`earthsciencekit.glaciology` -- Glen's flow law, plastic and shallow-ice profiles, sliding, crevasses, mass balance, sea-ice growth, postglacial rebound, ice-core dating, firn densification, degree-day melt and marine grounding lines.

  .. image:: _static/images/readme_glaciology.png
     :alt: Glen's-law velocity profiles, Vialov ice-sheet profiles, and Herron-Langway firn density with depth
     :width: 100%

- :mod:`earthsciencekit.hydrology` -- Catchment water balance, Darcy and Dupuit groundwater flow, Theis wells, infiltration, open channels, unit hydrographs, flood routing and evaporation.

  .. image:: _static/images/readme_hydrology.png
     :alt: Theis drawdown against time, Nash unit hydrographs, and a flood wave routed by the Muskingum method
     :width: 100%

- :mod:`earthsciencekit.oceanography` -- Tides, waves, seawater density, geostrophy, the Ekman spiral, wind-driven gyres and thermohaline boxes.

  .. image:: _static/images/readme_oceanography.png
     :alt: The Ekman spiral, Stommel's western-intensified gyre, and surface-wave phase speed against wavelength
     :width: 100%

- :mod:`earthsciencekit.plate_tectonics` -- Euler plate motion, spreading, subduction, isostasy, flexure and lithosphere cooling.

  .. image:: _static/images/readme_plate_tectonics.png
     :alt: Sea-floor depth and heat flow against age from half-space cooling, and lithospheric flexure under a line load
     :width: 100%

- :mod:`earthsciencekit.seismology` -- Elastic wave speeds, reflection, seismic moment, source pulses and damped response.

  .. image:: _static/images/readme_seismology.png
     :alt: Direct and head-wave refraction travel times, the Gutenberg-Richter law, and a damped oscillator
     :width: 100%

- :mod:`earthsciencekit.volcanology` -- Radiometric dating, isotope and trace-element geochemistry, magma viscosity and ascent, ground deformation and eruption columns.

  .. image:: _static/images/readme_volcanology.png
     :alt: The U-Pb concordia, trace elements in a crystallizing melt, and an eruption-column parcel overshooting its neutral level
     :width: 100%

Start here
----------

:doc:`quickstart` explains installation and basic use, and the
`quickstart notebook <https://colab.research.google.com/github/cpoli/earthsciencekit/blob/main/notebooks/quickstart.ipynb>`__
runs in Colab or `in JupyterLite <lite/lab/index.html?path=quickstart.ipynb>`__.
:doc:`numerical_methods` explains sampled histories and solver choices. The domain
guides describe equations, units, sign conventions and applicability limits;
:doc:`validation` documents the evidence behind this release, and
:doc:`limitations` identifies remaining scope limits. Version 0.3.0 includes nine
domains, with geomorphology to follow, and enforces 100% production statement and
branch coverage. That metric does not prove scientific correctness or model
applicability.

.. code-block:: python

   import numpy as np
   from earthsciencekit.climate import simulate_one_box
   from earthsciencekit.units import convert

   time = convert(np.linspace(0, 100, 201), "year", "s")
   result = simulate_one_box(4, 1.2, 2e8, time)
   print(result.summary())
   figure, axes = result.plot()

.. toctree::
   :maxdepth: 2
   :caption: Subpackages
   :hidden:

   subpackages/index

.. toctree::
   :maxdepth: 2
   :caption: History
   :hidden:

   history/index

.. toctree::
   :maxdepth: 2
   :caption: Examples
   :hidden:

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API
   :hidden:

   api/index

.. rubric:: Related toolkits

Sibling Python toolkits for learning and teaching, built the same way:

- `blockchainkit <https://cpoli.github.io/blockchainkit/>`__ -- cryptography, consensus and blockchain systems
- `chemistrykit <https://cpoli.github.io/chemistrykit/>`__ -- computational chemistry
- `engineeringkit <https://cpoli.github.io/engineeringkit/>`__ -- computational engineering
- `financialkit <https://cpoli.github.io/financialkit/>`__ -- quantitative finance
- `mathematicskit <https://cpoli.github.io/mathematicskit/>`__ -- computational mathematics
- `physicskit <https://cpoli.github.io/physicskit/>`__ -- computational physics
- `tbkit <https://cpoli.github.io/tbkit/>`__ -- tight-binding models, band topology and quantum transport
