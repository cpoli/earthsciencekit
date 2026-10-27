Breakthroughs in Glaciology
===========================

.. include:: /_generated/nav/glaciology.rst

Glaciology began with the recognition that ice once covered much of Europe
and that glaciers flow. The nineteenth century added the sinking of the land
under the ice and the physics of freezing seas. The twentieth century made
ice a material: Haskell measured the mantle by its rebound, Nye and Glen gave
ice a plastic and then a power-law rheology, and from them came crevasse
depths, sliding and the shallow-ice models that now simulate whole ice
sheets. Meier set the vocabulary of mass balance, Dansgaard turned ice cores
into thermometers and, with Johnsen, into clocks. Herron and Langway showed
how firn closes into bubbles, and Barnola read the atmosphere of past glacial
cycles from them. Reeh computed an ice sheet's melt from its air temperature,
and Weertman and Schoof
showed why marine ice sheets can collapse. This chronology traces the results
behind :mod:`earthsciencekit.glaciology`. Every example uses deterministic
synthetic inputs; none inverts real measurements.

.. contents:: Timeline
   :local:
   :depth: 1

1840 -- Agassiz's Ice Age
-------------------------

Louis Agassiz, building on Venetz and de Charpentier, read the erratic
boulders, polished and striated rock and moraines of the Swiss plain and the
Jura as the work of glaciers far larger than those of the Alps today, and
proposed that a great ice age had covered much of Europe. A cooling of a few
kelvin is enough: the equilibrium line, above which snow survives the summer,
falls by the cooling divided by the lapse rate.

*Implementation:* :func:`~earthsciencekit.glaciology.ela_shift` (*new*). The
example lowers the equilibrium line over a synthetic mountain range and
follows the area above it.

*References:* L. Agassiz, *Études sur les glaciers* (Neuchâtel: Jent et
Gassmann, 1840).

.. minigallery:: ../../examples/glaciology/plot_01_agassiz_ice_age.py

1843 -- Forbes's Viscous Glacier Flow
-------------------------------------

James David Forbes surveyed the Mer de Glace in 1842 and found that its
centre moves faster than its margins, its surface faster than its depths, and
all of it faster in summer. He concluded that a glacier is a viscous body
flowing under its own weight. For a slab of uniform viscosity on a slope the
velocity is a parabola in height above the bed.

*Implementation:* :func:`~earthsciencekit.glaciology.viscous_slab_velocity`
(*new*) and :func:`~earthsciencekit.glaciology.driving_stress` (*new*). The
example draws the parabolic profile for three viscosities and checks that its
mean is two-thirds of the surface speed.

*References:* J. D. Forbes, *Travels through the Alps of Savoy and Other
Parts of the Pennine Chain, with Observations on the Phenomena of Glaciers*
(Edinburgh: Adam and Charles Black, 1843).

.. minigallery:: ../../examples/glaciology/plot_02_forbes_viscous_glacier_flow.py

1865 -- Jamieson's Glacial Isostasy
-----------------------------------

Thomas Jamieson explained the raised beaches of Scotland, now well above the
sea, by the weight of the vanished ice sheet: the land had sunk under the ice
and risen again after it melted. In isostatic equilibrium the depression is
the ice thickness times the ratio of ice to mantle density, a little under a
third.

*Implementation:* :func:`~earthsciencekit.glaciology.isostatic_depression`
(*new*); compare :func:`~earthsciencekit.plate_tectonics.airy_root`. The
example depresses the bed under a synthetic ice sheet.

*References:* T. F. Jamieson, "On the history of the last geological changes
in Scotland," *Quarterly Journal of the Geological Society* 21 (1865):
161-203.

.. minigallery:: ../../examples/glaciology/plot_03_jamieson_glacial_isostasy.py

1891 -- Stefan's Sea-Ice Growth
-------------------------------

Josef Stefan treated the freezing of polar seas as a problem with a moving
boundary: the latent heat released at the base of the ice is conducted up
through the ice already formed. The thickness then grows as the square root
of the accumulated freezing degree-days, and Stefan compared the law with ice
thicknesses measured on Arctic expeditions.

*Implementation:* :func:`~earthsciencekit.glaciology.stefan_ice_thickness`
(*new*), :func:`~earthsciencekit.glaciology.stefan_growth_rate` (*new*) and
:func:`~earthsciencekit.glaciology.sea_ice_freeboard` (*new*). The example
grows ice through two synthetic winters.

*References:* J. Stefan, "Ueber die Theorie der Eisbildung, insbesondere
über die Eisbildung im Polarmeere," *Annalen der Physik* 278 (1891):
269-286.

.. minigallery:: ../../examples/glaciology/plot_04_stefan_sea_ice_growth.py

1935 -- Haskell's Mantle Viscosity from Rebound
-----------------------------------------------

Norman Haskell modelled the mantle under Fennoscandia as a viscous half-space
and computed how a depression of the size of the former ice sheet recovers
once the load is gone: exponentially, with a relaxation time proportional to
the viscosity and inversely proportional to the wavelength. Matching the
uplift recorded by raised shorelines gave about 1e21 Pa s, still the standard
value for the upper mantle.

*Implementation:* :func:`~earthsciencekit.glaciology.viscous_relaxation_time`
(*new*), :func:`~earthsciencekit.glaciology.haskell_viscosity` (*new*),
:func:`~earthsciencekit.glaciology.remaining_depression` (*new*) and
:func:`~earthsciencekit.glaciology.rebound_uplift_rate` (*new*). The example
recovers the viscosity from a synthetic uplift history.

*References:* N. A. Haskell, "The motion of a viscous fluid under a surface
load," *Physics* 6 (1935): 265-269.

.. minigallery:: ../../examples/glaciology/plot_05_haskell_mantle_viscosity.py

1951 -- Nye's Perfectly Plastic Ice Sheet
-----------------------------------------

John Nye treated ice as a perfectly plastic solid that yields at a fixed
stress of about 100 kPa, close to the basal stresses found under glaciers.
On a flat bed the ice sheet is then a parabola whose height depends only on
the distance from its margin, which reproduces the profiles of Greenland and
Antarctica remarkably well.

*Implementation:* :func:`~earthsciencekit.glaciology.nye_plastic_profile`
(*new*). The example draws the profile of an ice sheet 1000 km wide for three
yield stresses.

*References:* J. F. Nye, "The flow of glaciers and ice-sheets as a problem
in plasticity," *Proceedings of the Royal Society of London A* 207 (1951):
554-572.

.. minigallery:: ../../examples/glaciology/plot_06_nye_perfectly_plastic_ice_sheet.py

1955 -- Glen's Flow Law
-----------------------

John Glen measured the steady creep of polycrystalline ice under constant
compressive stress at several temperatures. The strain rate rose as roughly
the cube of the stress, with a rate factor that increases steeply towards the
melting point. Combined with Nye's (1952) laminar-flow analysis, Glen's law
gives the velocity profile of a glacier frozen to its bed.

*Implementation:* :func:`~earthsciencekit.glaciology.glen_strain_rate`
(*new*), :func:`~earthsciencekit.glaciology.glen_stress` (*new*),
:func:`~earthsciencekit.glaciology.glen_rate_factor` (*new*),
:func:`~earthsciencekit.glaciology.effective_viscosity` (*new*),
:func:`~earthsciencekit.glaciology.glen_velocity_profile` (*new*) and
:func:`~earthsciencekit.glaciology.glen_surface_velocity` (*new*). The
example plots the flow law at three temperatures and the velocity profile for
three exponents.

*References:* J. W. Glen, "The creep of polycrystalline ice," *Proceedings of
the Royal Society of London A* 228 (1955): 519-538.

.. minigallery:: ../../examples/glaciology/plot_07_glen_flow_law.py

1955 -- Nye's Crevasse Depth
----------------------------

John Nye argued that crevasses close where the pressure of the overlying ice
exceeds the tensile stress from stretching, which Glen's law gives from the
strain rate. Dry crevasses therefore reach only a few tens of metres, as
observed. Water standing in a crevasse adds its pressure and drives it
deeper, the mechanism later invoked for the break-up of ice shelves.

*Implementation:* :func:`~earthsciencekit.glaciology.nye_crevasse_depth`
(*new*) and :func:`~earthsciencekit.glaciology.water_filled_crevasse_depth`
(*new*). The example plots the depth of dry and water-filled crevasses
against the stretching rate.

*References:* J. F. Nye, "Comments on Dr. Loewe's letter and notes on
crevasses," *Journal of Glaciology* 2 (1955): 512-514.

.. minigallery:: ../../examples/glaciology/plot_08_nye_crevasse_depth.py

1957 -- Weertman's Sliding Law
------------------------------

Johannes Weertman gave the first theory of how temperate ice slides over a
rough bed: by regelation around small obstacles and enhanced creep around
large ones. The obstacle size at which the two are equally effective controls
the motion, and the sliding speed rises as a power of the basal stress, about
the square for Glen's exponent of 3.

*Implementation:* :func:`~earthsciencekit.glaciology.weertman_sliding_velocity`
(*new*), with :func:`~earthsciencekit.glaciology.glen_surface_velocity`. The
example adds sliding to internal deformation.

*References:* J. Weertman, "On the sliding of glaciers," *Journal of
Glaciology* 3 (1957): 33-38.

.. minigallery:: ../../examples/glaciology/plot_09_weertman_sliding_law.py

1962 -- Meier's Glacier Mass Balance
------------------------------------

Mark Meier proposed the definitions that made glacier mass budgets
comparable: the net balance at a point and over the glacier, the equilibrium
line where the balance is zero, and the accumulation-area ratio. Meier also
noted that the balance rises nearly linearly with altitude, which ties a
glacier's size to its equilibrium line.

*Implementation:* :func:`~earthsciencekit.glaciology.linear_mass_balance`
(*new*) and :func:`~earthsciencekit.glaciology.accumulation_area_ratio`
(*new*). The example finds the equilibrium line that balances a synthetic
glacier and its accumulation-area ratio.

*References:* M. F. Meier, "Proposed definitions for glacier mass budget
terms," *Journal of Glaciology* 4 (1962): 252-263.

.. minigallery:: ../../examples/glaciology/plot_10_meier_glacier_mass_balance.py

1964 -- Dansgaard's Isotope Thermometer
---------------------------------------

Willi Dansgaard analysed precipitation from stations around the world and
found that its oxygen-18 content falls by about 0.7 permil per kelvin of
mean annual temperature at mid and high latitudes, because the heavy isotopes
condense first as moist air cools. The relation turned the isotope profiles
of ice cores into temperature records.

*Implementation:* :func:`~earthsciencekit.glaciology.dansgaard_delta` (*new*)
and :func:`~earthsciencekit.glaciology.dansgaard_temperature` (*new*). The
example fits synthetic stations and reads the glacial cooling from the
isotope shift.

*References:* W. Dansgaard, "Stable isotopes in precipitation," *Tellus* 16
(1964): 436-468.

.. minigallery:: ../../examples/glaciology/plot_11_dansgaard_isotope_thermometer.py

1969 -- Dansgaard and Johnsen's Age-Depth Model
-----------------------------------------------

Willi Dansgaard and Sigfus Johnsen dated the Camp Century core, the first to
reach the bed of the Greenland ice sheet, with a flow model: the vertical
velocity falls linearly with depth down to a kink and then quadratically to
zero at the frozen bed. Integrating its reciprocal gave the first ice-core
time scale, older at depth than Nye's (1963) model of uniform thinning.

*Implementation:* :func:`~earthsciencekit.glaciology.dansgaard_johnsen_age`
(*new*), :func:`~earthsciencekit.glaciology.dansgaard_johnsen_downward_velocity`
(*new*) and :func:`~earthsciencekit.glaciology.nye_age` (*new*). The example
builds the Camp Century time scale and its annual-layer thickness.

*References:* W. Dansgaard and S. J. Johnsen, "A flow model and a time scale
for the ice core from Camp Century, Greenland," *Journal of Glaciology* 8
(1969): 215-223.

.. minigallery:: ../../examples/glaciology/plot_12_dansgaard_johnsen_age_depth_model.py

1974 -- Weertman's Marine Ice-Sheet Instability
-----------------------------------------------

Johannes Weertman analysed the junction between an ice sheet grounded below
sea level and its floating shelf. The ice at the junction has the flotation
thickness, set by the water depth, and the shelf spreads at a rate rising as
a power of that thickness (Weertman 1957). On a bed that deepens inland the
junction cannot be in stable equilibrium, so a marine ice sheet such as West
Antarctica could retreat irreversibly.

*Implementation:* :func:`~earthsciencekit.glaciology.flotation_thickness`
(*new*) and :func:`~earthsciencekit.glaciology.ice_shelf_strain_rate`
(*new*). The example compares the snowfall upstream with a discharge measure
on beds deepening seaward and inland.

*References:* J. Weertman, "Stability of the junction of an ice sheet and an
ice shelf," *Journal of Glaciology* 13 (1974): 3-11.

.. minigallery:: ../../examples/glaciology/plot_13_weertman_marine_ice_sheet_instability.py

1976 -- Mahaffy's Shallow-Ice Model
-----------------------------------

M. W. Mahaffy built the first three-dimensional numerical model of an ice
sheet and tested it on the Barnes Ice Cap. Because ice sheets are thin, the
flux at each point follows from the local thickness and slope, and the
thickness obeys a nonlinear diffusion equation, solved on a staggered grid.
Its steady state on a flat bed is Vialov's (1958) profile.

*Implementation:* :func:`~earthsciencekit.glaciology.shallow_ice_flux`
(*new*), :func:`~earthsciencekit.glaciology.shallow_ice_diffusivity` (*new*),
:func:`~earthsciencekit.glaciology.vialov_dome_height` (*new*),
:func:`~earthsciencekit.glaciology.vialov_profile` (*new*),
:func:`~earthsciencekit.glaciology.glacier_response_time` (*new*) and the
dynamic :func:`~earthsciencekit.glaciology.simulate_flowline_glacier`
(*new*). The example compares Vialov's profile with Nye's and follows a
mountain glacier after its equilibrium line rises by 100 m.

*References:* M. W. Mahaffy, "A three-dimensional numerical model of ice
sheets: tests on the Barnes Ice Cap, Northwest Territories," *Journal of
Geophysical Research* 81 (1976): 1059-1066.

.. minigallery:: ../../examples/glaciology/plot_14_mahaffy_shallow_ice_model.py

1980 -- Herron and Langway's Firn Densification
-----------------------------------------------

Michael Herron and Chester Langway fitted the density-depth profiles of
firn cores from Greenland and Antarctica with an empirical model in two
stages. Above 550 kg/m^3 grains settle and pack at a rate set by
temperature alone; below it they sinter more slowly, at a rate that also
grows with the square root of the accumulation. In each stage
``ln(density / (ice density - density))`` is linear in depth. The model puts
the pore close-off, near 830 kg/m^3, at 70 to 100 m and turns the gap
between the age of the ice and the age of its air into a computed quantity.

*Implementation:* :func:`~earthsciencekit.glaciology.herron_langway_density`
(*new*), :func:`~earthsciencekit.glaciology.herron_langway_depth` (*new*) and
:func:`~earthsciencekit.glaciology.herron_langway_mean_density` (*new*),
with :func:`~earthsciencekit.glaciology.gas_age_difference`. The example
compares the profiles at four core sites and maps the ice-age minus gas-age
difference against temperature and accumulation.

*References:* M. M. Herron and C. C. Langway, "Firn densification: an
empirical model," *Journal of Glaciology* 25 (1980): 373-385.

.. minigallery:: ../../examples/glaciology/plot_15_herron_langway_firn_densification.py

1987 -- Barnola's Vostok CO2 Record
-----------------------------------

Jean-Marc Barnola, Dominique Raynaud, Yevgeniy Korotkevich and Claude Lorius
measured the carbon dioxide in air bubbles of the Vostok core through the
last glacial cycle: about 190 ppm in glacial periods and 280 ppm in
interglacials, closely following the isotope temperature. The air is younger
than the ice around it because the firn stays open until the bubbles close,
a gap of several thousand years at Vostok.

*Implementation:* :func:`~earthsciencekit.glaciology.gas_age_difference`
(*new*), with :func:`~earthsciencekit.climate.co2_forcing`. The example
places a synthetic CO2 record at its ice depth, with the gap computed from
Herron and Langway's firn at each climate, and computes its forcing.

*References:* J. M. Barnola, D. Raynaud, Y. S. Korotkevich and C. Lorius,
"Vostok ice core provides 160,000-year record of atmospheric CO2," *Nature*
329 (1987): 408-414.

.. minigallery:: ../../examples/glaciology/plot_16_barnola_vostok_co2_record.py

1991 -- Reeh's Positive-Degree-Day Model
----------------------------------------

Niels Reeh estimated the melt of the whole Greenland ice sheet from air
temperature alone. He treated daily temperature as normally distributed
about a sinusoidal annual cycle, so that the expected positive degree-days
stay positive where the mean never thaws. The degree-days melt the year's
snow first, at about 3 mm of water per degree-day, then the ice beneath at
about 8, after Braithwaite's measurements. Meltwater refreezes in the snow
up to 60% of the snowfall before any runs off. Calov and Greve (2005) later
gave the expected degree-days in closed form. Reeh's volume of
*Polarforschung* is dated 1989 but appeared in 1991.

*Implementation:* :func:`~earthsciencekit.glaciology.expected_positive_degree_days`
(*new*) and :func:`~earthsciencekit.glaciology.pdd_surface_melt` (*new*). The
example finds the equilibrium line along a synthetic Greenland flow line.

*References:* N. Reeh, "Parameterization of melt rate and surface
temperature on the Greenland ice sheet," *Polarforschung* 59 (1989, printed
1991): 113-128. R. Calov and R. Greve, "A semi-analytical solution for the
positive degree-day model with stochastic temperature variations," *Journal
of Glaciology* 51 (2005): 173-175.

.. minigallery:: ../../examples/glaciology/plot_17_reeh_positive_degree_day_melt.py

2007 -- Schoof's Grounding-Line Flux
------------------------------------

Christian Schoof resolved the boundary layer at the grounding line of a
sliding marine ice sheet and derived the flux across it as a power of the
flotation thickness, about 4.75 for Glen's n = 3 and a sliding exponent of
1/3. Steady grounding lines are stable only where the bed rises inland,
confirming Weertman's instability and predicting hysteresis on overdeepened
beds.

*Implementation:* :func:`~earthsciencekit.glaciology.schoof_grounding_line_flux`
(*new*). The example finds the stable and unstable steady states on an
overdeepened bed.

*References:* C. Schoof, "Ice sheet grounding line dynamics: steady states,
stability, and hysteresis," *Journal of Geophysical Research* 112 (2007):
F03S28.

.. minigallery:: ../../examples/glaciology/plot_18_schoof_grounding_line_flux.py

See also
--------

- :doc:`/guides/glaciology` for conventions, equations and checks.
- `National Snow and Ice Data Center <https://nsidc.org/learn>`_
- `AntarcticGlaciers.org <https://www.antarcticglaciers.org/>`_
