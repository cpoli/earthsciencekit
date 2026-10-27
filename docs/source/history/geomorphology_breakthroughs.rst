Breakthroughs in Geomorphology
==============================

.. include:: /_generated/nav/geomorphology.rst

Geomorphology began when Hutton's followers read valleys as the work of the
rivers in them. In the nineteenth century Stokes gave the physics of a
settling grain, and Gilbert and Davis turned rivers and landscapes into
systems that grade themselves and age. The twentieth century made
sediment transport quantitative, through Hjulström's diagram, Shields's
threshold, Rouse's suspension profile, Bagnold's blown sand and the bedload
law of Meyer-Peter and Müller. It also measured whole landscapes: Strahler
ordered their networks, Leopold and Maddock their channels, Hack and Flint
their basins and profiles, Culling their hillslopes, and Wolman and Miller
asked which floods do the work. Howard and Kerby's stream-power law and
Lal's cosmogenic nuclides finally gave rates of erosion in years. This
chronology traces the results behind :mod:`earthsciencekit.geomorphology`.
Every example uses deterministic synthetic inputs; none inverts real
measurements. Horton's laws of stream numbers (1945), the other classic of
drainage networks, are in :doc:`/history/hydrology_breakthroughs`.

.. contents:: Timeline
   :local:
   :depth: 1

1802 -- Playfair's Law of Accordant Junctions
---------------------------------------------

John Playfair, in his exposition of Hutton's theory of the Earth, observed
that every river runs in a valley proportioned to its size, and that its
tributaries join it at exactly its own level. Neither would happen by chance:
the rivers must have cut the valleys themselves. In a steady landscape
eroded by stream power the elevation is proportional to the chi coordinate,
so a tributary and its trunk lie on one line in chi and meet at the same
height.

*Implementation:* :func:`~earthsciencekit.geomorphology.chi_coordinate`
(*new*) and :func:`~earthsciencekit.geomorphology.steady_state_slope`
(*new*). The example grades a trunk and a tributary and plots both against
chi.

*References:* J. Playfair, *Illustrations of the Huttonian Theory of the
Earth* (Edinburgh: William Creech, 1802), §99; J. T. Perron and L. Royden,
"An integral approach to bedrock river profile analysis," *Earth Surface
Processes and Landforms* 38 (2013): 570-576.

.. minigallery:: ../../examples/geomorphology/plot_01_playfair_accordant_junctions.py

1851 -- Stokes Settling
-----------------------

George Gabriel Stokes computed the viscous drag on a slowly moving sphere,
in a paper on the internal friction of fluids and its effect on pendulums.
Balanced against the submerged weight, the drag gives a settling velocity
proportional to the square of the diameter, which holds for silt and clay.
For sand and gravel the drag grows faster, and the Ferguson-Church formula
joins Stokes's law to a square-root law.

*Implementation:* :func:`~earthsciencekit.geomorphology.stokes_settling_velocity`
(*new*) and :func:`~earthsciencekit.geomorphology.ferguson_church_settling_velocity`
(*new*). The example compares both from 1 µm to 10 mm.

*References:* G. G. Stokes, "On the effect of the internal friction of
fluids on the motion of pendulums," *Transactions of the Cambridge
Philosophical Society* 9 (1851): 8-106; R. I. Ferguson and M. Church, "A
simple universal equation for grain settling velocity," *Journal of
Sedimentary Research* 74 (2004): 933-937.

.. minigallery:: ../../examples/geomorphology/plot_02_stokes_settling_velocity.py

1877 -- Gilbert's Graded River
------------------------------

Grove Karl Gilbert's report on the Henry Mountains of Utah set out the laws
of land sculpture. A river adjusts its slope until its capacity matches the
work it must do: steep in the headwaters, where little water flows, and
gentle downstream. The equilibrium long profile is concave. Under uniform
uplift and stream-power incision, a river grades itself in this way and
then cuts down exactly as fast as the land rises.

*Implementation:* :func:`~earthsciencekit.geomorphology.simulate_river_profile`
(*new*) and :func:`~earthsciencekit.geomorphology.stream_power_erosion`
(*new*). The example grows a graded profile from an uplifting plain.

*References:* G. K. Gilbert, *Report on the Geology of the Henry Mountains*
(Washington: Government Printing Office, 1877).

.. minigallery:: ../../examples/geomorphology/plot_03_gilbert_graded_river.py

1899 -- Davis's Geographical Cycle
----------------------------------

William Morris Davis described landscapes as passing through youth,
maturity and old age after a pulse of uplift, their relief falling quickly
at first and then ever more slowly toward a low peneplain. The cycle is
qualitative, but its shape depends on how erosion scales with slope: with
stream power linear in slope the erosion sweeps up the river as a wave,
while a steeper dependence gives Davis's slow decline. Isostatic rebound
means far more rock is removed than the mean elevation loses.

*Implementation:* :func:`~earthsciencekit.geomorphology.simulate_river_profile`
(*new*) with :func:`~earthsciencekit.plate_tectonics.airy_root`. The example
follows the relief of a river after uplift stops, for slope exponents of 1
and 2.

*References:* W. M. Davis, "The geographical cycle," *Geographical Journal*
14 (1899): 481-504.

.. minigallery:: ../../examples/geomorphology/plot_04_davis_geographical_cycle.py

1935 -- Hjulström's Erosion-Deposition Diagram
----------------------------------------------

Filip Hjulström, in his study of the River Fyris, gathered measurements into
one diagram of mean velocity against grain size, with fields of erosion,
transport and deposition. Once lifted, silt stays in suspension at
velocities far below those that eroded it, while gravel drops almost as soon
as it stops rolling. Hjulström's erosion curve also rises again for the
finest grains, from the cohesion of mud.

*Implementation:* :func:`~earthsciencekit.geomorphology.keulegan_mean_velocity`
(*new*) and :func:`~earthsciencekit.geomorphology.shear_velocity` (*new*).
The example rebuilds the diagram, without cohesion, from the Shields
threshold and the settling velocity.

*References:* F. Hjulström, "Studies of the morphological activity of rivers
as illustrated by the River Fyris," *Bulletin of the Geological Institution
of the University of Upsala* 25 (1935): 221-527; G. H. Keulegan, "Laws of
turbulent flow in open channels," *Journal of Research of the National
Bureau of Standards* 21 (1938): 707-741.

.. minigallery:: ../../examples/geomorphology/plot_05_hjulstrom_curve.py

1936 -- Shields's Threshold of Motion
-------------------------------------

Albert Shields measured in a Berlin flume the bed stress at which grains of
different sizes and densities begin to move. By similarity he collapsed the
results onto one curve: the critical ratio of bed stress to submerged grain
weight depends only on a grain Reynolds number. It falls through the sand
sizes and levels off near 0.055 for gravel.

*Implementation:* :func:`~earthsciencekit.geomorphology.shields_number`
(*new*), :func:`~earthsciencekit.geomorphology.dimensionless_grain_size`
(*new*) and :func:`~earthsciencekit.geomorphology.critical_shields_number`
(*new*, the Soulsby-Whitehouse fit). The example draws the curve with quartz
grains marked.

*References:* A. Shields, "Anwendung der Ähnlichkeitsmechanik und der
Turbulenzforschung auf die Geschiebebewegung," *Mitteilungen der Preußischen
Versuchsanstalt für Wasserbau und Schiffbau* 26 (Berlin, 1936); R. L.
Soulsby and R. J. S. W. Whitehouse, "Threshold of sediment motion in coastal
environments," *Proceedings of Pacific Coasts and Ports '97* (Christchurch,
1997): 149-154.

.. minigallery:: ../../examples/geomorphology/plot_06_shields_threshold_of_motion.py

1937 -- Rouse's Suspended-Sediment Profile
------------------------------------------

Hunter Rouse balanced the settling of suspended grains against turbulent
mixing with the eddy diffusivity of a logarithmic boundary layer, and found
the concentration profile as a power of height. Its exponent, the Rouse
number, compares the settling velocity with the shear velocity: fine sand
fills the whole depth of a river, coarse sand stays near the bed.

*Implementation:* :func:`~earthsciencekit.geomorphology.rouse_number` (*new*)
and :func:`~earthsciencekit.geomorphology.rouse_concentration` (*new*). The
example draws profiles for four sand sizes.

*References:* H. Rouse, "Modern conceptions of the mechanics of fluid
turbulence," *Transactions of the American Society of Civil Engineers* 102
(1937): 463-543.

.. minigallery:: ../../examples/geomorphology/plot_07_rouse_suspended_sediment_profile.py

1941 -- Bagnold's Physics of Blown Sand
---------------------------------------

Ralph Bagnold, an explorer of the Libyan desert, measured in a wind tunnel
how wind moves sand. Grains start to move at a threshold shear velocity
growing as the square root of their size, then hop in saltation, and the
mass carried grows as the cube of the shear velocity. A dune that keeps its
shape advances by the flux over its crest divided by its height.

*Implementation:* :func:`~earthsciencekit.geomorphology.threshold_shear_velocity`
(*new*), :func:`~earthsciencekit.geomorphology.bagnold_sand_flux` (*new*),
:func:`~earthsciencekit.geomorphology.kawamura_sand_flux` (*new*) and
:func:`~earthsciencekit.geomorphology.dune_migration_rate` (*new*). The
example compares the cube law with Kawamura's threshold form and finds the
speed of dunes of different heights.

*References:* R. A. Bagnold, *The Physics of Blown Sand and Desert Dunes*
(London: Methuen, 1941); R. Kawamura, "Study of sand movement by wind,"
*Reports of the Institute of Science and Technology, University of Tokyo* 5
(1951): 95-112.

.. minigallery:: ../../examples/geomorphology/plot_08_bagnold_blown_sand.py

1948 -- Meyer-Peter and Müller Bedload
--------------------------------------

Eugen Meyer-Peter and Robert Müller fitted flume experiments at ETH Zürich
with a bedload law still used in river engineering: the dimensionless
transport rate grows as the 3/2 power of the excess of the Shields stress
over 0.047. Wong and Parker (2006) found the coefficient about twice too
high for a plane bed.

*Implementation:* :func:`~earthsciencekit.geomorphology.meyer_peter_muller_transport`
(*new*), :func:`~earthsciencekit.geomorphology.bedload_flux` (*new*) and
:func:`~earthsciencekit.geomorphology.bed_shear_stress` (*new*), with
:func:`~earthsciencekit.hydrology.hydraulic_radius`. The example applies the
law to a gravel river in flood.

*References:* E. Meyer-Peter and R. Müller, "Formulas for bed-load
transport," *Proceedings of the 2nd Meeting of the International Association
for Hydraulic Structures Research* (Stockholm, 1948): 39-64; M. Wong and G.
Parker, "Reanalysis and correction of bed-load relation of Meyer-Peter and
Müller using their own database," *Journal of Hydraulic Engineering* 132
(2006): 1159-1168.

.. minigallery:: ../../examples/geomorphology/plot_09_meyer_peter_muller_bedload.py

1952 -- Strahler Stream Order
-----------------------------

Arthur Strahler simplified Horton's ordering of streams: every unbranched
headwater channel is first order, and the order rises by one only where two
streams of equal order meet. Unlike Horton's scheme it needs no choice of
which branch continues upstream, so it follows from the links of a network
alone.

*Implementation:* :func:`~earthsciencekit.geomorphology.strahler_order`
(*new*). The example orders a random network and draws it.

*References:* A. N. Strahler, "Hypsometric (area-altitude) analysis of
erosional topography," *Geological Society of America Bulletin* 63 (1952):
1117-1142; A. N. Strahler, "Quantitative analysis of watershed
geomorphology," *Transactions of the American Geophysical Union* 38 (1957):
913-920.

.. minigallery:: ../../examples/geomorphology/plot_10_strahler_stream_order.py

1953 -- Leopold and Maddock's Hydraulic Geometry
------------------------------------------------

Luna Leopold and Thomas Maddock plotted width, depth and velocity against
discharge at gauging stations and found power laws: downstream, width grows
about as the square root of discharge, depth as its 0.4 power and velocity
barely at all. Continuity requires the exponents to sum to one.

*Implementation:* :func:`~earthsciencekit.geomorphology.downstream_hydraulic_geometry`
(*new*), with :func:`~earthsciencekit.hydrology.manning_velocity`. The
example fits the laws to synthetic stations and finds the slope the channels
need.

*References:* L. B. Leopold and T. Maddock, *The Hydraulic Geometry of Stream
Channels and Some Physiographic Implications*, U.S. Geological Survey
Professional Paper 252 (1953).

.. minigallery:: ../../examples/geomorphology/plot_11_leopold_maddock_hydraulic_geometry.py

1957 -- Hack's Law
------------------

John Hack measured streams in Virginia and Maryland and found the length of
the main stream growing as about the 0.6 power of drainage area, with a
coefficient of 1.4 in miles. An exponent above one-half means basins grow
longer as they grow larger.

*Implementation:* :func:`~earthsciencekit.geomorphology.hack_length` (*new*).
The example fits the law to synthetic basins and converts Hack's coefficient
to SI.

*References:* J. T. Hack, *Studies of Longitudinal Stream Profiles in
Virginia and Maryland*, U.S. Geological Survey Professional Paper 294-B
(1957): 45-97.

.. minigallery:: ../../examples/geomorphology/plot_12_hack_law.py

1960 -- Culling's Hillslope Diffusion
-------------------------------------

W. E. H. Culling showed that if soil creeps downhill at a rate proportional
to the slope, the land surface obeys the diffusion equation and the
mathematics of heat conduction applies to hillslopes. A fault scarp rounds
into an error function whose steepest slope falls as one over the square
root of its age, which Hanks and others (1984) turned into a morphologic
age.

*Implementation:* :func:`~earthsciencekit.geomorphology.scarp_elevation`
(*new*), :func:`~earthsciencekit.geomorphology.scarp_maximum_gradient`
(*new*) and :func:`~earthsciencekit.geomorphology.morphologic_age` (*new*).
The example degrades a scarp and recovers its age.

*References:* W. E. H. Culling, "Analytical theory of erosion," *Journal of
Geology* 68 (1960): 336-344; T. C. Hanks, R. C. Bucknam, K. R. Lajoie and R.
E. Wallace, "Modification of wave-cut and faulting-controlled landforms,"
*Journal of Geophysical Research* 89 (1984): 5771-5790.

.. minigallery:: ../../examples/geomorphology/plot_13_culling_hillslope_diffusion.py

1960 -- Wolman and Miller's Magnitude and Frequency
---------------------------------------------------

M. Gordon Wolman and John Miller asked which events do most of the work of
shaping landscapes. Sediment transport rises faster than linearly with
discharge but large flows are rare, so the product of the two peaks at a
moderate, frequent flood, the effective discharge.

*Implementation:* :func:`~earthsciencekit.geomorphology.effective_discharge`
(*new*). The example bins a synthetic daily-flow record and compares the
peak with the closed form for a lognormal frequency.

*References:* M. G. Wolman and J. P. Miller, "Magnitude and frequency of
forces in geomorphic processes," *Journal of Geology* 68 (1960): 54-74.

.. minigallery:: ../../examples/geomorphology/plot_14_wolman_miller_effective_discharge.py

1974 -- Flint's Slope-Area Relation
-----------------------------------

J. J. Flint showed that channel slope falls as a power of drainage area
along the streams of a basin, with an exponent, the concavity index, near
one-half. Its coefficient, the steepness index, now maps rates of uplift.

*Implementation:* :func:`~earthsciencekit.geomorphology.flint_slope` (*new*).
The example fits the concavity and steepness of a synthetic slope-area plot.

*References:* J. J. Flint, "Stream gradient as a function of order,
magnitude, and discharge," *Water Resources Research* 10 (1974): 969-973.

.. minigallery:: ../../examples/geomorphology/plot_15_flint_slope_area_relation.py

1983 -- Howard and Kerby's Stream-Power Incision
------------------------------------------------

Alan Howard and Gordon Kerby measured channel lowering in Virginia badlands
and found it proportional to powers of drainage area and slope, the
stream-power law behind most models of bedrock rivers. A fall of base level
travels upstream as a knickpoint at a speed set by the drainage area.

*Implementation:* :func:`~earthsciencekit.geomorphology.knickpoint_celerity`
(*new*) and :func:`~earthsciencekit.geomorphology.simulate_river_profile`
(*new*). The example doubles the uplift rate and tracks the knickpoint.

*References:* A. D. Howard and G. Kerby, "Channel changes in badlands,"
*Geological Society of America Bulletin* 94 (1983): 739-752; K. X. Whipple
and G. E. Tucker, "Dynamics of the stream-power river incision model,"
*Journal of Geophysical Research* 104 (1999): 17661-17674.

.. minigallery:: ../../examples/geomorphology/plot_16_howard_kerby_knickpoint_migration.py

1991 -- Lal's Cosmogenic-Nuclide Erosion Rates
----------------------------------------------

Devendra Lal worked out how cosmic rays produce nuclides such as
beryllium-10 in rock near the ground, and how erosion and decay remove
them. A steadily eroding surface reaches a concentration inversely
proportional to its erosion rate, which measured in river sand gives the
erosion rate of a whole catchment.

*Implementation:* :func:`~earthsciencekit.geomorphology.cosmogenic_concentration`
(*new*), :func:`~earthsciencekit.geomorphology.cosmogenic_production_at_depth`
(*new*) and :func:`~earthsciencekit.geomorphology.steady_state_erosion_rate`
(*new*), with ``constants.LAMBDA_BE10``. The example follows accumulation to
steady state and inverts it.

*References:* D. Lal, "Cosmic ray labeling of erosion surfaces: in situ
nuclide production rates and erosion models," *Earth and Planetary Science
Letters* 104 (1991): 424-439.

.. minigallery:: ../../examples/geomorphology/plot_17_lal_cosmogenic_erosion_rates.py

See also
--------

- :doc:`/guides/geomorphology` for conventions, equations and checks.
- `Community Surface Dynamics Modeling System <https://csdms.colorado.edu/>`_
