Breakthroughs in Seismology
===========================

.. include:: /_generated/nav/seismology.rst

Seismology began as the theory of waves in elastic solids, decades before
there were instruments to record them. Once seismographs spread around 1900,
arrival times revealed the crust, the core and the inner core within less
than 30 years. The same records gave the first earthquake scales and the
statistics of aftershocks and catalogues. After 1960, wave theory and
earthquake statistics met in the seismic moment, a physical measure of
earthquake size, and travel times were fitted by a single reference Earth.
This chronology traces the results behind :mod:`earthsciencekit.seismology`.
Every example uses deterministic synthetic inputs or rounded published
values; none inverts real records.

.. contents:: Timeline
   :local:
   :depth: 1

1829 -- Poisson's Two Elastic Waves
-----------------------------------

Siméon Denis Poisson showed that the equations of an elastic solid allow two
kinds of wave, a faster compressional wave and a slower one transverse to the
direction of travel. For a solid whose two Lamé constants are equal, the
speed ratio is :math:`\sqrt{3}` and Poisson's ratio is 1/4, close to many
crustal rocks. Without rigidity the shear wave disappears, which later
showed that the outer core is liquid.

*Implementation:* :func:`~earthsciencekit.seismology.p_wave_speed`,
:func:`~earthsciencekit.seismology.s_wave_speed` and
:func:`~earthsciencekit.seismology.poisson_ratio`.

*References:* S.-D. Poisson, "Mémoire sur l'équilibre et le mouvement des
corps élastiques," *Mémoires de l'Académie Royale des Sciences de l'Institut
de France* 8 (1829): 357-570.

.. minigallery:: ../../examples/seismology/plot_01_poisson_two_elastic_waves.py

1885 -- Rayleigh's Surface Waves
--------------------------------

Lord Rayleigh found a third wave, confined to the free surface of an elastic
half-space and decaying with depth. Its speed is the root of a cubic that
depends only on Poisson's ratio, about 0.92 of the shear speed in a Poisson
solid. He predicted that, because it spreads in two dimensions rather than
three, it would dominate distant earthquake records, as it does.

*Implementation:* :func:`~earthsciencekit.seismology.rayleigh_wave_speed`
(*new*). The tests check the root against Rayleigh's secular equation and
Viktorov's approximation.

*References:* Lord Rayleigh, "On waves propagated along the plane surface of
an elastic solid," *Proceedings of the London Mathematical Society* 17
(1885): 4-11.

.. minigallery:: ../../examples/seismology/plot_02_rayleigh_surface_waves.py

1894 -- Omori's Aftershock Decay
--------------------------------

Fusakichi Omori counted the aftershocks of the 1891 Nobi earthquake and found
that their rate fell as :math:`K/(t + c)`. The total number therefore grows
only logarithmically with time, and a sequence can last for decades. Tokuji
Utsu generalized the exponent to :math:`p`, usually between 0.9 and 1.5, in
1961.

*Implementation:* :func:`~earthsciencekit.seismology.omori_rate`, in Utsu's
modified form. The example checks the late-time slope and the logarithmic
count.

*References:* F. Omori, "On the after-shocks of earthquakes," *Journal of the
College of Science, Imperial University of Tokyo* 7 (1894): 111-200.

.. minigallery:: ../../examples/seismology/plot_03_omori_aftershock_decay.py

1909 -- Mohorovičić's Discontinuity
-----------------------------------

On records of the Kupa Valley earthquake of 8 October 1909, Andrija
Mohorovičić identified two P phases. Close to the source the direct wave came
first, but farther out a faster phase overtook it. He explained it as a wave
refracted along a boundary where the speed jumps, at a depth of about 54 km.
The crossover distance fixes that depth, and the boundary under every
continent and ocean now bears his name.

*Implementation:* :func:`~earthsciencekit.seismology.travel_time`,
:func:`~earthsciencekit.seismology.head_wave_travel_time` (*new*) and
:func:`~earthsciencekit.seismology.crossover_distance` (*new*).

*References:* A. Mohorovičić, "Potres od 8. X 1909," *Godišnje izvješće
Zagrebačkog meteorološkog opservatorija za godinu 1909* 9, no. 4 (1910):
1-56.

.. minigallery:: ../../examples/seismology/plot_04_mohorovicic_discontinuity.py

1910 -- Reid's Elastic Rebound
------------------------------

Harry Fielding Reid compared triangulation surveys across the San Andreas
fault from before and after the 1906 San Francisco earthquake. He concluded
that the crust on either side bends elastically while the fault is locked,
storing the far-field motion as strain. The earthquake releases it in a
sudden slip of several metres, and the cycle starts again.

*Implementation:* :func:`~earthsciencekit.seismology.simulate_oscillator`,
used in the example for the coseismic release of the stored offset.

*References:* H. F. Reid, *The Mechanics of the Earthquake*, vol. 2 of *The
California Earthquake of April 18, 1906: Report of the State Earthquake
Investigation Commission* (Washington: Carnegie Institution of Washington,
1910).

.. minigallery:: ../../examples/seismology/plot_05_reid_elastic_rebound.py

1911 -- Love Waves
------------------

Seismograms showed large horizontal shear motions that Rayleigh waves could
not explain. Augustus Edward Hough Love showed that shear waves polarized
parallel to the surface can be trapped in a slow layer over a faster
half-space. Their speed depends on the period: short waves travel at the
speed of the layer and long waves at that of the half-space. That dispersion
became a standard way to measure the thickness of the crust.

*Implementation:* :func:`~earthsciencekit.seismology.love_wave_angular_frequency`
(*new*), the fundamental mode. The tests substitute it back into the
dispersion relation.

*References:* A. E. H. Love, *Some Problems of Geodynamics* (Cambridge:
Cambridge University Press, 1911).

.. minigallery:: ../../examples/seismology/plot_06_love_waves.py

1914 -- Gutenberg's Core Shadow
-------------------------------

Beno Gutenberg located the boundary of the Earth's core at a depth of about
2900 km, within a few kilometres of the modern value. Waves that enter a
core slower than the mantle bend towards the centre. Direct P waves therefore
fade near 103 degrees from the source, and core waves return only beyond
about 143 degrees, leaving a shadow zone between. On a sphere Snell's law
becomes the conservation of the ray parameter :math:`r \sin i / v`.

*Implementation:* :func:`~earthsciencekit.seismology.ray_parameter` (*new*).
The example traces straight rays through a two-shell toy Earth.

*References:* B. Gutenberg, "Über Erdbebenwellen. VII A. Beobachtungen an
Registrierungen von Fernbeben in Göttingen und Folgerungen über die
Konstitution des Erdkörpers," *Nachrichten von der Gesellschaft der
Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse* (1914):
125-176.

.. minigallery:: ../../examples/seismology/plot_07_gutenberg_core_shadow.py

1919 -- Zoeppritz's Reflection Coefficients
-------------------------------------------

Karl Zoeppritz derived how a plane wave divides into reflected and
transmitted P and S waves at a boundary between two elastic media. His
equations (published after his death) underlie reflection seismology. At
normal incidence they reduce to the contrast in impedance, density times
speed. The pressure transmission can then exceed one, but the energy
fractions still sum to one.

*Implementation:* :func:`~earthsciencekit.seismology.seismic_impedance`,
:func:`~earthsciencekit.seismology.reflection_coefficient` and
:func:`~earthsciencekit.seismology.transmission_coefficient`, at normal
incidence only.

*References:* K. Zoeppritz, "Erdbebenwellen VII. VIIb. Über Reflexion und
Durchgang seismischer Wellen durch Unstetigkeitsflächen," *Nachrichten von
der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische
Klasse* (1919): 66-84.

.. minigallery:: ../../examples/seismology/plot_08_zoeppritz_reflection.py

1935 -- Richter's Magnitude Scale
---------------------------------

Charles Richter ranked southern California earthquakes by the logarithm of
the peak amplitude on a standard Wood-Anderson seismometer, corrected for
distance. Magnitude 3 is a trace of 1 mm at 100 km, and each unit is ten
times the amplitude. It was the first instrumental measure of earthquake
size, and every later magnitude was calibrated against it.

*Implementation:* :func:`~earthsciencekit.seismology.local_magnitude`
(*new*), with Hutton and Boore's 1987 distance correction in place of
Richter's table.

*References:* C. F. Richter, "An instrumental earthquake magnitude scale,"
*Bulletin of the Seismological Society of America* 25 (1935): 1-32.

.. minigallery:: ../../examples/seismology/plot_09_richter_magnitude.py

1936 -- Lehmann's Inner Core
----------------------------

Inge Lehmann noticed weak P arrivals inside the core shadow zone, where none
should be. She explained them with an inner core in which waves travel
faster than in the outer core. Rays that meet it beyond the critical angle
are totally reflected and emerge in the shadow. Its radius, about 1220 km,
was later fixed by reflections and free oscillations.

*Implementation:* :func:`~earthsciencekit.seismology.critical_angle`. The
example shows that post-critical reflections from a toy inner core cover the
shadow zone.

*References:* I. Lehmann, "P'," *Publications du Bureau Central
Séismologique International, Série A, Travaux Scientifiques* 14 (1936):
87-115.

.. minigallery:: ../../examples/seismology/plot_10_lehmann_inner_core.py

1940 -- Jeffreys and Bullen's Travel-Time Tables
------------------------------------------------

Harold Jeffreys and Keith Bullen fitted the arrival times of thousands of
earthquakes to a spherically symmetric Earth. They published travel times of
P, S and the core phases out to 180 degrees, and their tables located
earthquakes for half a century. A ray keeps its ray parameter
:math:`r \sin i / v` all the way down and back, so tracing rays through a
layered sphere gives the travel-time curve, and the slope of that curve is
the ray parameter.

*Implementation:* :func:`~earthsciencekit.seismology.layered_ray_distance`
(*new*) and :func:`~earthsciencekit.seismology.layered_ray_travel_time`
(*new*), exact for speeds that follow the Mohorovičić law ``v = a r^b``
between tabulated depths. The example traces P and core phases through
rounded PREM-like speeds and finds the shadow zone.

*References:* H. Jeffreys and K. E. Bullen, *Seismological Tables* (London:
British Association for the Advancement of Science, Gray-Milne Trust, 1940).

.. minigallery:: ../../examples/seismology/plot_11_jeffreys_bullen_travel_times.py

1944 -- Gutenberg and Richter's Frequency-Magnitude Law
-------------------------------------------------------

Beno Gutenberg and Charles Richter counted California earthquakes by
magnitude and found that the number at or above magnitude :math:`M` falls as
:math:`10^{a - bM}`, with :math:`b` close to 1. Each magnitude unit brings
about ten times fewer events. The law is still the starting point of
seismic hazard estimates.

*Implementation:* :func:`~earthsciencekit.seismology.gutenberg_richter`. The
example recovers :math:`b` from a synthetic catalogue with Aki's 1965
maximum-likelihood estimator.

*References:* B. Gutenberg and C. F. Richter, "Frequency of earthquakes in
California," *Bulletin of the Seismological Society of America* 34 (1944):
185-188.

.. minigallery:: ../../examples/seismology/plot_12_gutenberg_richter_law.py

1953 -- Ricker's Wavelet
------------------------

Norman Ricker studied how a sharp pulse changes shape as it travels through
absorbing rock and proposed a symmetric wavelet, the second derivative of a
Gaussian, as its stable form. It has zero mean, a unit peak and a spectrum
that peaks at one frequency. It remains the default source pulse in
synthetic seismograms.

*Implementation:* :func:`~earthsciencekit.seismology.ricker_wavelet`. The
example checks the zero crossings, the zero mean and the spectral peak.

*References:* N. Ricker, "The form and laws of propagation of seismic
wavelets," *Geophysics* 18 (1953): 10-40.

.. minigallery:: ../../examples/seismology/plot_13_ricker_wavelet.py

1964 -- Knopoff's Q
-------------------

Leon Knopoff reviewed the evidence that rocks lose a nearly constant
fraction of wave energy per cycle over a wide band of frequencies. With a
constant quality factor :math:`Q`, amplitude decays as
:math:`e^{-\pi f t / Q}`, so high frequencies die first. The logarithm of a
ratio of two spectra is then linear in frequency, and its slope measures
:math:`Q`.

*Implementation:* :func:`~earthsciencekit.seismology.attenuated_amplitude`.
The example recovers :math:`Q` by the spectral-ratio method.

*References:* L. Knopoff, "Q," *Reviews of Geophysics* 2 (1964): 625-660.

.. minigallery:: ../../examples/seismology/plot_14_knopoff_q.py

1966 -- Aki's Seismic Moment
----------------------------

Keiiti Aki measured the long-period waves of the 1964 Niigata earthquake and
showed that they fix one number, the seismic moment: rigidity times fault
area times average slip. Unlike a magnitude read from one instrument, the
moment is a physical quantity that does not saturate for great earthquakes.
Aki's Niigata value was about :math:`3 \times 10^{20}` N m.

*Implementation:* :func:`~earthsciencekit.seismology.seismic_moment`.

*References:* K. Aki, "Generation and propagation of G waves from the
Niigata earthquake of June 16, 1964. Part 2. Estimation of earthquake moment,
released energy, and stress-strain drop from the G wave spectrum,"
*Bulletin of the Earthquake Research Institute* 44 (1966): 73-88.

.. minigallery:: ../../examples/seismology/plot_15_aki_seismic_moment.py

1979 -- Hanks and Kanamori's Moment Magnitude
---------------------------------------------

Thomas Hanks and Hiroo Kanamori defined a magnitude from the seismic moment,
:math:`M_w = \tfrac{2}{3}\log_{10} M_0 - 10.7` with the moment in dyne cm.
It agrees with the older scales where they work and keeps growing for great
earthquakes. Each magnitude unit multiplies the moment by :math:`10^{1.5}`,
about 32. In SI units the IASPEI standard writes it as
:math:`\tfrac{2}{3}(\log_{10} M_0 - 9.1)`.

*Implementation:* :func:`~earthsciencekit.seismology.moment_magnitude` and
:func:`~earthsciencekit.seismology.moment_from_magnitude`.

*References:* T. C. Hanks and H. Kanamori, "A moment magnitude scale,"
*Journal of Geophysical Research* 84 (1979): 2348-2350.

.. minigallery:: ../../examples/seismology/plot_16_hanks_kanamori_moment_magnitude.py

1981 -- Dziewonski and Anderson's PREM
--------------------------------------

Adam Dziewonski and Don Anderson fitted a single spherically symmetric model
of the Earth, the Preliminary Reference Earth Model, to body-wave travel
times, surface-wave dispersion, free oscillations and the Earth's mass and
moment of inertia. Its tables of speed and density against depth are still
the standard reference. Integrating slowness through them predicts arrival
times.

*Implementation:* :func:`~earthsciencekit.seismology.vertical_travel_time`
(*new*), exact for speeds linear between tabulated depths. The example
integrates a rounded PREM P-wave profile from the surface to the centre.

*References:* A. M. Dziewonski and D. L. Anderson, "Preliminary reference
Earth model," *Physics of the Earth and Planetary Interiors* 25 (1981):
297-356.

.. minigallery:: ../../examples/seismology/plot_17_prem_travel_time.py

See also
--------

- :doc:`/guides/seismology` for conventions, equations and checks.
- `USGS: Magnitude, energy release and intensity <https://www.usgs.gov/programs/earthquake-hazards/earthquake-magnitude-energy-release-and-shaking-intensity>`_
- `MIT: Essentials of Geophysics lecture notes <https://ocw.mit.edu/courses/12-201-essentials-of-geophysics-fall-2004/pages/lecture-notes/>`_
