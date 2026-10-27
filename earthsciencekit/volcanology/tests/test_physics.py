"""Physical identities, limits and dynamic checks for earthsciencekit.volcanology."""

import numpy as np
import pytest
from scipy.integrate import quad

from earthsciencekit.constants import LAMBDA_U235, LAMBDA_U238, SECONDS_PER_YEAR, U238_U235
from earthsciencekit.volcanology.systems.eruptions import (
    column_height,
    plume_rise_height,
    simulate_eruption_column,
    volcanic_explosivity_index,
)
from earthsciencekit.volcanology.systems.geochemistry import (
    delta_value,
    equilibrium_crystallization_melt,
    fractional_crystallization_melt,
    isotope_temperature,
    lattice_strain_partition,
    rayleigh_fractionation,
)
from earthsciencekit.volcanology.systems.geochronology import (
    concordia_ratios,
    decay_constant,
    isochron_fit,
    lead_lead_age,
    parent_fraction,
    radiocarbon_age,
    radiogenic_daughter,
    radiogenic_lead_ratio,
    radiometric_age,
)
from earthsciencekit.volcanology.systems.magma import (
    arrhenius_viscosity,
    dike_ascent_velocity,
    mogi_displacement,
    shaw_viscosity,
    water_solubility,
)

GA = 1e9 * SECONDS_PER_YEAR


def test_decay_conserves_atoms_and_age_inverts_accumulation():
    rate = decay_constant(1.25 * GA)
    time = np.array([0.1, 1, 4.5]) * GA
    remaining = parent_fraction(rate, time)
    np.testing.assert_allclose(parent_fraction(rate, 1.25 * GA), 0.5)
    np.testing.assert_allclose(remaining + radiogenic_daughter(remaining, rate, time), 1)
    np.testing.assert_allclose(radiometric_age(radiogenic_daughter(2.0, rate, time) / 2.0, rate), time)


def test_radiocarbon_half_life_is_libby():
    np.testing.assert_allclose(radiocarbon_age(0.5), 5568 * SECONDS_PER_YEAR, rtol=1e-4)
    assert radiocarbon_age(1.2) < 0


def test_lead_ratio_limits_and_inversion():
    np.testing.assert_allclose(radiogenic_lead_ratio(1e-3 * SECONDS_PER_YEAR), LAMBDA_U235 / (U238_U235 * LAMBDA_U238), rtol=1e-9)
    ages = np.array([0.5, 2.0, 4.55]) * GA
    np.testing.assert_allclose(lead_lead_age(radiogenic_lead_ratio(ages)), ages, rtol=1e-10)
    assert np.all(np.diff(radiogenic_lead_ratio(ages)) > 0)


def test_concordia_clocks_agree_and_discordia_is_a_chord():
    time = np.array([0.3, 1.8]) * GA
    points = concordia_ratios(time)
    np.testing.assert_allclose(radiometric_age(points[:, 0], LAMBDA_U235), time)
    np.testing.assert_allclose(radiometric_age(points[:, 1], LAMBDA_U238), time)
    # A partly reset zircon lies on the straight chord between formation and lead loss.
    mixed = 0.6 * points[1] + 0.4 * points[0]
    chord, offset = points[1] - points[0], mixed - points[0]
    np.testing.assert_allclose(chord[0] * offset[1] - chord[1] * offset[0], 0, atol=1e-14)


def test_isochron_recovers_age_and_initial_ratio():
    rate = 1.42e-11 / SECONDS_PER_YEAR
    parent = np.array([0.1, 0.8, 2.5, 6.0])
    daughter = 0.7030 + parent * np.expm1(rate * 2.0 * GA)
    np.testing.assert_allclose(isochron_fit(parent, daughter, rate), [2.0 * GA, 0.7030], rtol=1e-10)


@pytest.mark.parametrize(
    "args",
    [
        [[1.0], [1.0], 1.0],
        [[1, 2], [1, 2, 3], 1.0],
        [[1, 2], [1, 2], [1.0, 2.0]],
        [[1, 1], [1, 2], 1.0],
        [[1, 2], [2, 1], 1.0],
    ],
)
def test_isochron_rejections(args):
    with pytest.raises(ValueError):
        isochron_fit(*args)


def test_rayleigh_limits_and_delta():
    np.testing.assert_allclose(rayleigh_fractionation(0.002, 1, 0.99), 0.002)
    np.testing.assert_allclose(rayleigh_fractionation(0.002, [0.1, 0.5], 1), 0.002)
    shifted = delta_value(rayleigh_fractionation(0.002, 0.5, 1.01), 0.002)
    np.testing.assert_allclose(shifted, 1000 * np.expm1(np.log(0.5) * 0.01))
    assert shifted < 0


def test_crystallization_mass_balance_and_shared_tangent():
    partition, initial = 0.2, 100.0
    melt = 0.3
    removed, _ = quad(lambda f: partition * fractional_crystallization_melt(initial, f, partition), melt, 1)
    np.testing.assert_allclose(removed + melt * fractional_crystallization_melt(initial, melt, partition), initial, rtol=1e-10)
    np.testing.assert_allclose(equilibrium_crystallization_melt(initial, [1, 0], partition), [initial, initial / partition])
    step = 1e-6
    for function in (fractional_crystallization_melt, equilibrium_crystallization_melt):
        slope = (function(initial, 1, partition) - function(initial, 1 - step, partition)) / step
        np.testing.assert_allclose(slope, (partition - 1) * initial, rtol=1e-4)
    assert fractional_crystallization_melt(initial, 0.05, partition) > equilibrium_crystallization_melt(initial, 0.05, partition)


def test_lattice_strain_peaks_at_optimal_radius():
    radius = np.linspace(0.8e-10, 1.3e-10, 2001)
    partition = lattice_strain_partition(radius, 1.05e-10, 4.0, 250e9, 1600)
    np.testing.assert_allclose(radius[np.argmax(partition)], 1.05e-10, atol=1e-13)
    np.testing.assert_allclose(partition.max(), 4.0)
    # For the same misfit, an oversized ion strains the lattice more than an undersized one.
    assert lattice_strain_partition(1.2e-10, 1.05e-10, 4.0, 250e9, 1600) < lattice_strain_partition(0.9e-10, 1.05e-10, 4.0, 250e9, 1600)


def test_isotope_temperature_reference_and_sense():
    np.testing.assert_allclose(isotope_temperature(-1.0, -1.0), 289.65)
    assert isotope_temperature(1.0, 0.0) < isotope_temperature(0.0, 0.0)


def test_viscosity_pivot_and_reference():
    np.testing.assert_allclose(shaw_viscosity(1e4 / 1.5, [1.0, 2.5, 4.0]), 0.1 * np.exp(-6.4))
    assert shaw_viscosity(1100, 4.0) > shaw_viscosity(1100, 2.5)
    np.testing.assert_allclose(arrhenius_viscosity(1300, 50, 3e5, 1300), 50)


def test_solubility_dike_and_mogi_scaling():
    np.testing.assert_allclose(water_solubility(4e8, 4e-6) / water_solubility(1e8, 4e-6), 2)
    np.testing.assert_allclose(dike_ascent_velocity(4, 300, 100, 9.8) / dike_ascent_velocity(2, 300, 100, 9.8), 4)
    depth, volume, nu = 3000.0, 2e6, 0.25
    half = mogi_displacement(depth * np.sqrt(2 ** (2 / 3) - 1), depth, volume, nu)[1]
    np.testing.assert_allclose(half, mogi_displacement(0, depth, volume, nu)[1] / 2)
    uplift, _ = quad(lambda r: 2 * np.pi * r * mogi_displacement(r, depth, volume, nu)[1], 0, np.inf)
    np.testing.assert_allclose(uplift, 2 * (1 - nu) * volume, rtol=1e-8)


def test_plume_scaling_and_vei_bins():
    np.testing.assert_allclose(plume_rise_height(2e9, 0.01) / plume_rise_height(1e9, 0.01), 2**0.25)
    np.testing.assert_allclose(column_height(16e12) / column_height(1e12), 2)
    volume = [0, 9.9e3, 1e4, 9.9e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e15]
    np.testing.assert_array_equal(volcanic_explosivity_index(volume), [0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 8])


def test_eruption_column_follows_pure_plume_similarity():
    alpha, buoyancy_flux, origin = 0.12, 1.0, 1.0
    speed = 5 / (6 * alpha) * (0.9 * alpha * buoyancy_flux) ** (1 / 3)
    radius, velocity = 1.2 * alpha * origin, speed * origin ** (-1 / 3)
    time = np.linspace(0, 30, 61)
    result = simulate_eruption_column(radius, velocity, buoyancy_flux / (radius**2 * velocity), 0, time, entrainment=alpha)
    exact = (origin ** (4 / 3) + 4 / 3 * speed * time) ** 0.75
    np.testing.assert_allclose(result.values[:, 0], exact - origin, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(result.values[:, 1], speed * exact ** (-1 / 3), rtol=1e-6)


def test_eruption_column_reaches_mtt_height_and_settles():
    frequency = 0.02
    time = np.linspace(5, 3005, 1001)
    result = simulate_eruption_column(1, 1, 1, frequency, time)
    top = result.values[:, 0].max()
    np.testing.assert_allclose(top, plume_rise_height(1, frequency), rtol=0.01)
    assert 0.7 < result.values[-1, 0] / top < 0.9
    assert result.values[0, 0] == 0 and result.time[0] == 5
    with pytest.raises(ValueError):
        simulate_eruption_column(1, [1, 2], 1, frequency, time)


def test_scriven_growth_constant_limits_and_closed_form():
    from scipy.special import erfcx

    from earthsciencekit.volcanology import scriven_growth_constant

    # Without advection (gas as dense as the melt) the integral has a closed form.
    beta = np.array([0.05, 0.2, 0.7, 2.0, 6.0])
    supersaturation = 2 * beta**2 - 2 * np.sqrt(np.pi) * beta**3 * erfcx(beta)
    np.testing.assert_allclose(scriven_growth_constant(supersaturation, 1), beta, rtol=1e-9)
    np.testing.assert_allclose(scriven_growth_constant(1e-6, [0, 0.5, 1]), np.sqrt(0.5e-6), rtol=1e-3)
    # A light gas at large supersaturation follows Plesset and Zwick's sqrt(3 / pi) law.
    np.testing.assert_allclose(scriven_growth_constant(300, 0), np.sqrt(3 / np.pi) * 300, rtol=2e-3)
    assert np.all(np.diff(scriven_growth_constant(np.linspace(0.1, 5, 12), 0.05)) > 0)


def test_bubble_cell_grows_like_scriven_then_saturates():
    from earthsciencekit.constants import R_VAPOR
    from earthsciencekit.volcanology import (
        diffusive_supersaturation,
        scriven_bubble_radius,
        scriven_growth_constant,
        simulate_bubble_growth,
        water_solubility,
    )

    pressure, temperature, water = 5e7, 1100.0, 0.058
    surface = water_solubility(pressure, 4.1e-6)
    gas = pressure / (R_VAPOR * temperature)
    supersaturation = diffusive_supersaturation(water, surface, 2300, gas)
    time = np.linspace(0, 20, 11)
    # A huge melt shell is never depleted, so the quasi-static growth follows R**2 = R0**2 + 2 phi D t.
    unbounded = simulate_bubble_growth(1e-9, 1.0, water, pressure, temperature, time)
    np.testing.assert_allclose(unbounded.values[:, 0] ** 2, 1e-18 + 2 * supersaturation * 1e-11 * time, rtol=1e-6)
    # Scriven's full solution approaches it only when the supersaturation is small; here it grows 70 % faster.
    beta = scriven_growth_constant(supersaturation, gas / 2300)
    np.testing.assert_allclose(scriven_bubble_radius(time[1:], 1e-11, beta) / unbounded.values[1:, 0], 1.70, rtol=0.01)
    slight = simulate_bubble_growth(1e-9, 1.0, surface + 5e-4, pressure, temperature, time)
    beta = scriven_growth_constant(diffusive_supersaturation(surface + 5e-4, surface, 2300, gas), gas / 2300)
    np.testing.assert_allclose(scriven_bubble_radius(time[1:], 1e-11, beta), slight.values[1:, 0], rtol=0.08)
    cell = simulate_bubble_growth(1e-6, 6.2e-5, water, pressure, temperature, np.linspace(0, 5000, 51))
    radius, dissolved, vesicularity = cell.values.T
    assert np.all(np.diff(radius) > -1e-12) and np.isclose(dissolved[-1], surface, rtol=1e-6)
    # Water lost from the melt is gained by the bubble.
    melt_mass = 2300 * (6.2e-5**3 - 1e-6**3)
    np.testing.assert_allclose(melt_mass * (water - dissolved), gas * (radius**3 - 1e-6**3), rtol=1e-10)
    final = (1e-6**3 + melt_mass * (water - surface) / gas) ** (1 / 3)
    np.testing.assert_allclose(radius[-1], final, rtol=1e-6)
    assert 0 < vesicularity[0] < vesicularity[-1] < 1 and cell.names == ("radius", "dissolved_water", "vesicularity")


@pytest.mark.parametrize(
    "args",
    [
        (1e-6, 1e-7, 0.05, 5e7, 1100, [0, 1]),
        (1e-6, 1e-4, 0.01, 5e7, 1100, [0, 1]),
        ([1e-6, 2e-6], 1e-4, 0.05, 5e7, 1100, [0, 1]),
        (1e-6, 1e-4, 1.0, 5e7, 1100, [0, 1]),
    ],
)
def test_bubble_growth_rejects_bad_setup(args):
    from earthsciencekit.volcanology import simulate_bubble_growth

    with pytest.raises(ValueError):
        simulate_bubble_growth(*args)
