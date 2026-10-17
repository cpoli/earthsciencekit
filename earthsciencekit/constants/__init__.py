"""Explicit reference constants, not an evolving observational dataset.

SI units apply. Earth angular speed and material properties are educational
reference values. EPSILON_WATER is the ratio of dry-air to vapor gas constants.
"""

G = 6.67430e-11  # m^3 kg^-1 s^-2
MU0 = 1.2566370614359173e-6  # H/m, conventional nonmagnetic approximation
R_GAS = 8.31446261815324  # J/(mol K)
AVOGADRO = 6.02214076e23  # 1/mol
R_DRY_AIR = 287.05  # J/(kg K)
CP_DRY_AIR = 1004.0  # J/(kg K)
R_VAPOR = 461.5  # J/(kg K)
EPSILON_WATER = R_DRY_AIR / R_VAPOR
OMEGA_EARTH = 7.2921150e-5  # rad/s
SIGMA_SB = 5.670374419e-8  # W/(m^2 K^4)
GRAVITY = 9.80665  # m/s^2, standard gravity
EARTH_RADIUS = 6371008.8  # m, mean radius reference
SECONDS_PER_YEAR = 365.25 * 86400  # s, Julian year
LAMBDA_U238 = 1.55125e-10 / SECONDS_PER_YEAR  # 1/s, Jaffey et al. (1971)
LAMBDA_U235 = 9.8485e-10 / SECONDS_PER_YEAR  # 1/s, Jaffey et al. (1971)
U238_U235 = 137.88  # 1, present-day atomic ratio, Steiger and Jaeger (1977)
LIBBY_MEAN_LIFE = 8033 * SECONDS_PER_YEAR  # s, conventional radiocarbon mean life (5568-year half-life)
