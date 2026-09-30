"""
units.py -- laboratory units to SI.

Parameters are kept in the units they are measured in (rpm, mm, mPa.s) and
converted here, in one place, so a case file never carries a hand-converted
number.
"""

from math import pi


MM = 1e-3                    # mm -> m
UM = 1e-6                    # um -> m
RPM = 2.0 * pi / 60.0        # rpm -> rad/s
MPA_S = 1e-3                 # mPa.s (= cP) -> Pa.s


def rpm_to_rad_s(rpm):
    """Rotational speed, rpm -> rad/s."""
    return rpm * RPM


def kinematic_viscosity(mu_mPa_s, rho_kg_m3):
    """nu in m^2/s from dynamic viscosity in mPa.s and density in kg/m^3."""
    if rho_kg_m3 <= 0.0:
        raise ValueError(f"density must be positive; got {rho_kg_m3}")
    return mu_mPa_s * MPA_S / rho_kg_m3
