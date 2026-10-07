"""
blood.py -- viscosity of blood as a function of hematocrit and shear rate.

Quemada model with hematocrit-dependent coefficients, as tabulated by Hund,
Kameneva and Antaki (2017), "A quasi-mechanistic mathematical representation
for blood viscosity", Fluids 2, 10, Table 2:

    mu = mu_pl * (1 - K*hct/2)^-2
    K  = (k0 + k_inf*sqrt(g/g_c)) / (1 + sqrt(g/g_c))

with k0, k_inf and g_c functions of hematocrit. The Cokelet version gives all
three as exp(cubic in hct). The Das version replaces k0 by a0 + 2/(a1 + hct);
it rises with hematocrit at every shear rate, while the Cokelet version falls
with hematocrit between about 0.26 and 0.38 at some shear rates.

Hematocrit is a fraction (0.45, not 45). Viscosities are in mPa s.
"""

import numpy as np


MU_PLASMA_MPA_S = 1.23      # nominal plasma viscosity (Hund et al. 2017)

COKELET_K0 = (3.874, -10.41, 13.80, -6.738)
K_INF = (1.3435, -2.803, 2.711, -0.6479)
GAMMA_C = (-6.1508, 27.923, -25.60, 3.697)          # 1/s
DAS_K0 = (0.275363, 0.100158)


def _cubic(coef, x):
    return coef[0] + coef[1] * x + coef[2] * x**2 + coef[3] * x**3


def quemada(hct, gamma_dot, variant="das", mu_plasma=MU_PLASMA_MPA_S):
    """Quemada viscosity in mPa s.

    hct       : hematocrit as a fraction
    gamma_dot : shear rate, 1/s
    variant   : "das" or "cokelet" (the form of k0)
    Returns NaN where the model is undefined (K*hct/2 >= 1). Accepts scalars
    or NumPy arrays.
    """
    hct = np.asarray(hct, dtype=float)
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    if np.any(hct > 1.0):
        raise ValueError(f"hematocrit must be a fraction, not a percentage; got {hct}")

    if variant == "das":
        k0 = DAS_K0[0] + 2.0 / (DAS_K0[1] + hct)
    elif variant == "cokelet":
        k0 = np.exp(_cubic(COKELET_K0, hct))
    else:
        raise ValueError(f"variant must be 'das' or 'cokelet'; got {variant!r}")
    k_inf = np.exp(_cubic(K_INF, hct))
    gamma_c = np.exp(_cubic(GAMMA_C, hct))

    s = np.sqrt(gamma_dot / gamma_c)
    K = (k0 + k_inf * s) / (1.0 + s)
    packing = 1.0 - 0.5 * K * hct
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = mu_plasma * packing ** -2.0
    return np.where(packing > 0.0, mu, np.nan)


def couette_wall_shear_rate(omega, Ri, Ro):
    """Shear rate at the inner (rotating) wall of a Couette gap, 1/s.
    Newtonian: 2*omega*Ro^2 / (Ro^2 - Ri^2)."""
    return 2.0 * omega * Ro**2 / (Ro**2 - Ri**2)


def openfoam_nu_table(hct, rho, gamma_min=1e-3, gamma_max=1e5, per_decade=20, **quemada_kw):
    """Kinematic viscosity table for OpenFOAM's strainRateFunction model.

    Returns a one-line OpenFOAM list "((g1 nu1) (g2 nu2) ...)" of shear rate
    (1/s) and kinematic viscosity (m^2/s), log-spaced. OpenFOAM interpolates
    linearly between points and clamps outside the range, which is harmless
    here because the Quemada viscosity is finite as the shear rate tends to 0.
    """
    decades = np.log10(gamma_max) - np.log10(gamma_min)
    gamma = np.logspace(np.log10(gamma_min), np.log10(gamma_max), int(round(decades * per_decade)) + 1)
    mu = quemada(hct, gamma, **quemada_kw)
    if np.any(np.isnan(mu)):
        raise ValueError(f"the Quemada model is undefined for hematocrit {hct}")
    nu = mu * 1e-3 / rho
    return "(" + " ".join(f"({g:.6g} {n:.6g})" for g, n in zip(gamma, nu)) + ")"
