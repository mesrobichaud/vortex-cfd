"""
blood.py -- viscosity of blood as a function of hematocrit and shear rate.

Models:

- Quemada, with hematocrit-dependent coefficients, as tabulated by Hund,
  Kameneva and Antaki (2017), "A quasi-mechanistic mathematical representation
  for blood viscosity", Fluids 2, 10, Table 2:

      mu = mu_pl * (1 - K*hct/2)^-2
      K  = (k0 + k_inf*sqrt(g/g_c)) / (1 + sqrt(g/g_c))

  with k0, k_inf and g_c functions of hematocrit. The Cokelet version gives all
  three as exp(cubic in hct). The Das version replaces k0 by a0 + 2/(a1 + hct);
  it rises with hematocrit at every shear rate, while the Cokelet version falls
  with hematocrit between about 0.26 and 0.38 at some shear rates.

- MKM5, the 5-parameter modified Krieger model (Hund et al. 2017):

      mu   = mu_pl * (1 - hct/hct*)^-n
      n    = a + b*exp(-c*hct) + n_st      (n_st only above the threshold hct_st)
      n_st = beta * (1 + (lam*g)^nu_g)^-nu

- The single-phase fits compiled by Cho and Kensey (1991), Biorheology 28.

Each model's coefficients are a named parameter set (a dict) below; a set is
chosen with the `variant` (Quemada) or `params` (MKM5) argument, by name or as a
dict. Hematocrit is a fraction (0.45, not 45). Viscosities are in mPa s.
"""

import numpy as np


MU_PLASMA_MPA_S = 1.23      # nominal plasma viscosity (Hund et al. 2017)


def _check_fraction(hct):
    if np.any(hct > 1.0):
        raise ValueError(f"hematocrit must be a fraction, not a percentage; got {hct}")


# --- Quemada -------------------------------------------------------------------------------------------------------

COKELET_K0 = (3.874, -10.41, 13.80, -6.738)
K_INF = (1.3435, -2.803, 2.711, -0.6479)
GAMMA_C = (-6.1508, 27.923, -25.60, 3.697)          # 1/s
DAS_K0 = (0.275363, 0.100158)

# Parameter sets. k0_form: "cokelet" for k0 = exp(cubic), "das" for k0 = a0 + 2/(a1 + hct).
QUEMADA_SETS = {
    "das":     {"k0_form": "das",     "k0": DAS_K0,     "k_inf": K_INF, "gamma_c": GAMMA_C},
    "cokelet": {"k0_form": "cokelet", "k0": COKELET_K0, "k_inf": K_INF, "gamma_c": GAMMA_C},
}


def _cubic(coef, x):
    return coef[0] + coef[1] * x + coef[2] * x**2 + coef[3] * x**3


def quemada(hct, gamma_dot, variant="das", mu_plasma=MU_PLASMA_MPA_S):
    """Quemada viscosity in mPa s.

    hct       : hematocrit as a fraction
    gamma_dot : shear rate, 1/s
    variant   : a key of QUEMADA_SETS ("das" or "cokelet"), or a parameter set (dict)
    Returns NaN where the model is undefined (K*hct/2 >= 1). Accepts scalars
    or NumPy arrays.
    """
    hct = np.asarray(hct, dtype=float)
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    _check_fraction(hct)

    if isinstance(variant, str):
        if variant not in QUEMADA_SETS:
            raise ValueError(f"variant must be one of {sorted(QUEMADA_SETS)}; got {variant!r}")
        params = QUEMADA_SETS[variant]
    else:
        params = variant

    if params["k0_form"] == "das":
        k0 = params["k0"][0] + 2.0 / (params["k0"][1] + hct)
    elif params["k0_form"] == "cokelet":
        k0 = np.exp(_cubic(params["k0"], hct))
    else:
        raise ValueError(f"k0_form must be 'das' or 'cokelet'; got {params['k0_form']!r}")
    k_inf = np.exp(_cubic(params["k_inf"], hct))
    gamma_c = np.exp(_cubic(params["gamma_c"], hct))

    s = np.sqrt(gamma_dot / gamma_c)
    K = (k0 + k_inf * s) / (1.0 + s)
    packing = 1.0 - 0.5 * K * hct
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = mu_plasma * packing ** -2.0
    return np.where(packing > 0.0, mu, np.nan)


# --- MKM5 ----------------------------------------------------------------------------------------------------------

# Parameter sets. Both are fits to the data of Brooks, Goodwin and Seaman (1970), measured at 25 C.
# hct_star: critical hematocrit; hct_st: threshold below which there is no shear thinning.
MKM5_HUND2017 = {   # Hund et al. (2017), Table 3. Hund et al. used hct_st = 0.20; 0.0 is the value used in this project.
    "a": 0.0, "b": 8.71, "c": 2.87, "beta": 8.23, "lam": 108.0, "nu": 0.134, "nu_g": 2,
    "hct_star": 0.99, "hct_st": 0.0,
}
MKM5_OWEN2020 = {   # Schenkel and Halliday's (2020) refit, as printed in Owen et al. (2020, 2021), Table 1
    "a": 0.0, "b": 8.78084, "c": 2.82354, "beta": 16.27775, "lam": 1252.64407, "nu": 0.14275, "nu_g": 2,
    "hct_star": 0.95, "hct_st": 0.20,
}
MKM5_SETS = {"hund2017": MKM5_HUND2017, "owen2020": MKM5_OWEN2020}


def temp_correction(T, b_t=5.95, T_ref=25.0):
    """Factor on the plasma viscosity for temperature T (degrees C), relative to T_ref (Hund et al. 2017)."""
    return ((T + 273.15) / (T_ref + 273.15)) ** (-b_t)


def mkm5(hct, gamma_dot, T, mu_plasma=MU_PLASMA_MPA_S, ramp=None, params="hund2017"):
    """Apparent blood viscosity in mPa s (MKM5, Hund et al. 2017).

    hct       : hematocrit as a fraction (0.45, not 45)
    gamma_dot : shear rate in 1/s
    T         : temperature in degrees C; the plasma viscosity is scaled from 25 C
                with temp_correction. None: no correction.
    mu_plasma : plasma viscosity, mPa s
    ramp      : None (default) for the published model, where the shear-thinning
                term switches on abruptly at hct_st. Or a pair
                (hct_low, hct_high), e.g. (0.15, 0.25): the term then rises
                smoothly from 0 at hct_low to its full value at hct_high, which
                removes the step in viscosity at hct_st. This is a modification
                of the published model.
    params    : a key of MKM5_SETS, or a parameter set (dict)
    Returns NaN at and above hct_star, where the viscosity is infinite.
    Accepts scalars or NumPy arrays.
    """
    hct = np.asarray(hct, dtype=float)
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    _check_fraction(hct)
    p = MKM5_SETS[params] if isinstance(params, str) else params

    n_inf = p["a"] + p["b"] * np.exp(-p["c"] * hct)
    n_st = p["beta"] * (1.0 + (p["lam"] * gamma_dot) ** p["nu_g"]) ** (-p["nu"])

    # Fraction of the shear-thinning term that applies, 0 to 1
    if ramp is None:
        weight = np.where(hct > p["hct_st"], 1.0, 0.0)  # shear thinning only above hct_st
    else:
        hct_low, hct_high = ramp
        s = np.clip((hct - hct_low) / (hct_high - hct_low), 0.0, 1.0)
        weight = s * s * (3.0 - 2.0 * s)                # smooth step: zero slope at both ends
    n = n_inf + weight * n_st

    correction = 1.0 if T is None else temp_correction(T)
    with np.errstate(divide="ignore", invalid="ignore"):   # hct >= hct_star is masked below
        mu = mu_plasma * correction * (1.0 - hct / p["hct_star"]) ** (-n)
    valid = hct < p["hct_star"]                         # mu is infinite at hct_star itself
    return np.where(valid, mu, np.nan)


# --- Cho and Kensey (1991) -----------------------------------------------------------------------------------------
# Single-phase fits compiled by Cho and Kensey (1991), Biorheology 28, p. 243; viscosity in mPa s, shear rate g in 1/s.

_POISE = 100.0                                          # poise -> mPa s
CHO_KENSEY_MU0, CHO_KENSEY_MU_INF = 0.56 * _POISE, 0.0345 * _POISE


def powell_eyring(g, lam=5.383):
    x = lam * g
    return CHO_KENSEY_MU_INF + (CHO_KENSEY_MU0 - CHO_KENSEY_MU_INF) * np.arcsinh(x) / x


def modified_powell_eyring(g, lam=2.415, m=1.089):
    x = lam * g
    return CHO_KENSEY_MU_INF + (CHO_KENSEY_MU0 - CHO_KENSEY_MU_INF) * np.log1p(x) / x ** m


def cross(g, lam=1.007, m=1.028, mu0=CHO_KENSEY_MU0, mu_inf=CHO_KENSEY_MU_INF):
    return mu_inf + (mu0 - mu_inf) / (1.0 + (lam * g) ** m)


def simplified_cross(g):
    """Steffan et al."""
    return cross(g, lam=8.0, m=1.0, mu0=1.30 * _POISE, mu_inf=0.05 * _POISE)


def modified_cross(g, lam=3.736, m=2.406, a=0.254):
    return CHO_KENSEY_MU_INF + (CHO_KENSEY_MU0 - CHO_KENSEY_MU_INF) / (1.0 + (lam * g) ** m) ** a


def carreau(g, lam=3.313, n=0.3568):
    return CHO_KENSEY_MU_INF + (CHO_KENSEY_MU0 - CHO_KENSEY_MU_INF) * (1.0 + (lam * g) ** 2) ** ((n - 1.0) / 2.0)


def carreau_yasuda(g, lam=1.902, n=0.22, a=1.25):
    return CHO_KENSEY_MU_INF + (CHO_KENSEY_MU0 - CHO_KENSEY_MU_INF) * (1.0 + (lam * g) ** a) ** ((n - 1.0) / a)


def casson(g, k0=0.05 * 0.1, k1=0.04 * 0.1):
    """Nakamura and Sawada. k0: yield stress, Pa (0.05 dyne/cm^2); k1: Casson viscosity, Pa s (0.04 poise)."""
    return 1e3 * (np.sqrt(k0) + np.sqrt(k1 * g)) ** 2 / g


def power_law(g, m=0.42, n=0.61):
    """Liepsch and Moravec. m in poise s^(n-1)."""
    return m * _POISE * g ** (n - 1.0)


def power_law_davies(g):
    """Davies et al.: power law up to 226.5 1/s, constant 0.04 poise above."""
    return np.where(g <= 226.5, 0.35 * _POISE * g ** (0.60 - 1.0), 0.04 * _POISE)


# --- OpenFOAM and geometry helpers ---------------------------------------------------------------------------------

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
