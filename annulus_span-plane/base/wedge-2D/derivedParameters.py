# Run by pyFoamPrepareCase.py after default.parameters is read. Every
# parameter is already a variable here; any new variable becomes a parameter
# the templates can use. Temporary names are deleted at the end so they don't
# leak into the parameter record.

from cfdtools import units, blood

Ri = Ri_mm * units.MM
Ro = Ro_mm * units.MM
H = H_mm * units.MM
Hb = Hb_mm * units.MM
edge_radius = edge_radius_mm * units.MM

omega = units.rpm_to_rad_s(speed_rpm)

# Viscosity. nu is the constant kinematic viscosity for the Newtonian models;
# for the shear-thinning model it is the value at the bob-wall shear rate, kept
# as a reference (e.g. for Taylor-number estimates). nu_table is only used by
# the shear-thinning model.
gamma_ref = blood.couette_wall_shear_rate(omega, Ri, Ro)
nu_table = ""
if viscosity_model == "newtonian":
    nu = units.kinematic_viscosity(mu_mPa_s, rho_kg_m3)
elif viscosity_model in ("quemada_newtonian", "quemada"):
    mu_ref_mPa_s = float(blood.quemada(hct, gamma_ref))
    nu = units.kinematic_viscosity(mu_ref_mPa_s, rho_kg_m3)
    if viscosity_model == "quemada":
        nu_table = blood.openfoam_nu_table(hct, rho_kg_m3)
else:
    raise ValueError(f"unknown viscosity_model {viscosity_model!r}; "
                     "use newtonian, quemada_newtonian or quemada")

# The mesh is one wedge of angle_deg; multiply its forces and torques by this
# to get the full 360-degree device.
wedge_to_full = 360.0 / angle_deg

del units, blood
