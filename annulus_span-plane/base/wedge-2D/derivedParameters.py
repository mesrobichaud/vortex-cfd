# Run by pyFoamPrepareCase.py after default.parameters is read. Every
# parameter is already a variable here; any new variable becomes a parameter
# the templates can use. Temporary names are deleted at the end so they don't
# leak into the parameter record.

from cfdtools import units

Ri = Ri_mm * units.MM
Ro = Ro_mm * units.MM
H = H_mm * units.MM
Hb = Hb_mm * units.MM
edge_radius = edge_radius_mm * units.MM

omega = units.rpm_to_rad_s(speed_rpm)
nu = units.kinematic_viscosity(mu_mPa_s, rho_kg_m3)

# The mesh is one wedge of angle_deg; multiply its forces and torques by this
# to get the full 360-degree device.
wedge_to_full = 360.0 / angle_deg

del units
