import numpy as np
from physics import N_of_E, region_geometry

def ssa_factor(Egrid, t, B, P, D, h=1e-3, nu_sub=600):
    dE = Egrid * h
    Nm = N_of_E(Egrid - dE, t, P, D, nu=nu_sub)
    Np = N_of_E(Egrid + dE, t, P, D, nu=nu_sub)
    dg_dE = (Np / (Egrid + dE)**2 - Nm / (Egrid - dE)**2) / (2 * dE)
    nu = D["alpha"] * B * Egrid**2
    Psyn = D["beta"] * B**2 * Egrid**2
    dEdnu = 1.0 / (2 * D["alpha"] * B * Egrid)
    Rs, path, V = region_geometry(t, P, D)
    alpha_nu = np.maximum(-(2.99792458e10**2 / (8 * np.pi * nu**2)) * Psyn * Egrid**2 * dg_dE * dEdnu / V, 0.0)
    tau_ssa = alpha_nu * path
    f = np.where(tau_ssa < 1e-6, 1.0 - tau_ssa / 2, -np.expm1(-tau_ssa) / np.maximum(tau_ssa, 1e-300))
    return f, tau_ssa