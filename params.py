# This is the dictionary to get the parameters.
# If you want the Pulsar Period you call P["P_NS"]
#
#
me_c2_erg = 9.1093837e-28 * (2.99792458e10)**2 # Convert energy of the electron.

P = dict(
    # --- ejecta geometry ---
    n=1.0, m=9, Mej=2.0e30, kappa=0.2, x0=0.3,
    v_star0=3.0e9, t_star=1.0e3,

    # --- central neutron star ---
    P_NS=1.0e-3,
    B_NS=1.0e12,
    I_NS=1.0e45,
    R_NS=1.0e6,
    mu=1.0,
    sin2chi=2.0 / 3.0,

    # --- particle injection ---
    gamma_min=100,      
    gamma_max=1e7,     
    eta_nth=1.0, gamma=2,
    # Emin=8.0e-6, Emax=0.8,

    # --- numerical ---
    theta_max=1.0e4,
)

P["Emin"] = P["gamma_min"] * me_c2_erg
P["Emax"] = P["gamma_max"] * me_c2_erg

# observational bands: name, E1[eV], E2[eV], band-specific opacity [cm2/g] (None -> use P["kappa"])
# We should choose specific opacities for different bands.
BANDS = [
    ["Radio (1-10 GHz)",               4.1e-6, 4.1e-5, None],
    ["Optical (1.8-3.1 eV)",           1.8,   3.1,   None],
    ["X-ray (0.3-10 keV)",             3.0e2, 1.0e4, None],
    ["Gamma (0.1-10 MeV)",             1.0e5, 1.0e7, None],
    ["High-energy gamma (0.1-10 GeV)", 1.0e8, 1.0e10, None],
]


# NTH   : number of time steps (log-spaced in theta = t/t_star) used to build the light curves. 

# NBAND : number of energy points used to integrate the synchrotron emission INSIDE each observational band (radio, optical, X-ray...)

#   NU    : number of points in the INNER integral of N_of_E (the integral along the electron characteristic).  T

NTH = 400
NBAND = 50
NU = 800