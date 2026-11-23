# These are quantities that where derived from the "free" parameters. 
import numpy as np

def derived(P):
    """Quantities fixed by the parameters."""
    n, mu, m, g = P["n"], P["mu"], P["m"], P["gamma"]
    D = {}
    D["beta"] = 2 * 4.80320471e-10**4 * P["sin2chi"] / (3 * 9.1093837e-28**4 * 2.99792458e10**7)
    D["alpha"] = 3 * 4.80320471e-10 / (4 * np.pi * 9.1093837e-28**3 * 2.99792458e10**5)
    D["a"] = n * (1 + 2 * mu) - 1
    D["Rs0"] = P["v_star0"] * P["t_star"] / n
    D["bst"] = P["v_star0"] / 2.99792458e10

    D["k"] = 2.0
    D["L0"] = 8 * np.pi**4 * P["B_NS"]**2 * P["R_NS"]**6 / (3 * 2.99792458e10**3 * P["P_NS"]**4)
    D["tq"] = 3 * P["I_NS"] * 2.99792458e10**3 * P["P_NS"]**2 / (4 * np.pi**2 * P["B_NS"]**2 * P["R_NS"]**6)

    D["B0"] = P["B_NS"] * (P["R_NS"] / D["Rs0"])**mu
    D["B0_poynting_xcheck"] = np.sqrt(D["L0"] / (2.99792458e10 * D["Rs0"]**2))

    D["c0"] = D["beta"] * D["B0"]**2 * P["t_star"]

    x0 = P["x0"]
    D["tau1"] = ((3 - m) / (1 - m)) * P["kappa"] * P["Mej"] / (4 * np.pi * D["Rs0"]**2) \
        * (1 - x0**(1 - m)) / (1 - x0**(3 - m))
    D["q0"] = (2 - g) * P["eta_nth"] * D["L0"] / (P["Emax"]**(2 - g) - P["Emin"]**(2 - g))
    return D


def print_params(P, D):
    print("%-10s %-14s %s" % ("symbol", "value", "description"))
    for k, v in P.items():
        print("%-10s %-14s" % (k, v))
    print("\nderived:")
    print("  L0  = %.4e erg/s" % D["L0"])
    print("  tq  = %.4e s = %.3f days = %.3f yr" % (D["tq"], D["tq"] / 86400.0, D["tq"] / 86400.0 / 365.25))
    print("  R*,0 = %.3e cm" % D["Rs0"])
    print("  B*,0 = %.4e G" % D["B0"])
    print("  B*,0 (Poynting-equipartition cross-check, NOT used) = %.4e G" % D["B0_poynting_xcheck"])
    print("  tau(theta=1) = %.4e" % D["tau1"])
    print("  a = %.3g   q0 = %.3e" % (D["a"], D["q0"]))
    print("  total spin-down energy budget = %.3e erg" %
          (P["eta_nth"] * D["L0"] * D["tq"] / (D["k"] - 1)))