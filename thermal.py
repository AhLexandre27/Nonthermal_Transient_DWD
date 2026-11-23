# thermal.py
import numpy as np
from scipy.integrate import solve_ivp, quad
from params import BANDS

def thermal_and_check(th, Lnth, P, D):
    n = P["n"]; ts = P["t_star"]
    tau = lambda x: D["tau1"] * x**(-2 * n)
    logL = lambda x: np.interp(np.log(x), np.log(th), Lnth)
    Escale = ts * D["L0"]

    def rhs(x, y):
        Lx = logL(x)
        tu = tau(x)
        fab = -np.expm1(-tu)
        dth = ts * fab * Lx / Escale - (n / x) * y[0] * (1 + x**(1 - n) / (D["bst"] * (1 + tu)))
        return [dth]

    sol = solve_ivp(rhs, (th[0], th[-1]), [0.0], t_eval=th, method="LSODA",
                    rtol=1e-8, atol=1e-16)
    Eth = sol.y[0] * Escale
    Rs = D["Rs0"] * th**n
    Lth = 2.99792458e10 * Eth / (Rs * (1 + tau(th)))
    return Eth, Lth, tau(th)


def thermal_temperature(Lth_out, Rs):
    Lth_out = np.maximum(Lth_out, 0.0)
    return (Lth_out / (4 * np.pi * Rs**2 * 5.6703744e-5 + 1e-300))**0.25


def planck_band_fraction(x1, x2):
    x2c = min(x2, 700.0)
    if x2c <= x1 or x1 < 0:
        return 0.0
    val, _ = quad(lambda x: x**3 / np.expm1(x), x1, x2c, limit=200)
    return val / (np.pi**4 / 15.0)


def thermal_band_luminosities(th, Lth, Rs, P):
    Teff = thermal_temperature(Lth, Rs)
    Lth_band = np.zeros((len(BANDS), len(th)))
    for j in range(len(th)):
        if Teff[j] <= 0:
            continue
        kT = 1.380649e-16 * Teff[j]
        for b, (_, e1, e2, _) in enumerate(BANDS):
            x1 = e1 * 1.602176634e-12 / kT
            x2 = e2 * 1.602176634e-12 / kT
            Lth_band[b, j] = Lth[j] * planck_band_fraction(x1, x2)
    return Lth_band, Teff