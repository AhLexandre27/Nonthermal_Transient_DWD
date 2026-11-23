import sys
import numpy as np
from params import BANDS, NTH, NBAND, NU
from physics import N_of_E, B_ejecta, band_E
from ssa import ssa_factor

def lightcurves(P, D, nth=NTH, nband=NBAND, nu=NU, verbose=True):
    th = np.geomspace(1.0, P["theta_max"], nth)
    Eb0 = D["a"] / (D["beta"] * D["B0"]**2 * P["t_star"])
    Ebol = np.geomspace(0.3 * min(P["Emin"], Eb0) * P["theta_max"]**(-P["n"]), P["Emax"], 600)
    lnE = np.log(Ebol)

    Lint = np.zeros((len(BANDS), nth))
    Lssa = np.zeros((len(BANDS), nth))
    TauSSA = np.zeros((len(BANDS), nth))
    Lnth = np.zeros(nth)

    for j, t in enumerate(th):
        B = B_ejecta(t, P, D)
        Nb = N_of_E(Ebol, t, P, D, nu=nu)
        Lnth[j] = np.trapz(D["beta"] * B**2 * Ebol**3 * Nb, lnE)
        grids = []
        for (_, e1, e2, _) in BANDS:
            nu1, nu2 = e1 * 1.602176634e-12 / 6.62607015e-27, e2 * 1.602176634e-12 / 6.62607015e-27
            E1, E2 = band_E(nu1, nu2, B, D, P["Emax"])
            grids.append(np.geomspace(E1, E2, nband) if E1 < E2 else np.array([]))
        Eall = np.concatenate(grids) if grids else np.array([])
        Nall = N_of_E(Eall, t, P, D, nu=nu) if Eall.size else Eall
        i0 = 0
        for b, g_ in enumerate(grids):
            if len(g_) == 0:
                i0 += len(g_); continue
            Nj = Nall[i0:i0 + len(g_)]; i0 += len(g_)
            integrand = D["beta"] * B**2 * g_**3 * Nj
            Lint[b, j] = np.trapz(integrand, np.log(g_))
            f, tau_ssa = ssa_factor(g_, t, B, P, D)
            Lssa[b, j] = np.trapz(integrand * f, np.log(g_))
            TauSSA[b, j] = tau_ssa[len(tau_ssa) // 2]
        if verbose and (j % max(nth // 10, 1) == 0):
            print("  theta=%9.3e  (%d/%d)" % (t, j + 1, nth), file=sys.stderr)
    return th, Lint, Lssa, TauSSA, Lnth


def check_energy(th, Lnth, P, D):
    n, ts = P["n"], P["t_star"]
    t = th * ts
    Linj = P["eta_nth"] * D["L0"] * (1 + t / D["tq"])**(-D["k"])
    inj = np.trapz(Linj, t)
    rad = np.trapz(Lnth, t)
    print("\nEnergy balance: integrated L_syn / integrated L_injected = %.4f" % (rad / inj))


def check_ssa(th, TauSSA, P):
    print("\nSynchrotron self-absorption (SSA): peak tau at band center, over time")
    for b, (name, *_rest) in enumerate(BANDS):
        jmax = TauSSA[b].argmax()
        tmax = th[jmax] * P["t_star"] / 86400.0
        peak = TauSSA[b, jmax]
        flag = "  <-- SSA relevant here" if peak > 0.1 else ""
        print("  %-26s tau_ssa,max = %.3e (at t=%.2f d)%s" % (name, peak, tmax, flag))