# selftest.py
import numpy as np
from params import P
from physics import N_of_E

def selftest():
    Pt = dict(P); Pt.update(n=1.0, mu=1.0, gamma=1.5, Emin=1e-12, Emax=1e8, t_star=1.0)
    Dt = {}
    Dt["beta"] = 1.0; Dt["alpha"] = 1.0; Dt["a"] = Pt["n"] * (1 + 2 * Pt["mu"]) - 1
    Dt["B0"] = 1.0; Dt["c0"] = 1.0; Dt["q0"] = 1.0; Dt["k"] = 0.0; Dt["tq"] = 1e300
    th_ref = 200.0
    a_ref = Dt["a"]
    Eb = a_ref * th_ref**(2 * Pt["mu"] * Pt["n"] - 1)
    exact = {20: 1.3087e-5, 100: 2.4219e-7, 1000: 7.3972e-10}
    print("E/Eb      N_numeric      N_reference     ratio")
    for f, ex in exact.items():
        Nn = N_of_E(np.array([f * Eb]), th_ref, Pt, Dt, nu=3000)[0]
        print("%-6d  %.5e   %.5e   %.4f" % (f, Nn, ex, Nn / ex))