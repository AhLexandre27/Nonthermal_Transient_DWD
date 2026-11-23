
import numpy as np

def B_ejecta(th, P, D):
    return D["B0"] * th**(-P["mu"] * P["n"])

def L_pulsar(t, D):
    return D["L0"] * (1.0 + t / D["tq"])**(-D["k"])

def N_of_E(E, th, P, D, nu=800):
    n, g, a = P["n"], P["gamma"], D["a"]
    ts, tq, k, c0 = P["t_star"], D["tq"], D["k"], D["c0"]
    E = np.atleast_1d(np.asarray(E, float))
    N = np.zeros_like(E)
    eps = E * th**n
    uup = 1.0 / eps
    Smax = (c0 / a) * (1.0 - th**(-a))
    ok = (E < P["Emax"]) & (Smax > 0)
    if not ok.any():
        return N
    uu, ep = uup[ok], eps[ok]
    lo0 = np.maximum(uu - Smax, 0.2 * th**(-n) / P["Emax"])

    def thi_of(u, uu_):
        S = np.maximum(uu_ - u, 0.0)
        return (th**(-a) + a * S / c0)**(-1.0 / a)

    def Ei_of(u, uu_):
        return thi_of(u, uu_)**(-n) / u

    def root(target, lo, hi, uu_):
        a_, b_ = np.log(lo), np.log(hi)
        for _ in range(70):
            mid = 0.5 * (a_ + b_)
            big = Ei_of(np.exp(mid), uu_) > target
            a_ = np.where(big, mid, a_)
            b_ = np.where(big, b_, mid)
        r = np.exp(0.5 * (a_ + b_))
        r = np.where(Ei_of(lo, uu_) <= target, lo, r)
        r = np.where(Ei_of(hi, uu_) >= target, hi, r)
        return r

    u_lo = root(P["Emax"], lo0, uu, uu)
    u_hi = root(P["Emin"], lo0, uu, uu)
    good = u_lo < u_hi
    val = np.zeros_like(uu)
    if good.any():
        ul, uh, u_, e_ = u_lo[good], u_hi[good], uu[good], ep[good]
        x = np.linspace(0.0, 1.0, nu)
        u = ul[:, None] * (uh / ul)[:, None]**x[None, :]
        thi = thi_of(u, u_[:, None])
        Ei = thi**(-n) / u
        Q0 = D["q0"] * (1.0 + ts * thi / tq)**(-k)
        f = (thi**(a + 1) / c0) * ts * Q0 * Ei**(-g) * thi**(-n) * (1.0 / (u * e_[:, None]))**2
        val[good] = np.trapz(f * u, np.log(u), axis=1)
    N[ok] = val * th**n
    return N


def band_E(nu1, nu2, B, D, Emax):
    E1 = np.sqrt(nu1 / (D["alpha"] * B))
    E2 = min(np.sqrt(nu2 / (D["alpha"] * B)), Emax)
    return E1, E2


def region_geometry(t, P, D):
    Rs = D["Rs0"] * t**P["n"]
    path = Rs * (1 - P["x0"])
    V = (4 * np.pi / 3) * Rs**3 * (1 - P["x0"]**3)
    return Rs, path, V