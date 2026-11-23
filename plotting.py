# plotting.py
import numpy as np
import matplotlib.pyplot as plt
from params import BANDS
from constants import DAY

def make_plots(th, P, D, Lssa, Lth_band, Teff, tau):
    t = th * P["t_star"]
    td = t / DAY
    colors = ["#577590", "#2a9d8f", "#e76f51", "#6a4c93", "#999999"]
    figs = []

    Lnth_out = np.zeros_like(Lssa)
    for b, (name, e1, e2, kband) in enumerate(BANDS):
        kf = (kband if kband is not None else P["kappa"]) / P["kappa"]
        Lnth_out[b] = np.exp(-tau * kf) * Lssa[b]

    def new_fig(title):
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("t [days]"); ax.set_ylabel("L [erg/s]")
        ax.set_title(title)
        ax.grid(alpha=0.25, which="both")
        return fig, ax

    fig1, ax1 = new_fig("Non-thermal (synchrotron) luminosity per band")
    ymax1 = max(Lnth_out.max(), 1.0)
    for b, (name, *_r) in enumerate(BANDS):
        if Lnth_out[b].max() > 0:
            ax1.plot(td, Lnth_out[b], color=colors[b], lw=2, label=name)
    ax1.set_ylim(ymax1 * 1e-9, ymax1 * 5)
    ax1.legend(fontsize=8)
    figs.append(fig1)

    fig2, ax2 = new_fig("Thermal luminosity per band (blackbody split)")
    ymax2 = max(Lth_band.max(), 1.0)
    for b, (name, *_r) in enumerate(BANDS):
        if Lth_band[b].max() > 0:
            ax2.plot(td, Lth_band[b], color=colors[b], lw=2, label=name)
    ax2.set_ylim(ymax2 * 1e-9, ymax2 * 5)
    ax2.legend(fontsize=8)
    figs.append(fig2)

    Ltot = Lnth_out + Lth_band
    fig3, ax3 = new_fig("Total (non-thermal + thermal) luminosity per band")
    ymax3 = max(Ltot.max(), 1.0)
    for b, (name, *_r) in enumerate(BANDS):
        if Ltot[b].max() > 0:
            ax3.plot(td, Ltot[b], color=colors[b], lw=2, label=name)
    ax3.set_ylim(ymax3 * 1e-9, ymax3 * 5)
    ax3.legend(fontsize=8)
    figs.append(fig3)

    names = ["lightcurve_nonthermal.png", "lightcurve_thermal.png", "lightcurve_total.png"]
    for fig, name in zip(figs, names):
        fig.tight_layout()
        fig.savefig(name, dpi=150)
    plt.show()
    return figs, Lnth_out, Ltot