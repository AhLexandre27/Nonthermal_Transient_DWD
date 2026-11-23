# DWD Pulsar Transient Model

Light curves of the ejecta of a pulsar-powered double white dwarf (DWD) merger,
including non-thermal (synchrotron + self-absorption) and thermal emission.

---

## Physics (brief)

A central neutron star (pulsar) injects spin-down luminosity
`L(t) = L0 / (1 + t/tq)^2` into the expanding ejecta. The ejecta expands
self-similarly (`R* ∝ t^n`), diluting the magnetic field (`B* ∝ t^(-mu*n)`).
Relativistic electrons are injected with a power-law spectrum `E^-gamma` and
lose energy through adiabatic expansion and synchrotron radiation only (no
inverse Compton). The exact electron distribution `N(E,t)` is computed by the
method of characteristics. Synchrotron self-absorption (SSA) is applied to the
non-thermal light. The fraction of synchrotron light absorbed by the ejecta is
reprocessed into a thermal reservoir, which re-emits as a blackbody at the
effective temperature `T_eff(t)`.

---

## Files

### Core modules

| File | What it does |
|---|---|
| `constants.py` | CGS physical constants (`c`, `m_e`, `h`, `eV`, `DAY`, `kB`, `sigma_SB`, ...). |
| `params.py` | Free parameters `P` (ejecta, neutron star, injection), observational bands `BANDS`, and numerical resolution `NTH`, `NBAND`, `NU`. **Edit here.** |
| `derived.py` | `derived(P)` computes all derived quantities: `L0`, `tq`, `B0`, `tau1`, `q0`, `alpha`, `beta`, etc. `print_params` prints them. |
| `physics.py` | Pure physics: `B_ejecta` (field at the ejecta), `L_pulsar`, `N_of_E` (exact electron distribution), `band_E`, `region_geometry`. |
| `ssa.py` | `ssa_factor`: synchrotron self-absorption coefficient from the exact `N(E,t)`. |
| `thermal.py` | Thermal reservoir ODE (`thermal_and_check`), effective temperature, blackbody band split. |
| `lightcurves.py` | Main loop `lightcurves()`: computes non-thermal and SSA-corrected light curves per band. Also `check_energy` and `check_ssa` diagnostics. |
| `plotting.py` | `make_plots`: produces the three figures (non-thermal, thermal, total) per band. |
| `selftest.py` | Compares `N_of_E` against an independent reference solution. |

### Notebook

| File | What it does |
|---|---|
| `main.ipynb` | Runs the full pipeline: `derived` → `lightcurves` → `thermal` → `check` → `make_plots`. Produces the figures and prints peak luminosities. |

---

## How it connects


1. `derived(P)` turns the free parameters into physical quantities (`L0`, `tq`, `B0`, ...).
2. `lightcurves(P, D)` loops over time, computing `N(E,t)` and integrating the synchrotron emission per band, with SSA.
3. `thermal_and_check(th, Lnth, P, D)` solves the thermal reservoir ODE sourced by the absorbed synchrotron light.
4. `thermal_band_luminosities` splits the bolometric thermal output into bands via a blackbody at `T_eff(t)`.
5. `make_plots` produces the three figures.

---

## How to run

Open `main.ipynb` and run all cells. The three figures are saved as
`lightcurve_nonthermal.png`, `lightcurve_thermal.png`, `lightcurve_total.png`.
