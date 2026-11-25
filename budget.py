"""v0.8: uncertainty budget of the radially resolved model (GUM-style).

For each input x_i the script computes the sensitivity coefficient
    c_i = d ln(y) / d ln(x_i)
by central differences over +-5 %, for the outputs y = photons per flash, centre
temperature, collapse radius, flash width.  With an assumed relative standard uncertainty
u_i the contribution to u(ln y) is |c_i| u_i; contributions are combined in quadrature.

Usage:  python budget.py run     (about 15 minutes on two cores; writes budget.json)
        python budget.py report  (prints the tables from budget.json)
"""
import json
import sys
from multiprocessing import Pool

import numpy as np

STEP = 0.05
BAR = 1 / 1.01325
CASES = {
    "water": dict(run=dict(R0=4.5e-6, Pa_atm=1.35, f=26.5e3, pv=2337.0, Tw=293.15, mu_w=1.0e-3, sigma_w=0.0725,
                           rho_l=998.0, c_l=1481.0), ana=dict()),
    "acid": dict(run=dict(R0=13e-6, Pa_atm=1.4 * BAR, f=30e3, pv=5.3, Tw=298.0, mu_w=25e-3, sigma_w=0.055,
                          rho_l=1800.0, c_l=1470.0), ana=dict(half=25e-9, npts=1001)),
}
# input -> (kind, key).  kind: "run" = argument of profile.run, "glob" = module constant
INPUTS = {
    "drive pressure": ("run", "Pa_atm"),
    "ambient radius": ("run", "R0"),
    "frequency": ("run", "f"),
    "vapour pressure": ("run", "pv"),
    "viscosity": ("run", "mu_w"),
    "surface tension": ("run", "sigma_w"),
    "liquid density": ("run", "rho_l"),
    "sound speed": ("run", "c_l"),
    "vapour diffusivity": ("run", "d_scale"),
    "gas thermal conductivity": ("run", "k_scale"),
    "co-volume": ("glob", "B_AR"),
    "e-neutral cross-section": ("glob", "SIG_SCALE"),
}
# assumed relative standard uncertainties (1 sigma), per case, with the reason
U = {
    "water": {
        "drive pressure": (0.10, "fits tolerate about 10 % (Hopkins et al. 2005)"),
        "ambient radius": (0.05, "assumed, light-scattering fit"),
        "frequency": (0.001, "assumed, electronic"),
        "vapour pressure": (0.06, "liquid temperature known to 1 K"),
        "viscosity": (0.025, "liquid temperature known to 1 K"),
        "surface tension": (0.02, "assumed"),
        "liquid density": (0.001, "assumed"),
        "sound speed": (0.005, "assumed"),
        "vapour diffusivity": (0.40, "published simulations differ by a factor 2; rectangular, ln 2 / sqrt 3"),
        "gas thermal conductivity": (0.20, "assumed, dense-gas corrections"),
        "co-volume": (0.10, "assumed"),
        "e-neutral cross-section": (0.30, "approximate table"),
    },
}
U["acid"] = dict(U["water"])
U["acid"].update({
    "frequency": (0.066, "not reported; 30 to 37.8 kHz, rectangular"),
    "viscosity": (0.20, "strong dependence on concentration and temperature"),
    "vapour pressure": (0.50, "assumed"),
})
OUTS = ["photons", "Tc", "Rmin", "fwhm_uv", "p_max"]


def one(task):
    case, name, factor, extra = task
    import sbsl
    import profile as pr
    kw = dict(CASES[case]["run"]); ana = dict(CASES[case]["ana"]); ana.update(extra)
    if name is not None:
        kind, key = INPUTS[name]
        if kind == "run":
            kw[key] = kw.get(key, 1.0) * factor
        elif key == "B_AR":
            pr.B_AR = pr.B_AR * factor
        elif key == "SIG_SCALE":
            sbsl.SIG_SCALE = factor
    o = pr.analyse(pr.run(**kw), **ana)
    o.pop("profile", None)
    return dict(case=case, name=name, factor=factor, extra=extra, **{k: float(v) for k, v in o.items()})


def tasks():
    t = []
    for case in CASES:
        t.append((case, None, 1.0, {}))
        t.append((case, None, 1.0, {"lowering": False}))
        t.append((case, None, 1.0, {"chemistry": False}))
        for name in INPUTS:
            for fac in (1 - STEP, 1 + STEP):
                t.append((case, name, fac, {}))
        for fac in (0.90, 0.98, 0.99, 1.01, 1.02, 1.10):      # response to drive pressure in detail
            t.append((case, "drive pressure", fac, {}))
    return t


def report(res):
    for case in CASES:
        rows = [r for r in res if r["case"] == case]
        base = next(r for r in rows if r["name"] is None and not r["extra"])
        print(f"== {case}: baseline photons {base['photons']:.3g}, centre T {base['Tc']:.0f} K, Rmin {base['Rmin']*1e6:.3f} um, "
              f"flash {base['fwhm_uv']*1e9:.3f} ns, vapour {base['vapour_frac']*100:.1f} %")
        print(f"{'input':26s} {'c(photons)':>10s} {'c(T)':>7s} {'c(Rmin)':>8s} {'c(width)':>8s} {'u(input)':>9s} {'-> u(ln photons)':>17s} {'-> u(T)/T':>10s}")
        tot_p = tot_t = 0.0; contrib = []
        for name in INPUTS:
            lo = next(r for r in rows if r["name"] == name and abs(r["factor"] - (1 - STEP)) < 1e-9)
            hi = next(r for r in rows if r["name"] == name and abs(r["factor"] - (1 + STEP)) < 1e-9)
            d = np.log(1 + STEP) - np.log(1 - STEP)
            c = {k: (np.log(hi[k]) - np.log(lo[k])) / d for k in OUTS}
            u = U[case][name][0]
            up, ut = abs(c["photons"]) * u, abs(c["Tc"]) * u
            tot_p += up**2; tot_t += ut**2; contrib.append((name, up))
            print(f"{name:26s} {c['photons']:10.1f} {c['Tc']:7.2f} {c['Rmin']:8.2f} {c['fwhm_uv']:8.2f} {u*100:8.1f}% {up:17.2f} {ut*100:9.1f}%")
        tp = np.sqrt(tot_p)
        print(f"{'combined':26s} {'':10s} {'':7s} {'':8s} {'':8s} {'':9s} {tp:17.2f} {np.sqrt(tot_t)*100:9.1f}%")
        print(f"   -> photons uncertain by a factor exp({tp:.2f}) = {np.exp(tp):.3g} (1 sigma, linearised)")
        contrib.sort(key=lambda z: -z[1])
        print("   share of variance: " + ", ".join(f"{n} {100*v**2/tot_p:.0f} %" for n, v in contrib[:5]))
        cPa = next(v for n, v in contrib if n == "drive pressure") / U[case]["drive pressure"][0]
        rest = np.sqrt(max(tot_p - (cPa * U[case]["drive pressure"][0]) ** 2, 0))
        print(f"   without drive pressure the rest gives a factor {np.exp(rest):.3g}; to keep the pressure term below ln 2 "
              f"the drive must be known to {100*np.log(2)/cPa:.1f} %")
        print("   response to drive pressure: " + ", ".join(
            f"{(r['factor']-1)*100:+.0f} %: x{r['photons']/base['photons']:.3g}" for r in sorted(
                (r for r in rows if r["name"] == "drive pressure"), key=lambda r: r["factor"])))
        for tag, ex in [("ionisation lowering off", {"lowering": False}), ("chemistry off", {"chemistry": False})]:
            r = next(r for r in rows if r["name"] is None and r["extra"] == ex)
            print(f"   model switch, {tag}: photons x{r['photons']/base['photons']:.3g}, centre T {r['Tc']:.0f} K")
        print()


if __name__ == "__main__":
    if sys.argv[1:] == ["run"]:
        with Pool(2, maxtasksperchild=1) as pool:
            res = []
            for r in pool.imap_unordered(one, tasks()):
                res.append(r); json.dump(res, open("budget.json", "w"))
                print(len(res), r["case"], r["name"], r["factor"], r["extra"], f"{r['photons']:.3g}", flush=True)
    else:
        report(json.load(open("budget.json")))
