"""v0.11: Monte Carlo uncertainty budget (JCGM 101 style) for the radially resolved model.

All twelve inputs of budget.py are varied at once, with the uncertainties assumed there:
  - log-normal factor exp(u z) for inputs given as a relative standard uncertainty u
    (keeps conductivities, viscosities, cross-sections positive at u = 20-50 %);
  - vapour diffusivity: factor uniform in log between 1/2 and 2 (published simulations
    differ by a factor of two);
  - acid, acoustic frequency: uniform between 30 and 37.8 kHz (not reported).
Latin-hypercube sampling, 128 samples per case, one acoustic cycle per run (differs from
two cycles by 0.4 % in light output).

Usage:  python mc_budget.py run      (about an hour on two cores; writes mc.json)
        python mc_budget.py report
"""
import json
import math
import sys
from multiprocessing import Pool

import numpy as np
from scipy import stats
from scipy.stats import qmc

import budget

NS = 128
NAMES = list(budget.INPUTS)


def draws(case, seed):
    u01 = qmc.LatinHypercube(d=len(NAMES), seed=seed).random(NS)
    rows = []
    for row in u01:
        fac = {}
        for j, name in enumerate(NAMES):
            q = row[j]
            if name == "vapour diffusivity":
                fac[name] = math.exp((2 * q - 1) * math.log(2))
            elif case == "acid" and name == "frequency":
                fac[name] = (30e3 + q * 7.8e3) / 30e3
            else:
                fac[name] = math.exp(budget.U[case][name][0] * stats.norm.ppf(q))
        rows.append(fac)
    return rows


def one(task):
    case, i, fac = task
    import sbsl
    import profile as pr
    kw = dict(budget.CASES[case]["run"]); ana = dict(budget.CASES[case]["ana"])
    for name, f in fac.items():
        kind, key = budget.INPUTS[name]
        if kind == "run":
            kw[key] = kw.get(key, 1.0) * f
        elif key == "B_AR":
            pr.B_AR = pr.B_AR * f
        elif key == "SIG_SCALE":
            sbsl.SIG_SCALE = f
    o = pr.analyse(pr.run(cycles=1, **kw), **ana)
    o.pop("profile", None)
    return dict(case=case, i=i, fac=fac, **{k: float(v) for k, v in o.items()})


def tasks():
    t = []
    for case, seed in (("water", 11), ("acid", 12)):
        t.append((case, -1, {n: 1.0 for n in NAMES}))         # baseline, same settings
        t += [(case, i, f) for i, f in enumerate(draws(case, seed))]
    return t


def report(res):
    lin = {"water": 4.18, "acid": 7.07}           # linearised u(ln photons) from v0.8 (budget.py)
    rng = np.random.default_rng(0)
    for case in ("water", "acid"):
        rows = [r for r in res if r["case"] == case and r["i"] >= 0]
        base = next(r for r in res if r["case"] == case and r["i"] < 0)
        ph = np.array([r["photons"] for r in rows]); T = np.array([r["Tc"] for r in rows])
        lp = np.log(np.maximum(ph, 1e-3))
        q = np.percentile(ph, [2.5, 16, 50, 84, 97.5]); qT = np.percentile(T, [2.5, 50, 97.5])
        bs = np.array([np.percentile(rng.choice(ph, ph.size), [2.5, 97.5]) for _ in range(2000)])
        print(f"== {case}: {len(rows)} samples; baseline photons {base['photons']:.3g}, centre T {base['Tc']:.0f} K")
        print(f"   photons per flash: median {q[2]:.3g}; 68 % interval {q[1]:.3g} .. {q[3]:.3g}; 95 % interval {q[0]:.3g} .. {q[4]:.3g}")
        print(f"   95 % bounds, bootstrap spread: lower {np.percentile(bs[:,0],5):.3g}..{np.percentile(bs[:,0],95):.3g}, "
              f"upper {np.percentile(bs[:,1],5):.3g}..{np.percentile(bs[:,1],95):.3g}")
        print(f"   width of the 95 % interval: factor {q[4]/max(q[0],1e-3):.3g}; linearised GUM gives factor {math.exp(2*1.96*lin[case]):.3g}")
        print(f"   68 % interval as a factor either side of the median: down {q[2]/max(q[1],1e-3):.3g}, up {q[3]/q[2]:.3g}; "
              f"linearised: {math.exp(lin[case]):.3g} each way")
        print(f"   mean photons {ph.mean():.3g} = {ph.mean()/q[2]:.3g} x median; baseline is at the {stats.percentileofscore(ph, base['photons']):.0f}th percentile")
        print(f"   centre temperature: median {qT[1]:.0f} K, 95 % interval {qT[0]:.0f} .. {qT[2]:.0f} K")
        print(f"   samples with fewer than 1 % of baseline light: {np.mean(ph < 0.01*base['photons'])*100:.0f} %")
        # which inputs drive the spread: rank correlation of ln(photons) with ln(factor)
        cors = []
        for n in NAMES:
            lf = np.log([r["fac"][n] for r in rows])
            if np.ptp(lf) > 0:
                cors.append((n, stats.spearmanr(lf, lp)[0]))
        cors.sort(key=lambda z: -abs(z[1]))
        print("   rank correlation with ln(photons): " + ", ".join(f"{n} {c:+.2f}" for n, c in cors[:5]))
        print()


if __name__ == "__main__":
    if sys.argv[1:] == ["run"]:
        import os
        res = json.load(open("mc.json")) if os.path.exists("mc.json") else []
        done = {(r["case"], r["i"]) for r in res}
        todo = [t for t in tasks() if (t[0], t[1]) not in done]
        with Pool(2, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(one, todo):
                res.append(r); json.dump(res, open("mc.json", "w"))
                print(len(res), r["case"], r["i"], f"{r['photons']:.3g}", flush=True)
    else:
        report(json.load(open("mc.json")))
