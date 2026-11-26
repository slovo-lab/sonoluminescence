"""Conservation checks for the radially resolved model (profile.py).

1. Argon in the whole bubble: exact by construction (pressure follows from it).
2. Vapour in the whole bubble: changes only by evaporation and condensation at the wall.
3. Local check, not built in: for the inner part of the bubble (cells 0..k) the change of
   argon content must equal the argon flux through its outer face.
4. Pressure: the value fixed by argon conservation must agree with the rate dp/dt that
   the wall equation uses.
Usage: python conservation_check.py [N] [Pa_atm]
"""
import sys
import numpy as np
import profile as pr
from sbsl import kB
from vapor import B_AR


def check(N=40, Pa=1.4, window=3e-6, npts=3001):
    r = pr.run(Pa_atm=Pa, N=N)
    vy = r["vol_y"]
    tm, _ = pr.collapse_time(r)
    # the collapse lasts about 100 ps: sample it at 0.5 ps, the rest of the window more coarsely
    t = np.unique(np.concatenate([np.linspace(tm - window, tm + window, npts),
                                  np.linspace(tm - 3e-9, tm + 3e-9, 12001)]))
    R, U, p, T, x, n = pr.unpack(r, t)
    V = vy[None, :] * R[:, None] ** 3
    Nar_cell = n * (1 - x) * V
    Nv_cell = n * x * V
    out = {"argon_total": float(Nar_cell.sum(1).max() / Nar_cell.sum(1).min() - 1)}
    st = r["sol"].sol(t)
    D = [r["rhs"](ti, st[:, i], diag=True) for i, ti in enumerate(t)]
    F = np.array([d["F"] for d in D]); G = np.array([d["G"] for d in D]); dp = np.array([d["dp"] for d in D])
    # 3. local argon balance for inner regions, integrated over the window
    worst = 0.0
    for k in (N // 4, N // 2, 3 * N // 4, N - 2):
        inner = Nar_cell[:, :k + 1].sum(1)
        flux_out = F[:, k] - G[:, k]                                     # argon leaving through face k
        integ = np.concatenate([[0], np.cumsum(0.5 * (flux_out[1:] + flux_out[:-1]) * np.diff(t))])
        err = np.abs((inner - inner[0]) + integ).max() / inner.max()
        worst = max(worst, err)
    out["argon_local"] = worst
    # 2. vapour: interior content changes only through the last interior face
    vin = Nv_cell[:, :-1].sum(1)
    fl = -G[:, -1]
    integ = np.concatenate([[0], np.cumsum(0.5 * (fl[1:] + fl[:-1]) * np.diff(t))])
    out["vapour_wall"] = float(np.abs((vin - vin[0]) - integ).max() / np.abs(vin - vin[0]).max())
    # 4. pressure consistency
    dpn = np.gradient(p, t)
    sel = slice(10, -10)
    out["pressure_rate"] = float(np.abs(dpn[sel] - dp[sel]).max() / np.abs(dp[sel]).max())
    return out


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    Pa = float(sys.argv[2]) if len(sys.argv) > 2 else 1.4
    o = check(N, Pa)
    print(f"N={N} Pa={Pa}: argon total {o['argon_total']:.1e} | argon, inner regions {o['argon_local']:.2%} | "
          f"vapour vs wall flux {o['vapour_wall']:.2%} | dp/dt mismatch {o['pressure_rate']:.2%}")
