"""v0.6: radially resolved bubble interior (uniform pressure, non-uniform temperature
and composition).

Dynamics : Keller-Miksis wall equation.  Inside the bubble the pressure is uniform
           (low-Mach approximation); temperature T(r,t) and vapour mole fraction x(r,t)
           obey the energy and diffusion equations with heat conduction, vapour
           diffusion, interdiffusion enthalpy flux and Stefan flow at the wall.  Gas law
           p (1/n - b) = kT with one co-volume b for both species.
Chemistry: post-hoc in every shell at fixed volume and energy (as in v0.4).
Light    : formal solution of radiative transfer along straight rays through the
           shells (Kramers-Unsoeld continuum of Ar, H, O + electron-neutral
           bremsstrahlung, as in v0.3/v0.4).
"""
import numpy as np
from scipy.integrate import solve_ivp
import sbsl
import vapor
from sbsl import kB, h_pl, c0, qe, e2, rho, cw, P0, planck, kappa_en_per_pair
from vapor import (B_AR, D0_DIFF, N_REF, KAPPA_G, THETA_V, D_ATOM, E_I, e_w, dissociation,
                   ion_state)

T_REF = 300.0


def c_w_vec(T):
    th = THETA_V / np.asarray(T, float)[..., None]
    return kB * (3 + np.sum(th**2 * np.exp(th) / np.expm1(th) ** 2, axis=-1))


def grid(N=64, stretch=3.0):
    """Nodes in y = r/R, clustered towards the wall (y = 1)."""
    s = np.linspace(0, 1, N + 1)
    return np.tanh(stretch * s) / np.tanh(stretch)


def run(R0=4.5e-6, Pa_atm=1.3, f=26.5e3, cycles=2, pv=2337.0, Tw=293.15, mu_w=1.0e-3,
        sigma_w=0.0725, N=40, rtol=1e-5, d_scale=1.0, k_scale=1.0, rho_l=rho, c_l=cw):
    Pa = Pa_atm * 101325.0
    w = 2 * np.pi * f
    y = grid(N)
    yf = 0.5 * (y[1:] + y[:-1])                 # faces i+1/2, i = 0..N-1
    dy = np.diff(y)
    fl = np.concatenate([[0.0], yf])            # lower face of cell i (i = 0..N)
    fu = np.concatenate([yf, [1.0]])            # upper face of cell i
    vol_y = (fu**3 - fl**3) / 3                 # cell volumes / (4 pi R^3)
    p0 = P0 + 2 * sigma_w / R0
    cpw_minus = None
    # 3-point one-sided derivative weights at the wall
    h1, h2 = y[-1] - y[-2], y[-1] - y[-3]
    wN = (h1 + h2) / (h1 * h2); wN1 = -h2 / (h1 * (h2 - h1)); wN2 = h1 / (h2 * (h2 - h1))

    def fields(T, x, p):
        cv = (1 - x) * 1.5 * kB + x * c_w_vec(T)
        cp = cv + kB
        kTb = kB * T / (kB * T + B_AR * p)
        n = p / (kB * T + B_AR * p)
        kap = k_scale * KAPPA_G * (T / T_REF) ** 0.7
        nD = d_scale * D0_DIFF * N_REF * (T / Tw) ** 0.7
        return cv, cp, kTb, n, kap, nD

    def rhs(t, s):
        T = np.empty(N + 1); x = np.empty(N + 1)
        T[:N] = s[:N]; x[:N] = s[N:2 * N]
        p, R, U = s[2 * N], s[2 * N + 1], s[2 * N + 2]
        T[N] = Tw
        x[N] = min(pv / p, 0.999)
        x = np.clip(x, 0.0, 0.9999)
        cv, cp, kTb, n, kap, nD = fields(T, x, p)
        # face fluxes (per 4 pi), r^2 * flux
        gT = np.diff(T) / (R * dy); gx = np.diff(x) / (R * dy)
        rf2 = (R * yf) ** 2
        FT = rf2 * 0.5 * (kap[1:] + kap[:-1]) * gT
        Fx = rf2 * 0.5 * (nD[1:] + nD[:-1]) * gx
        gTw = (wN * T[N] + wN1 * T[N - 1] + wN2 * T[N - 2]) / R
        gxw = (wN * x[N] + wN1 * x[N - 1] + wN2 * x[N - 2]) / R
        FTw = R**2 * kap[N] * gTw; Fxw = R**2 * nD[N] * gxw
        FT_lo = np.concatenate([[0.0], FT]); FT_up = np.concatenate([FT, [FTw]])
        Fx_lo = np.concatenate([[0.0], Fx]); Fx_up = np.concatenate([Fx, [Fxw]])
        vol = vol_y * R**3
        # node gradients for the interdiffusion term
        dTdy = np.gradient(T, y); dxdy = np.gradient(x, y)
        dcp = c_w_vec(T) + kB - 2.5 * kB
        inter = nD * dcp * dxdy * dTdy / R**2
        Qp = (FT_up - FT_lo) / vol + inter
        # pressure from the wall condition
        A = np.sum(kB / (cp * p) * Qp * vol)
        Bc = np.sum(cv / cp * kTb * vol) / p
        s_wall = -(nD[N] / n[N]) * gxw / (1 - x[N])
        dp = (A - R**2 * (U + s_wall)) / Bc
        # velocity at faces and nodes
        r2u_f = np.cumsum(kB / (cp * p) * Qp * vol - dp / p * cv / cp * kTb * vol)   # upper faces
        r2u_lo = np.concatenate([[0.0], r2u_f[:-1]])
        u = np.zeros(N + 1)
        u[1:] = 0.5 * (r2u_lo[1:] + r2u_f[1:]) / (R * y[1:]) ** 2
        DT = (kB * T * dp + (kB * T + B_AR * p) * Qp) / (p * cp)
        Dx = (Fx_up - Fx_lo) / (vol * n)
        conv = (u - y * U) / R
        dTdt = DT - conv * dTdy
        dxdt = Dx - conv * dxdy
        Pac = -Pa * np.sin(w * t)
        M = U / c_l
        num = ((1 + M) * (p - P0 - Pac) / rho_l + R * dp / (rho_l * c_l)
               - 4 * mu_w * U / (rho_l * R) - 2 * sigma_w / (rho_l * R) - 1.5 * U**2 * (1 - M / 3))
        out = np.empty_like(s)
        out[:N] = dTdt[:N]; out[N:2 * N] = dxdt[:N]
        out[2 * N] = dp; out[2 * N + 1] = U; out[2 * N + 2] = num / ((1 - M) * R)
        return out

    s0 = np.concatenate([np.full(N, Tw), np.full(N, pv / p0), [p0, R0, 0.0]])
    atol = np.concatenate([np.full(N, 1e-3), np.full(N, 1e-9), [1e-3, 1e-13, 1e-8]])
    sol = solve_ivp(rhs, [0, cycles / f], s0, method="LSODA", rtol=rtol, atol=atol,
                    max_step=1 / f / 400, dense_output=True)
    return dict(sol=sol, y=y, vol_y=vol_y, N=N, f=f, R0=R0, Tw=Tw, pv=pv, Pa_atm=Pa_atm, c_l=c_l)


def unpack(res, t):
    """Fields at times t (array).  Returns R, U, p, T[t, node], x[t, node], n[t, node]."""
    N = res["N"]
    s = res["sol"].sol(np.atleast_1d(t))
    p, R, U = s[2 * N], s[2 * N + 1], s[2 * N + 2]
    T = np.vstack([s[:N], np.full((1, s.shape[1]), res["Tw"])]).T
    x = np.vstack([s[N:2 * N], np.minimum(res["pv"] / p, 0.999)[None, :]]).T
    x = np.clip(x, 0.0, 0.9999)
    n = p[:, None] / (kB * T + B_AR * p[:, None])
    return R, U, p, T, x, n


def collapse_time(res):
    f = res["f"]; t_end = res["sol"].t[-1]
    tc = np.linspace(t_end - 1 / f, t_end, 40001)
    R = res["sol"].sol(tc)[2 * res["N"] + 1]
    i0 = R.argmax(); seg = R[i0:]
    j = np.where((seg[1:-1] < seg[:-2]) & (seg[1:-1] <= seg[2:]))[0][0] + 1
    return tc[i0 + j], R.max()


def analyse(res, lam=np.linspace(200e-9, 800e-9, 41), half=0.6e-9, npts=401, chemistry=True,
            nray=32, emission=True):
    N = res["N"]; vol_y = res["vol_y"]; y = res["y"]
    tm, Rmax = collapse_time(res)
    tf = np.linspace(tm - half, tm + half, npts)
    R, U, p, T, x, n = unpack(res, tf)
    Npart = n * vol_y[None, :] * R[:, None] ** 3          # particles per cell / (4 pi)
    im = R.argmin()
    out = dict(Rmax=Rmax, Rmin=R.min(), p_max=float(p.max()), n_centre=float(n[im, 0]), n_mean=float(Npart[im].sum() / (R[im] ** 3 / 3)),
               vapour_frac=float((Npart * x)[im].sum() / Npart[im].sum()),
               Tc_dyn=T[:, 0].max(),
               Tavg_dyn=float(((Npart * T).sum(1) / Npart.sum(1)).max()))
    # conservation check over the whole run: argon number
    tt = np.linspace(res["sol"].t[-1] - 1 / res["f"], res["sol"].t[-1], 2001)
    R2, _, _, T2, x2, n2 = unpack(res, tt)
    Nar = (n2 * (1 - x2) * vol_y[None, :]).sum(1) * R2**3
    out["argon_drift"] = float(Nar.max() / Nar.min() - 1)
    out["Umax"] = float(abs(res["sol"].sol(np.linspace(tm - 20e-9, tm + 1e-9, 4001))[2 * N + 2]).max())
    out["Mach"] = out["Umax"] / res.get("c_l", cw)
    if not emission:
        return out
    # ---- post-hoc chemistry and ionisation in hot shells (fixed volume, fixed energy)
    hot = T > 3000.0
    Th, xh, nh = T[hot], x[hot], n[hot]
    budget = (1 - xh) * 1.5 * kB * Th + xh * e_w(Th)
    free = 1 - nh * B_AR

    def state(Tq):
        xd = dissociation(Tq, nh * xh / free) if chemistry else np.zeros_like(Tq)
        n_ar, n_h, n_o, n_m = nh * (1 - xh), 2 * xd * nh * xh, xd * nh * xh, (1 - xd) * nh * xh
        ne, dE, ions = ion_state(Tq, n_ar, n_h, n_o, n_ar + n_h + n_o + n_m)
        E = (1.5 * kB * Tq * ((1 - xh) + ne / nh) + xh * ((1 - xd) * e_w(Tq) + xd * (4.5 * kB * Tq + D_ATOM))
             + sum(ions[s] * (E_I[s] - dE) for s in ions) / nh)
        return E, xd, ne, dE, ions, dict(Ar=n_ar, H=n_h, O=n_o), n_m

    lo = np.full_like(Th, 250.0); hi = Th + 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        over = state(mid)[0] > budget
        hi = np.where(over, mid, hi); lo = np.where(over, lo, mid)
    Tq = 0.5 * (lo + hi)
    _, xd, ne, dE, ions, dens, n_m = state(Tq)
    Tfin = T.copy(); Tfin[hot] = Tq
    out["Tc"] = Tfin[:, 0].max()
    # ---- absorption coefficient per shell and wavelength
    nu = c0 / lam
    kap = np.zeros(T.shape + (len(lam),))
    xx = h_pl * nu[None, :] / (kB * Tq[:, None])
    pref = (16 * np.pi**2 / (3 * np.sqrt(3))) * e2**3 * kB * Tq[:, None] / (h_pl**4 * c0 * nu[None, :] ** 3)
    k_ion = sum(pref * (dens[s] - ions[s])[:, None] * np.exp(-(E_I[s] - dE)[:, None] / (kB * Tq[:, None]) + xx)
                for s in dens) * (1 - np.exp(-xx))
    Tg = np.geomspace(2900.0, Tq.max() * 1.001, 90)
    tab = np.log(np.array([kappa_en_per_pair(lam, Tt) for Tt in Tg]) + 1e-300)
    kp = np.exp(np.array([np.interp(np.log(Tq), np.log(Tg), tab[:, j]) for j in range(len(lam))]).T)
    n_neut = sum(dens[s] - ions[s] for s in dens) + n_m
    kap[hot] = k_ion + kp * (ne * n_neut)[:, None]
    B = np.zeros_like(kap)
    B[hot] = planck(lam[None, :], Tq[:, None])
    # ---- ray tracing: shells bounded by faces; node N (wall half-cell) is cold, dropped
    yf = 0.5 * (y[1:] + y[:-1])
    edges = np.concatenate([[0.0], yf[:-1]])            # inner edges of shells 0..N-1, outer = yf
    outer = yf                                          # outer radius of shell i (i = 0..N-1)
    pb = (np.arange(nray) + 0.5) / nray * outer[-1]     # impact parameters / R
    wts = 2 * np.pi * pb * (outer[-1] / nray)
    half_len = np.sqrt(np.maximum(outer[None, :] ** 2 - pb[:, None] ** 2, 0)) - \
        np.sqrt(np.maximum(edges[None, :] ** 2 - pb[:, None] ** 2, 0))          # (ray, shell)
    Pl = np.zeros((npts, len(lam)))
    for it in range(npts):
        if not hot[it].any():
            continue
        k_s = kap[it, :N, :]; B_s = B[it, :N, :]                       # (shell, lam)
        dtau = half_len[:, :, None] * R[it] * k_s[None, :, :]          # (ray, shell, lam), one half
        # far half: light generated in shell i passes shells < i (far side, inward), then the
        # whole near half.  Order of traversal: far outer -> centre -> near outer.
        tau_in = np.cumsum(dtau, axis=1) - dtau                        # depth from centre side to shell i (exclusive)
        tau_half = dtau.sum(axis=1, keepdims=True)
        tau_out = tau_half - tau_in - dtau                             # depth from shell i outwards (exclusive)
        src = B_s[None] * (-np.expm1(-dtau))
        near = (src * np.exp(-tau_out)).sum(axis=1)
        far = (src * np.exp(-(tau_in + tau_half))).sum(axis=1)
        I = near + far                                                  # (ray, lam)
        Pl[it] = 4 * np.pi * R[it] ** 2 * (wts[:, None] * I).sum(axis=0)
    rate = np.trapezoid(Pl * lam[None, :] / (h_pl * c0), lam, axis=1)
    Ptot = np.trapezoid(Pl, lam, axis=1)

    def fwhm(sig):
        idx = np.where(sig >= sig.max() / 2)[0]
        return tf[idx[-1]] - tf[idx[0]]

    uv = (lam >= 300e-9) & (lam <= 400e-9); red = (lam >= 590e-9) & (lam <= 650e-9)
    spec = np.trapezoid(Pl, tf, axis=0)
    Ts = np.linspace(3000, 60000, 1141)
    ok = spec > 0
    Tcol = min(Ts, key=lambda Tt: np.sum((np.log(spec[ok]) - np.log(planck(lam[ok], Tt))
                                          - np.mean(np.log(spec[ok]) - np.log(planck(lam[ok], Tt)))) ** 2))
    ip = Ptot.argmax()
    tau_c = (kap[ip, :N, :] * (np.diff(np.concatenate([[0.0], yf])) * R[ip])[:, None]).sum(0).max()
    out.update(photons=np.trapezoid(rate, tf), Ppeak=Ptot.max(), tau=float(tau_c),
               fwhm_uv=fwhm(np.trapezoid(Pl[:, uv], lam[uv], axis=1)),
               fwhm_red=fwhm(np.trapezoid(Pl[:, red], lam[red], axis=1)), T_colour=Tcol,
               x_centre=float(x[im, 0]), t_peak_minus_tRmin=float(tf[ip] - tf[im]))
    ne_full = np.zeros_like(T); ne_full[hot] = ne
    out["ne_centre"] = float(ne_full[ip, 0]); out["ne_max"] = float(ne_full.max())
    out["T_lightweighted"] = float(Tfin[ip, 0])
    out["profile"] = (y, Tfin[ip], x[ip], R[ip])
    return out
