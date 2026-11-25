"""v0.4: water vapour in the bubble (issue #2).

Dynamics : Keller-Miksis wall equation; heat and vapour exchange through diffusive
           boundary layers (Toegel, Gompf, Pecha & Lohse, PRL 85, 3165 (2000)).
           Vapour is carried as inert H2O in the dynamics.
Chemistry: post-hoc, at fixed volume and fixed energy: equilibrium H2O <-> 2H + O from
           partition functions, plus Saha ionisation of Ar, H and O with the v0.2
           lowering of the ionisation potential.
Light    : Kramers-Unsoeld continuum of Ar, H and O, plus electron-neutral
           bremsstrahlung (v0.3; the argon cross-section is used for every neutral).
"""
import numpy as np
from scipy.integrate import solve_ivp
import sbsl
from sbsl import (kB, h_pl, c0, me, qe, eps0, e2, rho, cw, sigma, mu, P0, T0,
                  emissivity, planck, kappa_en_per_pair, N_LIQ, DE_LIQ)

AMU = 1.66053907e-27
B_AR = 4.37e-29                     # co-volume per argon atom, m^3: reproduces the hard
                                    # core R0/8.86 used in v0.1-v0.3
B_W = B_AR * 3.05 / 3.20            # water, scaled by the ratio of van der Waals b
D0_DIFF = 23.55e-6                  # diffusivity of vapour in argon at n = 2.5e25 m^-3
N_REF = 2.5e25
KAPPA_G = 0.018                     # thermal conductivity of argon, W/(m K)
THETA_V = np.array([2295.0, 5255.0, 5400.0])   # vibrational temperatures of H2O, K
THETA_R = (40.1, 20.9, 13.4)                   # rotational temperatures of H2O, K
D_ATOM = 9.51 * qe                  # H2O -> 2H + O at 0 K
E_I = {"Ar": 15.76 * qe, "H": 13.598 * qe, "O": 13.618 * qe}
G_I = {"Ar": 12.0, "H": 1.0, "O": 8.0 / 9.0}   # 2 g+/g0
PV_20C = 2337.0                     # saturated vapour pressure of water at 20 C, Pa


def e_w(T):
    """Internal energy of one H2O molecule (J), translation + rotation + vibration."""
    T = np.asarray(T, float)[..., None]
    return (kB * T * (3 + np.sum((THETA_V / T) / np.expm1(THETA_V / T), axis=-1, keepdims=True)))[..., 0]


def c_w(T):
    th = THETA_V / T
    return kB * (3 + np.sum(th**2 * np.exp(th) / np.expm1(th) ** 2))


def run(R0=4.5e-6, Pa_atm=1.3, f=26.5e3, cycles=3, pv=PV_20C, rtol=1e-9):
    Pa = Pa_atm * 101325.0
    w = 2 * np.pi * f
    n_sat = pv / (kB * T0)
    V0 = 4 * np.pi / 3 * R0**3
    Nw0 = n_sat * V0
    # argon amount: total pressure balances P0 + 2 sigma / R0
    N_ar = (P0 + 2 * sigma / R0 - pv) * V0 / (kB * T0)
    N_ar = N_ar / (1 + (P0 + 2 * sigma / R0 - pv) * B_AR / (kB * T0))   # co-volume correction

    def rhs(t, y):
        R, U, T, Nw = y
        Nw = max(Nw, 0.0)
        V = 4 * np.pi / 3 * R**3
        dV = 4 * np.pi * R**2 * U
        Ntot = N_ar + Nw
        free = V - N_ar * B_AR - Nw * B_W
        D = D0_DIFF * N_REF * V / Ntot
        chi = KAPPA_G * V / (2.5 * kB * N_ar + 4 * kB * Nw)
        aU = abs(U) + 1e-30
        l_d = min(np.sqrt(R * D / aU), R / np.pi)
        l_t = min(np.sqrt(R * chi / aU), R / np.pi)
        dNw = 4 * np.pi * R**2 * D * (n_sat - Nw / V) / l_d
        Q = 4 * np.pi * R**2 * KAPPA_G * (T0 - T) / l_t
        p = Ntot * kB * T / free
        Cv = 1.5 * kB * N_ar + c_w(T) * Nw
        dT = (Q - p * dV + (4 * kB * T0 - float(e_w(T))) * dNw) / Cv
        dp = p * (dNw / Ntot + dT / T - (dV - B_W * dNw) / free)
        Pac = -Pa * np.sin(w * t)
        M = U / cw
        num = ((1 + M) * (p - P0 - Pac) / rho + R * dp / (rho * cw)
               - 4 * mu * U / (rho * R) - 2 * sigma / (rho * R)
               - 1.5 * U**2 * (1 - M / 3))
        return [U, num / ((1 - M) * R), dT, dNw]

    sol = solve_ivp(rhs, [0, cycles / f], [R0, 0.0, T0, Nw0], method="LSODA", rtol=rtol,
                    atol=[1e-13, 1e-9, 1e-6, max(Nw0, 1.0) * 1e-9], max_step=1 / f / 400,
                    dense_output=True)
    return dict(sol=sol, R0=R0, N_ar=N_ar, f=f, Pa_atm=Pa_atm, pv=pv)


# ---------------------------------------------------------------- equilibrium state
def _q_trans(m, T):
    return (2 * np.pi * m * kB * T / h_pl**2) ** 1.5


def K_diss(T):
    """Equilibrium constant n_H^2 n_O / n_H2O (m^-6) from partition functions."""
    q_h = _q_trans(1.008 * AMU, T) * 2
    q_o = _q_trans(15.999 * AMU, T) * 9
    rot = np.sqrt(np.pi * T**3 / np.prod(THETA_R)) / 2
    vib = np.prod(1 / (-np.expm1(-THETA_V / np.asarray(T, float)[..., None])), axis=-1)
    q_w = _q_trans(18.015 * AMU, T) * rot * vib
    return q_h**2 * q_o / q_w * np.exp(-D_ATOM / (kB * T))


def dissociation(T, nw, eta=None):
    """Fraction x of H2O split into 2H + O; nw = vapour number density (free volume).

    eta = packing fraction (1 - free volume / volume).  If given, the reaction, which
    makes three particles out of one, is additionally suppressed by the hard-core
    factor exp(-2 eta / (1 - eta)) (equal co-volumes assumed for all species).
    """
    T = np.asarray(T, float); nw = np.asarray(nw, float)
    if np.all(nw <= 0):
        return np.zeros(np.broadcast(T, nw).shape)
    c = K_diss(T) / np.maximum(nw, 1e-300) ** 2        # 4 x^3 / (1 - x) = c
    if eta is not None:
        c = c * np.exp(-2 * eta / (1 - eta))
    lo = np.zeros(np.broadcast(T, nw).shape); hi = np.ones_like(lo)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        big = 4 * mid**3 > c * (1 - mid)
        hi = np.where(big, mid, hi); lo = np.where(big, lo, mid)
    return 0.5 * (lo + hi)


def ion_state(T, n_ar, n_h, n_o, n_all, chemistry=True, lowering=True):
    """Electron density and lowering dE for a mixture of Ar, H and O (all per m^3)."""
    T = np.asarray(T, float)
    base = DE_LIQ * n_all / N_LIQ if lowering else np.zeros_like(T)
    dE = base
    ne = np.full_like(T, 1.0)
    dens = {"Ar": n_ar, "H": n_h, "O": n_o}
    for _ in range(80):
        S = {s: G_I[s] * _q_trans(me, T) * np.exp(-(E_I[s] - dE) / (kB * T)) for s in dens}
        lo = np.full_like(T, -10.0); hi = np.log(np.maximum(n_ar + n_h + n_o, 1.0))
        for _ in range(50):                       # solve sum_s n_s S_s/(ne + S_s) = ne
            mid = 0.5 * (lo + hi); x = np.exp(mid)
            g = sum(dens[s] * S[s] / (x + S[s]) for s in dens) - x
            lo = np.where(g > 0, mid, lo); hi = np.where(g > 0, hi, mid)
        ne_new = np.exp(0.5 * (lo + hi))
        if not lowering:
            ne = ne_new
            break
        lamD = np.sqrt(eps0 * kB * T / (2 * qe**2 * np.maximum(ne_new, 1.0)))
        dE_new = np.minimum(base + kB * T / 4 * ((6 * e2 / (lamD * kB * T) + 1) ** (2 / 3) - 1), 12 * qe)
        done = np.max(np.abs(ne_new - ne) / (ne + 1e-300)) < 1e-8
        ne = 0.5 * ne + 0.5 * ne_new if not done else ne_new
        dE = 0.5 * dE + 0.5 * dE_new
        if done:
            break
    S = {s: G_I[s] * _q_trans(me, T) * np.exp(-(E_I[s] - dE) / (kB * T)) for s in dens}
    ions = {s: dens[s] * S[s] / (ne + S[s]) for s in dens}
    return ne, dE, ions


def analyse(res, lam=np.linspace(200e-9, 800e-9, 121), chemistry=True, npts=2001, excluded=False):
    sol, f, N_ar = res["sol"], res["f"], res["N_ar"]
    t_end = sol.t[-1]
    tc = np.linspace(t_end - 1 / f, t_end, 40001)
    R = sol.sol(tc)[0]
    i0 = R.argmax(); seg = R[i0:]
    j = np.where((seg[1:-1] < seg[:-2]) & (seg[1:-1] <= seg[2:]))[0][0] + 1
    tm = tc[i0 + j]
    tf = np.linspace(tm - 1e-9, tm + 1e-9, npts)
    Rf, Uf, Tm, Nw = sol.sol(tf)
    Nw = np.maximum(Nw, 0.0)
    V = 4 * np.pi / 3 * Rf**3
    free = V - N_ar * B_AR - Nw * B_W
    n_ar = N_ar / V; n_w = Nw / V
    budget = 1.5 * kB * N_ar * Tm + Nw * e_w(Tm)       # thermal energy from the dynamics

    def state(T):
        x = dissociation(T, Nw / free, (1 - free / V) if excluded else None) if chemistry else np.zeros_like(T)
        n_h, n_o, n_m = 2 * x * n_w, x * n_w, (1 - x) * n_w
        ne, dE, ions = ion_state(T, n_ar, n_h, n_o, n_ar + n_h + n_o + n_m)
        E = (1.5 * kB * T * (N_ar + ne * V) + Nw * ((1 - x) * e_w(T) + x * (4.5 * kB * T + D_ATOM))
             + V * sum(ions[s] * (E_I[s] - dE) for s in ions))
        return E, x, ne, dE, ions, (n_h, n_o, n_m)

    lo = np.full_like(Tm, 250.0); hi = Tm.copy() + 1.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        over = state(mid)[0] > budget
        hi = np.where(over, mid, hi); lo = np.where(over, lo, mid)
    Tf = 0.5 * (lo + hi)
    _, x, ne, dE, ions, (n_h, n_o, n_m) = state(Tf)
    dens = {"Ar": n_ar, "H": n_h, "O": n_o}
    L = lam[None, :]; nu = c0 / L; xx = h_pl * nu / (kB * Tf[:, None])
    pref = (16 * np.pi**2 / (3 * np.sqrt(3))) * e2**3 * kB * Tf[:, None] / (h_pl**4 * c0 * nu**3)
    k_ion = sum(pref * (dens[s] - ions[s])[:, None] * np.exp(-(E_I[s] - dE)[:, None] / (kB * Tf[:, None]) + xx)
                for s in dens) * (1 - np.exp(-xx))
    Tg = np.geomspace(max(Tf.min(), 300.0), Tf.max() * 1.001, 90)
    tab = np.log(np.array([kappa_en_per_pair(lam, T) for T in Tg]) + 1e-300)
    kp = np.exp(np.array([np.interp(np.log(Tf), np.log(Tg), tab[:, j]) for j in range(len(lam))]).T)
    n_neut = sum(dens[s] - ions[s] for s in dens) + n_m
    k = k_ion + kp * (ne * n_neut)[:, None]
    Pl = 4 * np.pi**2 * Rf[:, None] ** 2 * planck(L, Tf[:, None]) * emissivity(k * Rf[:, None])
    rate = np.trapezoid(Pl * L / (h_pl * c0), lam, axis=1)
    Ptot = np.trapezoid(Pl, lam, axis=1)
    ip = Ptot.argmax()

    def fwhm(sig):
        idx = np.where(sig >= sig.max() / 2)[0]
        return tf[idx[-1]] - tf[idx[0]]

    uv = (lam >= 300e-9) & (lam <= 400e-9); red = (lam >= 590e-9) & (lam <= 650e-9)
    spec = np.trapezoid(Pl, tf, axis=0)
    Ts = np.linspace(3000, 60000, 1141)
    Tcol = min(Ts, key=lambda T: np.sum((np.log(spec) - np.log(planck(lam, T))
                                         - np.mean(np.log(spec) - np.log(planck(lam, T)))) ** 2))
    return dict(Rmax=R.max(), Rmin=Rf.min(), Mach=abs(Uf).max() / cw, T_dyn=Tm.max(), Tmax=Tf.max(),
                vapour_frac=(Nw / (N_ar + Nw))[Rf.argmin()], packing=(1 - free / V).max(), x=x[ip], ne_frac=(ne / (n_ar + n_w))[ip],
                e_from_HO=((ions["H"] + ions["O"]) / np.maximum(ne, 1e-300))[ip], dE_eV=dE[ip] / qe,
                photons=np.trapezoid(rate, tf), Ppeak=Ptot.max(), tau=(k * Rf[:, None]).max(),
                fwhm_uv=fwhm(np.trapezoid(Pl[:, uv], lam[uv], axis=1)),
                fwhm_red=fwhm(np.trapezoid(Pl[:, red], lam[red], axis=1)), T_colour=Tcol)
