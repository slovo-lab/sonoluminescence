"""Minimal single-bubble sonoluminescence model (argon bubble in water).

Dynamics : Keller-Miksis equation, van der Waals hard core.
Heating  : uniform bubble, Peclet-dependent polytropic exponent + conductive loss
           (Hilgenfeldt-Grossmann-Lohse 1999 style).
Light    : Saha ionisation + Kramers-Unsoeld continuum absorption (ff + bf,
           hydrogenic), Kirchhoff's law for a homogeneous sphere.
"""
import numpy as np
from scipy.integrate import solve_ivp

kB = 1.380649e-23; h_pl = 6.62607015e-34; c0 = 2.99792458e8
me = 9.1093837e-31; qe = 1.602176634e-19; eps0 = 8.8541878128e-12
e2 = qe**2 / (4 * np.pi * eps0)  # Gaussian e^2 in SI units (J*m)

# water, 20 C
rho = 998.0; cw = 1481.0; sigma = 0.0725; mu = 1.0e-3; P0 = 101325.0; T0 = 293.15
GAMMA = 5.0 / 3.0
CHI0 = 2.0e-5          # thermal diffusivity of argon at ambient, m^2/s
E_ION = 15.76 * qe     # argon


def run(R0=4.5e-6, Pa_atm=1.3, f=26.5e3, cycles=3, dE=0.0,
        rtol=1e-10, ipd="none"):
    Pa = Pa_atm * 101325.0
    w = 2 * np.pi * f
    hc = R0 / 8.86                     # van der Waals hard core radius (argon)
    pg0 = P0 + 2 * sigma / R0
    N = pg0 * (4 * np.pi / 3) * (R0**3 - hc**3) / (kB * T0)

    def rhs(t, y):
        R, U, T = y
        d = R**3 - hc**3
        chi = CHI0 * d / (R0**3 - hc**3)          # diffusivity ~ 1/density
        Pe = R * abs(U) / chi + 1e-30
        kap = 1 + (GAMMA - 1) * np.exp(-5.8 / Pe**0.6)
        dT = -(kap - 1) * 3 * R**2 * U * T / d - chi * (T - T0) / R**2
        pg = pg0 * (R0**3 - hc**3) / d * T / T0
        dpg = pg * (dT / T - 3 * R**2 * U / d)
        Pac = -Pa * np.sin(w * t)
        M = U / cw
        num = ((1 + M) * (pg - P0 - Pac) / rho + R * dpg / (rho * cw)
               - 4 * mu * U / (rho * R) - 2 * sigma / (rho * R)
               - 1.5 * U**2 * (1 - M / 3))
        return [U, num / ((1 - M) * R), dT]

    sol = solve_ivp(rhs, [0, cycles / f], [R0, 0.0, T0], method="LSODA",
                    rtol=rtol, atol=[1e-13, 1e-9, 1e-6], max_step=1 / f / 400,
                    dense_output=True)
    return dict(sol=sol, R0=R0, hc=hc, N=N, f=f, Pa_atm=Pa_atm, dE=dE, ipd=ipd)


def ionisation(T, n, dE=0.0):
    """Saha degree of ionisation for argon (g+/g0 ~ 6)."""
    S = 2 * 6 * (2 * np.pi * me * kB * T / h_pl**2) ** 1.5 * np.exp(-(E_ION - dE) / (kB * T)) / n
    return (-S + np.sqrt(S * S + 4 * S)) / 2


N_LIQ = 2.11e28                 # number density of liquid argon, m^-3
DE_LIQ = (1.75 + 0.17) * qe     # lowering in liquid argon: |P+| + |V0| (15.76 -> 13.84 eV)


def lowering(T, n, mode="none", dE0=0.0, start=None, iters=400):
    """Density-dependent lowering of the ionisation potential, in joules.

    mode "none"    : fixed shift dE0 only
    mode "neutral" : dE0 + polarisation of the dense neutral gas, scaled linearly
                     from the measured liquid-argon value
    mode "plasma"  : dE0 + Stewart-Pyatt screening by free charges (self-consistent)
    mode "both"    : neutral + plasma
    Returns (dE, alpha).  `start` = initial guess for alpha (None -> ideal Saha).
    """
    T = np.asarray(T, float); n = np.asarray(n, float)
    base = dE0 + (DE_LIQ * n / N_LIQ if mode in ("neutral", "both") else 0.0)
    base = np.broadcast_to(base, np.broadcast(T, n).shape).astype(float)
    a = ionisation(T, n, base) if start is None else np.full(base.shape, float(start))
    if mode not in ("plasma", "both"):
        return base, ionisation(T, n, base)
    dE = base
    for _ in range(iters):
        ne = np.maximum(a * n, 1.0)
        lamD = np.sqrt(eps0 * kB * T / (2 * qe**2 * ne))      # electrons + singly charged ions
        K = e2 / (lamD * kB * T)
        dE = base + kB * T / 4 * ((6 * K + 1) ** (2 / 3) - 1)   # Stewart-Pyatt, z = z* = 1
        dE = np.minimum(dE, E_ION)
        a_new = ionisation(T, n, dE)
        if np.max(np.abs(a_new - a) / (a + 1e-300)) < 1e-10:
            a = a_new
            break
        a = 0.5 * a + 0.5 * a_new
    return dE, a


def kappa(lam, T, n, dE=0.0):
    """Kramers-Unsoeld absorption coefficient (1/m), incl. stimulated emission."""
    nu = c0 / lam
    x = h_pl * nu / (kB * T)
    a = ionisation(T, n, dE)
    n0 = n * (1 - a)
    k = (16 * np.pi**2 / (3 * np.sqrt(3))) * e2**3 * kB * T * n0 / (h_pl**4 * c0 * nu**3)
    return k * np.exp(-(E_ION - dE) / (kB * T) + x) * (1 - np.exp(-x))


# --- electron-neutral bremsstrahlung ------------------------------------------------
# Momentum-transfer cross-section of argon for electrons.  APPROXIMATE values typed in
# from memory of the standard swarm-derived sets (Milloy et al. 1977 below 4 eV, Phelps-
# type compilations above); expect errors of a few tens of per cent.  To be replaced by
# an LXCat data set.  Units: eV, 1e-20 m^2.
_SIG_E = np.array([1e-3, 0.01, 0.05, 0.1, 0.23, 0.5, 1, 2, 3, 5, 7, 10, 12, 15, 20, 30, 50, 100])
_SIG_Q = np.array([7.5, 4.6, 1.6, 0.59, 0.095, 0.35, 1.38, 3.0, 4.1, 6.7, 10.0, 15.0, 15.2, 14.1, 11.0, 6.5, 3.5, 1.9])
ALPHA_FS = e2 / (h_pl / (2 * np.pi) * c0)
SIG_SCALE = 1.0


def sigma_tr(E):
    """Momentum-transfer cross-section (m^2) at electron energy E (J)."""
    x = np.log(np.clip(E / qe, _SIG_E[0], _SIG_E[-1]))
    return SIG_SCALE * 1e-20 * np.exp(np.interp(x, np.log(_SIG_E), np.log(_SIG_Q)))


def kappa_en_per_pair(lam, T, sigma=sigma_tr):
    """Electron-neutral free-free absorption coefficient divided by n_e*n_0 (m^5).

    Emission from the Dalgarno-Lane / Kas'yanov-Starostin cross-section
        dsigma/domega = 4 alpha/(3 pi omega m c^2) sqrt(E'/E) [E sigma(E') + E' sigma(E)],
    E' = E - hbar*omega, averaged over a Maxwellian; absorption from Kirchhoff's law.
    lam: (L,) array, T: scalar.
    """
    hbar = h_pl / (2 * np.pi)
    w = 2 * np.pi * c0 / np.asarray(lam, float)[:, None]
    s = np.linspace(0, 7, 561)[None, :]
    Ep = kB * T * s**2                       # E' (energy after emission)
    E = Ep + hbar * w
    dsdw = 4 * ALPHA_FS / (3 * np.pi * w * me * c0**2) * np.sqrt(Ep / E) * (E * sigma(Ep) + Ep * sigma(E))
    f = 2 / np.sqrt(np.pi) * (kB * T) ** -1.5 * np.sqrt(E) * np.exp(-E / (kB * T))
    v = np.sqrt(2 * E / me)
    integrand = f * v * dsdw * (2 * kB * T * s)          # dE = 2 kT s ds
    eps_w = hbar * w[:, 0] * np.trapezoid(integrand, s[0], axis=1)   # W m^3 per (rad/s)
    nu = c0 / np.asarray(lam, float)
    Bnu = 2 * h_pl * nu**3 / c0**2 / np.expm1(h_pl * nu / (kB * T))
    return eps_w / (2 * Bnu)


def emissivity(tau):
    tau = np.asarray(tau, float)
    big = 1 + np.exp(-2 * tau) / np.maximum(tau, 1e-300) + np.expm1(-2 * tau) / (2 * np.maximum(tau, 1e-300) ** 2)
    return np.where(tau < 1e-3, 4 * tau / 3, big)


def planck(lam, T):
    return 2 * h_pl * c0**2 / lam**5 / np.expm1(h_pl * c0 / (lam * kB * T))


def analyse(res, lam=np.linspace(200e-9, 800e-9, 121), energy=False, eneutral=False):
    sol, f, hc, N = res["sol"], res["f"], res["hc"], res["N"]
    t_end = sol.t[-1]
    # last full cycle
    tc = np.linspace(t_end - 1 / f, t_end, 40001)
    R = sol.sol(tc)[0]
    Rmax = R.max()
    # main collapse = first minimum after the maximum
    i0 = R.argmax()
    seg = R[i0:]
    # first local minimum
    j = np.where((seg[1:-1] < seg[:-2]) & (seg[1:-1] <= seg[2:]))[0][0] + 1
    tm = tc[i0 + j]
    # refine on fine grid (2 ns window, 0.5 ps steps)
    tf = np.linspace(tm - 1e-9, tm + 1e-9, 4001)
    Rf, Uf, Tf = sol.sol(tf)
    n = N / (4 * np.pi / 3 * Rf**3)
    L = lam[None, :]
    mode = res.get("ipd", "none")
    if energy:
        # Post-hoc energy bookkeeping: the dynamics give 1.5 kT per atom as if nothing
        # ionised.  Find the temperature T' at which thermal energy of atoms + electrons
        # plus the energy spent on ionisation equals that budget (bisection in T').
        lo = np.full_like(Tf, 300.0); hi = Tf.copy()
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            d_, a_ = lowering(mid, n, mode, res["dE"])
            u = 1.5 * kB * mid * (1 + a_) + a_ * (E_ION - d_)
            hi = np.where(u > 1.5 * kB * Tf, mid, hi); lo = np.where(u > 1.5 * kB * Tf, lo, mid)
        Tf = 0.5 * (lo + hi)
    dEt, at = lowering(Tf, n, mode, res["dE"])
    k = kappa(L, Tf[:, None], n[:, None], dEt[:, None])
    k_ion = k
    if eneutral:
        Tg = np.geomspace(max(Tf.min(), 300.0), Tf.max() * 1.001, 90)
        tab = np.log(np.array([kappa_en_per_pair(lam, T) for T in Tg]) + 1e-300)   # (T, L)
        kp = np.exp(np.array([np.interp(np.log(Tf), np.log(Tg), tab[:, j]) for j in range(len(lam))]).T)
        k = k + kp * (at * n)[:, None] * ((1 - at) * n)[:, None]
    Pl = 4 * np.pi**2 * Rf[:, None] ** 2 * planck(L, Tf[:, None]) * emissivity(k * Rf[:, None])  # W/m
    Nph_rate = np.trapezoid(Pl * L / (h_pl * c0), lam, axis=1)
    out = dict(Rmax=Rmax, Rmin=Rf.min(), Tmax=Tf.max(), Umax=abs(Uf).max(),
               alpha=at.max(), dE_eV=dEt[at.argmax()] / qe,
               ion_energy_ratio=(at * E_ION / (1.5 * kB * Tf)).max(),
               photons=np.trapezoid(Nph_rate, tf), Ppeak=np.trapezoid(Pl, lam, axis=1).max(),
               tau_peak=(k * Rf[:, None]).max())

    def fwhm(sig):
        hm = sig.max() / 2
        idx = np.where(sig >= hm)[0]
        return tf[idx[-1]] - tf[idx[0]]

    uv = (lam >= 300e-9) & (lam <= 400e-9)
    red = (lam >= 590e-9) & (lam <= 650e-9)
    out["fwhm_uv"] = fwhm(np.trapezoid(Pl[:, uv], lam[uv], axis=1))
    out["fwhm_red"] = fwhm(np.trapezoid(Pl[:, red], lam[red], axis=1))
    out["fwhm_all"] = fwhm(Nph_rate)
    ip = np.trapezoid(Pl, lam, axis=1).argmax()
    out["frac_ion_uv"] = float(k_ion[ip, np.abs(lam - 300e-9).argmin()] / k[ip, np.abs(lam - 300e-9).argmin()])
    out["frac_ion_red"] = float(k_ion[ip, np.abs(lam - 700e-9).argmin()] / k[ip, np.abs(lam - 700e-9).argmin()])
    # time-integrated spectrum -> best-fit blackbody colour temperature (shape only)
    spec = np.trapezoid(Pl, tf, axis=0)
    Ts = np.linspace(5000, 60000, 1101)
    best = min(Ts, key=lambda T: np.sum((np.log(spec) - np.log(planck(lam, T))
                                         - np.mean(np.log(spec) - np.log(planck(lam, T)))) ** 2))
    out["T_colour"] = best
    out["curve"] = (tc, R)
    out["fine"] = (tf, Rf, Tf, Nph_rate)
    return out


if __name__ == "__main__":
    import sys
    print("Pa  Rmax/um Rmin/um Tmax/K  alpha   photons  FWHM_uv FWHM_red/ps  Tcol  tau")
    for Pa in [1.2, 1.25, 1.3, 1.35, 1.4]:
        o = analyse(run(Pa_atm=Pa))
        print(f"{Pa:4.2f} {o['Rmax']*1e6:6.1f} {o['Rmin']*1e6:6.3f} {o['Tmax']:8.0f} {o['alpha']:.2e} "
              f"{o['photons']:.2e} {o['fwhm_uv']*1e12:6.0f} {o['fwhm_red']*1e12:6.0f} {o['T_colour']:7.0f} {o['tau_peak']:.2e}")
