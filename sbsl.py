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
        rtol=1e-10):
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
    return dict(sol=sol, R0=R0, hc=hc, N=N, f=f, Pa_atm=Pa_atm, dE=dE)


def ionisation(T, n, dE=0.0):
    """Saha degree of ionisation for argon (g+/g0 ~ 6)."""
    S = 2 * 6 * (2 * np.pi * me * kB * T / h_pl**2) ** 1.5 * np.exp(-(E_ION - dE) / (kB * T)) / n
    return (-S + np.sqrt(S * S + 4 * S)) / 2


def kappa(lam, T, n, dE=0.0):
    """Kramers-Unsoeld absorption coefficient (1/m), incl. stimulated emission."""
    nu = c0 / lam
    x = h_pl * nu / (kB * T)
    a = ionisation(T, n, dE)
    n0 = n * (1 - a)
    k = (16 * np.pi**2 / (3 * np.sqrt(3))) * e2**3 * kB * T * n0 / (h_pl**4 * c0 * nu**3)
    return k * np.exp(-(E_ION - dE) / (kB * T) + x) * (1 - np.exp(-x))


def emissivity(tau):
    tau = np.asarray(tau, float)
    big = 1 + np.exp(-2 * tau) / np.maximum(tau, 1e-300) + np.expm1(-2 * tau) / (2 * np.maximum(tau, 1e-300) ** 2)
    return np.where(tau < 1e-3, 4 * tau / 3, big)


def planck(lam, T):
    return 2 * h_pl * c0**2 / lam**5 / np.expm1(h_pl * c0 / (lam * kB * T))


def analyse(res, lam=np.linspace(200e-9, 800e-9, 121)):
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
    k = kappa(L, Tf[:, None], n[:, None], res["dE"])
    Pl = 4 * np.pi**2 * Rf[:, None] ** 2 * planck(L, Tf[:, None]) * emissivity(k * Rf[:, None])  # W/m
    Nph_rate = np.trapezoid(Pl * L / (h_pl * c0), lam, axis=1)
    out = dict(Rmax=Rmax, Rmin=Rf.min(), Tmax=Tf.max(), Umax=abs(Uf).max(),
               alpha=ionisation(Tf.max(), n[Tf.argmax()], res["dE"]),
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
