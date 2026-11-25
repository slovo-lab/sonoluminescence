"""Issue #3: electron-neutral bremsstrahlung.  Baseline = v0.2 (lowering 'both' + energy budget)."""
import numpy as np
import sbsl
from sbsl import run, analyse

# --- check against the classical (Drude) limit for a constant cross-section, h*nu << kT
T = 15000.0; sig0 = 1e-19; lam = np.array([200e-6, 50e-6])
k_num = sbsl.kappa_en_per_pair(lam, T, sigma=lambda E: sig0 + 0 * E)
w = 2 * np.pi * sbsl.c0 / lam
vbar = np.sqrt(8 * sbsl.kB * T / (np.pi * sbsl.me))
k_drude = sbsl.qe**2 / (sbsl.eps0 * sbsl.me * sbsl.c0 * w**2) * (4 / 3) * sig0 * vbar
print("Drude-limit check, numerical/analytic:", np.round(k_num / k_drude, 3), " (h nu / kT =", np.round(sbsl.h_pl * sbsl.c0 / lam / (sbsl.kB * T), 3), ")")

print("\ncase        Pa   Tmax   photons   peak/mW  tau      FWHM uv/red  T_colour  ion share of absorption uv/red")
out = {}
for tag, en in [("v0.2", False), ("v0.3 +e-n", True)]:
    for Pa in [1.2, 1.25, 1.3, 1.35, 1.4]:
        o = analyse(run(Pa_atm=Pa, ipd="both"), energy=True, eneutral=en); out[tag, Pa] = o
        print(f"{tag:10s} {Pa:4.2f} {o['Tmax']:6.0f}  {o['photons']:.2e}  {o['Ppeak']*1e3:7.3f}  {o['tau_peak']:.1e}  "
              f"{o['fwhm_uv']*1e12:4.0f}/{o['fwhm_red']*1e12:4.0f}     {o['T_colour']:6.0f}    {o['frac_ion_uv']:.2f}/{o['frac_ion_red']:.2f}")
print("\ngain from e-n bremsstrahlung (photons):", {Pa: round(out['v0.3 +e-n', Pa]['photons'] / out['v0.2', Pa]['photons'], 2) for Pa in [1.2, 1.3, 1.4]})
print("brightness ratio 1.4/1.3 atm:", round(out['v0.3 +e-n', 1.4]['photons'] / out['v0.3 +e-n', 1.3]['photons'], 1))
print("\nsensitivity to the cross-section (Pa = 1.35 atm):")
for sc in [0.7, 1.0, 1.3]:
    sbsl.SIG_SCALE = sc
    o = analyse(run(Pa_atm=1.35, ipd="both"), energy=True, eneutral=True)
    print(f"  sigma x{sc}: photons {o['photons']:.2e}, peak {o['Ppeak']*1e3:.3f} mW, T_colour {o['T_colour']:.0f} K")
sbsl.SIG_SCALE = 1.0
print("\nbubble size scan with e-n (Pa = 1.35 atm):")
for R0 in [3e-6, 4.5e-6, 6e-6]:
    o = analyse(run(R0=R0, Pa_atm=1.35, ipd="both"), energy=True, eneutral=True)
    print(f"  R0={R0*1e6:.1f} um: photons {o['photons']:.2e}, peak {o['Ppeak']*1e3:.3f} mW, FWHM uv/red {o['fwhm_uv']*1e12:.0f}/{o['fwhm_red']*1e12:.0f} ps")
