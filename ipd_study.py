"""Issue #1: does density-dependent ionisation-potential lowering close the brightness gap?"""
import numpy as np
import sbsl
from sbsl import run, analyse, lowering

print("mode     Pa   Tmax   alpha    dE/eV  photons   peak/mW  tau      FWHM uv/red  T_colour  E_ion*alpha/(1.5kT)")
out = {}
for mode in ["none", "neutral", "plasma", "both"]:
    for Pa in [1.2, 1.25, 1.3, 1.35, 1.4]:
        o = analyse(run(Pa_atm=Pa, ipd=mode)); out[mode, Pa] = o
        print(f"{mode:8s} {Pa:4.2f} {o['Tmax']:6.0f} {o['alpha']:.2e} {o['dE_eV']:5.2f}  {o['photons']:.2e}  "
              f"{o['Ppeak']*1e3:7.3f}  {o['tau_peak']:.1e}  {o['fwhm_uv']*1e12:4.0f}/{o['fwhm_red']*1e12:4.0f}     "
              f"{o['T_colour']:6.0f}    {o['ion_energy_ratio']:.2f}")
print("\nwith ionisation energy taken out of the heat budget (post-hoc):")
for mode in ["none", "both"]:
    for Pa in [1.2, 1.25, 1.3, 1.35, 1.4]:
        o = analyse(run(Pa_atm=Pa, ipd=mode), energy=True); out[mode + "+E", Pa] = o
        print(f"{mode+'+E':8s} {Pa:4.2f} {o['Tmax']:6.0f} {o['alpha']:.2e} {o['dE_eV']:5.2f}  {o['photons']:.2e}  "
              f"{o['Ppeak']*1e3:7.3f}  {o['tau_peak']:.1e}  {o['fwhm_uv']*1e12:4.0f}/{o['fwhm_red']*1e12:4.0f}     "
              f"{o['T_colour']:6.0f}    {o['ion_energy_ratio']:.2f}")
print("\nbrightness ratio 1.4 atm / 1.3 atm:")
for mode in ["none", "neutral", "plasma", "both", "none+E", "both+E"]:
    print(f"  {mode:8s} {out[mode,1.4]['photons']/out[mode,1.3]['photons']:.0f}x")
print("\nuniqueness of the self-consistent solution (start from alpha=1e-12 vs alpha=1):")
for T in [10000, 15000, 20000, 25000]:
    for n in [5e27, 1e28]:
        lo = lowering(T, n, "both", start=1e-12)[1]; hi = lowering(T, n, "both", start=1.0)[1]
        print(f"  T={T} n={n:.0e}  alpha_low={float(lo):.3e} alpha_high={float(hi):.3e}")
