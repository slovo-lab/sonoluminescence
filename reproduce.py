"""Reproduce every number quoted in README.md.  Usage: python reproduce.py"""
import numpy as np
import sbsl
from sbsl import run, analyse


def row(tag, res):
    o = analyse(res)
    n = res["N"] / (4 * np.pi / 3 * o["Rmin"] ** 3)
    print(f"{tag:22s} Rmax {o['Rmax']*1e6:5.1f} um | Rmin {o['Rmin']*1e6:.3f} um | Tmax {o['Tmax']:6.0f} K | "
          f"Mach {o['Umax']/sbsl.cw:.2f} | n {n:.2e} m^-3 | photons {o['photons']:.2e} | "
          f"FWHM uv/red {o['fwhm_uv']*1e12:4.0f}/{o['fwhm_red']*1e12:4.0f} ps | "
          f"peak {o['Ppeak']*1e3:.3f} mW | tau {o['tau_peak']:.1e} | T_colour {o['T_colour']:.0f} K")


print("== driving pressure scan (R0 = 4.5 um, 26.5 kHz)")
for Pa in [1.2, 1.25, 1.3, 1.35, 1.4]:
    row(f"Pa={Pa} atm", run(Pa_atm=Pa))
print("== ambient radius scan")
for R0 in [3e-6, 4.5e-6, 6e-6, 7.5e-6]:
    for Pa in [1.3, 1.4]:
        row(f"R0={R0*1e6:.1f} um Pa={Pa}", run(R0=R0, Pa_atm=Pa))
print("== numerical convergence (Pa = 1.35 atm)")
for rt in [1e-8, 1e-10, 1e-12]:
    row(f"rtol={rt:g}", run(Pa_atm=1.35, rtol=rt))
row("cycles=6", run(Pa_atm=1.35, cycles=6))
print("== ionisation-potential lowering (Pa = 1.3 atm)")
for dE in [0, 1, 2, 3]:
    row(f"dE={dE} eV", run(Pa_atm=1.3, dE=dE * sbsl.qe))
