"""Issue #2: how much does water vapour cool the collapse?"""
import time
import numpy as np
import vapor
from vapor import run, analyse, PV_20C

HDR = ("case              Pa   Rmax  Rmin   Mach  vap%  T_dyn  Tmax   diss  e-/atom  e-from-H,O  dE/eV  photons   peak/mW  tau    "
       "FWHM uv/red  T_col")


def row(tag, Pa, o):
    print(f"{tag:16s} {Pa:4.2f} {o['Rmax']*1e6:5.1f} {o['Rmin']*1e6:5.3f} {o['Mach']:5.2f} {o['vapour_frac']*100:5.1f} "
          f"{o['T_dyn']:6.0f} {o['Tmax']:6.0f} {o['x']:5.2f} {o['ne_frac']:.2e}   {o['e_from_HO']:5.2f}     {o['dE_eV']:4.2f}  "
          f"{o['photons']:.2e}  {o['Ppeak']*1e3:7.3f}  {o['tau']:.1e}  {o['fwhm_uv']*1e12:4.0f}/{o['fwhm_red']*1e12:4.0f}   {o['T_colour']:6.0f}")


if __name__ == "__main__":
    print(HDR); out = {}
    for tag, pv, chem in [("dry", 0.0, True), ("wet, inert H2O", PV_20C, False), ("wet, chemistry", PV_20C, True)]:
        for Pa in [1.25, 1.3, 1.35, 1.4, 1.45, 1.5]:
            o = analyse(run(Pa_atm=Pa, pv=pv), chemistry=chem); out[tag, Pa] = o; row(tag, Pa, o)
    print("\nbrightness wet(chemistry)/dry:", {Pa: round(out['wet, chemistry', Pa]['photons'] / out['dry', Pa]['photons'], 3) for Pa in [1.3, 1.4, 1.5]})
    print("brightness ratio per 0.1 atm (1.4/1.3): dry", round(out['dry', 1.4]['photons'] / out['dry', 1.3]['photons'], 1),
          " wet", round(out['wet, chemistry', 1.4]['photons'] / out['wet, chemistry', 1.3]['photons'], 1))
    print("\nwith hard-core suppression of dissociation (Toegel, Hilgenfeldt & Lohse 2002):")
    for Pa in [1.25, 1.3, 1.35, 1.4, 1.45, 1.5]:
        o = analyse(run(Pa_atm=Pa), excluded=True); row("wet, excl. vol.", Pa, o)
    print("  peak packing fraction at 1.4 atm:", round(o["packing"], 2))
    print("\nvapour pressure only (Pa = 1.4 atm; other water properties kept at 20 C):")
    for label, pv in [("0 C", 611.0), ("10 C", 1228.0), ("20 C", 2337.0), ("30 C", 4246.0)]:
        row("pv of " + label, 1.4, analyse(run(Pa_atm=1.4, pv=pv)))
    print("\nconvergence (wet, Pa = 1.4 atm):")
    for rt, npts in [(1e-8, 2001), (1e-10, 2001), (1e-9, 4001)]:
        row(f"rtol={rt:g} n={npts}", 1.4, analyse(run(Pa_atm=1.4, rtol=rt), npts=npts))
