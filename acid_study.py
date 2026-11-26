"""Test without tuning: argon bubbles in 85 % sulfuric acid.

Measurements (Flannigan, Hopkins, Camara, Putterman & Suslick, PRL 96, 204301 (2006)):
  dim bubble   : R0 = 13 um, effective Pa = 1.4 bar; collapse radius 3.1 um;
                 T = 10,000 K (Ar lines); n_Ar = 1.0-1.2e21 cm^-3; pressure 1100-1600 bar;
                 n_e ~ 3.8e17 cm^-3
  bright bubble: R0 = 17 um, effective Pa = 1.7 bar; collapse radius 3.5 um;
                 n_Ar = 1.8e21 cm^-3; pressure 3700 bar; wall speed up to 400 m/s;
                 about 1e4 times brighter than the dim bubble; flash 5-10 ns
Liquid (Hopkins et al., PRL 95, 254301 (2005)): density 1800 kg/m^3, sound speed 1470 m/s,
  viscosity 25 mPa s, surface tension 55 mN/m, vapour pressure 0.04 torr.
The acoustic frequency is not given in the 2006 paper; 30 kHz is assumed (37.8 kHz, used
in the 2005 paper, is run as a check).
"""
import numpy as np
import profile as pr

SA = dict(rho_l=1800.0, c_l=1470.0, mu_w=25e-3, sigma_w=0.055, pv=5.3, Tw=298.0)
BAR = 1 / 1.01325


def line(tag, o):
    print(f"{tag:30s} Rmax/R0 {o['Rmax']/o['R0']:4.1f}  Rmin {o['Rmin']*1e6:4.2f} um  wall {o['Umax']:4.0f} m/s  p {o['p_max']/1e5:5.0f} bar  "
          f"n centre {o['n_centre']/1e27:4.2f}e21 cm^-3  T centre {o['Tc']:6.0f} K (mass-avg {o['Tavg_dyn']:5.0f})  "
          f"n_e {o['ne_centre']/1e6:.1e} cm^-3  photons {o['photons']:.2e}  FWHM {o['fwhm_uv']*1e9:5.2f} ns  T_colour {o['T_colour']:5.0f} K")


def case(tag, R0, Pa_bar, f=30e3, **over):
    kw = {**SA, **over}
    o = pr.analyse(pr.run(R0=R0, Pa_atm=Pa_bar * BAR, f=f, **kw), half=25e-9, npts=1001)
    o["R0"] = R0
    line(tag, o)
    return o


if __name__ == "__main__":
    dim = case("dim    (13 um, 1.4 bar)", 13e-6, 1.4)
    bright = case("bright (17 um, 1.7 bar)", 17e-6, 1.7)
    print(f"brightness bright/dim: {bright['photons']/dim['photons']:.0f}x")
    print("-- sensitivity of the dim bubble")
    case("dim, 37.8 kHz", 13e-6, 1.4, f=37.8e3)
    case("dim, viscosity 10 mPa s", 13e-6, 1.4, mu_w=10e-3)
    case("dim, Pa +10 %", 13e-6, 1.54)
    case("dim, Pa -10 %", 13e-6, 1.26)
    print("-- water reference from the same model (R0 = 4.5 um, 26.5 kHz, 20 C)")
    for Pa in [1.3, 1.4]:
        o = pr.analyse(pr.run(Pa_atm=Pa)); o["R0"] = 4.5e-6
        line(f"water, {Pa} atm", o)
        print(f"   bright acid / water: {bright['photons']/o['photons']:.0f}x")
    print("-- what if the acid had the vapour pressure of water at 20 C?")
    wet = case("bright, pv = 2337 Pa", 17e-6, 1.7, pv=2337.0)
    print(f"   photons relative to the real acid: {wet['photons']/bright['photons']:.3f}; vapour share at collapse {wet['vapour_frac']*100:.1f} %")
    print("-- size and duration bookkeeping, bright acid vs water at 1.4 atm")
    print(f"   atoms ratio (R0^3 at equal ambient pressure): {(17/4.5)**3:.0f}; flash duration ratio: {bright['fwhm_uv']/o['fwhm_uv']:.0f}")
