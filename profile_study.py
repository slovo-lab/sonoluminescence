"""v0.6 (preliminary): radially resolved interior.  Benchmark case A and the standard case."""
import numpy as np
import profile as pr

BAR = 1e5 / 101325.0
A = dict(R0=4.5e-6, Pa_atm=1.2 * BAR, f=26.5e3, Tw=298.0, mu_w=0.89e-3, sigma_w=0.072)


def line(tag, o):
    print(f"{tag:22s} Rmax {o['Rmax']*1e6:5.1f} Rmin {o['Rmin']*1e6:5.3f} | vapour total {o['vapour_frac']*100:4.1f} % centre {o['x_centre']*100:4.1f} % | "
          f"T centre inert {o['Tc_dyn']:6.0f} K, mass-avg {o['Tavg_dyn']:6.0f} K, centre after chemistry {o['Tc']:6.0f} K | "
          f"photons {o['photons']:.2e} peak {o['Ppeak']*1e3:.4f} mW FWHM uv/red {o['fwhm_uv']*1e12:3.0f}/{o['fwhm_red']*1e12:3.0f} ps | argon drift {o['argon_drift']*100:.1f} %")


if __name__ == "__main__":
    print("== benchmark case A (references: centre T dry 20,900 / 17,900 K; inert vapour 9,700 / 17,000 K; chemistry 7,000 / 8,900 K; vapour 14 % / 7.7 %)")
    line("A dry", pr.analyse(pr.run(**A, pv=0.0)))
    line("A wet", pr.analyse(pr.run(**A, pv=3169.0)))
    print("== standard case (R0 = 4.5 um, 26.5 kHz, water 20 C)")
    for Pa in [1.3, 1.4]:
        line(f"wet Pa={Pa}", pr.analyse(pr.run(Pa_atm=Pa)))
