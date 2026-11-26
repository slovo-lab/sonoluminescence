"""Issue #9: trapped vapour share of v0.4 against published simulations that resolve the
bubble interior.

Reference numbers as collected in the comparison tables of arXiv:physics/0501041
(Sec. IV.A, Table II), which also quotes Storey & Szeri (2000):
  case A  argon, R0 = 4.5 um, Pa = 1.2 bar, 26.5 kHz, water at 298 K
          vapour share  : 14 % (Storey & Szeri, moment not specified); 7.7 % at Rmin (0501041)
          Rmax / Rmin   : 31.3 / 0.70 um wet, 28.0 / 0.80 um dry (Storey & Szeri)
          peak T, dry / inert vapour / chemistry :
                          20,900 / 9,700 / 7,000 K (Storey & Szeri)
                          17,900 / 17,000 / 8,900 K (0501041)
  case B  argon, R0 = 6.0 um, Pa = 1.4 bar, 20.6 kHz, water at 293.15 K
          vapour share  : 33 % (Storey & Szeri); 23 % at Rmin (0501041)
The published temperatures are peak values inside the bubble; this model has one
uniform temperature.
"""
import numpy as np
from vapor import run, analyse

BAR = 1e5 / 101325.0
CASES = {
    "A": dict(R0=4.5e-6, Pa_atm=1.2 * BAR, f=26.5e3, Tw=298.0, pv=3169.0, mu_w=0.89e-3, sigma_w=0.0720),
    "B": dict(R0=6.0e-6, Pa_atm=1.4 * BAR, f=20.6e3, Tw=293.15, pv=2337.0, mu_w=1.0e-3, sigma_w=0.0725),
}


def line(tag, o):
    print(f"{tag:26s} Rmax {o['Rmax']*1e6:5.1f}  Rmin {o['Rmin']*1e6:5.3f}  vapour {o['vapour_frac']*100:5.1f} %  "
          f"T inert {o['T_dyn']:6.0f} K  T chem {o['Tmax']:6.0f} K  diss {o['x']:.2f}")


for name, kw in CASES.items():
    print(f"== case {name}")
    dry = analyse(run(**{**kw, "pv": 0.0}))
    print(f"{'dry':26s} Rmax {dry['Rmax']*1e6:5.1f}  Rmin {dry['Rmin']*1e6:5.3f}  T before/after ionisation budget "
          f"{dry['T_dyn']:6.0f} / {dry['Tmax']:6.0f} K")
    line("wet", analyse(run(**kw)))
    for sc in [0.25, 0.5, 2.0, 4.0]:
        line(f"wet, diffusivity x{sc}", analyse(run(**kw, d_scale=sc)))

print("== case A with Pa read as 1.2 atm instead of 1.2 bar")
kw = {**CASES["A"], "Pa_atm": 1.2}
dry = analyse(run(**{**kw, "pv": 0.0}))
print(f"{'dry':26s} Rmax {dry['Rmax']*1e6:5.1f}  Rmin {dry['Rmin']*1e6:5.3f}  T before/after ionisation budget "
      f"{dry['T_dyn']:6.0f} / {dry['Tmax']:6.0f} K")
line("wet", analyse(run(**kw)))

print("== v0.4 standard case (R0 = 4.5 um, 26.5 kHz, 20 C) with doubled diffusivity")
for Pa in [1.3, 1.35, 1.4, 1.5]:
    for sc in [1.0, 2.0]:
        o = analyse(run(Pa_atm=Pa, d_scale=sc))
        print(f"Pa {Pa:4.2f} atm, D x{sc:3.1f}:  vapour {o['vapour_frac']*100:5.1f} %  T {o['Tmax']:6.0f} K  "
              f"photons {o['photons']:.2e}  peak {o['Ppeak']*1e3:.4f} mW  FWHM uv {o['fwhm_uv']*1e12:4.0f} ps")
