"""v0.12: radiative attachment of electrons to H and O atoms (issue #10).

Water vapour breaks up into H and O atoms; a free electron that attaches to one emits a
photon.  In equilibrium this is the absorption by H- and O- ions (Kirchhoff's law), the
process that makes the solar surface opaque.  Same runs, light counted with and without it.
"""
import numpy as np
import profile as pr
import budget

BAR = 1 / 1.01325


def line(tag, a, b):
    print(f"{tag:28s} photons {a['photons']:.3g} -> {b['photons']:.3g} (x{b['photons']/a['photons']:.2f}) | peak {a['Ppeak']*1e3:.3g} -> {b['Ppeak']*1e3:.3g} mW | "
          f"tau {a['tau']:.2g} -> {b['tau']:.2g} | FWHM uv/red {b['fwhm_uv']*1e12:.0f}/{b['fwhm_red']*1e12:.0f} ps | T_colour {a['T_colour']:.0f} -> {b['T_colour']:.0f} K")


if __name__ == "__main__":
    print("== water, R0 = 4.5 um, 26.5 kHz, 20 C")
    for Pa in (1.3, 1.35, 1.4):
        r = pr.run(Pa_atm=Pa)
        line(f"water {Pa} atm", pr.analyse(r), pr.analyse(r, attachment=True))
    r = pr.run(Pa_atm=1.35)
    for sc in (0.5, 2.0):
        line(f"water 1.35 atm, O- x{sc}", pr.analyse(r), pr.analyse(r, attachment=True, o_scale=sc))
    print("== water at 0 C vapour pressure (611 Pa), 1.4 atm")
    r = pr.run(Pa_atm=1.4, pv=611.0)
    line("cold water 1.4 atm", pr.analyse(r), pr.analyse(r, attachment=True))
    print("== sulfuric acid (almost no vapour), dim bubble")
    kw = dict(budget.CASES["acid"]["run"]); ana = budget.CASES["acid"]["ana"]
    r = pr.run(**kw)
    line("acid dim", pr.analyse(r, **ana), pr.analyse(r, attachment=True, **ana))
