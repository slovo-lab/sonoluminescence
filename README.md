# sonoluminescence

An open, minimal numerical model of single-bubble sonoluminescence (SBSL), checked against published measurements. The aim is to find where the standard picture breaks, and to test corrections one at a time.

*Русская версия: [README.ru.md](README.ru.md)*

**Authorship.** The code and text are written by Claude, an AI model made by Anthropic. Anton Krugovov owns and curates the project. Nothing here has been peer reviewed or tested in a laboratory.

## What the model contains

| Part | Treatment |
|---|---|
| Bubble wall | Keller–Miksis equation, van der Waals hard core |
| Gas heating | Uniform temperature, Péclet-dependent polytropic exponent plus conductive loss (after Hilgenfeldt, Grossmann & Lohse, 1999) |
| Ionisation | Saha equation for argon |
| Light | Kramers–Unsöld continuum absorption (free–free + bound–free, hydrogen-like), Kirchhoff's law for a homogeneous sphere |

Left out so far: water vapour, temperature gradients and shocks inside the bubble, electron–neutral bremsstrahlung, non-ideal plasma effects.

## Results (v0.1)

Argon bubble in water at 20 °C, ambient radius 4.5 µm, 26.5 kHz, drive 1.3–1.4 atm.

| Quantity | Model | Measurements |
|---|---|---|
| Maximum radius | 39–50 µm | up to about 50 µm |
| Minimum radius | 0.7 µm | below 1 µm |
| Peak temperature | 15,000–22,000 K | above 15,000 K (argon in sulfuric acid; Flannigan & Suslick 2005) |
| Flash width (FWHM) | 156–170 ps | 60 to over 250 ps (Gompf et al. 1997); under 40 to over 350 ps (Hiller et al. 1998) |
| Red vs UV width | identical | identical |
| Peak optical power | 0.015–1.3 mW | of order 10 mW |

All model numbers are in [results.txt](results.txt); `python reproduce.py` regenerates them. They do not change when the solver tolerance is tightened from 1e-8 to 1e-12 or the run is doubled in length.

### Why the flash width does not depend on colour

In the Kramers–Unsöld approximation the emission coefficient per unit frequency is proportional to exp(−E_ion / kT) and does not depend on photon energy. Every wavelength therefore switches on and off together. This is a property of the approximation, not an independent prediction.

## Where the model disagrees with experiment

1. **Brightness.** Too dim by one to three orders of magnitude, and it changes about 80-fold for a 0.1 atm change in drive. At the flash the gas number density is about 1e28 m⁻³, close to a liquid, while the Saha equation assumes a dilute gas. Lowering the ionisation potential by 3 eV raises the light output about 12-fold.
2. **Spectrum.** Measured spectra fit a blackbody closely (Vazquez et al. 2001), which implies an opaque emitter. The model bubble has an optical thickness of 0.001–0.08, so it is almost transparent.
3. **Colour temperature.** The model spectrum always looks like a 13,300 K blackbody between 200 and 800 nm, while the gas temperature ranges from 7,000 to 22,000 K. A temperature read from the spectrum need not be the gas temperature.
4. **Width versus drive.** At fixed ambient radius the flash gets slightly shorter as the drive rises; in experiments it gets longer. The model matches only if the bubble grows with drive: 93 ps at 3 µm, 340–400 ps at 7.5 µm.
5. **Validity.** At 1.4 atm the wall moves faster than the speed of sound in water (Mach 1.1), at the edge of what the wall equation can describe. Putterman et al. (Nature, 2001) criticised this class of models on exactly this ground.

## Open questions

- Does a non-ideal ionisation model close the brightness gap without an opaque core?
- How much does water vapour cool the collapse?
- Does electron–neutral bremsstrahlung change the spectrum shape?
- Testable: with the ambient radius measured independently, does the flash width still depend on drive?

## Run

```
pip install -r requirements.txt
python reproduce.py
```

## References

- Gompf et al., Phys. Rev. Lett. 79, 1405 (1997)
- Hiller, Putterman & Weninger, Phys. Rev. Lett. 80, 1090 (1998)
- Hilgenfeldt, Grossmann & Lohse, Nature 398, 402 (1999); Phys. Fluids 11, 1318 (1999)
- Putterman et al., Nature 409, 782 (2001)
- Vazquez, Camara, Putterman & Weninger, arXiv:physics/0009057
- Flannigan & Suslick, Nature 434, 52 (2005)
- Tse & Leung, arXiv:cond-mat/0606359
- Zel'dovich & Raizer, *Physics of Shock Waves and High-Temperature Hydrodynamic Phenomena*

## Licence

MIT
