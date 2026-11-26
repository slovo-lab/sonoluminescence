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

Left out of v0.1: water vapour, temperature gradients and shocks inside the bubble, electron–neutral bremsstrahlung, non-ideal plasma effects. The last two are added in v0.2 and v0.3 below.

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

## v0.2: ionisation in a dense gas (issue #1)

Question: does a density-dependent lowering of the ionisation potential close the brightness gap?

Two effects were added, both applied to the Saha equation and to the absorption coefficient:

- **Neutral gas.** Polarisation of the surrounding atoms, scaled linearly in density from the measured liquid-argon value (15.76 → 13.84 eV at 2.11e28 m⁻³).
- **Free charges.** Stewart–Pyatt screening, solved self-consistently with the ionisation degree.

A third change turned out to matter as much: taking the energy spent on ionisation out of the heat budget (post-hoc, not yet fed back into the wall dynamics).

| At 1.4 atm | v0.1 | + energy budget | + energy budget + lowering |
|---|---|---|---|
| Peak temperature | 22,200 K | 19,000 K | 17,600 K |
| Lowering | 0 | 0 | 3.1 eV |
| Ionisation degree | 5.0% | 2.3% | 4.0% |
| Peak optical power | 1.27 mW | 0.29 mW | 0.91 mW |
| Optical thickness | 0.08 | 0.02 | 0.08 |
| Flash width | 156 ps | 188 ps | 184 ps |

**Answer: no.** The lowering is 2–3 eV and raises the light output about 3-fold, but v0.1 overstated the brightness about 4-fold by ignoring the energy cost of ionisation. Net, the model stays about 10 times too dim at the strongest drive and more at weaker drive, and the bubble stays transparent. No opaque core forms, and the self-consistent solution is unique (no abrupt transition) between 10,000 and 25,000 K.

Side results:

- The energy budget halves the sensitivity to drive (1.4 atm / 1.3 atm brightness ratio falls from 83 to 31) and brings the peak temperature closer to the measured 15,000 K.
- The flash width becomes nearly independent of drive (180–186 ps for 1.3–1.4 atm).

Caveats: the energy correction is post-hoc; the neutral-gas term is a linear interpolation; Stewart–Pyatt is used outside the fully ionised plasmas it was built for. Kappus, Bataller & Putterman (PRL 111, 234301, 2013) argue for a much larger lowering, which this model does not reproduce.

All numbers: [results_ipd.txt](results_ipd.txt), regenerated by `python ipd_study.py`.

## v0.3: electron–neutral bremsstrahlung (issue #3)

At ionisation degrees of a few per cent there are 25–100 neutral atoms per ion, so electrons radiate mostly while scattering on neutrals. v0.3 adds this channel (Dalgarno–Lane cross-section, Maxwell average, Kirchhoff's law) on top of v0.2. In the long-wavelength limit the code reproduces the classical Drude absorption coefficient to within 0.1%.

| Drive | Peak power v0.2 | Peak power v0.3 | Optical thickness v0.3 | Spectrum colour v0.3 | Flash UV / red v0.3 |
|---|---|---|---|---|---|
| 1.30 atm | 0.03 mW | 0.13 mW | 0.04 | 8,350 K | 194 / 207 ps |
| 1.35 atm | 0.25 mW | 0.75 mW | 0.15 | 9,100 K | 196 / 210 ps |
| 1.40 atm | 0.91 mW | 2.27 mW | 0.39 | 9,750 K | 202 / 225 ps |

What changed:

- **Brightness.** Up 3-fold at 1.4 atm and 6-fold at 1.3 atm. The gap to the measured order of 10 mW shrinks to about 4-fold at the strongest drive. At the flash peak neutrals supply 55–67% of the absorption in the UV and 78–89% in the red.
- **Drive sensitivity.** The 1.4 atm / 1.3 atm brightness ratio falls from 31 to 17.
- **Opacity.** The optical thickness reaches 0.39, so the bubble is no longer transparent at strong drive.
- **Spectrum.** The fixed 13,300 K look of v0.1 is gone; the colour temperature now moves from 8,350 to 9,750 K with drive.

New disagreement: the red flash is now 7–11% longer than the UV flash, while experiments report identical widths. Whether 15–20 ps is within the experimental resolution needs checking.

Caveats: the argon cross-section table is approximate (scaling it by ±30% changes brightness by about ±22%); binary electron–atom collisions are a rough picture at near-liquid density.

All numbers: [results_en.txt](results_en.txt), regenerated by `python en_study.py`.

## v0.4: water vapour (issue #2)

The wall dynamics now exchange heat and vapour with the liquid through diffusive boundary layers (Toegel, Gompf, Pecha & Lohse 2000). Vapour is carried as inert H₂O during the collapse; afterwards, at fixed volume and energy, the code solves the equilibrium H₂O ⇌ 2H + O (from partition functions) together with the ionisation of Ar, H and O. Light is computed as in v0.3. Code: `vapor.py`, numbers: [results_vapor.txt](results_vapor.txt).

| Drive | Vapour share at collapse | Peak temperature, dry | Peak temperature, wet | Photons per flash, dry | Photons per flash, wet |
|---|---|---|---|---|---|
| 1.30 atm | 22% | 17,600 K | 8,000 K | 1.4 million | 59 |
| 1.35 atm | 29% | 19,700 K | 8,200 K | 4.0 million | 99 |
| 1.40 atm | 37% | 21,300 K | 8,100 K | 7.9 million | 86 |
| 1.50 atm | 57% | 24,100 K | 7,200 K | 19 million | 15 |

**Answer: vapour cools the collapse by more than half and makes the model tens of thousands of times too dim at 20 °C.** The 4-fold gap of v0.3 was that small only because the bubble was dry. This is a known problem, not a new one: Toegel et al. (2000) asked whether vapour prevents upscaling of sonoluminescence, and Storey & Szeri (2002) found temperatures of about 7,000 K, which they called too low for bremsstrahlung.

What the runs show:

- **Heat capacity matters more than chemistry.** With inert vapour the gas reaches 10,400–11,700 K and gives 2,000–9,000 photons. Dissociation removes another 2,300–3,600 K.
- **The proposed fix is not enough here.** Hard-core suppression of dissociation (Toegel, Hilgenfeldt & Lohse 2002) raises the temperature to about 9,100 K and the light 4-fold, to about 400 photons.
- **Stronger drive stops helping.** Brightness peaks near 1.35 atm and then falls, because a harder-driven bubble traps more vapour.
- **Cold water rescues the model.** With the vapour pressure of 0 °C water (all else unchanged) the 1.4 atm bubble holds 7% vapour, reaches 18,100 K and emits 2.3 million photons at 3.9 mW peak, the measured order of magnitude. Experiments do show much brighter flashes in cold water, but the model's 26,000-fold change between 0 and 20 °C is far too steep.
- **Flash width.** 123–181 ps wet (1.3–1.4 atm), inside the measured range; 384–451 ps dry, outside it.
- **Electrons.** In the wet bubble 37% of the free electrons come from hydrogen and oxygen atoms.

A second finding concerns the dry model itself: the boundary-layer heat model gives 29,600 K at 1.35 atm before the ionisation budget, where the v0.1 heat model gave 18,400 K. The choice of heat-loss model changes the answer as much as any physics added since.

Caveats: chemistry and ionisation are post-hoc (issue #6); only atoms are produced (no OH, H₂, O₂), which understates cooling at intermediate temperatures; the argon cross-section stands in for all neutrals; no light from water or its fragments.

## v0.5: is the trapped vapour share right? (issue #9)

The v0.4 model was run at two cases for which simulations resolving the bubble interior are published. Reference numbers are taken from the comparison tables in arXiv:physics/0501041, which also quotes Storey & Szeri (2000); they have not yet been checked against the original papers. Script: `vapour_benchmark.py`, numbers: [results_benchmark.txt](results_benchmark.txt).

| | This model | Storey & Szeri | arXiv:physics/0501041 |
|---|---|---|---|
| **Case A** (4.5 µm, 1.2 bar, 26.5 kHz, 298 K): vapour share | 14.9% | 14% | 7.7% |
| Case A: temperature, inert vapour | 8,000 K | 9,700 K | 17,000 K |
| Case A: temperature, with chemistry | 6,100 K | 7,000 K | 8,900 K |
| Case A: temperature, dry | 9,100 K | 20,900 K | 17,900 K |
| **Case B** (6.0 µm, 1.4 bar, 20.6 kHz, 293 K): vapour share | 50% | 33% | 23% |

Reading:

- **Moderate drive: agreement.** The vapour share matches Storey & Szeri within one percentage point, and the wet temperatures are 13–17% below theirs.
- **Strong drive: too much vapour.** The boundary-layer model traps 1.5 to 2.2 times more than the resolved simulations. Doubling the vapour diffusivity brings case B to 23%.
- **The references disagree with each other** by about a factor of two on the vapour share, and by 7,300 K on the inert-vapour temperature.
- **Dry bubble.** The uniform model is 2 times colder than the resolved peak temperatures. Published values are peaks inside the bubble; a uniform model can only give an average.

Effect on v0.4, standard case with doubled diffusivity (vapour share 11–17% instead of 22–37% at 1.3–1.4 atm): temperature about 11,400–12,000 K, 7,500–22,000 photons per flash, peak power 0.016–0.068 mW. The gap to the measured order of 10 mW shrinks from tens of thousands to roughly 150–600 times. The conclusion of v0.4 stands in weaker form: a uniform bubble with a realistic amount of vapour is too dim at room temperature.

## v0.6: temperature and vapour resolved along the radius

`profile.py` drops the uniform-bubble assumption: pressure stays uniform, while temperature and vapour share vary from centre to wall (heat conduction, vapour diffusion, Stefan flow at the wall). Chemistry and ionisation are post-hoc per shell; light is traced along rays through the shells. Numbers: [results_profile.txt](results_profile.txt).

Numerical quality: argon and vapour are conserved (see v0.10); successive acoustic cycles repeat to four digits; going from 40 to 80 nodes changes the light output by 2.6% and the centre temperature by 0.4%.

- Benchmark case A: vapour 5.5% in total and 14% at the centre; centre temperature 7,400 K with chemistry. The references give 7.7% and 14% vapour, 7,000 and 8,900 K.
- Standard case, 1.3–1.4 atm: vapour 7–11% in total but 19–28% at the centre; centre temperature 10,000–10,300 K with chemistry; 8,000–48,000 photons; peak 0.023–0.18 mW.
- Vapour collects in the hot centre, so the centre is not the dry, bright core one might hope for. Against a peak power of order 10 mW the model is 57–440 times too dim.

## v0.7: a test without tuning — argon bubbles in sulfuric acid

The v0.6 code was run for 85% sulfuric acid with only the liquid properties changed (density, sound speed, viscosity, surface tension, vapour pressure 0.04 torr) and the bubble size and drive taken from the experiment. Nothing was fitted. Script: `acid_study.py`, numbers: [results_acid.txt](results_acid.txt). Measurements: Flannigan, Hopkins, Camara, Putterman & Suslick (2006); liquid data and brightness ratios: Hopkins et al. (2005).

| | Model | Measured |
|---|---|---|
| **Dim bubble** (R0 13 µm, 1.4 bar): collapse radius | 2.95 µm | 3.1 µm |
| Gas density at the centre | 1.18 × 10²¹ cm⁻³ | 1.0–1.2 × 10²¹ cm⁻³ |
| Temperature at the centre | 12,100 K | 10,000 K |
| Pressure | 2,170 bar | 1,100–1,600 bar |
| Electron density | 4.4 × 10¹⁸ cm⁻³ | 3.8 × 10¹⁷ cm⁻³ |
| **Bright bubble** (R0 17 µm, 1.7 bar): collapse radius | 3.51 µm | 3.5 µm |
| Gas density at the centre | 1.56 × 10²¹ cm⁻³ | 1.8 × 10²¹ cm⁻³ |
| Pressure | 4,980 bar | 3,700 bar |
| Flash duration | 3.8 ns | 5–10 ns |
| Brightness, bright / dim | 330 times | about 10,000 times |
| Brightness, bright acid / water | 380–2,300 times | 1,500–2,700 times |

Reading:

- **Mechanics and thermodynamics agree.** Collapse radius within 5%, gas density within 13%, temperature within 21%, pressure within a factor of 1.3 to 2. The model also gives nanosecond flashes in acid against 0.1 ns in water, as observed.
- **The acid/water brightness ratio comes out at the measured order.** Bookkeeping: the acid bubble holds 54 times more atoms and shines 38 times longer.
- **Vapour is still decisive.** The same acid bubble given the vapour pressure of water traps 7.8% vapour, cools from 16,500 K to 6,600 K and emits 98,000 times less light.
- **Brightness and electron density are not predicted, only bracketed** (see v0.8). The electron density is 12 times too high at 30 kHz and 0.6 times the measured value at 37.8 kHz.

Caveats: the acoustic frequency is not given in the 2006 paper and 30 kHz is assumed; at 37.8 kHz the dim bubble gives 9,100 K, 1,180 bar and 0.90 × 10²¹ cm⁻³. Lowering the viscosity from 25 to 10 mPa·s raises the light 67-fold.

## v0.8: uncertainty budget

`budget.py` computes, for each input, the sensitivity coefficient c = d ln(output) / d ln(input) by central differences over ±5%, multiplies it by an assumed relative standard uncertainty of that input, and combines the contributions in quadrature. Two reference cases: water at 20 °C (R0 4.5 µm, 1.35 atm, 26.5 kHz) and the dim argon bubble in sulfuric acid. Numbers: [results_budget.txt](results_budget.txt); the assumed uncertainties and their reasons are listed in the script.

| Input | Assumed uncertainty | Water: c for brightness | Water: share of variance | Acid: c for brightness | Acid: share of variance |
|---|---|---|---|---|---|
| Drive pressure | 4% / 10% | 25 | 6% | 68 | 92% |
| Vapour diffusivity | 40% | 9.6 | 84% | 0.0 | 0% |
| Gas thermal conductivity | 20% | −6.1 | 9% | −1.6 | 0% |
| Viscosity | 2.5% / 20% | −0.2 | 0% | −8.5 | 6% |
| Frequency | 0.1% / 6.6% | −4.6 | 0% | −15 | 2% |
| Vapour pressure | 6% / 50% | −7.6 | 1% | 0.0 | 0% |
| Ambient radius | 5% | −2.1 | 0% | 1.5 | 0% |
| e–neutral cross-section | 30% | 0.7 | 0% | 0.9 | 0% |

| | Water | Acid |
|---|---|---|
| Brightness, combined uncertainty (1σ, linearised) | factor 65 | factor 1,200 |
| Centre temperature, combined uncertainty | 24% | 45% |
| Same without the drive-pressure term | factor 58 | factor 7.5 |
| Drive-pressure accuracy needed to keep its term below a factor of 2 | 2.8% | 1.0% |

Reading:

- **Brightness scales as the 25th power of the drive pressure in water and the 68th power in acid.** In acid, with the 10% tolerance quoted for experimental fits, the pressure term alone makes brightness uncertain by a factor of several hundred. In water, with a hydrophone precision of 4%, the pressure term is a factor of 2.7.
- **In water the largest term is a model parameter**, the vapour diffusivity, for which published simulations differ by a factor of two; it carries 84% of the variance. In acid the vapour terms vanish and the drive pressure carries 92%.
- **What would make brightness a test:** the pressure at the bubble to 1–3%, and for water an independent fix of the trapped vapour.
- **Robust outputs.** The collapse radius has coefficients near 1 for every input, which is why the mechanics agree with the measurements at the 5% level.
- The response is not symmetric: in water −10% in drive gives 170 times less light, +10% only 4 times more, because a harder-driven bubble traps more vapour. The linearised figure overstates the upward range.
- Model switches at the baseline: without chemistry the water bubble is 15 times brighter; without ionisation-potential lowering it is 2.6 times dimmer.

Two errors in the first version of `profile.py` were found while preparing this budget and are fixed: the argon count grew by about 4% per cycle, and the search for the moment of collapse could miss it by a fraction of a nanosecond, which made the light output jump by a factor of 10 between otherwise identical runs. The numbers in v0.6 and v0.7 above are from the corrected code.

## v0.10: mass conservation closed

The earlier fix (v0.8) made the argon count exact by computing the pressure from it. That hid the remaining error instead of removing it: vapour was carried as a mole fraction at the grid nodes, in non-conservative form, so the inconsistency moved into the vapour. Over one cycle the vapour content of the bubble differed from the vapour that crossed the wall by 3.3%, 3.9% and 7.3% on grids of 20, 40 and 80 nodes. Growing with grid refinement, it was an error in the formulation, not in the grid.

`profile.py` now carries the number of vapour molecules in each cell and changes it only by fluxes through the cell faces (bulk flow relative to the moving face plus diffusion). `conservation_check.py` tests four things, sampling the collapse every 0.5 ps:

| Check | 20 nodes | 40 nodes | 80 nodes |
|---|---|---|---|
| Argon in the whole bubble | 2.5 × 10⁻¹⁴ | 1.8 × 10⁻¹⁴ | 1.8 × 10⁻¹⁴ |
| Vapour change against wall flux | 0.15% | 0.82% | 1.91% |
| Argon in inner regions against flux through their boundary | 15.7% | 6.6% | 3.3% |
| Pressure from conservation against the rate in the wall equation | 0.22% | 0.17% | 0.12% |

Reading:

- **Vapour.** The cell changes add up exactly to the flux through the wall; what remains is the time-integration tolerance. At 40 nodes the residual is 0.82% with the default tolerance, 0.077% at ten times tighter and 0.005% at a hundred times tighter.
- **Argon, locally.** The error halves with every doubling of the grid: a first-order discretisation error, not a leak. It comes from moving temperature along the nodes while the fluxes are computed at the faces.
- **Pressure.** The value fixed by conservation and the rate used to move the wall agree to 0.2%.

Effect on results: photon counts change by 1–2% and temperatures by under 0.3% in every reference case. No conclusion above changes. The tables in v0.6–v0.8 are recomputed with the corrected code.

## v0.9: literature check

What was already known:

- **Sulfuric acid.** Moshaii, Imani & Silatani (Phys. Rev. E 80, 046325, 2009) simulated an argon bubble in water and in several concentrations of sulfuric acid to quantify the effects of vapour pressure and viscosity. Kim, Byun & Kwak report calculated wall velocity, peak temperature and pressure in agreement with the acid experiments. The comparison in v0.7 therefore confirms published work and is not new.
- **Sensitivity.** The review by Brenner, Hilgenfeldt & Lohse (Rev. Mod. Phys. 74, 425, 2002) already states that neither the ambient radius nor the drive pressure is known a priori, that fitted values depend on model details, that the hydrophone precision is about 0.05 bar, and that errors in the interior model and in the emission model may compensate each other. The qualitative point of v0.8 is in that review.

What the searches did not find: a budget that assigns each input a sensitivity coefficient and an uncertainty and adds them up, as in v0.8. A few searches are not proof that none exists.

One correction followed from this check: v0.8 first assumed 10% for the drive pressure in water; the review's 0.05 bar is 4%, and the tables above now use it. For the acid the 10% tolerance of Hopkins et al. is kept.

## v0.11: Monte Carlo uncertainty budget

`mc_budget.py` varies all twelve inputs of v0.8 at once, with the same assumed uncertainties, and runs the full model for each draw (Latin-hypercube sampling, 128 draws per case, one acoustic cycle per run). Numbers: [results_mc.txt](results_mc.txt).

| | Water (R0 4.5 µm, 1.35 atm) | Acid (dim bubble) |
|---|---|---|
| Photons per flash, median | 20,000 | 8,500 |
| 68% interval | 320 to 480,000 | 0.2 to 1.5 million |
| 95% interval | 23 to 1.7 million | 0 to 15 million |
| 68% interval as a factor below / above the median | 62 / 24 | 35,000 / 180 |
| Linearised budget (v0.8), each way | 65 | 1,200 |
| Mean as a multiple of the median | 13 | 210 |
| Draws with under 1% of the baseline light | 12% | 38% |
| Centre temperature, 95% interval | 7,600 to 18,900 K | 2,900 to 17,600 K |
| Input most correlated with brightness (rank correlation) | vapour diffusivity, +0.90 | drive pressure, +0.94 |

Reading:

- **Water.** The linearised budget was about right on the low side and too wide on the high side: the light cannot rise far, because a harder collapse traps more vapour. How fast vapour escapes the bubble decides the answer, as in v0.8.
- **Acid.** The linearised budget was wrong in kind, not in degree. With the 10% drive-pressure tolerance quoted for the experimental fits, the model cannot say whether the dim bubble glows at all: in 38% of the draws it gives under 1% of the baseline light, and the 95% interval starts at zero. The distribution has a threshold, not a spread.
- **Mean and median differ by a factor of 13 and 210.** An average over draws, or over flashes of a jittering bubble, is dominated by the rare hot collapses. Comparing a model's typical flash with a measured average brightness is biased upward by that factor.
- The rank correlations of frequency in water (−0.25) and of sound speed in acid (+0.16) are probably chance: their assumed uncertainties are 0.1% and 0.5% and their linear effect is below 1%. With twelve inputs and 128 draws, one or two such values are expected.

Caveats: one acoustic cycle per run differs from two cycles by 0.4% in water but 19% in acid (baseline 66,000 against 55,500 photons), which is small next to the spread. 128 draws fix the 95% bounds to within a factor of about 2 (bootstrap in the results file).

## v0.12: radiative attachment of electrons to H and O (issue #10)

When water breaks up, a free electron can attach to a hydrogen or oxygen atom and emit a photon; in equilibrium this is the absorption by H⁻ and O⁻ ions, the process that makes the surface of the Sun opaque. This is not new for sonoluminescence: Hammer & Frommhold (2002) included it, and An (2006) found it dominant for argon bubbles in cold water. It was missing here. `attach_study.py`, numbers: [results_attach.txt](results_attach.txt). The H⁻ cross-section is the Wishart (1979) fit; the O⁻ cross-section is an approximate plateau, varied ×0.5 and ×2.

| Case | Photons without | Photons with | Factor |
|---|---|---|---|
| Water 20 °C, 1.30 atm | 8,000 | 29,000 | 3.7 |
| Water 20 °C, 1.35 atm | 24,000 | 73,000 (67,000–85,000 for O⁻ ×0.5–×2) | 3.0 |
| Water 20 °C, 1.40 atm | 48,000 | 131,000 | 2.8 |
| Water, vapour pressure of 0 °C, 1.40 atm | 3.1 million | 3.6 million | 1.15 |
| Sulfuric acid, dim bubble | 55,500 | 56,200 | 1.01 |

- In water at 20 °C the light rises about threefold; peak power reaches 0.07–0.43 mW, so the gap to the measured order of 10 mW shrinks from 57–440 to about 23–130 times. It does not close.
- In acid the effect is nil, as it should be with almost no vapour. In cold water, where the bubble is hotter and the fragments are ionised, it adds 15%.
- The flash becomes slightly longer and its colour temperature rises by 400–600 K.

## Open questions

- ~~Does a non-ideal ionisation model close the brightness gap without an opaque core?~~ No, see v0.2.
- Feed the ionisation energy into the wall dynamics instead of correcting afterwards.
- ~~How much does water vapour cool the collapse?~~ By more than half, see v0.4.
- ~~Is the trapped vapour share (22–37%) right?~~ Right at moderate drive, 1.5–2.2 times too high at strong drive, see v0.5.
- Check the benchmark numbers against the original papers.
- ~~Does a hot centre (non-uniform temperature) supply the missing light?~~ No: vapour collects at the centre, see v0.6.
- What acoustic frequency was used for the argon bubbles in sulfuric acid?
- Add xenon, to use the xenon measurements in acid and water.
- ~~Replace the linearised budget by a Monte Carlo one.~~ Done, see v0.11.
- Add OH band emission (the remaining part of issue #10).
- ~~Has an uncertainty budget of this kind been published before?~~ Not found in a first search; the qualitative point is in the 2002 review, see v0.9.
- Do water and its fragments (OH, H, O) emit enough light to fill the gap?
- Which heat-loss model is right for the dry bubble?
- ~~Does electron–neutral bremsstrahlung change the spectrum shape?~~ Yes, see v0.3.
- Is a 7–11% longer red flash compatible with the measurements?
- Replace the approximate argon cross-section table with an LXCat data set.
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
- Dalgarno & Lane, Astrophys. J. 145, 623 (1966)
- Milloy, Crompton, Rees & Robertson, Aust. J. Phys. 30, 61 (1977)
- Kappus, Bataller & Putterman, Phys. Rev. Lett. 111, 234301 (2013)
- Stewart & Pyatt, Astrophys. J. 144, 1203 (1966)
- Toegel, Gompf, Pecha & Lohse, Phys. Rev. Lett. 85, 3165 (2000)
- Toegel, Hilgenfeldt & Lohse, Phys. Rev. Lett. 88, 034301 (2002)
- Hopkins, Putterman, Kappus, Suslick & Camara, Phys. Rev. Lett. 95, 254301 (2005)
- Flannigan, Hopkins, Camara, Putterman & Suslick, Phys. Rev. Lett. 96, 204301 (2006)
- Storey & Szeri, Proc. R. Soc. Lond. A 456, 1685 (2000)
- Storey & Szeri, Phys. Rev. Lett. 88, 074301 (2002)
- arXiv:physics/0501041
- Brenner, Hilgenfeldt & Lohse, Rev. Mod. Phys. 74, 425 (2002)
- Moshaii, Imani & Silatani, Phys. Rev. E 80, 046325 (2009)
- Hammer & Frommhold, Phys. Rev. E 65, 046309 (2002)
- An, Phys. Rev. E 74, 026304 (2006)
- Wishart, Mon. Not. R. Astron. Soc. 187, 59P (1979)
- JCGM 101:2008, Propagation of distributions using a Monte Carlo method
- Tse & Leung, arXiv:cond-mat/0606359
- Zel'dovich & Raizer, *Physics of Shock Waves and High-Temperature Hydrodynamic Phenomena*

## Licence

MIT
