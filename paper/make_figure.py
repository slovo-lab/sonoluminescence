"""Figure 1 of the preprint: distribution of photons per flash from the Monte Carlo budget."""
import json, math, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = pathlib.Path(__file__).resolve().parent
res = json.load(open(HERE.parent / "mc.json"))
LIN = {"water": 4.18, "acid": 7.07}          # linearised u(ln photons), v0.8
TITLE = {"en": {"water": "Argon bubble in water, 20 °C", "acid": "Argon bubble in 85 % sulfuric acid"},
         "ru": {"water": "Пузырёк аргона в воде, 20 °C", "acid": "Пузырёк аргона в 85%-й серной кислоте"}}
LAB = {"en": dict(x="Photons per flash", y="Monte Carlo draws", base="baseline", band="linearised budget, ±1σ",
                  mid="68 % of draws", zero="<0.1"),
       "ru": dict(x="Фотонов за вспышку", y="Выборок Монте-Карло", base="базовый расчёт", band="линейный бюджет, ±1σ",
                  mid="68 % выборок", zero="<0,1")}
BLUE, INK, INK2, GRID = "#257de0", "#1f1f1f", "#6b6b6b", "#e6e6e6"

for lang in ("en", "ru"):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=False)
    for ax, case in zip(axes, ("water", "acid")):
        ph = np.array([r["photons"] for r in res if r["case"] == case and r["i"] >= 0])
        base = next(r["photons"] for r in res if r["case"] == case and r["i"] < 0)
        floor = 1e-1
        lp = np.log10(np.maximum(ph, floor))
        bins = np.arange(-1, 8.01, 0.5)
        ax.hist(lp, bins=bins, color=BLUE, edgecolor="white", linewidth=2)
        q16, q50, q84 = np.percentile(ph, [16, 50, 84])
        ymax = ax.get_ylim()[1]
        lo, hi = math.log10(base) - LIN[case] / math.log(10), math.log10(base) + LIN[case] / math.log(10)
        ax.plot([max(lo, -1), min(hi, 8)], [ymax * 1.30] * 2, color=INK2, lw=2, ls=(0, (4, 3)))
        ax.text(max(lo, -1), ymax * 1.34, LAB[lang]["band"], ha="left", va="bottom", color=INK2, fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", pad=1))
        ax.plot([math.log10(max(q16, floor)), math.log10(q84)], [ymax * 1.15] * 2, color=INK, lw=2)
        ax.text(math.log10(max(q16, floor)), ymax * 1.19, LAB[lang]["mid"], ha="left", va="bottom", color=INK, fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", pad=1))
        ax.plot([math.log10(base)] * 2, [0, ymax * 0.99], color=INK, lw=1)
        ax.text(math.log10(base), ymax * 1.00, LAB[lang]["base"], ha="center", va="bottom", color=INK, fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
        ax.set_ylim(0, ymax * 1.45)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.set_xlim(-1, 8)
        ticks = list(range(-1, 9))
        ax.set_xticks(ticks)
        ax.set_xticklabels([LAB[lang]["zero"] if t == -1 else ("1" if t == 0 else f"$10^{t}$") for t in ticks], fontsize=8)
        ax.set_title(TITLE[lang][case], fontsize=10, color=INK, loc="left")
        ax.set_xlabel(LAB[lang]["x"], fontsize=9, color=INK2)
        ax.set_ylabel(LAB[lang]["y"], fontsize=9, color=INK2)
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2)
    fig.tight_layout()
    fig.savefig(HERE / f"fig1_{lang}.png", dpi=200)
print("ok")
