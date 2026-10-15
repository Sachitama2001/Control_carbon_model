"""Conceptual comparison at matched origin slope and saturation."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts/spac_carbon/rational_20261015/conceptual_comparison"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    font = FontProperties(fname="/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf")
    x = np.linspace(0, 8, 801)
    exponential = -np.expm1(-x)
    rational = x/(1+x)
    assert np.all(exponential >= rational)
    assert exponential[0] == rational[0] == 0
    with plt.rc_context({"font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False}):
        fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
        ax.fill_between(x, rational, exponential, color="#dddddd", alpha=.6)
        ax.plot(x, exponential, color="#0072B2", lw=2.8,
                label=r"指数関数：$1-e^{-qC_L}$")
        ax.plot(x, rational, color="#D55E00", lw=2.8, ls="--",
                label=r"有理関数：$qC_L/(1+qC_L)$")
        ax.axhline(1, color=".45", ls=":", lw=1.3)
        ax.set(xlim=(0, 8), ylim=(0, 1.1))
        ax.set_xlabel(r"Dimensionless leaf carbon, $qC_L$ (–)", fontsize=13)
        ax.set_ylabel(r"Normalized photosynthetic response, $A/(ph)$ (–)", fontsize=13)
        ax.legend(prop=font, loc="lower right", frameon=True)
        ax.grid(alpha=.15)
        for extension in ["png", "pdf", "svg"]:
            fig.savefig(OUTPUT / f"canopy_response_comparison.{extension}", dpi=200)
        plt.close(fig)
    print(OUTPUT / "canopy_response_comparison.png")


if __name__ == "__main__":
    main()
