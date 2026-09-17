"""Generate figures/results for the capacity-attractor exploration note."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from control_carbon.capacity_feedback import (
    CubicCapacityFeedback,
    critical_ramp_duration,
    periodic_attractor,
    periodic_capacity,
    simulate_capacity_ramp,
)


def main():
    model = CubicCapacityFeedback()
    output = Path(__file__).resolve().parents[1]/"artifacts/capacity_feedback"
    output.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    stocks = np.linspace(0, 11, 500)
    axes[0].plot(stocks, stocks, "k--", label="identity x")
    for driver in (0, .5, 1):
        axes[0].plot(stocks, model.capacity(stocks, driver),
                     label=f"capacity, lambda={driver:g}")
    axes[0].set(xlabel="Current carbon stock x",
                ylabel="Pointwise capacity c-hat(x, lambda)", ylim=(0, 11))
    axes[0].legend(fontsize=8)

    records = {}
    for duration, color in ((20, "tab:red"), (120, "tab:blue")):
        result = simulate_capacity_ramp(duration, model)
        t, x, z = result["time"], result["stock"], result["driver"]
        low = np.ones_like(t)*model.low
        threshold = model.threshold_initial+model.threshold_shift*z
        high = model.high_initial+model.high_shift*z
        axes[1].plot(t, x, color=color, lw=2,
                     label=f"stock, T={duration} yr ({result['outcome']})")
        axes[1].plot(t, high, color=color, ls=":", alpha=.7)
        if duration == 20:
            axes[1].plot(t, threshold, "k--", lw=1, label="basin boundary")
            axes[1].plot(t, low, "k:", lw=1, label="low capacity branch")
        records[str(duration)] = {
            "outcome": result["outcome"],
            "final_stock": float(x[-1]),
            "minimum_stock": float(np.min(x)),
        }
    axes[1].set(xlabel="Time [yr]", ylabel="Carbon stock / branch",
                xlim=(0, 180), ylim=(0, 11))
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output/"capacity_map_and_ramp.pdf")
    plt.close(fig)

    t = np.linspace(0, 2, 1000)
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    ax.plot(t, periodic_capacity(t), label="instantaneous QSE/capacity")
    ax.plot(t, periodic_attractor(t), label="attracting periodic trajectory")
    ax.set(xlabel="Time [yr]", ylabel="Carbon stock", title="Periodic forcing baseline")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output/"periodic_capacity_lag.pdf")
    plt.close(fig)

    bracket = critical_ramp_duration(model, iterations=15)
    summary = {
        "classification": "model-agnostic mathematical counterexample; not VISIT",
        "parameters": model.__dict__,
        "duration_bracket_years": bracket,
        "outcomes": records,
        "frozen_matrix_eigenvalue": -model.relaxation,
    }
    (output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
