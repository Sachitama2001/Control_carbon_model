"""Reproduce the uncalibrated C--N benchmark and TeX figures."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from control_carbon.minimal_nitrogen_tipping import (
    Parameters, growth, equilibria, jacobian, simulate, classify,
    critical_duration, fast_step_certificate,
)


def main():
    out = Path(__file__).resolve().parents[1]/"artifacts/minimal_nitrogen_tipping"
    out.mkdir(parents=True, exist_ok=True)
    results = {"provenance": "Model-agnostic hypothesis, not a forest fit",
               "certificate": fast_step_certificate(), "ramps": {}, "critical": {}}
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    n = np.linspace(0, 12, 500)
    axes[0].plot(n, growth(n), label="N excess inhibition")
    axes[0].plot(n, growth(n, Parameters(inhibition=False)), label="Monod control")
    axes[0].axhline(.1, color="k", ls=":", label="turnover m")
    axes[0].set(xlabel="Mineral N / N0", ylabel="NPP / live C [1/yr]")
    axes[0].legend(fontsize=7)
    for duration in (10., 40., 100.):
        # Same calendar center, final time AND cumulative nitrogen input.
        t, y = simulate(duration, hold=1000-100-duration/2)
        calendar = t+100-duration/2
        for ax, j in zip(axes[1:], (0, 1)):
            ax.plot(calendar, y[:, j], label=f"T={duration:g} yr")
            ax.set_xlim(40, 350)
            ax.set_xlabel("Calendar time [yr]")
        results["ramps"][str(duration)] = {"final": y[-1].tolist(),
            "outcome": classify(y[-1]), "integrated_influx": 9300.}
    axes[1].set_ylabel("Live carbon / (q N0)")
    axes[2].set_ylabel("Mineral N / N0")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out/"nitrogen_response.pdf")
    plt.close(fig)
    for label, p, kwargs in [
        ("baseline", Parameters(), {}),
        ("refined", Parameters(), {"rtol": 1e-10, "max_step": 1.}),
        ("smoothstep", Parameters(), {"smooth": True}),
        ("recycle_k0.02", Parameters(decomposition=.02), {}),
        ("recycle_k0.1", Parameters(decomposition=.1), {}),
        ("recycle_k1", Parameters(decomposition=1.), {}),
        ("recycle_k10", Parameters(decomposition=10.), {}),
    ]:
        bracket = critical_duration(p, iterations=13, **kwargs)
        spectral = max(np.max(np.linalg.eigvals(jacobian(equilibria(i,p)["forest"],p)).real)
                       for i in np.linspace(3,10,101))
        results["critical"][label] = {"duration_bracket_years": bracket,
                                      "max_forest_real_eigenvalue": float(spectral)}
        print(label, bracket, flush=True)
    for duration in (10., 40.):
        p = Parameters(inhibition=False)
        results[f"monod_{duration}"] = classify(simulate(duration,p)[1][-1],p=p)
    (out/"summary.json").write_text(json.dumps(results, indent=2)+"\n")
    print(out/"summary.json")


if __name__ == "__main__":
    main()
