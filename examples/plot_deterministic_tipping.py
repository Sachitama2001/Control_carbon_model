"""Reproduce deterministic mechanism checks and a temperature-only soil control."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from control_carbon import (
    SaturatingFeedback,
    SmoothRamp,
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
    VISITSoilParameters,
    run_feedback_ramp,
    visit_decomposition_scalars,
    visit_soil_fixed_point,
)
from control_carbon.visit_source_map import (
    VISIT_PARAMETER_WORKBOOK_SHA256,
    VISIT_SOURCE_COMMIT,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/deterministic_tipping.png"))
    args = parser.parse_args()

    model = SaturatingFeedback()
    configurations = {
        "B": ("gain", SmoothRamp(3, 1.5, 0.02)),
        "R_slow": ("half_saturation", SmoothRamp(1, 10, 0.01)),
        "R_fast": ("half_saturation", SmoothRamp(1, 10, 10)),
    }
    experiments = {
        name: run_feedback_ramp(model, driver, protocol)
        for name, (driver, protocol) in configurations.items()
    }
    data = {}
    reports = {}
    convergence_errors = []
    for name, result in experiments.items():
        driver, protocol = configurations[name]
        refined = run_feedback_ramp(
            model, driver, protocol, rtol=1e-11, atol=1e-13, steps_per_ramp=200,
        )
        error = float(max(
            np.max(np.abs(result.ramp.states - refined.ramp.states[::2])),
            np.max(np.abs(result.hold.states - refined.hold.states)),
        ))
        if error > 1e-6 or result.final_basin != refined.final_basin:
            raise RuntimeError(f"Refinement changed the result for {name}")
        convergence_errors.append(error)
        data[f"{name}_time"] = np.concatenate((result.ramp.times, result.hold.times[1:]))
        data[f"{name}_stock"] = np.concatenate((result.ramp.states[:, 0], result.hold.states[1:, 0]))
        data[f"{name}_driver"] = np.concatenate((result.ramp.drivers, result.hold.drivers[1:]))
        reports[name] = {
            "driver": driver, "protocol": asdict(protocol),
            "max_driver_speed": protocol.max_driver_speed,
            "initial_stock": float(result.ramp.states[0, 0]),
            "hold_time": float(result.hold.times[-1] - result.hold.times[0]),
            "final_stock": float(result.hold.states[-1, 0]),
            "final_basin": result.final_basin,
            "distance_to_attractor": result.distance_to_attractor,
            "final_residual": result.final_residual,
            "refinement_max_stock_difference": error,
            "rtol": result.ramp.rtol, "atol": result.ramp.atol,
            "ramp_max_step": result.ramp.max_step, "hold_max_step": result.hold.max_step,
            "method": result.ramp.method,
        }

    rates = np.logspace(-2, 1, 19)
    scan = [run_feedback_ramp(model, "half_saturation", SmoothRamp(1, 10, rate)) for rate in rates]
    data["scan_rates"] = rates
    data["scan_final_stock"] = np.array([result.hold.states[-1, 0] for result in scan])
    data["scan_ramp_end_stock"] = np.array([result.ramp.states[-1, 0] for result in scan])
    data["scan_final_basin"] = np.array([result.final_basin for result in scan])

    soil = VISITSoilParameters(1.1, 0.3, 0.75, 0.18, 0.08, 0.025, 0.5, 0.5, 0.5, 0.6, 0.3, 0.1)
    decomposition = VISITDecompositionParameters(0.22, 0.1, 0.21, 0.07)
    environment = VISITDecompositionEnvironment(10, 8, 50, 100, 0.8, 0.9, 120, 300)
    litter = np.full(6, 0.01)
    temperatures = np.linspace(0, 25, 51)
    soil_totals = []
    soil_residuals = []
    for temperature in temperatures:
        scalars = visit_decomposition_scalars(
            replace(environment, soil_temperature_10cm=temperature), decomposition,
        )
        equilibrium = visit_soil_fixed_point(soil, litter, scalars.litter, scalars.humus)
        if not equilibrium.converged or not equilibrium.stable:
            raise RuntimeError("Temperature control lacks a stable unique frozen fixed point")
        soil_totals.append(equilibrium.state.sum())
        soil_residuals.append(equilibrium.residual_norm)
    data["soil_temperature_C"] = temperatures
    data["soil_total_carbon"] = np.array(soil_totals)
    data["soil_fixed_point_residual"] = np.array(soil_residuals)

    gain_values = np.linspace(2, 4, 201)
    branches = np.array([model.equilibria("gain", gain) if gain > 2 else [0, 1, 1] for gain in gain_values])
    data["branch_gain"] = gain_values
    data["branch_stock"] = branches
    half_saturations = np.linspace(1, 10, 201)
    moving_roots = np.array([model.equilibria("half_saturation", value) for value in half_saturations])
    data["moving_half_saturation"] = half_saturations
    data["moving_roots"] = moving_roots

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    stable_color, fast_color, boundary_color = "#166970", "#bf423f", "#70634b"
    axis = axes[0, 0]
    reference = np.interp(10, temperatures, soil_totals)
    axis.plot(temperatures, np.array(soil_totals) / reference, color=stable_color)
    axis.set(title="A. VISIT soil: temperature-only control", xlabel="Upper-soil temperature (C)", ylabel="Frozen total C / total C at 10 C")
    axis.text(0.04, 0.08, "9 pools; native daily algebra\nOther environment and litter fixed", transform=axis.transAxes, fontsize=9)

    axis = axes[0, 1]
    axis.plot([1, 4], [0, 0], color=stable_color)
    axis.plot(gain_values, branches[:, 2], color=stable_color, label="Stable feedback roots")
    axis.plot(gain_values, branches[:, 1], "--", color=boundary_color, label="Unstable threshold")
    axis.plot([1, 4], [1, 4], ":", color="#555555", label="No-feedback equilibrium")
    axis.scatter([2], [1], color=fast_color, zorder=5, label="Fold: g = 2")
    axis.set(title="B. Analytic branches (H = m = 1)", xlabel="Feedback gain g = P / (m H)", ylabel="Normalized carbon")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[0, 2]
    axis.plot(data["B_time"], data["B_stock"], color=stable_color, label="Simulated carbon")
    qse = np.array([model.equilibria("gain", value)[-1] if value >= 2 else np.nan for value in data["B_driver"]])
    axis.plot(data["B_time"], qse, "--", color=boundary_color, label="High frozen branch")
    axis.axvline(experiments["B"].protocol.end_time, color="#aaaaaa", linestyle=":", label="Final hold begins")
    axis.set(title="C. B experiment: branch disappears", xlabel="Model time", ylabel="Normalized carbon")
    forcing_axis = axis.twinx()
    forcing_axis.plot(data["B_time"], data["B_driver"], color=fast_color, alpha=0.5)
    forcing_axis.axhline(2, color=fast_color, linewidth=0.7, linestyle=":")
    forcing_axis.set_ylabel("Gain g(t)", color=fast_color)
    axis.legend(frameon=False, fontsize=8, loc="center right")

    axis = axes[1, 0]
    axis.plot(half_saturations, moving_roots[:, 2], "--", color=boundary_color, label="High frozen branch")
    axis.plot(half_saturations, moving_roots[:, 1], ":", color=boundary_color, label="Basin boundary")
    for name, color, label in [("R_slow", stable_color, "Slow: r = 0.01"), ("R_fast", fast_color, "Fast: r = 10")]:
        axis.plot(data[f"{name}_driver"], data[f"{name}_stock"], color=color, linewidth=2, label=label)
        axis.scatter([10], [data[f"{name}_stock"][-1]], color=color, marker="s", s=25)
    axis.scatter([1], [model.equilibria("half_saturation", 1)[-1]], color="#222222", marker="*", s=90)
    axis.set(title="D. R experiment: same path, two rates", xlabel="Effective driver H (g = 3 fixed)", ylabel="Normalized carbon")
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[1, 1]
    axis.semilogx(rates, data["scan_final_stock"], "o-", color=stable_color, label="After 40-time-unit hold")
    axis.semilogx(rates, data["scan_ramp_end_stock"], ".--", color=fast_color, label="At ramp end")
    axis.axhline(moving_roots[-1, 1], color=boundary_color, linestyle=":", label="Final basin boundary")
    axis.set(title="E. Rate scan (not a critical-rate estimate)", xlabel="Ramp rate r (inverse model time)", ylabel="Normalized carbon")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[1, 2]
    axis.bar(list(experiments), np.maximum(convergence_errors, 1e-16), color=[boundary_color, stable_color, fast_color])
    axis.set_yscale("log")
    axis.set(title="F. Solver refinement check", ylabel="Maximum stock difference", xlabel="Tighter tolerances; half internal max-step")
    axis.text(0.04, 0.94, "Default: rtol 1e-9, atol 1e-11\nRefined: rtol 1e-11, atol 1e-13", transform=axis.transAxes, va="top", fontsize=9)

    for axis in axes.flat:
        axis.grid(alpha=0.18)
    fig.suptitle("Deterministic checks: source-grounded soil control and hypothetical feedback", fontsize=15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    np.savez_compressed(args.output.with_suffix(".npz"), **data)
    manifest = {
        "scope": "Deterministic synthetic mechanisms, not a VISIT/TKY tipping prediction",
        "feedback_model": asdict(model), "experiments": reports,
        "rate_scan": {"rates": rates.tolist(), "hold_time": 40, "rtol": 1e-9, "atol": 1e-11},
        "soil_control": {
            "source_commit": VISIT_SOURCE_COMMIT,
            "source_functions": ["visit_local/soil_proc.c::f_cycle_soil", "visit_local/decomposition.c::frl", "visit_local/decomposition.c::frh"],
            "approximation": "source-grounded reduced compartment model; native daily carbon algebra",
            "workbook_sha256": VISIT_PARAMETER_WORKBOOK_SHA256, "land_cover_class": 4,
            "varying": "soil_temperature_10cm", "temperature_range_C": [0, 25],
            "reference_environment": asdict(environment),
            "soil_parameters": asdict(soil), "decomposition_parameters": asdict(decomposition),
            "fixed_litter_inputs_Mg_C_ha_day": litter.tolist(),
            "max_fixed_point_residual": float(max(soil_residuals)),
            "warning": "Controlled synthetic environment, not measured TKY weather or hydrological equilibrium",
        },
    }
    args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(reports, indent=2))
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())