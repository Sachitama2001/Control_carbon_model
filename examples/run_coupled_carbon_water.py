"""Run the first reproducible idealized carbon--water experiment.

The experiment starts from a coupled frozen-forcing equilibrium, applies a
100-day precipitation reduction, and then restores the baseline climate.  It
writes numerical arrays and a provenance-rich JSON manifest; it does not claim
that the illustrative parameter set is calibrated to VISITc or a field site.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from control_carbon.coupled_carbon_water import (
    COUPLED_STATE_NAMES,
    COUPLED_STATE_UNITS,
    CoupledCarbonWaterForcing,
    CoupledCarbonWaterParameters,
    coupled_equilibrium,
    coupled_model_provenance_manifest,
    evaluate_coupled_carbon_water,
    instantaneous_carbon_capacity,
    simulate_coupled_carbon_water,
)


def protocol(time: float) -> CoupledCarbonWaterForcing:
    # A smooth window avoids making adaptive-solver differences at a forcing
    # discontinuity look like state-equation error.
    dry_window = 0.5 * (
        np.tanh((time - 40.0) / 0.75) - np.tanh((time - 140.0) / 0.75)
    )
    precipitation = 1.5 - 1.35 * dry_window
    return CoupledCarbonWaterForcing(precipitation=precipitation)


def run(output_directory: Path) -> dict[str, object]:
    parameters = CoupledCarbonWaterParameters()
    baseline = CoupledCarbonWaterForcing(precipitation=1.5)
    equilibrium = coupled_equilibrium(
        baseline,
        [2.0, 35.0, 6.0, 120.0, 0.3, 7.0, 3.0, 155.0],
        parameters,
    )
    if not equilibrium.converged or equilibrium.stability != "stable":
        raise RuntimeError("illustrative baseline equilibrium was not stable")

    times = np.linspace(0.0, 240.0, 481)
    trajectory = simulate_coupled_carbon_water(
        equilibrium.state, times, protocol, parameters, max_step=0.25
    )
    refined = simulate_coupled_carbon_water(
        equilibrium.state, times, protocol, parameters, max_step=0.125
    )
    refinement_error = float(np.max(np.abs(trajectory.states - refined.states)))

    carbon_capacity = np.empty((times.size, 4))
    carbon_capacity_total = np.empty(times.size)
    carbon_potential = np.empty(times.size)
    carbon_budget_residual = np.empty(times.size)
    water_budget_residual = np.empty(times.size)
    precipitation = np.empty(times.size)
    for index, (time, state) in enumerate(zip(times, trajectory.states, strict=True)):
        forcing = protocol(float(time))
        capacity = instantaneous_carbon_capacity(state, forcing, parameters)
        evaluation = evaluate_coupled_carbon_water(state, forcing, parameters)
        carbon_capacity[index] = capacity.state
        carbon_capacity_total[index] = capacity.total
        carbon_potential[index] = capacity.total - state[:4].sum()
        carbon_budget_residual[index] = evaluation.carbon_budget_residual
        water_budget_residual[index] = evaluation.water_budget_residual
        precipitation[index] = forcing.precipitation

    output_directory.mkdir(parents=True, exist_ok=True)
    arrays_path = output_directory / "coupled_carbon_water_drydown.npz"
    manifest_path = output_directory / "coupled_carbon_water_drydown.json"
    np.savez_compressed(
        arrays_path,
        time_days=times,
        states=trajectory.states,
        precipitation=precipitation,
        equilibrium=equilibrium.state,
        carbon_capacity=carbon_capacity,
        carbon_capacity_total=carbon_capacity_total,
        signed_carbon_potential=carbon_potential,
        carbon_budget_residual=carbon_budget_residual,
        water_budget_residual=water_budget_residual,
    )
    manifest: dict[str, object] = {
        "experiment": "idealized precipitation dry-down and re-wetting",
        "approximation": "source-grounded reduced continuous carbon-water model",
        "state_names": list(COUPLED_STATE_NAMES),
        "state_units": list(COUPLED_STATE_UNITS),
        "time_unit": "day",
        "forcing_protocol": {
            "baseline_precipitation_mm_day": 1.5,
            "dry_precipitation_mm_day": 0.15,
            "dry_interval_days": [40.0, 140.0],
            "tanh_transition_scale_day": 0.75,
        },
        "equilibrium": {
            "state": equilibrium.state.tolist(),
            "residual_infinity_norm": equilibrium.residual_norm,
            "stability": equilibrium.stability,
            "maximum_real_eigenvalue": float(np.max(equilibrium.eigenvalues.real)),
        },
        "numerics": {
            "method": "DOP853",
            "reported_max_step_day": trajectory.max_step,
            "refined_max_step_day": refined.max_step,
            "maximum_state_refinement_difference": refinement_error,
            "maximum_absolute_carbon_budget_residual": float(
                np.max(np.abs(carbon_budget_residual))
            ),
            "maximum_absolute_water_budget_residual": float(
                np.max(np.abs(water_budget_residual))
            ),
        },
        "files": {"arrays": arrays_path.name},
        "provenance": coupled_model_provenance_manifest(),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "arrays": arrays_path,
        "manifest": manifest_path,
        "equilibrium_stability": equilibrium.stability,
        "refinement_error": refinement_error,
        "carbon_budget_residual": manifest["numerics"][
            "maximum_absolute_carbon_budget_residual"
        ],
        "water_budget_residual": manifest["numerics"][
            "maximum_absolute_water_budget_residual"
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/coupled_carbon_water"),
    )
    result = run(parser.parse_args().output_dir)
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
