"""Plot native VISIT and tangent responses to a soil-temperature pulse."""

from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from control_carbon import (
    DECOMPOSITION_ENVIRONMENT_NAMES,
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
    VISITSoilParameters,
    build_native_visit_soil_bridge,
    compare_native_visit_soil_environment_perturbation,
    visit_decomposition_scalars,
    visit_soil_fixed_point,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--visit-source",
        type=Path,
        default=Path("/mnt/d/CT/VISIT-matrix/visit_local"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/visit_soil_temperature_pulse.png"),
    )
    args = parser.parse_args()

    soil = VISITSoilParameters(
        1.1, 0.3, 0.75, 0.18, 0.08, 0.025,
        0.5, 0.5, 0.5, 0.6, 0.3, 0.1,
    )
    decomposition = VISITDecompositionParameters(0.22, 0.1, 0.21, 0.07)
    environment = VISITDecompositionEnvironment(
        10.0, 8.0, 50.0, 100.0, 0.8, 0.9, 120.0, 300.0
    )
    litter_inputs = np.ones(6)
    scalars = visit_decomposition_scalars(environment, decomposition)
    fixed_point = visit_soil_fixed_point(
        soil, litter_inputs, scalars.litter, scalars.humus
    ).state

    n_steps = 120
    pulse_day = 5
    pulse_kelvin = 2.0
    inputs = np.repeat(litter_inputs[np.newaxis, :], n_steps, axis=0)
    environments = (environment,) * n_steps
    perturbations = np.zeros((n_steps, len(DECOMPOSITION_ENVIRONMENT_NAMES)))
    temperature_column = DECOMPOSITION_ENVIRONMENT_NAMES.index(
        "soil_temperature_10cm"
    )
    perturbations[pulse_day, temperature_column] = pulse_kelvin

    with TemporaryDirectory() as directory:
        executable = build_native_visit_soil_bridge(
            args.visit_source, Path(directory) / "visit_soil_bridge"
        )
        comparison = compare_native_visit_soil_environment_perturbation(
            executable,
            args.visit_source,
            fixed_point,
            inputs,
            soil,
            decomposition,
            environments,
            perturbations,
        )

    days = np.arange(n_steps)
    state_days = np.arange(n_steps + 1)
    native_state_difference = comparison.native_state_difference
    python_state_difference = (
        comparison.perturbed_python_states - comparison.baseline_python_states
    )
    tangent_state_difference = comparison.tangent_state_difference

    def project(states: np.ndarray) -> np.ndarray:
        return np.column_stack((states[:, :6].sum(axis=1), states[:, 6:].sum(axis=1)))

    native_projection = project(native_state_difference)
    tangent_projection = project(tangent_state_difference)
    native_rh = comparison.native_heterotrophic_respiration_difference
    python_rh = (
        comparison.perturbed_python_heterotrophic_respiration
        - comparison.baseline_python_heterotrophic_respiration
    )
    tangent_rh = comparison.tangent_heterotrophic_respiration_difference

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    axes[0, 0].step(days, perturbations[:, temperature_column], where="post", color="#b23a2b")
    axes[0, 0].set(title="Forcing", xlabel="Day", ylabel="Upper-soil temperature perturbation (K)")

    axes[0, 1].plot(days, native_rh, label="Native VISIT C", color="#15616d", linewidth=2.2)
    axes[0, 1].plot(days, python_rh, label="Python nonlinear", color="#ff7d00", linestyle="--")
    axes[0, 1].plot(days, tangent_rh, label="Tangent prediction", color="#6a4c93", linestyle=":", linewidth=2.0)
    axes[0, 1].set(title="Heterotrophic respiration response", xlabel="Day", ylabel="Delta Rh (Mg C ha$^{-1}$ day$^{-1}$)")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].plot(native_projection[:, 0], native_projection[:, 1], color="#15616d", linewidth=2.2, label="Native VISIT C")
    axes[1, 0].plot(tangent_projection[:, 0], tangent_projection[:, 1], color="#6a4c93", linestyle="--", label="Tangent prediction")
    axes[1, 0].scatter([0.0], [0.0], marker="*", s=90, color="#222222", label="Reference fixed point")
    marker_days = (pulse_day + 1, 30, 60, 120)
    marker_offsets = ((4, 8), (6, -8), (6, 5), (8, 5))
    for marker_day, marker_offset in zip(marker_days, marker_offsets, strict=True):
        axes[1, 0].scatter(*native_projection[marker_day], s=24, color="#15616d")
        axes[1, 0].annotate(
            str(marker_day),
            native_projection[marker_day],
            xytext=marker_offset,
            textcoords="offset points",
            fontsize=8,
        )
    axes[1, 0].set(
        title="Projected state perturbation trajectory",
        xlabel="Delta total litter C (Mg C ha$^{-1}$)",
        ylabel="Delta total humus C (Mg C ha$^{-1}$)",
    )
    axes[1, 0].legend(frameon=False)

    state_error = np.max(np.abs(native_state_difference - tangent_state_difference), axis=1)
    rh_error = np.abs(native_rh - tangent_rh)
    error_floor = 1e-16
    state_line = axes[1, 1].semilogy(
        state_days,
        np.maximum(state_error, error_floor),
        label="Max pool-state error",
        color="#2a9d8f",
    )[0]
    axes[1, 1].set(
        title="Tangent approximation error",
        xlabel="Day",
        ylabel="State error (Mg C ha$^{-1}$)",
    )
    rh_axis = axes[1, 1].twinx()
    rh_line = rh_axis.semilogy(
        days,
        np.maximum(rh_error, error_floor),
        label="Rh error",
        color="#d1495b",
    )[0]
    rh_axis.set_ylabel("Rh error (Mg C ha$^{-1}$ day$^{-1}$)")
    axes[1, 1].legend(handles=(state_line, rh_line), frameon=False, loc="lower right")

    fig.suptitle("TKY class 4 soil response to a 2 K one-day upper-soil pulse")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())