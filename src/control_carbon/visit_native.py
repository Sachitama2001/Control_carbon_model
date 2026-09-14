"""Minimal execution bridge to the authoritative VISIT soil C routines."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
from pathlib import Path
import subprocess
from typing import Final, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .visit_decomposition import (
    DECOMPOSITION_ENVIRONMENT_NAMES,
    DecompositionTemperatureMode,
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
    visit_decomposition_regime,
    visit_decomposition_scalars,
)
from .visit_soil import VISITSoilParameters, visit_soil_daily_step
from .visit_soil import visit_soil_discrete_system, visit_soil_fixed_point
from .visit_plant import (
    VISITPlantAllocationParameters,
    VISITPlantAllocationResult,
    VISITPlantRespirationParameters,
    VISITPlantRespirationResult,
    VISITPlantStructuralState,
    VISITPlantTurnoverContext,
    VISITPlantTurnoverParameters,
    VISITPlantTurnoverResult,
)
from .visit_source_map import VISIT_SOURCE_COMMIT


NATIVE_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2] / "native" / "visit_soil_bridge.c"
)
NATIVE_PLANT_TURNOVER_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2]
    / "native"
    / "visit_plant_turnover_bridge.c"
)
NATIVE_PLANT_RESPIRATION_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2]
    / "native"
    / "visit_plant_respiration_bridge.c"
)
NATIVE_PLANT_ALLOCATION_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2]
    / "native"
    / "visit_plant_allocation_bridge.c"
)
@dataclass(frozen=True)
class NativeVISITSoilTrajectory:
    """Outputs from direct execution of VISIT's C soil routines."""

    states: np.ndarray
    heterotrophic_respiration: np.ndarray
    litter_scalars: np.ndarray
    humus_scalars: np.ndarray
    environments: np.ndarray
    executable: Path
    source_commit: str
    source_sha256: dict[str, str]
    excluded_processes: tuple[str, ...]


@dataclass(frozen=True)
class NativeVISITSoilComparison:
    """Native C and Python soil trajectories under identical inputs."""

    native: NativeVISITSoilTrajectory
    python_states: np.ndarray
    python_heterotrophic_respiration: np.ndarray

    @property
    def state_error(self) -> np.ndarray:
        return self.native.states - self.python_states

    @property
    def heterotrophic_respiration_error(self) -> np.ndarray:
        return (
            self.native.heterotrophic_respiration
            - self.python_heterotrophic_respiration
        )


@dataclass(frozen=True)
class NativeVISITSoilPerturbationComparison:
    """Native baseline/perturbation differences and the soil IRF prediction."""

    baseline: NativeVISITSoilTrajectory
    perturbed: NativeVISITSoilTrajectory
    irf_state_difference: np.ndarray
    irf_heterotrophic_respiration_difference: np.ndarray

    @property
    def native_state_difference(self) -> np.ndarray:
        return self.perturbed.states - self.baseline.states

    @property
    def native_heterotrophic_respiration_difference(self) -> np.ndarray:
        return (
            self.perturbed.heterotrophic_respiration
            - self.baseline.heterotrophic_respiration
        )

    @property
    def state_error(self) -> np.ndarray:
        return self.native_state_difference - self.irf_state_difference

    @property
    def heterotrophic_respiration_error(self) -> np.ndarray:
        return (
            self.native_heterotrophic_respiration_difference
            - self.irf_heterotrophic_respiration_difference
        )


@dataclass(frozen=True)
class VISITSoilEnvironmentalLinearization:
    """Daily tangent system along a baseline soil/environment trajectory."""

    baseline_states: np.ndarray
    baseline_heterotrophic_respiration: np.ndarray
    state_jacobians: np.ndarray
    environment_state_jacobians: np.ndarray
    output_state_jacobians: np.ndarray
    environment_output_jacobians: np.ndarray
    environment_names: tuple[str, ...]

    def predict(self, environment_perturbations: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        perturbations = np.asarray(environment_perturbations, dtype=float)
        expected_shape = (self.state_jacobians.shape[0], len(self.environment_names))
        if perturbations.shape != expected_shape or not np.all(np.isfinite(perturbations)):
            raise ValueError(
                f"environment_perturbations must be finite with shape {expected_shape}"
            )
        state_difference = np.zeros_like(self.baseline_states)
        output_difference = np.empty(self.state_jacobians.shape[0], dtype=float)
        for step, environment_difference in enumerate(perturbations):
            output_difference[step] = (
                self.output_state_jacobians[step] @ state_difference[step]
                + self.environment_output_jacobians[step] @ environment_difference
            )
            state_difference[step + 1] = (
                self.state_jacobians[step] @ state_difference[step]
                + self.environment_state_jacobians[step] @ environment_difference
            )
        return state_difference, output_difference


@dataclass(frozen=True)
class NativeVISITSoilEnvironmentPerturbationComparison:
    """Native/Python nonlinear responses and tangent environmental prediction."""

    baseline_native: NativeVISITSoilTrajectory
    perturbed_native: NativeVISITSoilTrajectory
    baseline_python_states: np.ndarray
    perturbed_python_states: np.ndarray
    baseline_python_heterotrophic_respiration: np.ndarray
    perturbed_python_heterotrophic_respiration: np.ndarray
    tangent_state_difference: np.ndarray
    tangent_heterotrophic_respiration_difference: np.ndarray

    @property
    def native_state_difference(self) -> np.ndarray:
        return self.perturbed_native.states - self.baseline_native.states

    @property
    def native_heterotrophic_respiration_difference(self) -> np.ndarray:
        return (
            self.perturbed_native.heterotrophic_respiration
            - self.baseline_native.heterotrophic_respiration
        )

    @property
    def native_python_state_error(self) -> np.ndarray:
        return self.perturbed_native.states - self.perturbed_python_states

    @property
    def native_python_heterotrophic_respiration_error(self) -> np.ndarray:
        return (
            self.perturbed_native.heterotrophic_respiration
            - self.perturbed_python_heterotrophic_respiration
        )

    @property
    def tangent_state_error(self) -> np.ndarray:
        return self.native_state_difference - self.tangent_state_difference

    @property
    def tangent_heterotrophic_respiration_error(self) -> np.ndarray:
        return (
            self.native_heterotrophic_respiration_difference
            - self.tangent_heterotrophic_respiration_difference
        )


def build_native_visit_soil_bridge(
    visit_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile the bridge against unmodified VISIT soil/decomposition sources."""
    source_root = Path(visit_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        (
            "soil_proc.c",
            "decomposition.c",
            "structure.h",
            "prototype.h",
            "definition.h",
            "setting.h",
        ),
    )
    if not NATIVE_BRIDGE_SOURCE.is_file():
        raise FileNotFoundError(f"missing native bridge source: {NATIVE_BRIDGE_SOURCE}")
    command = [
        compiler,
        "-std=c11",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wno-unused-parameter",
        "-I",
        str(source_root),
        str(NATIVE_BRIDGE_SOURCE),
        str(source_root / "soil_proc.c"),
        str(source_root / "decomposition.c"),
        "-lm",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def build_native_visit_plant_turnover_bridge(
    visit_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile a harness around the litterfall block in native plant_process."""
    source_root = Path(visit_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        (
            "plant_proc.c",
            "litterfall.c",
            "structure.h",
            "prototype.h",
            "definition.h",
            "setting.h",
        ),
    )
    command = [
        compiler,
        "-std=c11",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wno-unused-parameter",
        "-I",
        str(source_root),
        str(NATIVE_PLANT_TURNOVER_BRIDGE_SOURCE),
        str(source_root / "plant_proc.c"),
        str(source_root / "litterfall.c"),
        "-lm",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def run_native_visit_plant_turnover(
    executable: str | Path,
    state: VISITPlantStructuralState,
    parameters: VISITPlantTurnoverParameters,
    context: VISITPlantTurnoverContext,
    *,
    vegetation_type: int = 4,
) -> VISITPlantTurnoverResult:
    """Execute native plant_process with non-turnover processes stubbed out."""
    values = (
        *state.as_array(),
        parameters.leaf_rate,
        parameters.stem_rate,
        parameters.root_rate,
        parameters.deciduous_shedding_fraction,
        context.season,
        context.crop_stage,
        context.lai,
        context.deep_soil_water_potential,
        context.air_temperature,
        context.photosynthesis_type,
        vegetation_type,
    )
    completed = subprocess.run(
        [str(Path(executable).resolve())],
        input=" ".join(str(value) for value in values) + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    output = np.asarray(
        [float(value) for value in completed.stdout.strip().split(",")], dtype=float
    )
    if output.shape != (6,):
        raise ValueError(f"native plant turnover bridge returned shape {output.shape}")
    regime = "native-plant-process"
    return VISITPlantTurnoverResult(
        state=VISITPlantStructuralState(*output[:3]),
        litterfall=output[3:],
        regime=regime,
    )


def build_native_visit_plant_respiration_bridge(
    visit_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile native coefficient-update and plant-respiration functions."""
    source_root = Path(visit_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        (
            "ecophysiology.c",
            "respiration.c",
            "structure.h",
            "prototype.h",
            "definition.h",
            "setting.h",
        ),
    )
    command = [
        compiler,
        "-std=c11",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wno-unused-parameter",
        "-ffunction-sections",
        "-fdata-sections",
        "-I",
        str(source_root),
        str(NATIVE_PLANT_RESPIRATION_BRIDGE_SOURCE),
        str(source_root / "ecophysiology.c"),
        str(source_root / "respiration.c"),
        "-Wl,--gc-sections",
        "-lm",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def run_native_visit_plant_respiration(
    executable: str | Path,
    state: VISITPlantStructuralState,
    parameters: VISITPlantRespirationParameters,
    *,
    surface_temperature: float,
    upper_soil_temperature: float,
    allocation_fluxes: ArrayLike,
    positive_epp: bool = True,
) -> VISITPlantRespirationResult:
    """Execute native Q10, size-specific, maintenance, and growth respiration."""
    allocations = _finite_array(allocation_fluxes, (3,), "allocation_fluxes")
    values = (
        *state.as_array(),
        parameters.growth_foliage,
        parameters.growth_stem,
        parameters.growth_root,
        parameters.maintenance_foliage_15c,
        parameters.maintenance_stem_sapwood,
        parameters.maintenance_stem_heartwood,
        parameters.maintenance_root_fine,
        parameters.maintenance_root_coarse,
        parameters.q10_foliage_base,
        parameters.q10_stem_base,
        parameters.q10_root_base,
        parameters.size_stem,
        parameters.size_root,
        surface_temperature,
        upper_soil_temperature,
        *allocations,
        int(positive_epp),
    )
    completed = subprocess.run(
        [str(Path(executable).resolve())],
        input=" ".join(str(value) for value in values) + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    output = np.asarray(
        [float(value) for value in completed.stdout.strip().split(",")], dtype=float
    )
    if output.shape != (16,):
        raise ValueError(f"native plant respiration bridge returned shape {output.shape}")
    return VISITPlantRespirationResult(
        q10=output[0:3],
        specific_maintenance=output[3:6],
        maintenance=output[6:9],
        growth=output[9:12],
        sapwood_fine_root=output[12:14],
        heartwood_coarse_root=output[14:16],
    )


def build_native_visit_plant_allocation_bridge(
    visit_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile a minimal harness around native ``f_allocation``."""
    source_root = Path(visit_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        ("allocation.c", "structure.h", "prototype.h", "definition.h", "setting.h"),
    )
    subprocess.run(
        [
            compiler,
            "-std=c11",
            "-O2",
            "-Wall",
            "-Wextra",
            "-Wno-unused-parameter",
            "-ffunction-sections",
            "-fdata-sections",
            "-I",
            str(source_root),
            str(NATIVE_PLANT_ALLOCATION_BRIDGE_SOURCE),
            str(source_root / "allocation.c"),
            "-Wl,--gc-sections",
            "-lm",
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output


def run_native_visit_plant_allocation(
    executable: str | Path,
    parameters: VISITPlantAllocationParameters,
    *,
    lai: float,
    optimum_lai: float,
    season: int,
    crop_stage: int,
    epp: float,
    gpp: float,
    maintenance_respiration: ArrayLike,
) -> VISITPlantAllocationResult:
    """Execute native ``allocation.c::f_allocation`` for one plant layer."""
    maintenance = _finite_array(
        maintenance_respiration, (3,), "maintenance_respiration"
    )
    values = (
        parameters.foliage_fraction,
        parameters.aboveground_woody_fraction,
        parameters.specific_leaf_area,
        lai,
        optimum_lai,
        season,
        crop_stage,
        epp,
        gpp,
        *maintenance,
    )
    completed = subprocess.run(
        [str(Path(executable).resolve())],
        input=" ".join(str(value) for value in values) + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    output = np.asarray(
        [float(value) for value in completed.stdout.strip().split(",")], dtype=float
    )
    if output.shape != (9,):
        raise ValueError(f"native plant allocation bridge returned shape {output.shape}")
    return VISITPlantAllocationResult(
        fractions=output[0:4],
        fluxes=output[4:8],
        total_translocation=float(output[8]),
        regime="native-f-allocation",
    )


def run_native_visit_soil(
    executable: str | Path,
    visit_source_root: str | Path,
    initial_state: ArrayLike,
    litter_inputs: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environment: VISITDecompositionEnvironment | Sequence[VISITDecompositionEnvironment],
) -> NativeVISITSoilTrajectory:
    """Run baseline ``EX_DECTMP=0`` carbon dynamics with daily environments."""
    initial = _finite_array(initial_state, (9,), "initial_state")
    inputs = np.asarray(litter_inputs, dtype=float)
    if inputs.ndim != 2 or inputs.shape[1] != 6 or not np.all(np.isfinite(inputs)):
        raise ValueError("litter_inputs must be finite with shape (n_steps, 6)")
    if np.any(initial < 0) or np.any(inputs < 0):
        raise ValueError("native soil states and litter inputs must be non-negative")

    parameter_values = (
        soil_parameters.sr_lf,
        soil_parameters.sr_lc,
        soil_parameters.sr_lr,
        soil_parameters.sr_ha,
        soil_parameters.sr_hi,
        soil_parameters.sr_hp,
        soil_parameters.f_co2_lf,
        soil_parameters.f_co2_lc,
        soil_parameters.f_co2_lr,
        soil_parameters.f_hm_a,
        soil_parameters.f_hm_i,
        soil_parameters.f_hm_p,
        decomposition_parameters.kml,
        decomposition_parameters.kmh,
        decomposition_parameters.kmsl,
        decomposition_parameters.kmsh,
    )
    environments = _environment_sequence(environment, inputs.shape[0])
    lines = [
        " ".join(str(value) for value in (*initial, *parameter_values)),
        str(inputs.shape[0]),
        *(
            " ".join(
                str(value)
                for value in (*input_values, *_environment_values(environment_values))
            )
            for input_values, environment_values in zip(inputs, environments, strict=True)
        ),
    ]
    completed = subprocess.run(
        [str(Path(executable).resolve())],
        input="\n".join(lines) + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    rows = np.asarray(
        [[float(value) for value in line.split(",")] for line in completed.stdout.splitlines()],
        dtype=float,
    )
    if rows.shape != (inputs.shape[0], 13):
        raise ValueError(f"native bridge returned unexpected shape {rows.shape}")
    expected_steps = np.arange(inputs.shape[0], dtype=float)
    if not np.array_equal(rows[:, 0], expected_steps):
        raise ValueError("native bridge returned non-consecutive steps")
    states = np.vstack((initial, rows[:, 1:10]))
    source_root = Path(visit_source_root).resolve()
    source_files = {
        "soil_proc.c": source_root / "soil_proc.c",
        "decomposition.c": source_root / "decomposition.c",
        "bridge": NATIVE_BRIDGE_SOURCE,
    }
    return NativeVISITSoilTrajectory(
        states=states,
        heterotrophic_respiration=rows[:, 10],
        litter_scalars=rows[:, 11],
        humus_scalars=rows[:, 12],
        environments=np.asarray(
            [_environment_values(value) for value in environments], dtype=float
        ),
        executable=Path(executable).resolve(),
        source_commit=VISIT_SOURCE_COMMIT,
        source_sha256={name: _sha256(path) for name, path in source_files.items()},
        excluded_processes=(
            "nitrogen routines are linked as no-op stubs after the carbon update",
            "stable isotope branch is disabled by the pinned SCI_SCHEME=0",
        ),
    )


def compare_native_visit_soil_to_python(
    executable: str | Path,
    visit_source_root: str | Path,
    initial_state: ArrayLike,
    litter_inputs: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environment: VISITDecompositionEnvironment | Sequence[VISITDecompositionEnvironment],
) -> NativeVISITSoilComparison:
    """Compare native C and source-grounded Python baseline soil trajectories."""
    inputs = np.asarray(litter_inputs, dtype=float)
    native = run_native_visit_soil(
        executable,
        visit_source_root,
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        environment,
    )
    environments = _environment_sequence(environment, inputs.shape[0])
    python_states, python_respiration = _simulate_python_soil(
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        environments,
    )
    return NativeVISITSoilComparison(
        native=native,
        python_states=python_states,
        python_heterotrophic_respiration=python_respiration,
    )


def compare_native_visit_soil_perturbation_to_irf(
    executable: str | Path,
    visit_source_root: str | Path,
    baseline_litter_inputs: ArrayLike,
    input_perturbations: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environment: VISITDecompositionEnvironment,
) -> NativeVISITSoilPerturbationComparison:
    """Compare a native C perturbation experiment with the exact soil IRF."""
    baseline_inputs = _finite_array(
        baseline_litter_inputs, (6,), "baseline_litter_inputs"
    )
    perturbations = np.asarray(input_perturbations, dtype=float)
    if (
        perturbations.ndim != 2
        or perturbations.shape[1] != 6
        or not np.all(np.isfinite(perturbations))
    ):
        raise ValueError("input_perturbations must be finite with shape (n_steps, 6)")
    perturbed_inputs = baseline_inputs[np.newaxis, :] + perturbations
    if np.any(baseline_inputs < 0) or np.any(perturbed_inputs < 0):
        raise ValueError("baseline and perturbed litter inputs must be non-negative")
    scalars = visit_decomposition_scalars(
        environment,
        decomposition_parameters,
        DecompositionTemperatureMode.BASELINE,
    )
    fixed_point = visit_soil_fixed_point(
        soil_parameters, baseline_inputs, scalars.litter, scalars.humus
    )
    baseline_series = np.repeat(
        baseline_inputs[np.newaxis, :], perturbations.shape[0], axis=0
    )
    baseline = run_native_visit_soil(
        executable,
        visit_source_root,
        fixed_point.state,
        baseline_series,
        soil_parameters,
        decomposition_parameters,
        environment,
    )
    perturbed = run_native_visit_soil(
        executable,
        visit_source_root,
        fixed_point.state,
        perturbed_inputs,
        soil_parameters,
        decomposition_parameters,
        environment,
    )
    system = visit_soil_discrete_system(
        soil_parameters, scalars.litter, scalars.humus
    )
    state_difference, output_difference = system.forced_response(perturbations)
    return NativeVISITSoilPerturbationComparison(
        baseline=baseline,
        perturbed=perturbed,
        irf_state_difference=state_difference,
        irf_heterotrophic_respiration_difference=output_difference[:, 0],
    )


def linearize_visit_soil_environment_trajectory(
    initial_state: ArrayLike,
    litter_inputs: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environments: Sequence[VISITDecompositionEnvironment],
    *,
    relative_step: float = 1e-6,
    absolute_step: float = 1e-7,
) -> VISITSoilEnvironmentalLinearization:
    """Construct the daily tangent system for smooth baseline environments."""
    inputs = np.asarray(litter_inputs, dtype=float)
    if inputs.ndim != 2 or inputs.shape[1] != 6 or not np.all(np.isfinite(inputs)):
        raise ValueError("litter_inputs must be finite with shape (n_steps, 6)")
    daily_environments = _environment_sequence(environments, inputs.shape[0])
    if relative_step <= 0 or absolute_step <= 0:
        raise ValueError("finite-difference steps must be positive")
    baseline_states, baseline_respiration = _simulate_python_soil(
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        daily_environments,
    )
    n_steps = inputs.shape[0]
    state_jacobians = np.empty((n_steps, 9, 9), dtype=float)
    environment_state_jacobians = np.empty((n_steps, 9, 8), dtype=float)
    output_state_jacobians = np.empty((n_steps, 9), dtype=float)
    environment_output_jacobians = np.empty((n_steps, 8), dtype=float)
    for step, environment in enumerate(daily_environments):
        scalars = visit_decomposition_scalars(
            environment,
            decomposition_parameters,
            DecompositionTemperatureMode.BASELINE,
        )
        system = visit_soil_discrete_system(
            soil_parameters, scalars.litter, scalars.humus
        )
        state_jacobians[step] = system.A
        output_state_jacobians[step] = system.C[0]
        baseline_regime = visit_decomposition_regime(
            environment, decomposition_parameters
        )
        for column, name in enumerate(DECOMPOSITION_ENVIRONMENT_NAMES):
            value = float(getattr(environment, name))
            difference_step = max(absolute_step, relative_step * max(1.0, abs(value)))
            lower = replace(environment, **{name: value - difference_step})
            upper = replace(environment, **{name: value + difference_step})
            if (
                visit_decomposition_regime(lower, decomposition_parameters)
                != baseline_regime
                or visit_decomposition_regime(upper, decomposition_parameters)
                != baseline_regime
            ):
                raise ValueError(
                    f"finite difference for {name} crosses a nonsmooth decomposition branch"
                )
            lower_result = _python_soil_step(
                baseline_states[step],
                inputs[step],
                soil_parameters,
                decomposition_parameters,
                lower,
            )
            upper_result = _python_soil_step(
                baseline_states[step],
                inputs[step],
                soil_parameters,
                decomposition_parameters,
                upper,
            )
            environment_state_jacobians[step, :, column] = (
                upper_result.state - lower_result.state
            ) / (2.0 * difference_step)
            environment_output_jacobians[step, column] = (
                upper_result.fluxes.heterotrophic_respiration
                - lower_result.fluxes.heterotrophic_respiration
            ) / (2.0 * difference_step)
    return VISITSoilEnvironmentalLinearization(
        baseline_states=baseline_states,
        baseline_heterotrophic_respiration=baseline_respiration,
        state_jacobians=state_jacobians,
        environment_state_jacobians=environment_state_jacobians,
        output_state_jacobians=output_state_jacobians,
        environment_output_jacobians=environment_output_jacobians,
        environment_names=DECOMPOSITION_ENVIRONMENT_NAMES,
    )


def compare_native_visit_soil_environment_perturbation(
    executable: str | Path,
    visit_source_root: str | Path,
    initial_state: ArrayLike,
    litter_inputs: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    baseline_environments: Sequence[VISITDecompositionEnvironment],
    environment_perturbations: ArrayLike,
) -> NativeVISITSoilEnvironmentPerturbationComparison:
    """Compare native, nonlinear Python, and tangent environmental responses."""
    inputs = np.asarray(litter_inputs, dtype=float)
    environments = _environment_sequence(baseline_environments, inputs.shape[0])
    perturbations = np.asarray(environment_perturbations, dtype=float)
    expected_shape = (inputs.shape[0], len(DECOMPOSITION_ENVIRONMENT_NAMES))
    if perturbations.shape != expected_shape or not np.all(np.isfinite(perturbations)):
        raise ValueError(f"environment_perturbations must be finite with shape {expected_shape}")
    perturbed_environments = tuple(
        VISITDecompositionEnvironment(
            **{
                name: float(getattr(environment, name) + perturbation[index])
                for index, name in enumerate(DECOMPOSITION_ENVIRONMENT_NAMES)
            }
        )
        for environment, perturbation in zip(environments, perturbations, strict=True)
    )
    baseline_native = run_native_visit_soil(
        executable,
        visit_source_root,
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        environments,
    )
    perturbed_native = run_native_visit_soil(
        executable,
        visit_source_root,
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        perturbed_environments,
    )
    baseline_python_states, baseline_python_respiration = _simulate_python_soil(
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        environments,
    )
    perturbed_python_states, perturbed_python_respiration = _simulate_python_soil(
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        perturbed_environments,
    )
    linearization = linearize_visit_soil_environment_trajectory(
        initial_state,
        inputs,
        soil_parameters,
        decomposition_parameters,
        environments,
    )
    tangent_states, tangent_respiration = linearization.predict(perturbations)
    return NativeVISITSoilEnvironmentPerturbationComparison(
        baseline_native=baseline_native,
        perturbed_native=perturbed_native,
        baseline_python_states=baseline_python_states,
        perturbed_python_states=perturbed_python_states,
        baseline_python_heterotrophic_respiration=baseline_python_respiration,
        perturbed_python_heterotrophic_respiration=perturbed_python_respiration,
        tangent_state_difference=tangent_states,
        tangent_heterotrophic_respiration_difference=tangent_respiration,
    )


def _simulate_python_soil(
    initial_state: ArrayLike,
    inputs: np.ndarray,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environments: Sequence[VISITDecompositionEnvironment],
) -> tuple[np.ndarray, np.ndarray]:
    states = np.empty((inputs.shape[0] + 1, 9), dtype=float)
    respiration = np.empty(inputs.shape[0], dtype=float)
    states[0] = _finite_array(initial_state, (9,), "initial_state")
    for step, (input_values, environment) in enumerate(
        zip(inputs, environments, strict=True)
    ):
        result = _python_soil_step(
            states[step],
            input_values,
            soil_parameters,
            decomposition_parameters,
            environment,
        )
        states[step + 1] = result.state
        respiration[step] = result.fluxes.heterotrophic_respiration
    return states, respiration


def _python_soil_step(
    state: ArrayLike,
    litter_input: ArrayLike,
    soil_parameters: VISITSoilParameters,
    decomposition_parameters: VISITDecompositionParameters,
    environment: VISITDecompositionEnvironment,
):
    scalars = visit_decomposition_scalars(
        environment,
        decomposition_parameters,
        DecompositionTemperatureMode.BASELINE,
    )
    return visit_soil_daily_step(
        state,
        litter_input,
        soil_parameters,
        scalars.litter,
        scalars.humus,
    )




def _environment_sequence(
    environment: VISITDecompositionEnvironment | Sequence[VISITDecompositionEnvironment],
    n_steps: int,
) -> tuple[VISITDecompositionEnvironment, ...]:
    if isinstance(environment, VISITDecompositionEnvironment):
        return (environment,) * n_steps
    environments = tuple(environment)
    if len(environments) != n_steps:
        raise ValueError(f"environment sequence must contain {n_steps} daily values")
    if not all(isinstance(value, VISITDecompositionEnvironment) for value in environments):
        raise TypeError("environment sequence values must be VISITDecompositionEnvironment")
    return environments


def _environment_values(environment: VISITDecompositionEnvironment) -> tuple[float, ...]:
    return tuple(float(getattr(environment, name)) for name in DECOMPOSITION_ENVIRONMENT_NAMES)


def _finite_array(values: ArrayLike, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def _validate_native_revision(
    source_root: Path, relative_paths: tuple[str, ...]
) -> None:
    required = tuple(source_root / path for path in relative_paths)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing native VISIT sources: " + ", ".join(missing))
    revision = subprocess.run(
        ["git", "-C", str(source_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if revision != VISIT_SOURCE_COMMIT:
        raise ValueError(
            f"native VISIT source revision {revision} does not match {VISIT_SOURCE_COMMIT}"
        )
    source_diff = subprocess.run(
        ["git", "-C", str(source_root), "diff", "--quiet", "HEAD", "--", *relative_paths],
        check=False,
    )
    if source_diff.returncode != 0:
        raise ValueError("native VISIT bridge sources differ from the pinned commit")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()