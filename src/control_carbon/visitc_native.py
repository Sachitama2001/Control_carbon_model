"""Native-C validation bridges for current VISITc hydrology and PM algebra."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import subprocess
from typing import Final

import numpy as np

from .visitc_hydrology import (
    VISITCHydrologyForcing,
    VISITCHydrologyParameters,
    VISITCHydrologyState,
    VISITCHydrologyStep,
    VISITCPenmanMonteithEnvironment,
    VISITCPenmanMonteithFluxes,
    VISITCCanopyConductanceParameters,
    visitc_hydrology_step,
    visitc_penman_monteith_fluxes,
    visitc_canopy_conductance,
    visitc_irradiance_extinction,
    visitc_leaf_area_index,
)
from .visitc_source_map import VISITC_SOURCE_COMMIT


NATIVE_VISITC_HYDROLOGY_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2] / "native" / "visitc_hydrology_bridge.c"
)
NATIVE_VISITC_PM_BRIDGE_SOURCE: Final = (
    Path(__file__).resolve().parents[2] / "native" / "visitc_pm_bridge.c"
)


@dataclass(frozen=True)
class NativeVISITCHydrologyResult:
    state: VISITCHydrologyState
    diagnostics: np.ndarray
    executable: Path
    source_commit: str
    source_sha256: dict[str, str]

    @property
    def state_vector(self) -> np.ndarray:
        return np.asarray(
            (
                self.state.snow_water,
                self.state.upper_soil_water,
                self.state.whole_soil_water,
            )
        )


@dataclass(frozen=True)
class NativeVISITCHydrologyComparison:
    native: NativeVISITCHydrologyResult
    python: VISITCHydrologyStep
    native_vector: np.ndarray
    python_vector: np.ndarray

    @property
    def error(self) -> np.ndarray:
        return self.native_vector - self.python_vector


@dataclass(frozen=True)
class NativeVISITCPMComparison:
    native_pm_vector: np.ndarray
    python_pm_vector: np.ndarray
    native_radiation_vector: np.ndarray
    python_radiation_vector: np.ndarray
    python: VISITCPenmanMonteithFluxes
    executable: Path
    source_commit: str
    source_sha256: dict[str, str]

    @property
    def pm_error(self) -> np.ndarray:
        return self.native_pm_vector - self.python_pm_vector

    @property
    def radiation_error(self) -> np.ndarray:
        return self.native_radiation_vector - self.python_radiation_vector


def build_native_visitc_hydrology_bridge(
    visitc_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile ``f_hydrology`` with only its PM calls replaced by inputs."""
    source_root = Path(visitc_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        (
            "hydro_balance.c",
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
        str(NATIVE_VISITC_HYDROLOGY_BRIDGE_SOURCE),
        str(source_root / "hydro_balance.c"),
        "-lm",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def build_native_visitc_pm_bridge(
    visitc_source_root: str | Path,
    output_path: str | Path,
    *,
    compiler: str = "cc",
) -> Path:
    """Compile native atmospheric, radiation, and PM functions directly."""
    source_root = Path(visitc_source_root).resolve()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    _validate_native_revision(
        source_root,
        (
            "hydro_flows.c",
            "radiation.c",
            "ecophysiology.c",
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
        str(NATIVE_VISITC_PM_BRIDGE_SOURCE),
        str(source_root / "hydro_flows.c"),
        str(source_root / "radiation.c"),
        str(source_root / "ecophysiology.c"),
        "-Wl,--gc-sections",
        "-lm",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def run_native_visitc_hydrology(
    executable: str | Path,
    source_root: str | Path,
    state: VISITCHydrologyState,
    forcing: VISITCHydrologyForcing,
    parameters: VISITCHydrologyParameters,
) -> NativeVISITCHydrologyResult:
    """Execute one native ``f_hydrology`` day with effective PM inputs."""
    values: tuple[object, ...] = (
        state.snow_water,
        state.upper_soil_water,
        state.whole_soil_water,
        forcing.precipitation,
        forcing.air_temperature,
        forcing.deep_soil_temperature,
        forcing.potential_interception_tree,
        forcing.potential_interception_c3,
        forcing.potential_interception_c4,
        forcing.potential_soil_evaporation,
        forcing.potential_transpiration_tree,
        forcing.potential_transpiration_c3,
        forcing.potential_transpiration_c4,
        parameters.field_capacity_upper,
        parameters.field_capacity_whole,
        parameters.hydraulic_conductivity,
        parameters.tree_leaf_area_index,
        parameters.c3_leaf_area_index,
        parameters.c4_leaf_area_index,
        parameters.c3_understory_fraction,
        parameters.c4_understory_fraction,
        parameters.site_id or "_",
    )
    output = _run_csv(executable, (), values, 18)
    root = Path(source_root).resolve()
    return NativeVISITCHydrologyResult(
        state=VISITCHydrologyState(*output[:3]),
        diagnostics=output[3:],
        executable=Path(executable).resolve(),
        source_commit=VISITC_SOURCE_COMMIT,
        source_sha256={
            name: _sha256(root / name)
            for name in ("hydro_balance.c", "structure.h", "setting.h")
        },
    )


def compare_native_visitc_hydrology_to_python(
    executable: str | Path,
    source_root: str | Path,
    state: VISITCHydrologyState,
    forcing: VISITCHydrologyForcing,
    parameters: VISITCHydrologyParameters,
) -> NativeVISITCHydrologyComparison:
    """Compare all exposed one-day native stores and flux diagnostics."""
    native = run_native_visitc_hydrology(
        executable, source_root, state, forcing, parameters
    )
    python = visitc_hydrology_step(state, forcing, parameters)
    native_vector = np.concatenate((native.state_vector, native.diagnostics))
    python_vector = _python_hydrology_vector(python)
    return NativeVISITCHydrologyComparison(
        native=native,
        python=python,
        native_vector=native_vector,
        python_vector=python_vector,
    )


def compare_native_visitc_pm_to_python(
    executable: str | Path,
    source_root: str | Path,
    environment: VISITCPenmanMonteithEnvironment,
    *,
    upper_soil_water: float,
    field_capacity_upper: float,
) -> NativeVISITCPMComparison:
    """Compare native radiation/atmospheric/PM functions with Python."""
    python = visitc_penman_monteith_fluxes(
        environment,
        upper_soil_water=upper_soil_water,
        field_capacity_upper=field_capacity_upper,
    )
    radiation = python.net_radiation
    native_pm = _run_csv(
        executable,
        ("pm",),
        (
            environment.air_temperature,
            environment.air_pressure,
            environment.vapor_pressure,
            environment.vapor_pressure_deficit,
            environment.wind_speed,
            environment.day_length,
            radiation.tree,
            radiation.c3,
            radiation.c4,
            radiation.ground,
            environment.canopy_conductance_tree,
            environment.canopy_conductance_c3,
            environment.canopy_conductance_c4,
            upper_soil_water,
            field_capacity_upper,
        ),
        12,
    )
    p = environment.radiation
    native_radiation = _run_csv(
        executable,
        ("radiation",),
        (
            environment.air_temperature,
            environment.surface_temperature,
            environment.vapor_pressure,
            environment.cloud_fraction,
            environment.incoming_shortwave,
            *p.leaf_area_index,
            *p.extinction_initial,
            *p.extinction_radiation,
            *p.albedo,
            p.c3_understory_fraction,
            p.c4_understory_fraction,
        ),
        17,
    )
    python_pm = np.asarray(
        (
            python.saturated_vapor_pressure,
            python.saturation_vapor_pressure_slope,
            python.air_density,
            python.aerodynamic_resistance,
            python.potential_soil_evaporation,
            python.potential_interception_tree,
            python.potential_interception_c3,
            python.potential_interception_c4,
            python.potential_transpiration_tree,
            python.potential_transpiration_c3,
            python.potential_transpiration_c4,
            python.soil_resistance,
        )
    )
    python_radiation = np.asarray(
        (
            radiation.tree,
            radiation.c3,
            radiation.c4,
            radiation.ground,
            radiation.ecosystem,
            *radiation.cover_fractions,
            *radiation.absorbed_shortwave,
            *radiation.partitioned_longwave,
        )
    )
    root = Path(source_root).resolve()
    return NativeVISITCPMComparison(
        native_pm_vector=native_pm,
        python_pm_vector=python_pm,
        native_radiation_vector=native_radiation,
        python_radiation_vector=python_radiation,
        python=python,
        executable=Path(executable).resolve(),
        source_commit=VISITC_SOURCE_COMMIT,
        source_sha256={
            name: _sha256(root / name)
            for name in ("hydro_flows.c", "radiation.c", "definition.h")
        },
    )


def compare_native_visitc_canopy_conductance_to_python(
    executable: str | Path,
    parameters: VISITCCanopyConductanceParameters,
) -> float:
    """Return native-minus-Python error for ``f_canopy_cond``."""
    native = _run_csv(
        executable,
        ("conductance",),
        tuple(parameters.__dict__.values()),
        1,
    )[0]
    python = visitc_canopy_conductance(parameters).conductance
    return float(native - python)


def compare_native_visitc_lai_to_python(
    executable: str | Path, foliage_carbon: float, specific_leaf_area: float
) -> float:
    """Return native-minus-Python error for ``lai_mass``."""
    native = _run_csv(
        executable, ("lai",), (foliage_carbon, specific_leaf_area), 1
    )[0]
    return float(native - visitc_leaf_area_index(foliage_carbon, specific_leaf_area))


def compare_native_visitc_extinction_to_python(
    executable: str | Path, initial_extinction: float, solar_height_degrees: float
) -> float:
    """Return native-minus-Python error for ``irr_attn``."""
    native = _run_csv(
        executable,
        ("extinction",),
        (initial_extinction, solar_height_degrees),
        1,
    )[0]
    return float(
        native
        - visitc_irradiance_extinction(initial_extinction, solar_height_degrees)
    )


def _python_hydrology_vector(step: VISITCHydrologyStep) -> np.ndarray:
    f = step.fluxes
    return np.asarray(
        (
            step.state.snow_water,
            step.state.upper_soil_water,
            step.state.whole_soil_water,
            f.snow_fraction,
            f.thaw,
            f.interception_tree,
            f.interception_c3,
            f.interception_c4,
            f.interception,
            f.upper_runoff,
            f.soil_evaporation,
            f.transpiration_tree,
            f.transpiration_c3,
            f.transpiration_c4,
            f.transpiration,
            f.lower_runoff,
            f.actual_evapotranspiration,
            f.potential_evapotranspiration,
        )
    )


def _run_csv(
    executable: str | Path,
    arguments: tuple[str, ...],
    values: tuple[object, ...],
    expected_size: int,
) -> np.ndarray:
    completed = subprocess.run(
        [str(Path(executable).resolve()), *arguments],
        input=" ".join(str(value) for value in values) + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    result = np.asarray(
        [float(value) for value in completed.stdout.strip().split(",")], dtype=float
    )
    if result.shape != (expected_size,):
        raise ValueError(f"native VISITc bridge returned shape {result.shape}")
    return result


def _validate_native_revision(
    source_root: Path, relative_paths: tuple[str, ...]
) -> None:
    missing = [str(source_root / path) for path in relative_paths if not (source_root / path).is_file()]
    if missing:
        raise FileNotFoundError("missing native VISITc sources: " + ", ".join(missing))
    revision = subprocess.run(
        ["git", "-C", str(source_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if revision != VISITC_SOURCE_COMMIT:
        raise ValueError(
            f"native VISITc source revision {revision} does not match {VISITC_SOURCE_COMMIT}"
        )
    source_diff = subprocess.run(
        ["git", "-C", str(source_root), "diff", "--quiet", "HEAD", "--", *relative_paths],
        check=False,
    )
    if source_diff.returncode != 0:
        raise ValueError("native VISITc bridge sources differ from the pinned commit")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


__all__ = [
    "NativeVISITCHydrologyComparison",
    "NativeVISITCHydrologyResult",
    "NativeVISITCPMComparison",
    "build_native_visitc_hydrology_bridge",
    "build_native_visitc_pm_bridge",
    "compare_native_visitc_hydrology_to_python",
    "compare_native_visitc_canopy_conductance_to_python",
    "compare_native_visitc_extinction_to_python",
    "compare_native_visitc_lai_to_python",
    "compare_native_visitc_pm_to_python",
    "run_native_visitc_hydrology",
]
