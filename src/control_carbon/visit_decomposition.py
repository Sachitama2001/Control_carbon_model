"""Environmental decomposition scalars reconstructed from VISIT source.

The source globals ``SA_PARA`` and ``SA_PARA_VAR`` are represented here as
explicit sensitivity-multiplier arguments. Callers can therefore reproduce a
time-varying sensitivity experiment by supplying the applicable value each day,
without hidden mutable global state.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, IntEnum
from typing import Final

import numpy as np

from .provenance import SourceRef
from .visit_source_map import VISIT_SOURCE_COMMIT, VISIT_SOURCE_REPOSITORY


DECOMPOSITION_ENVIRONMENT_NAMES: Final[tuple[str, ...]] = (
    "soil_temperature_10cm",
    "soil_temperature_deep",
    "soil_water_upper",
    "soil_water_lower",
    "soil_aperture_upper",
    "soil_aperture_whole",
    "field_capacity_30cm",
    "field_capacity",
)


class DecompositionTemperatureMode(IntEnum):
    """Meanings of the source macro ``EX_DECTMP``."""

    BASELINE = 0
    WARMING = 1
    POOL_SPECIFIC = 2
    WARMING_POOL_SPECIFIC = 3

    @property
    def uses_pool_specific_response(self) -> bool:
        return self in (self.POOL_SPECIFIC, self.WARMING_POOL_SPECIFIC)


class HumusBranchPolicy(Enum):
    """Policy for the unreachable ``frh`` branches in the source snapshot."""

    SOURCE = "source"
    INFERRED_STYPE_FIX = "inferred-stype-fix"


class VISITSourceUndefinedError(RuntimeError):
    """Raised when the selected native source path has undefined behavior."""


@dataclass(frozen=True)
class VISITDecompositionEnvironment:
    """Daily environmental values consumed by ``frl`` and ``frh``."""

    soil_temperature_10cm: float
    soil_temperature_deep: float
    soil_water_upper: float
    soil_water_lower: float
    soil_aperture_upper: float
    soil_aperture_whole: float
    field_capacity_30cm: float
    field_capacity: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("decomposition environment values must be finite")
        if np.any(values[2:] < 0):
            raise ValueError("soil water, aperture, and field capacity must be non-negative")
        if self.field_capacity_30cm == 0 or self.field_capacity == 0:
            raise ValueError("field capacities must be positive")


@dataclass(frozen=True)
class VISITDecompositionParameters:
    """Moisture and aperture coefficients from ``struct Schar``."""

    kml: float
    kmh: float
    kmsl: float
    kmsh: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("decomposition response parameters must be finite and positive")


@dataclass(frozen=True)
class VISITDecompositionScalars:
    """Effective daily multipliers passed to the 9-pool soil update."""

    litter: float
    humus: np.ndarray
    mode: DecompositionTemperatureMode
    humus_branch_policy: HumusBranchPolicy


def visit_decomposition_regime(
    environment: VISITDecompositionEnvironment,
    parameters: VISITDecompositionParameters,
) -> tuple[bool, bool, str, str]:
    """Identify temperature-threshold and moisture-limiting source branches."""
    def limiting_branch(
        water: float,
        half_saturation: float,
        field_capacity: float,
        aperture: float,
        aperture_half_saturation: float,
    ) -> str:
        water_factor = 0.8 * water / (half_saturation * field_capacity + water) + 0.2
        aperture_factor = (
            0.4
            * aperture
            * aperture_half_saturation
            / (aperture_half_saturation + aperture)
            + 0.6
        )
        return "water" if water_factor < aperture_factor else "aperture"

    return (
        environment.soil_temperature_10cm > -20.0,
        environment.soil_temperature_deep > -20.0,
        limiting_branch(
            environment.soil_water_upper,
            parameters.kml,
            environment.field_capacity_30cm,
            environment.soil_aperture_upper,
            parameters.kmsl,
        ),
        limiting_branch(
            environment.soil_water_lower,
            parameters.kmh,
            environment.field_capacity,
            environment.soil_aperture_whole,
            parameters.kmsh,
        ),
    )


VISIT_STYPE_CALL_SOURCE = SourceRef(
    path="visit_local/soil_proc.c",
    symbol="frl(..., stype) / frh(..., stype)",
    role="select decomposition temperature-response formula by call-site integer",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="f_cycle_soil",
    assumptions=(
        "stype is a call-site selector, not a parameter-file or forcing value",
        "EX_DECTMP 0/1 passes 0 to both functions",
        "EX_DECTMP 2/3 passes 1 to frl and 0/1/2 to active/intermediate/passive frh",
    ),
    approximation_level="native-source",
)

VISIT_DECOMPOSITION_SOURCE = SourceRef(
    path="visit_local/decomposition.c",
    symbol="frl / frh",
    role="temperature, moisture, and soil-aperture decomposition multipliers",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="frl / frh",
    assumptions=(
        "SA_PARA temperature and water multipliers are supplied explicitly",
        "frh stype 1/2 is undefined in the source snapshot because all branches test stype == 0",
    ),
    approximation_level="native-source-with-explicit-undefined-branch",
)


def _response_multiplier(value: float, name: str) -> float:
    multiplier = float(value)
    if not np.isfinite(multiplier) or multiplier <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return multiplier


def _lloyd_taylor_factor(
    temperature: float,
    reference_offset: float,
    temperature_offset: float,
    temperature_sensitivity_multiplier: float,
) -> float:
    """Evaluate the Lloyd-Taylor form used by ``decomposition.c``."""
    if temperature <= -20.0:
        return 0.01
    return float(
        0.01
        + np.exp(
            308.56
            * (
                1.0 / reference_offset
                - 1.0
                / (
                    temperature
                    + temperature_offset * temperature_sensitivity_multiplier
                )
            )
        )
    )


def _moisture_aperture_factor(
    water: float,
    half_saturation: float,
    field_capacity: float,
    aperture: float,
    aperture_half_saturation: float,
    soil_water_sensitivity_multiplier: float,
) -> float:
    """Evaluate the source moisture limit and soil-aperture limit."""
    water_denominator = (
        half_saturation * soil_water_sensitivity_multiplier * field_capacity + water
    )
    water_factor = 0.8 * water / water_denominator + 0.2
    aperture_factor = (
        0.4
        * aperture
        * aperture_half_saturation
        / (aperture_half_saturation + aperture)
        + 0.6
    )
    return float(min(water_factor, aperture_factor))


def visit_litter_decomposition_scalar(
    environment: VISITDecompositionEnvironment,
    parameters: VISITDecompositionParameters,
    mode: DecompositionTemperatureMode,
    *,
    temperature_sensitivity_multiplier: float = 1.0,
    soil_water_sensitivity_multiplier: float = 1.0,
) -> float:
    """Reproduce ``frl`` with ``stype`` selected as in ``f_cycle_soil``."""
    temperature_multiplier = _response_multiplier(
        temperature_sensitivity_multiplier, "temperature_sensitivity_multiplier"
    )
    water_multiplier = _response_multiplier(
        soil_water_sensitivity_multiplier, "soil_water_sensitivity_multiplier"
    )
    mode = DecompositionTemperatureMode(mode)
    if mode.uses_pool_specific_response:
        reference_offset, temperature_offset = 66.02, 56.02
    else:
        reference_offset, temperature_offset = 56.02, 46.02
    temperature_factor = _lloyd_taylor_factor(
        environment.soil_temperature_10cm,
        reference_offset,
        temperature_offset,
        temperature_multiplier,
    )
    moisture_factor = _moisture_aperture_factor(
        environment.soil_water_upper,
        parameters.kml,
        environment.field_capacity_30cm,
        environment.soil_aperture_upper,
        parameters.kmsl,
        water_multiplier,
    )
    return temperature_factor * moisture_factor


def visit_humus_decomposition_scalars(
    environment: VISITDecompositionEnvironment,
    parameters: VISITDecompositionParameters,
    mode: DecompositionTemperatureMode,
    *,
    branch_policy: HumusBranchPolicy = HumusBranchPolicy.SOURCE,
    temperature_sensitivity_multiplier: float = 1.0,
    soil_water_sensitivity_multiplier: float = 1.0,
) -> np.ndarray:
    """Reproduce ``frh`` or apply the explicitly selected inferred branch fix.

    In source-faithful mode, ``EX_DECTMP`` 2/3 raises because calls with
    ``stype`` 1 and 2 leave the C local variable ``fth`` uninitialized.
    """
    temperature_multiplier = _response_multiplier(
        temperature_sensitivity_multiplier, "temperature_sensitivity_multiplier"
    )
    water_multiplier = _response_multiplier(
        soil_water_sensitivity_multiplier, "soil_water_sensitivity_multiplier"
    )
    mode = DecompositionTemperatureMode(mode)
    branch_policy = HumusBranchPolicy(branch_policy)
    if mode.uses_pool_specific_response and branch_policy is HumusBranchPolicy.SOURCE:
        raise VISITSourceUndefinedError(
            "VISIT decomposition.c::frh leaves fth uninitialized for stype 1/2 "
            "under EX_DECTMP 2/3"
        )

    moisture_factor = _moisture_aperture_factor(
        environment.soil_water_lower,
        parameters.kmh,
        environment.field_capacity,
        environment.soil_aperture_whole,
        parameters.kmsh,
        water_multiplier,
    )
    response_offsets = ((56.02, 46.02),) * 3
    if mode.uses_pool_specific_response:
        response_offsets = ((56.02, 46.02), (50.02, 40.02), (44.02, 34.02))
    return np.asarray(
        [
            _lloyd_taylor_factor(
                environment.soil_temperature_deep,
                reference_offset,
                temperature_offset,
                temperature_multiplier,
            )
            * moisture_factor
            for reference_offset, temperature_offset in response_offsets
        ],
        dtype=float,
    )


def visit_decomposition_scalars(
    environment: VISITDecompositionEnvironment,
    parameters: VISITDecompositionParameters,
    mode: DecompositionTemperatureMode = DecompositionTemperatureMode.BASELINE,
    *,
    humus_branch_policy: HumusBranchPolicy = HumusBranchPolicy.SOURCE,
    temperature_sensitivity_multiplier: float = 1.0,
    soil_water_sensitivity_multiplier: float = 1.0,
) -> VISITDecompositionScalars:
    """Calculate all effective decomposition multipliers for one day."""
    mode = DecompositionTemperatureMode(mode)
    humus_branch_policy = HumusBranchPolicy(humus_branch_policy)
    litter = visit_litter_decomposition_scalar(
        environment,
        parameters,
        mode,
        temperature_sensitivity_multiplier=temperature_sensitivity_multiplier,
        soil_water_sensitivity_multiplier=soil_water_sensitivity_multiplier,
    )
    humus = visit_humus_decomposition_scalars(
        environment,
        parameters,
        mode,
        branch_policy=humus_branch_policy,
        temperature_sensitivity_multiplier=temperature_sensitivity_multiplier,
        soil_water_sensitivity_multiplier=soil_water_sensitivity_multiplier,
    )
    return VISITDecompositionScalars(
        litter=litter,
        humus=humus,
        mode=mode,
        humus_branch_policy=humus_branch_policy,
    )