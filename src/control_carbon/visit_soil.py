"""Source-grounded daily carbon update for the VISIT 9-pool soil subsystem."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import ArrayLike

from .provenance import SourceRef
from .visit_source_map import VISIT_SOURCE_COMMIT, VISIT_SOURCE_REPOSITORY


SOIL_POOL_NAMES: Final[tuple[str, ...]] = (
    "ltr_tf",
    "ltr_tc",
    "ltr_tr",
    "ltr_gf",
    "ltr_gc",
    "ltr_gr",
    "msl_a",
    "msl_i",
    "msl_p",
)
LITTER_INPUT_NAMES: Final[tuple[str, ...]] = (
    "li_tf",
    "li_tc",
    "li_tr",
    "li_gf",
    "li_gc",
    "li_gr",
)

VISIT_SOIL_PROCESS_SOURCE = SourceRef(
    path="visit_local/soil_proc.c",
    symbol="f_cycle_soil",
    role="daily litter and humus carbon mass balance",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="f_cycle_soil",
    units="Mg C ha^-1 and Mg C ha^-1 day^-1",
    assumptions=(
        "carbon pools only",
        "decomposition environmental scalars supplied as effective inputs",
        "stable-isotope and nitrogen updates excluded",
    ),
    approximation_level="native-source-carbon-subsystem",
)


@dataclass(frozen=True)
class VISITSoilParameters:
    """Source parameters from ``struct Schar`` used by ``f_cycle_soil``."""

    sr_lf: float
    sr_lc: float
    sr_lr: float
    sr_ha: float
    sr_hi: float
    sr_hp: float
    f_co2_lf: float
    f_co2_lc: float
    f_co2_lr: float
    f_hm_a: float
    f_hm_i: float
    f_hm_p: float

    def __post_init__(self) -> None:
        values = np.asarray(
            (
                self.sr_lf,
                self.sr_lc,
                self.sr_lr,
                self.sr_ha,
                self.sr_hi,
                self.sr_hp,
                self.f_co2_lf,
                self.f_co2_lc,
                self.f_co2_lr,
                self.f_hm_a,
                self.f_hm_i,
                self.f_hm_p,
            ),
            dtype=float,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("soil parameters must be finite")
        if np.any(values[:6] < 0):
            raise ValueError("decomposition coefficients must be non-negative")
        if np.any((values[6:] < 0) | (values[6:] > 1)):
            raise ValueError("partition fractions must lie in [0, 1]")

    @property
    def humification_total(self) -> float:
        """Return the source's litter-to-humus partition sum."""
        return self.f_hm_a + self.f_hm_i + self.f_hm_p

    def base_decomposition_coefficients(self) -> np.ndarray:
        """Return source ``sr_*`` values in the 9-pool ordering."""
        return np.asarray(
            (
                self.sr_lf,
                self.sr_lc,
                self.sr_lr,
                self.sr_lf,
                self.sr_lc,
                self.sr_lr,
                self.sr_ha,
                self.sr_hi,
                self.sr_hp,
            ),
            dtype=float,
        )

    def litter_co2_fractions(self) -> np.ndarray:
        """Return respiration fractions in the six-litter ordering."""
        return np.asarray(
            (
                self.f_co2_lf,
                self.f_co2_lc,
                self.f_co2_lr,
                self.f_co2_lf,
                self.f_co2_lc,
                self.f_co2_lr,
            ),
            dtype=float,
        )

    def humification_fractions(self) -> np.ndarray:
        """Return active/intermediate/passive humification fractions."""
        return np.asarray((self.f_hm_a, self.f_hm_i, self.f_hm_p), dtype=float)


@dataclass(frozen=True)
class VISITSoilFluxes:
    """Carbon fluxes calculated from the pre-update soil state."""

    degradation: np.ndarray
    microbial_respiration: np.ndarray
    humus_formation: np.ndarray
    heterotrophic_respiration: float
    unaccounted_carbon: float


@dataclass(frozen=True)
class VISITSoilStepResult:
    """Result of one native-order soil carbon update."""

    state: np.ndarray
    fluxes: VISITSoilFluxes
    clipped_pools: np.ndarray


def _validated_vector(values: ArrayLike, size: int, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},)")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    if np.any(array < 0):
        raise ValueError(f"{name} must be non-negative")
    return array


def soil_decomposition_fractions(
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
) -> np.ndarray:
    """Return daily degraded fractions ``sr / 1000 * f_tm`` by pool."""
    if not np.isfinite(litter_scalar) or litter_scalar < 0:
        raise ValueError("litter_scalar must be finite and non-negative")
    humus_scalars_array = _validated_vector(humus_scalars, 3, "humus_scalars")
    environmental_scalars = np.concatenate(
        (np.full(6, litter_scalar, dtype=float), humus_scalars_array)
    )
    return parameters.base_decomposition_coefficients() / 1000.0 * environmental_scalars


def visit_soil_daily_matrices(
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the affine daily update before native negative-pool clipping.

    The returned matrices satisfy ``x[k+1] = A x[k] + B u[k]`` whenever
    ``f_cycle_soil`` would not activate its post-update negative-pool clipping.
    """
    decomposition = soil_decomposition_fractions(
        parameters, litter_scalar, humus_scalars
    )
    matrix_a = np.diag(1.0 - decomposition)
    co2_fractions = parameters.litter_co2_fractions()
    humification = parameters.humification_fractions()
    matrix_a[6:, :6] = np.outer(
        humification, decomposition[:6] * (1.0 - co2_fractions)
    )
    matrix_b = np.zeros((9, 6), dtype=float)
    matrix_b[:6] = np.eye(6)
    return matrix_a, matrix_b


def visit_soil_daily_step(
    state: ArrayLike,
    litter_inputs: ArrayLike,
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
    *,
    clip_negative: bool = True,
) -> VISITSoilStepResult:
    """Reproduce the carbon-pool algebra and update order of ``f_cycle_soil``."""
    state_array = _validated_vector(state, 9, "state")
    litter_inputs_array = _validated_vector(litter_inputs, 6, "litter_inputs")
    decomposition_fractions = soil_decomposition_fractions(
        parameters, litter_scalar, humus_scalars
    )
    degradation = state_array * decomposition_fractions
    co2_fractions = parameters.litter_co2_fractions()
    litter_respiration = degradation[:6] * co2_fractions
    humus_formation = np.outer(
        parameters.humification_fractions(),
        degradation[:6] * (1.0 - co2_fractions),
    )
    microbial_respiration = np.concatenate((litter_respiration, degradation[6:]))

    updated_state = np.empty(9, dtype=float)
    updated_state[:6] = state_array[:6] + litter_inputs_array - degradation[:6]
    updated_state[6:] = state_array[6:] + humus_formation.sum(axis=1) - degradation[6:]
    carbon_before_clipping = float(updated_state.sum())
    clipped_pools = updated_state < 0.0
    if clip_negative:
        updated_state = np.maximum(updated_state, 0.0)

    heterotrophic_respiration = float(microbial_respiration.sum())
    expected_carbon = float(
        state_array.sum() + litter_inputs_array.sum() - heterotrophic_respiration
    )
    fluxes = VISITSoilFluxes(
        degradation=degradation,
        microbial_respiration=microbial_respiration,
        humus_formation=humus_formation,
        heterotrophic_respiration=heterotrophic_respiration,
        unaccounted_carbon=expected_carbon - carbon_before_clipping,
    )
    return VISITSoilStepResult(
        state=updated_state,
        fluxes=fluxes,
        clipped_pools=clipped_pools,
    )