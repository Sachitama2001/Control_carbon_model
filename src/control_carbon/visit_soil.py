"""Source-grounded daily carbon update for the VISIT 9-pool soil subsystem.

Native daily timing is preserved: degradation and respiration use the stocks
at the start of the day, and the current day's litter inputs are added during
the subsequent mass-balance update. Consequently a litter pulse first affects
heterotrophic respiration on the following day.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Final

import numpy as np
from numpy.typing import ArrayLike

from .provenance import SourceRef
from .state_space import DiscreteFixedPointResult, DiscreteLTI, discrete_fixed_point
from .visit_decomposition import (
    DECOMPOSITION_ENVIRONMENT_NAMES,
    DecompositionTemperatureMode,
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
    visit_decomposition_regime,
    visit_decomposition_scalars,
)
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


@dataclass(frozen=True)
class VISITSoilTrajectoryComparison:
    """IRF prediction and direct daily simulation under the same perturbation."""

    days: np.ndarray
    baseline_states: np.ndarray
    perturbed_states: np.ndarray
    irf_predicted_states: np.ndarray
    baseline_heterotrophic_respiration: np.ndarray
    perturbed_heterotrophic_respiration: np.ndarray
    irf_predicted_heterotrophic_respiration: np.ndarray
    phase_axis_names: tuple[str, str]
    baseline_phase: np.ndarray
    perturbed_phase: np.ndarray
    irf_predicted_phase: np.ndarray

    @property
    def state_error(self) -> np.ndarray:
        """Direct simulation minus the IRF-based state prediction."""
        return self.perturbed_states - self.irf_predicted_states

    @property
    def heterotrophic_respiration_error(self) -> np.ndarray:
        """Direct simulation minus the IRF-based respiration prediction."""
        return (
            self.perturbed_heterotrophic_respiration
            - self.irf_predicted_heterotrophic_respiration
        )

    @property
    def output_days(self) -> np.ndarray:
        """Days corresponding to same-index daily respiration outputs."""
        return self.days[:-1]


@dataclass(frozen=True)
class VISITSoilEnvironmentalSensitivity:
    """Implicit and direct finite-difference fixed-point sensitivities."""

    fixed_point: DiscreteFixedPointResult
    environment_names: tuple[str, ...]
    state_sensitivity: np.ndarray
    heterotrophic_respiration_sensitivity: np.ndarray
    finite_difference_state_sensitivity: np.ndarray
    finite_difference_heterotrophic_respiration_sensitivity: np.ndarray
    finite_difference_steps: np.ndarray
    decomposition_regime: tuple[bool, bool, str, str]

    @property
    def state_sensitivity_error(self) -> np.ndarray:
        return self.state_sensitivity - self.finite_difference_state_sensitivity

    @property
    def heterotrophic_respiration_sensitivity_error(self) -> np.ndarray:
        return (
            self.heterotrophic_respiration_sensitivity
            - self.finite_difference_heterotrophic_respiration_sensitivity
        )


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
    Respiration at index ``k`` is calculated from ``x[k]`` before ``u[k]`` is
    added to the litter stocks.
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


def visit_soil_discrete_system(
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
) -> DiscreteLTI:
    """Build the daily soil system with heterotrophic respiration as output."""
    matrix_a, matrix_b = visit_soil_daily_matrices(
        parameters, litter_scalar, humus_scalars
    )
    decomposition = soil_decomposition_fractions(
        parameters, litter_scalar, humus_scalars
    )
    respiration_coefficients = np.concatenate(
        (
            decomposition[:6] * parameters.litter_co2_fractions(),
            decomposition[6:],
        )
    )
    return DiscreteLTI(
        A=matrix_a,
        B=matrix_b,
        C=respiration_coefficients[np.newaxis, :],
        D=np.zeros((1, 6), dtype=float),
        dt=1.0,
    )


def visit_soil_fixed_point(
    parameters: VISITSoilParameters,
    litter_inputs: ArrayLike,
    litter_scalar: float,
    humus_scalars: ArrayLike,
) -> DiscreteFixedPointResult:
    """Solve the constant-environment daily soil carbon fixed point."""
    litter_inputs_array = _validated_vector(litter_inputs, 6, "litter_inputs")
    system = visit_soil_discrete_system(parameters, litter_scalar, humus_scalars)
    result = discrete_fixed_point(system, litter_inputs_array)
    if np.any(result.state < 0):
        raise ValueError("soil fixed point contains negative carbon pools")
    return result


def visit_soil_environmental_sensitivity(
    parameters: VISITSoilParameters,
    litter_inputs: ArrayLike,
    decomposition_parameters: VISITDecompositionParameters,
    environment: VISITDecompositionEnvironment,
    *,
    relative_step: float = 1e-5,
    absolute_step: float = 1e-6,
) -> VISITSoilEnvironmentalSensitivity:
    """Differentiate the fixed point and equilibrium Rh by environment.

    The implicit discrete relation is
    ``(I-A) dx*/de = (dA/de) x*`` because the litter-input matrix is constant.
    Direct finite differences of independently solved fixed points are returned
    as a numerical validation reference.
    """
    litter_input_array = _validated_vector(litter_inputs, 6, "litter_inputs")
    if relative_step <= 0 or absolute_step <= 0:
        raise ValueError("finite-difference steps must be positive")

    def system_for(value: VISITDecompositionEnvironment) -> DiscreteLTI:
        scalars = visit_decomposition_scalars(
            value,
            decomposition_parameters,
            DecompositionTemperatureMode.BASELINE,
        )
        return visit_soil_discrete_system(
            parameters, scalars.litter, scalars.humus
        )

    system = system_for(environment)
    fixed_point = discrete_fixed_point(system, litter_input_array)
    fixed_point_matrix = np.eye(system.n_state) - system.A
    baseline_regime = visit_decomposition_regime(
        environment, decomposition_parameters
    )
    n_environment = len(DECOMPOSITION_ENVIRONMENT_NAMES)
    state_sensitivity = np.empty((system.n_state, n_environment), dtype=float)
    respiration_sensitivity = np.empty(n_environment, dtype=float)
    finite_state_sensitivity = np.empty_like(state_sensitivity)
    finite_respiration_sensitivity = np.empty_like(respiration_sensitivity)
    difference_steps = np.empty(n_environment, dtype=float)
    for column, name in enumerate(DECOMPOSITION_ENVIRONMENT_NAMES):
        value = float(getattr(environment, name))
        difference_step = max(absolute_step, relative_step * max(1.0, abs(value)))
        lower_environment = replace(
            environment, **{name: value - difference_step}
        )
        upper_environment = replace(
            environment, **{name: value + difference_step}
        )
        if (
            visit_decomposition_regime(lower_environment, decomposition_parameters)
            != baseline_regime
            or visit_decomposition_regime(
                upper_environment, decomposition_parameters
            )
            != baseline_regime
        ):
            raise ValueError(
                f"finite difference for {name} crosses a nonsmooth decomposition branch"
            )
        lower_system = system_for(lower_environment)
        upper_system = system_for(upper_environment)
        derivative_a = (upper_system.A - lower_system.A) / (2.0 * difference_step)
        derivative_c = (upper_system.C - lower_system.C) / (2.0 * difference_step)
        derivative_state = np.linalg.solve(
            fixed_point_matrix, derivative_a @ fixed_point.state
        )
        state_sensitivity[:, column] = derivative_state
        respiration_sensitivity[column] = (
            derivative_c @ fixed_point.state + system.C @ derivative_state
        )[0]

        lower_fixed_point = discrete_fixed_point(lower_system, litter_input_array)
        upper_fixed_point = discrete_fixed_point(upper_system, litter_input_array)
        finite_state_sensitivity[:, column] = (
            upper_fixed_point.state - lower_fixed_point.state
        ) / (2.0 * difference_step)
        lower_respiration = (lower_system.C @ lower_fixed_point.state)[0]
        upper_respiration = (upper_system.C @ upper_fixed_point.state)[0]
        finite_respiration_sensitivity[column] = (
            upper_respiration - lower_respiration
        ) / (2.0 * difference_step)
        difference_steps[column] = difference_step
    return VISITSoilEnvironmentalSensitivity(
        fixed_point=fixed_point,
        environment_names=DECOMPOSITION_ENVIRONMENT_NAMES,
        state_sensitivity=state_sensitivity,
        heterotrophic_respiration_sensitivity=respiration_sensitivity,
        finite_difference_state_sensitivity=finite_state_sensitivity,
        finite_difference_heterotrophic_respiration_sensitivity=(
            finite_respiration_sensitivity
        ),
        finite_difference_steps=difference_steps,
        decomposition_regime=baseline_regime,
    )


def visit_soil_heterotrophic_respiration_irf(
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
    n_steps: int,
) -> np.ndarray:
    """Return daily respiration responses to unit litter pulses.

    Columns follow ``LITTER_INPUT_NAMES``. Row zero is zero because native
    ``f_cycle_soil`` decomposes the pre-input litter stocks before adding the
    current day's litter inputs.
    """
    system = visit_soil_discrete_system(parameters, litter_scalar, humus_scalars)
    return system.impulse_response(n_steps)[:, 0, :]


def visit_soil_litter_humus_projection(states: ArrayLike) -> np.ndarray:
    """Project 9-pool states onto total litter and total humus carbon."""
    state_array = np.asarray(states, dtype=float)
    if state_array.ndim != 2 or state_array.shape[1] != 9:
        raise ValueError("states must have shape (n_steps, 9)")
    if not np.all(np.isfinite(state_array)):
        raise ValueError("states must be finite")
    return np.column_stack((state_array[:, :6].sum(axis=1), state_array[:, 6:].sum(axis=1)))


def compare_visit_soil_irf_to_daily_simulation(
    parameters: VISITSoilParameters,
    baseline_litter_inputs: ArrayLike,
    input_perturbations: ArrayLike,
    litter_scalar: float,
    humus_scalars: ArrayLike,
) -> VISITSoilTrajectoryComparison:
    """Compare IRF convolution with direct source-order daily soil updates.

    The baseline begins at its fixed point. This fixed-coefficient soil slice is
    linear before clipping, so disagreement exposes clipping or implementation
    errors. Later native VISIT comparisons can reuse the same result contract.
    """
    baseline_inputs = _validated_vector(
        baseline_litter_inputs, 6, "baseline_litter_inputs"
    )
    perturbations = np.asarray(input_perturbations, dtype=float)
    if perturbations.ndim != 2 or perturbations.shape[1] != 6:
        raise ValueError("input_perturbations must have shape (n_steps, 6)")
    if not np.all(np.isfinite(perturbations)):
        raise ValueError("input_perturbations must be finite")
    perturbed_inputs = baseline_inputs[np.newaxis, :] + perturbations
    if np.any(perturbed_inputs < 0):
        raise ValueError("perturbed litter inputs must be non-negative")

    system = visit_soil_discrete_system(parameters, litter_scalar, humus_scalars)
    fixed_point = visit_soil_fixed_point(
        parameters, baseline_inputs, litter_scalar, humus_scalars
    )
    n_steps = perturbations.shape[0]
    baseline_states = np.empty((n_steps + 1, 9), dtype=float)
    perturbed_states = np.empty((n_steps + 1, 9), dtype=float)
    baseline_respiration = np.empty(n_steps, dtype=float)
    perturbed_respiration = np.empty(n_steps, dtype=float)
    baseline_states[0] = fixed_point.state
    perturbed_states[0] = fixed_point.state

    for step in range(n_steps):
        baseline_result = visit_soil_daily_step(
            baseline_states[step],
            baseline_inputs,
            parameters,
            litter_scalar,
            humus_scalars,
        )
        perturbed_result = visit_soil_daily_step(
            perturbed_states[step],
            perturbed_inputs[step],
            parameters,
            litter_scalar,
            humus_scalars,
        )
        baseline_states[step + 1] = baseline_result.state
        perturbed_states[step + 1] = perturbed_result.state
        baseline_respiration[step] = (
            baseline_result.fluxes.heterotrophic_respiration
        )
        perturbed_respiration[step] = (
            perturbed_result.fluxes.heterotrophic_respiration
        )

    irf_state_difference, irf_respiration_difference = system.forced_response(
        perturbations
    )
    irf_predicted_states = baseline_states + irf_state_difference
    irf_predicted_respiration = (
        baseline_respiration + irf_respiration_difference[:, 0]
    )
    baseline_phase = visit_soil_litter_humus_projection(baseline_states)
    perturbed_phase = visit_soil_litter_humus_projection(perturbed_states)
    irf_predicted_phase = visit_soil_litter_humus_projection(irf_predicted_states)
    return VISITSoilTrajectoryComparison(
        days=np.arange(n_steps + 1, dtype=float),
        baseline_states=baseline_states,
        perturbed_states=perturbed_states,
        irf_predicted_states=irf_predicted_states,
        baseline_heterotrophic_respiration=baseline_respiration,
        perturbed_heterotrophic_respiration=perturbed_respiration,
        irf_predicted_heterotrophic_respiration=irf_predicted_respiration,
        phase_axis_names=("total litter carbon", "total humus carbon"),
        baseline_phase=baseline_phase,
        perturbed_phase=perturbed_phase,
        irf_predicted_phase=irf_predicted_phase,
    )


def visit_soil_daily_step(
    state: ArrayLike,
    litter_inputs: ArrayLike,
    parameters: VISITSoilParameters,
    litter_scalar: float,
    humus_scalars: ArrayLike,
    *,
    clip_negative: bool = True,
) -> VISITSoilStepResult:
    """Reproduce the carbon-pool algebra and pre-input update order of ``f_cycle_soil``."""
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