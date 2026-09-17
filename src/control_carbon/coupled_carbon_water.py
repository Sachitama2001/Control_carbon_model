"""Minimal carbon--water model for matrix and non-equilibrium analysis.

The model shares four conceptual components (leaf, stem, root, soil) between
carbon and water, yielding eight dynamic states.  It is a deliberately reduced
continuous-time synthesis, not a transcription of native VISITc.  Carbon uses
a donor-column compartment matrix.  Water uses an incidence matrix and
potential-gradient fluxes, preserving the different physical meanings of the
two systems.

All rates use days as the time unit. Carbon stocks are Mg C ha-1 and water
stores are mm. Parameter defaults are illustrative and must not be presented
as calibrated VISITc values.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .hydraulic_relations import (
    HYDRAULIC_RELATION_PROVENANCE,
    PlantPressureVolumeParameters,
    VISITCSoilTexture,
    plant_pressure_volume_potential,
    visitc_soil_total_potential,
)
from .provenance import provenance_manifest
from .visitc_source_map import (
    VISITC_HYDROLOGY_SOURCE,
    VISITC_LOCATION_SOURCE,
    VISITC_STRUCTURE_SOURCE,
)


CARBON_STATE_NAMES = ("leaf_carbon", "stem_carbon", "root_carbon", "soil_carbon")
WATER_STATE_NAMES = ("leaf_water", "stem_water", "root_water", "soil_water")
COUPLED_STATE_NAMES = CARBON_STATE_NAMES + WATER_STATE_NAMES
COUPLED_STATE_UNITS = ("Mg C ha-1",) * 4 + ("mm",) * 4
WATER_FLUX_NAMES = ("soil_to_root", "root_to_stem", "stem_to_leaf")

# A column is one internal flow; -1 marks its donor and +1 its receiver.
WATER_INCIDENCE_MATRIX = np.array(
    [
        [0.0, 0.0, 1.0],
        [0.0, 1.0, -1.0],
        [1.0, -1.0, 0.0],
        [-1.0, 0.0, 0.0],
    ]
)

COUPLED_MODEL_LITERATURE = {
    "matrix_review": "https://doi.org/10.1029/2022MS003008",
    "disequilibrium_decomposition": "https://doi.org/10.1029/2021JG006764",
    "FETCH2": "https://doi.org/10.1002/2016JG003467",
    "SurEau-Ecos": "https://doi.org/10.5194/gmd-15-5593-2022",
    "rate_tipping_framework": "https://doi.org/10.1098/rsta.2011.0306",
}


def coupled_model_provenance_manifest() -> dict[str, object]:
    """Provenance for the VISITc concepts used in this reduced synthesis."""
    manifest = provenance_manifest(
        model="four-component carbon-water synthesis",
        approximation_level="source-grounded reduced continuous model",
        sources=(
            VISITC_STRUCTURE_SOURCE,
            VISITC_HYDROLOGY_SOURCE,
            VISITC_LOCATION_SOURCE,
        ),
    )
    manifest["literature"] = dict(COUPLED_MODEL_LITERATURE)
    manifest["hydraulic_relations"] = HYDRAULIC_RELATION_PROVENANCE
    manifest["assumptions"] = [
        "single vegetation type; no snow, fire, harvest, nitrogen, or phosphorus states",
        "four aggregated carbon pools and four finite water stores",
        "VISITc soil retention plus aggregated SurEau-Ecos symplasmic plant pressure-volume curves",
        "remaining carbon and external water flux defaults are illustrative, not calibrated VISITc parameters",
    ]
    return manifest


@dataclass(frozen=True)
class CoupledCarbonWaterForcing:
    """External/effective forcing held constant at an instant."""

    precipitation: float = 2.5
    photosynthesis_scale: float = 1.0
    decomposition_scale: float = 1.0
    potential_transpiration: float = 1.8
    potential_soil_evaporation: float = 0.35

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("forcing values must be finite and nonnegative")


@dataclass(frozen=True)
class CoupledCarbonWaterParameters:
    """Illustrative parameters for the eight-state reduced model."""

    carbon_allocation: tuple[float, float, float] = (0.35, 0.25, 0.40)
    plant_turnover_rates: tuple[float, float, float] = (
        1.0 / 365.0,
        1.0 / (30.0 * 365.0),
        1.0 / (3.0 * 365.0),
    )
    soil_respiration_rate: float = 1.0 / (18.0 * 365.0)
    maximum_gpp: float = 0.022
    photosynthesis_water_half_saturation: float = 0.25
    decomposition_water_half_saturation: float = 0.35
    water_capacities: tuple[float, float, float, float] = (1.5, 12.0, 4.0, 180.0)
    hydraulic_conductances: tuple[float, float, float] = (1.4, 1.0, 0.7)
    hydraulic_carbon_half_saturation: tuple[float, float, float] = (1.0, 8.0, 0.5)
    soil_texture: int = int(VISITCSoilTexture.MEDIUM)
    soil_gravitational_potential: float = -1.0
    plant_osmotic_potentials: tuple[float, float, float] = (-1.5, -1.3, -1.1)
    plant_bulk_moduli: tuple[float, float, float] = (12.0, 10.0, 8.0)
    plant_minimum_relative_water: float = 1.0e-3
    transpiration_water_half_saturation: float = 0.15
    transpiration_leaf_carbon_half_saturation: float = 0.5
    evaporation_water_half_saturation: float = 0.20
    drainage_rate: float = 0.012
    field_capacity_fraction: float = 0.72
    runoff_fraction: float = 0.04

    def __post_init__(self) -> None:
        allocation = _parameter_vector(self.carbon_allocation, 3, "carbon_allocation")
        turnover = _parameter_vector(
            self.plant_turnover_rates, 3, "plant_turnover_rates", positive=True
        )
        capacities = _parameter_vector(
            self.water_capacities, 4, "water_capacities", positive=True
        )
        conductances = _parameter_vector(
            self.hydraulic_conductances, 3, "hydraulic_conductances", positive=True
        )
        half_saturation = _parameter_vector(
            self.hydraulic_carbon_half_saturation,
            3,
            "hydraulic_carbon_half_saturation",
            positive=True,
        )
        osmotic = _parameter_vector(
            self.plant_osmotic_potentials, 3, "plant_osmotic_potentials"
        )
        bulk_moduli = _parameter_vector(
            self.plant_bulk_moduli, 3, "plant_bulk_moduli", positive=True
        )
        scalars = np.asarray(
            [
                self.soil_respiration_rate,
                self.maximum_gpp,
                self.photosynthesis_water_half_saturation,
                self.decomposition_water_half_saturation,
                self.transpiration_water_half_saturation,
                self.transpiration_leaf_carbon_half_saturation,
                self.evaporation_water_half_saturation,
                self.drainage_rate,
                self.plant_minimum_relative_water,
                self.soil_gravitational_potential,
            ]
        )
        if not np.all(np.isfinite(scalars)):
            raise ValueError("scalar parameters must be finite")
        if np.any(scalars[:-1] <= 0.0):
            raise ValueError("rate and half-saturation parameters must be positive")
        if np.any(osmotic >= 0.0):
            raise ValueError("plant_osmotic_potentials must be negative")
        if np.any(bulk_moduli <= np.abs(osmotic)):
            raise ValueError("plant bulk moduli must exceed osmotic magnitudes")
        try:
            VISITCSoilTexture(self.soil_texture)
        except ValueError as error:
            raise ValueError("soil_texture must be VISITc selector 0, 1, or 2") from error
        if not 0.0 < self.plant_minimum_relative_water < 1.0:
            raise ValueError("plant_minimum_relative_water must lie in (0, 1)")
        if not np.isclose(allocation.sum(), 1.0):
            raise ValueError("carbon_allocation must sum to one")
        if np.any(allocation < 0.0):
            raise ValueError("carbon_allocation must be nonnegative")
        if not np.isfinite(self.field_capacity_fraction) or not 0.0 < self.field_capacity_fraction <= 1.0:
            raise ValueError("field_capacity_fraction must lie in (0, 1]")
        if not np.isfinite(self.runoff_fraction) or not 0.0 <= self.runoff_fraction <= 1.0:
            raise ValueError("runoff_fraction must lie in [0, 1]")
        object.__setattr__(self, "carbon_allocation", tuple(allocation))
        object.__setattr__(self, "plant_turnover_rates", tuple(turnover))
        object.__setattr__(self, "water_capacities", tuple(capacities))
        object.__setattr__(self, "hydraulic_conductances", tuple(conductances))
        object.__setattr__(
            self, "hydraulic_carbon_half_saturation", tuple(half_saturation)
        )
        object.__setattr__(self, "plant_osmotic_potentials", tuple(osmotic))
        object.__setattr__(self, "plant_bulk_moduli", tuple(bulk_moduli))


@dataclass(frozen=True)
class CoupledCarbonWaterFluxes:
    """Fluxes used in one right-hand-side evaluation."""

    gross_primary_production: float
    plant_turnover: np.ndarray
    heterotrophic_respiration: float
    carbon_input: np.ndarray
    carbon_matrix: np.ndarray
    water_internal: np.ndarray
    water_potential: np.ndarray
    effective_hydraulic_conductance: np.ndarray
    transpiration: float
    soil_evaporation: float
    drainage: float
    runoff: float
    water_input: np.ndarray
    water_output: np.ndarray


@dataclass(frozen=True)
class CoupledCarbonWaterEvaluation:
    derivative: np.ndarray
    fluxes: CoupledCarbonWaterFluxes
    carbon_budget_residual: float
    water_budget_residual: float


@dataclass(frozen=True)
class CarbonCapacity:
    """Frozen-coefficient carbon capacity at an instantaneous water state."""

    state: np.ndarray
    total: float
    input_rate: float
    ecosystem_transit_time: float
    carbon_matrix: np.ndarray
    carbon_input: np.ndarray


@dataclass(frozen=True)
class CapacityChangeDecomposition:
    productivity_effect: float
    transit_time_effect: float
    interaction_effect: float
    reconstructed_change: float
    direct_change: float


@dataclass(frozen=True)
class CoupledJacobianBlocks:
    carbon_carbon: np.ndarray
    carbon_water: np.ndarray
    water_carbon: np.ndarray
    water_water: np.ndarray
    full: np.ndarray


@dataclass(frozen=True)
class CoupledEquilibrium:
    state: np.ndarray
    residual_norm: float
    converged: bool
    jacobian: np.ndarray
    eigenvalues: np.ndarray
    stability: str
    condition_number: float
    function_evaluations: int
    message: str


@dataclass(frozen=True)
class CoupledTrajectory:
    times: np.ndarray
    states: np.ndarray
    forcings: tuple[CoupledCarbonWaterForcing, ...]
    function_evaluations: int
    rtol: float
    atol: float
    max_step: float


def _parameter_vector(
    values: ArrayLike, size: int, name: str, *, positive: bool = False
) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (size,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape ({size},)")
    if positive and np.any(array <= 0.0):
        raise ValueError(f"{name} must be positive")
    return array


def _state(values: ArrayLike, *, allow_negative: bool = False) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (8,) or not np.all(np.isfinite(array)):
        raise ValueError("state must be finite with shape (8,)")
    if not allow_negative and np.any(array < 0.0):
        raise ValueError("state must be nonnegative")
    return array


def _saturation(value: float, half_saturation: float) -> float:
    value = max(float(value), 0.0)
    return value / (value + half_saturation)


def _carbon_and_water_fluxes(
    state: np.ndarray,
    forcing: CoupledCarbonWaterForcing,
    parameters: CoupledCarbonWaterParameters,
) -> CoupledCarbonWaterFluxes:
    carbon = np.maximum(state[:4], 0.0)
    water = np.maximum(state[4:], 0.0)
    capacity = np.asarray(parameters.water_capacities)
    relative_water = water / capacity

    photosynthesis_water = np.sqrt(
        _saturation(relative_water[0], parameters.photosynthesis_water_half_saturation)
        * _saturation(relative_water[3], parameters.photosynthesis_water_half_saturation)
    )
    gpp = parameters.maximum_gpp * forcing.photosynthesis_scale * photosynthesis_water

    turnover_rates = np.asarray(parameters.plant_turnover_rates)
    plant_turnover = turnover_rates * carbon[:3]
    decomposition_water = _saturation(
        relative_water[3], parameters.decomposition_water_half_saturation
    )
    soil_loss_rate = (
        parameters.soil_respiration_rate
        * forcing.decomposition_scale
        * decomposition_water
    )
    respiration = soil_loss_rate * carbon[3]

    carbon_matrix = np.zeros((4, 4), dtype=float)
    carbon_matrix[:3, :3] = -np.diag(turnover_rates)
    carbon_matrix[3, :3] = turnover_rates
    carbon_matrix[3, 3] = -soil_loss_rate
    carbon_input = np.zeros(4, dtype=float)
    carbon_input[:3] = np.asarray(parameters.carbon_allocation) * gpp

    half_carbon = np.asarray(parameters.hydraulic_carbon_half_saturation)
    carbon_controls = np.array(
        [
            _saturation(carbon[2], half_carbon[0]),
            np.sqrt(
                _saturation(carbon[2], half_carbon[0])
                * _saturation(carbon[1], half_carbon[1])
            ),
            np.sqrt(
                _saturation(carbon[1], half_carbon[1])
                * _saturation(carbon[0], half_carbon[2])
            ),
        ]
    )
    conductance = np.asarray(parameters.hydraulic_conductances) * carbon_controls
    plant_potential = np.asarray(
        [
            plant_pressure_volume_potential(
                water[index],
                PlantPressureVolumeParameters(
                    saturated_water=capacity[index],
                    osmotic_potential=parameters.plant_osmotic_potentials[index],
                    bulk_modulus=parameters.plant_bulk_moduli[index],
                    minimum_relative_water=parameters.plant_minimum_relative_water,
                ),
            )
            for index in range(3)
        ]
    )
    soil_potential = visitc_soil_total_potential(
        water[3],
        capacity[3],
        parameters.soil_texture,
        gravitational_potential=parameters.soil_gravitational_potential,
    )
    water_potential = np.concatenate((plant_potential, (soil_potential,)))
    water_internal = conductance * np.array(
        [
            water_potential[3] - water_potential[2],
            water_potential[2] - water_potential[1],
            water_potential[1] - water_potential[0],
        ]
    )

    transpiration = (
        forcing.potential_transpiration
        * _saturation(
            relative_water[0], parameters.transpiration_water_half_saturation
        )
        * _saturation(
            carbon[0], parameters.transpiration_leaf_carbon_half_saturation
        )
    )
    soil_evaporation = forcing.potential_soil_evaporation * _saturation(
        relative_water[3], parameters.evaporation_water_half_saturation
    )
    field_capacity = parameters.field_capacity_fraction * capacity[3]
    drainage = parameters.drainage_rate * max(water[3] - field_capacity, 0.0)
    runoff = (
        parameters.runoff_fraction
        * forcing.precipitation
        * min(max(relative_water[3], 0.0), 1.0) ** 2
    )
    water_input = np.array([0.0, 0.0, 0.0, forcing.precipitation])
    water_output = np.array(
        [transpiration, 0.0, 0.0, soil_evaporation + drainage + runoff]
    )
    return CoupledCarbonWaterFluxes(
        gross_primary_production=gpp,
        plant_turnover=plant_turnover,
        heterotrophic_respiration=respiration,
        carbon_input=carbon_input,
        carbon_matrix=carbon_matrix,
        water_internal=water_internal,
        water_potential=water_potential,
        effective_hydraulic_conductance=conductance,
        transpiration=transpiration,
        soil_evaporation=soil_evaporation,
        drainage=drainage,
        runoff=runoff,
        water_input=water_input,
        water_output=water_output,
    )


def evaluate_coupled_carbon_water(
    state: ArrayLike,
    forcing: CoupledCarbonWaterForcing,
    parameters: CoupledCarbonWaterParameters | None = None,
    *,
    allow_negative_trial_state: bool = False,
) -> CoupledCarbonWaterEvaluation:
    """Evaluate fluxes, derivatives, and exact numerical budget residuals."""
    parameters = parameters or CoupledCarbonWaterParameters()
    state_array = _state(state, allow_negative=allow_negative_trial_state)
    fluxes = _carbon_and_water_fluxes(state_array, forcing, parameters)
    carbon_derivative = fluxes.carbon_input + fluxes.carbon_matrix @ state_array[:4]
    water_derivative = (
        WATER_INCIDENCE_MATRIX @ fluxes.water_internal
        + fluxes.water_input
        - fluxes.water_output
    )
    derivative = np.concatenate((carbon_derivative, water_derivative))
    carbon_residual = float(
        carbon_derivative.sum()
        - (
            fluxes.gross_primary_production
            - fluxes.heterotrophic_respiration
        )
    )
    water_residual = float(
        water_derivative.sum()
        - (
            forcing.precipitation
            - fluxes.transpiration
            - fluxes.soil_evaporation
            - fluxes.drainage
            - fluxes.runoff
        )
    )
    return CoupledCarbonWaterEvaluation(
        derivative=derivative,
        fluxes=fluxes,
        carbon_budget_residual=carbon_residual,
        water_budget_residual=water_residual,
    )


def coupled_carbon_water_rhs(
    state: ArrayLike,
    forcing: CoupledCarbonWaterForcing,
    parameters: CoupledCarbonWaterParameters | None = None,
) -> np.ndarray:
    """Return the eight-state derivative for a single forcing value."""
    return evaluate_coupled_carbon_water(state, forcing, parameters).derivative


def instantaneous_carbon_capacity(
    state: ArrayLike,
    forcing: CoupledCarbonWaterForcing,
    parameters: CoupledCarbonWaterParameters | None = None,
) -> CarbonCapacity:
    """Freeze the current water-dependent carbon coefficients and solve capacity.

    This is a Luo-style diagnostic ``-M_C^-1 u_C``. It is not generally an
    equilibrium of the coupled system because water and hydraulic fluxes may
    subsequently change.
    """
    parameters = parameters or CoupledCarbonWaterParameters()
    state_array = _state(state)
    fluxes = _carbon_and_water_fluxes(state_array, forcing, parameters)
    try:
        capacity = np.linalg.solve(-fluxes.carbon_matrix, fluxes.carbon_input)
    except np.linalg.LinAlgError as error:
        raise ValueError("instantaneous carbon matrix has no finite capacity") from error
    input_rate = float(fluxes.carbon_input.sum())
    total = float(capacity.sum())
    transit_time = total / input_rate if input_rate > 0.0 else np.inf
    return CarbonCapacity(
        state=capacity,
        total=total,
        input_rate=input_rate,
        ecosystem_transit_time=transit_time,
        carbon_matrix=fluxes.carbon_matrix.copy(),
        carbon_input=fluxes.carbon_input.copy(),
    )


def decompose_capacity_change(
    baseline: CarbonCapacity, comparison: CarbonCapacity
) -> CapacityChangeDecomposition:
    """Decompose ``delta(mu*tau)`` into productivity, time, and interaction."""
    delta_input = comparison.input_rate - baseline.input_rate
    delta_time = (
        comparison.ecosystem_transit_time - baseline.ecosystem_transit_time
    )
    productivity = baseline.ecosystem_transit_time * delta_input
    transit = baseline.input_rate * delta_time
    interaction = delta_input * delta_time
    reconstructed = productivity + transit + interaction
    return CapacityChangeDecomposition(
        productivity_effect=float(productivity),
        transit_time_effect=float(transit),
        interaction_effect=float(interaction),
        reconstructed_change=float(reconstructed),
        direct_change=float(comparison.total - baseline.total),
    )


def coupled_jacobian_blocks(
    state: ArrayLike,
    forcing: CoupledCarbonWaterForcing,
    parameters: CoupledCarbonWaterParameters | None = None,
    *,
    relative_step: float = 1e-6,
) -> CoupledJacobianBlocks:
    """Central finite-difference Jacobian partitioned into C/W blocks."""
    parameters = parameters or CoupledCarbonWaterParameters()
    state_array = _state(state)
    if not np.isfinite(relative_step) or relative_step <= 0.0:
        raise ValueError("relative_step must be finite and positive")
    jacobian = np.empty((8, 8), dtype=float)
    for column in range(8):
        step = relative_step * max(1.0, abs(state_array[column]))
        lower = state_array.copy()
        upper = state_array.copy()
        lower[column] -= step
        upper[column] += step
        if lower[column] < 0.0:
            baseline = evaluate_coupled_carbon_water(
                state_array, forcing, parameters
            ).derivative
            high = evaluate_coupled_carbon_water(
                upper, forcing, parameters
            ).derivative
            jacobian[:, column] = (high - baseline) / step
        else:
            low = evaluate_coupled_carbon_water(
                lower, forcing, parameters
            ).derivative
            high = evaluate_coupled_carbon_water(
                upper, forcing, parameters
            ).derivative
            jacobian[:, column] = (high - low) / (2.0 * step)
    return CoupledJacobianBlocks(
        carbon_carbon=jacobian[:4, :4].copy(),
        carbon_water=jacobian[:4, 4:].copy(),
        water_carbon=jacobian[4:, :4].copy(),
        water_water=jacobian[4:, 4:].copy(),
        full=jacobian,
    )


def schur_complement_carbon_jacobian(
    blocks: CoupledJacobianBlocks,
) -> np.ndarray:
    """Eliminate fast water anomalies: J_CC - J_CW J_WW^-1 J_WC."""
    try:
        water_response = np.linalg.solve(
            blocks.water_water, blocks.water_carbon
        )
    except np.linalg.LinAlgError as error:
        raise ValueError("water Jacobian is singular; Schur reduction is undefined") from error
    return blocks.carbon_carbon - blocks.carbon_water @ water_response


def coupled_equilibrium(
    forcing: CoupledCarbonWaterForcing,
    initial_guess: ArrayLike,
    parameters: CoupledCarbonWaterParameters | None = None,
    *,
    residual_tolerance: float = 1e-9,
    stability_tolerance: float = 1e-8,
) -> CoupledEquilibrium:
    """Solve the frozen-forcing coupled equilibrium with nonnegative bounds."""
    parameters = parameters or CoupledCarbonWaterParameters()
    guess = _state(initial_guess)
    if residual_tolerance <= 0.0 or stability_tolerance <= 0.0:
        raise ValueError("tolerances must be positive")
    solution = least_squares(
        lambda trial: evaluate_coupled_carbon_water(
            trial, forcing, parameters
        ).derivative,
        guess,
        bounds=(np.zeros(8), np.full(8, np.inf)),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=5000,
    )
    residual = float(np.linalg.norm(solution.fun, ord=np.inf))
    blocks = coupled_jacobian_blocks(solution.x, forcing, parameters)
    eigenvalues = np.linalg.eigvals(blocks.full)
    converged = bool(solution.success and residual <= residual_tolerance)
    stability = "unresolved"
    if converged:
        if np.any(eigenvalues.real > stability_tolerance):
            stability = "unstable"
        elif np.all(eigenvalues.real < -stability_tolerance):
            stability = "stable"
        else:
            stability = "nonhyperbolic"
    return CoupledEquilibrium(
        state=solution.x.copy(),
        residual_norm=residual,
        converged=converged,
        jacobian=blocks.full,
        eigenvalues=eigenvalues,
        stability=stability,
        condition_number=float(np.linalg.cond(blocks.full)),
        function_evaluations=solution.nfev,
        message=str(solution.message),
    )


def simulate_coupled_carbon_water(
    initial_state: ArrayLike,
    times: ArrayLike,
    forcing: CoupledCarbonWaterForcing
    | Callable[[float], CoupledCarbonWaterForcing],
    parameters: CoupledCarbonWaterParameters | None = None,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_step: float = 0.25,
) -> CoupledTrajectory:
    """Integrate the reduced model under constant or time-varying forcing."""
    parameters = parameters or CoupledCarbonWaterParameters()
    initial = _state(initial_state)
    times_array = np.asarray(times, dtype=float)
    if (
        times_array.ndim != 1
        or times_array.size < 2
        or not np.all(np.isfinite(times_array))
        or np.any(np.diff(times_array) <= 0.0)
    ):
        raise ValueError("times must be finite and strictly increasing")
    if rtol <= 0.0 or atol <= 0.0 or max_step <= 0.0:
        raise ValueError("solver tolerances and max_step must be positive")
    forcing_function = forcing if callable(forcing) else lambda _time: forcing

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        forcing_value = forcing_function(time)
        if not isinstance(forcing_value, CoupledCarbonWaterForcing):
            raise TypeError("forcing callable must return CoupledCarbonWaterForcing")
        # Tiny negative solver stages are handled by the smooth boundary
        # extension used internally; accepted output states are checked below.
        return evaluate_coupled_carbon_water(
            state,
            forcing_value,
            parameters,
            allow_negative_trial_state=True,
        ).derivative

    solution = solve_ivp(
        rhs,
        (times_array[0], times_array[-1]),
        initial,
        t_eval=times_array,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(f"coupled integration failed: {solution.message}")
    if not np.all(np.isfinite(solution.y)):
        raise RuntimeError("coupled integration returned nonfinite states")
    if np.min(solution.y) < -10.0 * atol:
        raise RuntimeError("coupled integration left the nonnegative state domain")
    states = np.maximum(solution.y.T, 0.0)
    forcing_values = tuple(forcing_function(float(time)) for time in times_array)
    if not all(
        isinstance(value, CoupledCarbonWaterForcing) for value in forcing_values
    ):
        raise TypeError("forcing callable must return CoupledCarbonWaterForcing")
    return CoupledTrajectory(
        times=times_array.copy(),
        states=states,
        forcings=forcing_values,
        function_evaluations=solution.nfev,
        rtol=float(rtol),
        atol=float(atol),
        max_step=float(max_step),
    )


__all__ = [
    "CARBON_STATE_NAMES",
    "COUPLED_MODEL_LITERATURE",
    "COUPLED_STATE_NAMES",
    "COUPLED_STATE_UNITS",
    "WATER_FLUX_NAMES",
    "WATER_INCIDENCE_MATRIX",
    "WATER_STATE_NAMES",
    "CapacityChangeDecomposition",
    "CarbonCapacity",
    "CoupledCarbonWaterEvaluation",
    "CoupledCarbonWaterFluxes",
    "CoupledCarbonWaterForcing",
    "CoupledCarbonWaterParameters",
    "CoupledEquilibrium",
    "CoupledJacobianBlocks",
    "CoupledTrajectory",
    "coupled_carbon_water_rhs",
    "coupled_equilibrium",
    "coupled_jacobian_blocks",
    "coupled_model_provenance_manifest",
    "decompose_capacity_change",
    "evaluate_coupled_carbon_water",
    "instantaneous_carbon_capacity",
    "schur_complement_carbon_jacobian",
    "simulate_coupled_carbon_water",
]
