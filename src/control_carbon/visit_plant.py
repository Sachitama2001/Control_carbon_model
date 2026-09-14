"""Source-grounded structural plant-carbon slices from VISIT."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .provenance import SourceRef
from .visit_source_map import VISIT_SOURCE_COMMIT, VISIT_SOURCE_REPOSITORY


VISIT_LITTERFALL_SOURCE = SourceRef(
    path="visit_local/litterfall.c",
    symbol="f_lf / f_lc / f_lr",
    role="foliage, stem, and root litterfall fluxes",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="f_lf / f_lc / f_lr",
    units="Mg C ha^-1 day^-1",
    assumptions=(
        "turnover rates are supplied after f_ecophysiology::f_mortality",
        "state is the plant_process state after emergence and crop planting",
    ),
    approximation_level="native-source-plant-turnover-slice",
)

VISIT_PLANT_TURNOVER_SOURCE = SourceRef(
    path="visit_local/plant_proc.c",
    symbol="plant_process litterfall block",
    role="seasonal, drought, and crop-stage litterfall selection and pool update",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="plant_process",
    units="Mg C ha^-1 and Mg C ha^-1 day^-1",
    assumptions=(
        "photosynthesis, respiration, allocation, NSC, and survival reallocation excluded",
        "FIX_PHENOLOGY does not alter litterfall algebra after season is supplied",
    ),
    approximation_level="native-source-plant-turnover-slice",
)

VISIT_PLANT_RESPIRATION_SOURCE = SourceRef(
    path="visit_local/respiration.c / visit_local/ecophysiology.c",
    symbol="f_rfm/f_rcm/f_rrm/f_rfg/f_rcg/f_rrg; f_q10_ar/f_spcfc_resp",
    role="maintenance and growth respiration with daily coefficient updates",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="f_q10_ar / f_spcfc_resp / respiration functions",
    units="Mg C ha^-1 day^-1",
    assumptions=("nitrogen concentration multiplier is 1.0 in the source",),
    approximation_level="native-source-plant-respiration-slice",
)

VISIT_PLANT_ALLOCATION_SOURCE = SourceRef(
    path="visit_local/allocation.c",
    symbol="f_allocation",
    role="piecewise allocation of GPP/EPP to foliage, stem, root, and grain",
    repository=VISIT_SOURCE_REPOSITORY,
    commit=VISIT_SOURCE_COMMIT,
    function="f_allocation",
    units="fractions and Mg C ha^-1 day^-1",
    approximation_level="native-source-plant-allocation-slice",
)


@dataclass(frozen=True)
class VISITPlantStructuralState:
    """Foliage, stem/branch, and root carbon immediately before litterfall."""

    foliage: float
    stem: float
    root: float

    def __post_init__(self) -> None:
        values = np.asarray((self.foliage, self.stem, self.root), dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("plant structural carbon must be finite and non-negative")

    def as_array(self) -> np.ndarray:
        return np.asarray((self.foliage, self.stem, self.root), dtype=float)


@dataclass(frozen=True)
class VISITPlantTurnoverParameters:
    """Daily rates already updated by VISIT's ``f_mortality`` path."""

    leaf_rate: float
    stem_rate: float
    root_rate: float
    deciduous_shedding_fraction: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
            raise ValueError("turnover rates and shedding fraction must lie in [0, 1]")


@dataclass(frozen=True)
class VISITPlantTurnoverContext:
    """Regime variables read by the litterfall block of ``plant_process``."""

    season: int
    crop_stage: int
    lai: float
    deep_soil_water_potential: float
    air_temperature: float
    photosynthesis_type: int

    def __post_init__(self) -> None:
        values = np.asarray(
            (self.lai, self.deep_soil_water_potential, self.air_temperature),
            dtype=float,
        )
        if not np.all(np.isfinite(values)) or self.lai < 0.0:
            raise ValueError("turnover context values must be finite and LAI non-negative")


@dataclass(frozen=True)
class VISITPlantTurnoverResult:
    """Structural state and litter flux after the native litterfall block."""

    state: VISITPlantStructuralState
    litterfall: np.ndarray
    regime: str


@dataclass(frozen=True)
class VISITPlantRespirationParameters:
    """TKY-compatible source parameters used by plant respiration routines."""

    growth_foliage: float
    growth_stem: float
    growth_root: float
    maintenance_foliage_15c: float
    maintenance_stem_sapwood: float
    maintenance_stem_heartwood: float
    maintenance_root_fine: float
    maintenance_root_coarse: float
    q10_foliage_base: float
    q10_stem_base: float
    q10_root_base: float
    size_stem: float
    size_root: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("respiration parameters must be finite and non-negative")
        if np.any(values[8:11] <= 0.0) or np.any(values[11:] <= 0.0):
            raise ValueError("Q10 and size parameters must be positive")


@dataclass(frozen=True)
class VISITPlantRespirationResult:
    """Daily maintenance/growth respiration and source intermediate values."""

    maintenance: np.ndarray
    growth: np.ndarray
    q10: np.ndarray
    specific_maintenance: np.ndarray
    sapwood_fine_root: np.ndarray
    heartwood_coarse_root: np.ndarray

    @property
    def total(self) -> float:
        return float(self.maintenance.sum() + self.growth.sum())


@dataclass(frozen=True)
class VISITPlantAllocationParameters:
    foliage_fraction: float
    aboveground_woody_fraction: float
    specific_leaf_area: float

    def __post_init__(self) -> None:
        if not 0.0 < self.foliage_fraction <= 1.0:
            raise ValueError("foliage_fraction must lie in (0, 1]")
        if not 0.0 <= self.aboveground_woody_fraction <= 1.0:
            raise ValueError("aboveground_woody_fraction must lie in [0, 1]")
        if not np.isfinite(self.specific_leaf_area) or self.specific_leaf_area <= 0.0:
            raise ValueError("specific_leaf_area must be finite and positive")


@dataclass(frozen=True)
class VISITPlantAllocationResult:
    fractions: np.ndarray
    fluxes: np.ndarray
    total_translocation: float
    regime: str


def visit_plant_allocation(
    parameters: VISITPlantAllocationParameters,
    *,
    lai: float,
    optimum_lai: float,
    season: int,
    crop_stage: int,
    epp: float,
    gpp: float,
    maintenance_respiration: np.ndarray,
) -> VISITPlantAllocationResult:
    """Reproduce ``allocation.c::f_allocation`` without updating plant pools."""
    maintenance = np.asarray(maintenance_respiration, dtype=float)
    values = np.asarray((lai, optimum_lai, epp, gpp), dtype=float)
    if maintenance.shape != (3,) or not np.all(np.isfinite(maintenance)):
        raise ValueError("maintenance_respiration must be finite with shape (3,)")
    if not np.all(np.isfinite(values)) or lai < 0.0 or optimum_lai < 0.0:
        raise ValueError("allocation inputs must be finite and LAI non-negative")

    if epp > 0.0:
        if lai > optimum_lai:
            foliage_fraction = 0.0
            regime = "positive-epp-above-optimum-lai"
        else:
            foliage_demand = (
                (optimum_lai - lai)
                * 100.0
                * 2.0
                / 2.2
                / parameters.specific_leaf_area
            )
            baseline_foliage_flux = epp * parameters.foliage_fraction
            if foliage_demand <= baseline_foliage_flux:
                foliage_fraction = parameters.foliage_fraction
            else:
                foliage_fraction = min(foliage_demand / epp, 0.05)
            regime = "positive-epp-below-optimum-lai"
        season_factor = 0.7 if season == 0 else 1.0
        foliage_fraction *= season_factor
        stem_fraction = (
            (1.0 - foliage_fraction) * parameters.aboveground_woody_fraction
        )
        root_fraction = (
            (1.0 - foliage_fraction)
            * (1.0 - parameters.aboveground_woody_fraction)
        )
        grain_fraction = 0.7 if crop_stage == 4 else 0.0
        structural_scale = 1.0 - grain_fraction
        fractions = np.asarray(
            (foliage_fraction, stem_fraction, root_fraction, grain_fraction)
        )
        fluxes = np.asarray(
            (
                structural_scale * foliage_fraction * epp,
                structural_scale * stem_fraction * epp,
                structural_scale * root_fraction * epp,
                grain_fraction * epp,
            )
        )
        total_translocation = float(fluxes.sum())
    else:
        foliage_fraction = parameters.foliage_fraction
        stem_fraction = (
            (1.0 - foliage_fraction) * parameters.aboveground_woody_fraction
        )
        root_fraction = (
            (1.0 - foliage_fraction)
            * (1.0 - parameters.aboveground_woody_fraction)
        )
        fractions = np.asarray(
            (foliage_fraction, stem_fraction, root_fraction, 0.0)
        )
        fluxes = np.asarray(
            (
                foliage_fraction * gpp - maintenance[0],
                stem_fraction * gpp - maintenance[1],
                root_fraction * gpp - maintenance[2],
                0.0,
            )
        )
        total_translocation = float(
            (foliage_fraction + stem_fraction + root_fraction) * epp
        )
        regime = "nonpositive-epp"
    return VISITPlantAllocationResult(
        fractions=fractions,
        fluxes=fluxes,
        total_translocation=total_translocation,
        regime=regime,
    )


def visit_plant_turnover_step(
    state: VISITPlantStructuralState,
    parameters: VISITPlantTurnoverParameters,
    context: VISITPlantTurnoverContext,
) -> VISITPlantTurnoverResult:
    """Reproduce the litterfall and immediate pool subtraction in plant_process."""
    if context.season == 3:
        if context.lai > 0.05:
            leaf_litter = parameters.deciduous_shedding_fraction * state.foliage
            regime = "seasonal-shedding"
        else:
            leaf_litter = state.foliage
            regime = "seasonal-final-shedding"
    else:
        leaf_rate = 0.7 if context.crop_stage == 5 else parameters.leaf_rate
        leaf_litter = leaf_rate * state.foliage
        regime = "crop-harvest" if context.crop_stage == 5 else "baseline-turnover"
        if (
            context.deep_soil_water_potential < -148.075
            and context.lai > 0.05
            and context.air_temperature > -5.0
            and context.photosynthesis_type in (3, 4)
        ):
            leaf_litter += 0.012 * state.foliage
            regime = "drought-enhanced-turnover"

    stem_rate = 0.7 if context.crop_stage == 5 else parameters.stem_rate
    root_rate = 0.7 if context.crop_stage == 5 else parameters.root_rate
    litterfall = np.asarray(
        (leaf_litter, stem_rate * state.stem, root_rate * state.root), dtype=float
    )
    updated = state.as_array() - litterfall
    if np.any(updated < 0.0):
        raise ValueError("source litterfall would make a structural pool negative")
    return VISITPlantTurnoverResult(
        state=VISITPlantStructuralState(*updated),
        litterfall=litterfall,
        regime=regime,
    )


def visit_plant_respiration(
    state: VISITPlantStructuralState,
    parameters: VISITPlantRespirationParameters,
    *,
    surface_temperature: float,
    upper_soil_temperature: float,
    allocation_fluxes: np.ndarray,
    positive_epp: bool = True,
) -> VISITPlantRespirationResult:
    """Reproduce daily coefficient updates and plant respiration equations."""
    temperatures = np.asarray((surface_temperature, upper_soil_temperature))
    allocations = np.asarray(allocation_fluxes, dtype=float)
    if not np.all(np.isfinite(temperatures)):
        raise ValueError("respiration temperatures must be finite")
    if allocations.shape != (3,) or not np.all(np.isfinite(allocations)):
        raise ValueError("allocation_fluxes must be finite with shape (3,)")
    if positive_epp and np.any(allocations < 0.0):
        raise ValueError("positive-EPP allocation fluxes must be non-negative")

    acclimation = np.exp(-0.009 * (surface_temperature - 15.0))
    q10 = acclimation * np.asarray(
        (
            parameters.q10_foliage_base,
            parameters.q10_stem_base,
            parameters.q10_root_base,
        )
    )
    stem_power = 1.0 - 0.33334 * state.stem / (parameters.size_stem + state.stem)
    root_power = 1.0 - 0.33334 * state.root / (parameters.size_root + state.root)
    sapwood = min(state.stem**stem_power, state.stem)
    fine_root = min(state.root**root_power, state.root)
    heartwood = state.stem - sapwood
    coarse_root = state.root - fine_root
    stem_specific = (
        parameters.maintenance_stem_sapwood * sapwood
        + parameters.maintenance_stem_heartwood * heartwood
    ) / (state.stem + 0.00001)
    root_specific = (
        parameters.maintenance_root_fine * fine_root
        + parameters.maintenance_root_coarse * coarse_root
    ) / (state.root + 0.00001)
    specific = np.asarray(
        (parameters.maintenance_foliage_15c, stem_specific, root_specific)
    )
    temperature_factors = np.asarray(
        (
            np.exp(np.log(q10[0]) / 10.0 * (surface_temperature - 15.0)),
            np.exp(np.log(q10[1]) / 10.0 * (surface_temperature - 15.0)),
            np.exp(np.log(q10[2]) / 10.0 * (upper_soil_temperature - 15.0)),
        )
    )
    maintenance = state.as_array() * specific / 1000.0 * temperature_factors
    growth_rates = np.asarray(
        (parameters.growth_foliage, parameters.growth_stem, parameters.growth_root)
    )
    growth = growth_rates * allocations if positive_epp else np.zeros(3)
    return VISITPlantRespirationResult(
        maintenance=maintenance,
        growth=growth,
        q10=q10,
        specific_maintenance=specific,
        sapwood_fine_root=np.asarray((sapwood, fine_root)),
        heartwood_coarse_root=np.asarray((heartwood, coarse_root)),
    )