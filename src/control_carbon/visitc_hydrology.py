"""Auditable partial transcription of the current VISITc water balance.

This module follows ``point/hydro_balance.c::f_hydrology`` at commit
``5202debd...`` in source order.  Penman--Monteith helper routines are not
transcribed here: their seven potential fluxes are explicit effective inputs.
That boundary lets us test the store-update algebra before reconstructing
radiation, aerodynamic resistance, and canopy conductance.

The routine intentionally preserves a source-level accounting issue: baseflow
is subtracted from ``sww`` and then included in ``ro2``, which is subtracted
again.  ``water_budget_residual`` exposes the resulting loss instead of
silently correcting it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .visitc_source_map import VISITC_HYDROLOGY_SOURCE


@dataclass(frozen=True)
class VISITCHydrologyState:
    """Native water stores at the beginning of a daily update, in millimetres."""

    snow_water: float
    upper_soil_water: float
    whole_soil_water: float

    def __post_init__(self) -> None:
        _finite_nonnegative(self.__dict__, "state")

    @property
    def total_storage(self) -> float:
        """Diagnostic sum used to audit the native update algebra."""
        return self.snow_water + self.upper_soil_water + self.whole_soil_water


@dataclass(frozen=True)
class VISITCHydrologyForcing:
    """Meteorology plus potential Penman--Monteith fluxes for one day."""

    precipitation: float
    air_temperature: float
    deep_soil_temperature: float
    potential_interception_tree: float = 0.0
    potential_interception_c3: float = 0.0
    potential_interception_c4: float = 0.0
    potential_soil_evaporation: float = 0.0
    potential_transpiration_tree: float = 0.0
    potential_transpiration_c3: float = 0.0
    potential_transpiration_c4: float = 0.0

    def __post_init__(self) -> None:
        values = self.__dict__
        if not all(np.isfinite(value) for value in values.values()):
            raise ValueError("forcing must be finite")
        nonnegative = {key: value for key, value in values.items() if "temperature" not in key}
        _finite_nonnegative(nonnegative, "water forcing")


@dataclass(frozen=True)
class VISITCHydrologyParameters:
    """Parameters read by the transcribed portion of ``f_hydrology``."""

    field_capacity_upper: float
    field_capacity_whole: float
    hydraulic_conductivity: float
    tree_leaf_area_index: float = 0.0
    c3_leaf_area_index: float = 0.0
    c4_leaf_area_index: float = 0.0
    c3_understory_fraction: float = 0.0
    c4_understory_fraction: float = 0.0
    site_id: str = ""

    def __post_init__(self) -> None:
        numeric = {
            key: value for key, value in self.__dict__.items() if key != "site_id"
        }
        _finite_nonnegative(numeric, "hydrology parameter")
        if self.field_capacity_upper <= 0.0 or self.field_capacity_whole <= 0.0:
            raise ValueError("field capacities must be positive")
        for name in ("c3_understory_fraction", "c4_understory_fraction"):
            if getattr(self, name) > 1.0:
                raise ValueError("understory fractions must lie in [0, 1]")
        if self.c3_understory_fraction + self.c4_understory_fraction > 1.0:
            raise ValueError("C3 and C4 understory fractions must sum to at most one")


@dataclass(frozen=True)
class VISITCHydrologyFluxes:
    """Daily source-named flux diagnostics, all in mm day-1."""

    snow_fraction: float
    snowfall: float
    rainfall: float
    thaw: float
    interception_tree: float
    interception_c3: float
    interception_c4: float
    liquid_input_to_soil: float
    upper_runoff: float
    soil_evaporation: float
    transpiration_tree: float
    transpiration_c3: float
    transpiration_c4: float
    baseflow: float
    lower_bucket_runoff: float
    lower_runoff: float
    redistribution: float
    actual_evapotranspiration: float
    potential_evapotranspiration: float
    upper_clip_correction: float
    whole_clip_correction: float

    @property
    def interception(self) -> float:
        return self.interception_tree + self.interception_c3 + self.interception_c4

    @property
    def transpiration(self) -> float:
        return self.transpiration_tree + self.transpiration_c3 + self.transpiration_c4


@dataclass(frozen=True)
class VISITCHydrologyStep:
    """State and audit diagnostics returned by one native-order update."""

    state: VISITCHydrologyState
    fluxes: VISITCHydrologyFluxes
    water_budget_residual: float
    expected_source_residual: float
    source: object = VISITC_HYDROLOGY_SOURCE

    @property
    def exposes_baseflow_double_subtraction(self) -> bool:
        """Whether nonzero baseflow contributes the documented extra loss."""
        return self.fluxes.baseflow > 0.0 and np.isclose(
            self.expected_source_residual
            - self.fluxes.upper_clip_correction
            - self.fluxes.whole_clip_correction,
            -self.fluxes.baseflow,
        )


def _finite_nonnegative(values: dict[str, float], label: str) -> None:
    if not all(np.isfinite(value) and value >= 0.0 for value in values.values()):
        raise ValueError(f"{label} values must be finite and nonnegative")


def _supply_limited_flux(supply: float, potential: float) -> float:
    """VISITc's quadratic actual-flux relation with its native 0.85 factor."""
    if supply < 0.0 or potential < 0.0:
        raise ValueError("supply and potential flux must be nonnegative")
    coefficient = 0.85
    total = supply + potential
    discriminant = total * total - 4.0 * coefficient * supply * potential
    # Round-off can only make a mathematically nonnegative discriminant tiny negative.
    discriminant = max(discriminant, 0.0)
    return (total - math.sqrt(discriminant)) / (2.0 * coefficient)


def _bucket_runoff(gain: float, dry_index: float) -> float:
    cubic_sum = max(gain**3 + dry_index**3, 0.0)
    return max(cubic_sum**0.33333 - dry_index, 0.0)


def visitc_hydrology_step(
    state: VISITCHydrologyState,
    forcing: VISITCHydrologyForcing,
    parameters: VISITCHydrologyParameters,
) -> VISITCHydrologyStep:
    """Evaluate one source-ordered partial ``f_hydrology`` update.

    The result is a transcription, not a corrected hydrological model.  In
    particular, ``whole_soil_water`` retains the source variable name ``sww``
    even though the update algebra often behaves like a second bucket.
    """
    snow_fraction = 1.0 / (1.0 + math.exp(0.75 * (forcing.air_temperature - 2.0)))
    snowfall = snow_fraction * forcing.precipitation
    rainfall = (1.0 - snow_fraction) * forcing.precipitation

    if state.snow_water > 0.1:
        thaw_fraction = (1.0 / 11.0) / (
            1.0 + math.exp(-0.5 * (forcing.air_temperature - 4.0))
        )
        thaw_multiplier = 1.0 + 10.0 / (0.05 * state.snow_water + 1.0)
        thaw = thaw_fraction * thaw_multiplier * state.snow_water
    else:
        thaw = state.snow_water
    snow_water = state.snow_water + snowfall - thaw

    tree_capture = min(rainfall, parameters.tree_leaf_area_index * 0.25)
    interception_tree = _supply_limited_flux(
        tree_capture, forcing.potential_interception_tree
    )

    c3_available = parameters.c3_understory_fraction * (
        rainfall - interception_tree
    )
    c3_capture = min(
        c3_available,
        parameters.c3_understory_fraction * parameters.c3_leaf_area_index * 0.25,
    )
    interception_c3 = _supply_limited_flux(
        c3_capture, forcing.potential_interception_c3
    )

    # The native C4 expression also subtracts only tree interception.
    c4_available = parameters.c4_understory_fraction * (
        rainfall - interception_tree
    )
    c4_capture = min(
        c4_available,
        parameters.c4_understory_fraction * parameters.c4_leaf_area_index * 0.25,
    )
    interception_c4 = _supply_limited_flux(
        c4_capture, forcing.potential_interception_c4
    )
    interception = interception_tree + interception_c3 + interception_c4
    liquid_input = rainfall - interception + thaw

    upper_runoff = _bucket_runoff(
        liquid_input, parameters.field_capacity_upper - state.upper_soil_water
    )
    upper = state.upper_soil_water + liquid_input - upper_runoff

    soil_evaporation = max(
        _supply_limited_flux(upper, forcing.potential_soil_evaporation), 0.0
    )
    upper -= soil_evaporation
    transpiration_c3 = _supply_limited_flux(
        upper, forcing.potential_transpiration_c3
    )
    transpiration_c4 = _supply_limited_flux(
        upper, forcing.potential_transpiration_c4
    )
    upper -= transpiration_c3 + transpiration_c4

    whole = state.whole_soil_water
    transpiration_tree = _supply_limited_flux(
        whole, forcing.potential_transpiration_tree
    )

    if forcing.deep_soil_temperature > 0.0:
        baseflow_rate = 0.003 if parameters.site_id == "QHB" else 0.001
        baseflow = baseflow_rate * whole
    else:
        baseflow = 0.0
    whole -= baseflow

    lower_bucket_runoff = _bucket_runoff(
        upper_runoff, parameters.field_capacity_whole - whole
    )
    lower_runoff = lower_bucket_runoff + baseflow

    maximum_redistribution = (
        parameters.hydraulic_conductivity * 1000.0 * 3600.0 * 24.0
    )
    capacity_ratio = parameters.field_capacity_upper / parameters.field_capacity_whole
    redistribution = (whole * capacity_ratio - upper) / (1.0 + capacity_ratio)
    if redistribution > 0.0:
        redistribution = min(0.5 * redistribution, maximum_redistribution)
    else:
        redistribution = max(redistribution, -maximum_redistribution)
    upper += redistribution
    whole -= redistribution

    upper_clip_correction = max(-upper, 0.0)
    whole_clip_correction = max(-whole, 0.0)
    upper = max(upper, 0.0)
    whole = max(whole, 0.0)

    # Preserve lines 229--230: ro2 already contains baseflow, even though
    # baseflow was removed explicitly above.
    whole += upper_runoff - lower_runoff - transpiration_tree

    actual_evapotranspiration = (
        interception
        + soil_evaporation
        + transpiration_tree
        + transpiration_c3
        + transpiration_c4
    )
    potential_evapotranspiration = (
        forcing.potential_soil_evaporation
        + forcing.potential_interception_tree
        + forcing.potential_interception_c3
        + forcing.potential_interception_c4
        + forcing.potential_transpiration_tree
        + forcing.potential_transpiration_c3
        + forcing.potential_transpiration_c4
    )
    new_state = VISITCHydrologyState(snow_water, upper, whole)
    fluxes = VISITCHydrologyFluxes(
        snow_fraction=snow_fraction,
        snowfall=snowfall,
        rainfall=rainfall,
        thaw=thaw,
        interception_tree=interception_tree,
        interception_c3=interception_c3,
        interception_c4=interception_c4,
        liquid_input_to_soil=liquid_input,
        upper_runoff=upper_runoff,
        soil_evaporation=soil_evaporation,
        transpiration_tree=transpiration_tree,
        transpiration_c3=transpiration_c3,
        transpiration_c4=transpiration_c4,
        baseflow=baseflow,
        lower_bucket_runoff=lower_bucket_runoff,
        lower_runoff=lower_runoff,
        redistribution=redistribution,
        actual_evapotranspiration=actual_evapotranspiration,
        potential_evapotranspiration=potential_evapotranspiration,
        upper_clip_correction=upper_clip_correction,
        whole_clip_correction=whole_clip_correction,
    )
    storage_change = new_state.total_storage - state.total_storage
    external_balance = (
        forcing.precipitation - actual_evapotranspiration - lower_runoff
    )
    residual = storage_change - external_balance
    expected = -baseflow + upper_clip_correction + whole_clip_correction
    return VISITCHydrologyStep(
        state=new_state,
        fluxes=fluxes,
        water_budget_residual=residual,
        expected_source_residual=expected,
    )


__all__ = [
    "VISITCHydrologyFluxes",
    "VISITCHydrologyForcing",
    "VISITCHydrologyParameters",
    "VISITCHydrologyState",
    "VISITCHydrologyStep",
    "visitc_hydrology_step",
]
