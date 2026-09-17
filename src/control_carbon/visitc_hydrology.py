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

from dataclasses import dataclass, replace
import math

import numpy as np

from .visitc_source_map import (
    VISITC_HYDROLOGY_SOURCE,
    VISITC_HYDRO_FLOWS_SOURCE,
    VISITC_RADIATION_SOURCE,
    VISITC_CANOPY_CONDUCTANCE_SOURCE,
)


@dataclass(frozen=True)
class VISITCCanopyRadiationParameters:
    """Inputs to VISITc's LAI-dependent ``f_net_rad`` partition."""

    leaf_area_index: tuple[float, float, float]
    extinction_initial: tuple[float, float, float]
    extinction_radiation: tuple[float, float, float]
    albedo: tuple[float, float, float, float]
    c3_understory_fraction: float = 0.0
    c4_understory_fraction: float = 0.0

    def __post_init__(self) -> None:
        for name, size in (
            ("leaf_area_index", 3),
            ("extinction_initial", 3),
            ("extinction_radiation", 3),
            ("albedo", 4),
        ):
            values = np.asarray(getattr(self, name), dtype=float)
            if values.shape != (size,) or not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must be finite with shape ({size},)")
            if np.any(values < 0.0):
                raise ValueError(f"{name} must be nonnegative")
        fractions = np.asarray(
            (self.c3_understory_fraction, self.c4_understory_fraction), dtype=float
        )
        if (
            not np.all(np.isfinite(fractions))
            or np.any(fractions < 0.0)
            or fractions.sum() > 1.0
        ):
            raise ValueError("understory fractions must be nonnegative and sum to <= 1")
        if np.any(np.asarray(self.albedo) > 1.0):
            raise ValueError("albedo values must lie in [0, 1]")


@dataclass(frozen=True)
class VISITCNetRadiation:
    """Layer net radiation and source intermediate diagnostics, W m-2."""

    tree: float
    c3: float
    c4: float
    ground: float
    ecosystem: float
    cover_fractions: tuple[float, float, float, float]
    absorbed_shortwave: tuple[float, float, float, float]
    partitioned_longwave: tuple[float, float, float, float]
    source: object = VISITC_RADIATION_SOURCE


@dataclass(frozen=True)
class VISITCPenmanMonteithEnvironment:
    """Prepared daily atmosphere, radiation, and conductance inputs.

    Temperatures are degree Celsius, pressure and vapour quantities hPa,
    wind is m s-1, day length is hours, radiation is W m-2, and canopy
    conductance is mmol H2O m-2 s-1.
    """

    air_temperature: float
    surface_temperature: float
    air_pressure: float
    vapor_pressure: float
    vapor_pressure_deficit: float
    wind_speed: float
    day_length: float
    incoming_shortwave: float
    cloud_fraction: float
    canopy_conductance_tree: float
    canopy_conductance_c3: float
    canopy_conductance_c4: float
    radiation: VISITCCanopyRadiationParameters

    def __post_init__(self) -> None:
        numeric = np.asarray(
            tuple(
                value
                for name, value in self.__dict__.items()
                if name != "radiation"
            ),
            dtype=float,
        )
        if not np.all(np.isfinite(numeric)):
            raise ValueError("Penman-Monteith environment must be finite")
        if self.air_pressure <= 0.0:
            raise ValueError("air_pressure must be positive")
        for name in (
            "vapor_pressure",
            "vapor_pressure_deficit",
            "wind_speed",
            "day_length",
            "incoming_shortwave",
            "canopy_conductance_tree",
            "canopy_conductance_c3",
            "canopy_conductance_c4",
        ):
            if getattr(self, name) < 0.0:
                raise ValueError(f"{name} must be nonnegative")
        if not 0.0 <= self.cloud_fraction <= 1.0:
            raise ValueError("cloud_fraction must lie in [0, 1]")


@dataclass(frozen=True)
class VISITCPenmanMonteithFluxes:
    """Source-exact PM potentials and their prepared dependencies."""

    saturated_vapor_pressure: float
    saturation_vapor_pressure_slope: float
    air_density: float
    aerodynamic_resistance: float
    soil_resistance: float
    net_radiation: VISITCNetRadiation
    potential_interception_tree: float
    potential_interception_c3: float
    potential_interception_c4: float
    potential_soil_evaporation: float
    potential_transpiration_tree: float
    potential_transpiration_c3: float
    potential_transpiration_c4: float
    source: object = VISITC_HYDRO_FLOWS_SOURCE


@dataclass(frozen=True)
class VISITCCanopyConductanceParameters:
    """Direct inputs to ``ecophysiology.c::f_canopy_cond``."""

    photosynthetic_capacity: float
    radiation_extinction: float
    light_use_efficiency: float
    canopy_top_ppfd: float
    leaf_area_index: float
    atmospheric_co2: float
    co2_compensation_point: float
    vapor_pressure_deficit: float
    minimum_stomatal_conductance: float
    ball_berry_slope: float
    vpd_scale: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("canopy-conductance inputs must be finite")
        if self.radiation_extinction <= 0.0 or self.vpd_scale <= 0.0:
            raise ValueError("radiation_extinction and vpd_scale must be positive")
        if self.leaf_area_index < 0.0 or self.vapor_pressure_deficit < 0.0:
            raise ValueError("LAI and VPD must be nonnegative")
        if self.atmospheric_co2 == self.co2_compensation_point:
            raise ValueError("atmospheric CO2 must differ from compensation point")


@dataclass(frozen=True)
class VISITCCanopyConductance:
    gross_photosynthesis_proxy: float
    conductance: float
    co2_factor: float
    vpd_factor: float
    source: object = VISITC_CANOPY_CONDUCTANCE_SOURCE


def visitc_leaf_area_index(
    foliage_carbon: float, specific_leaf_area: float
) -> float:
    """Transcribe ``lai_mass`` from foliage C and SLA (cm2 gDM-1)."""
    values = np.asarray((foliage_carbon, specific_leaf_area), dtype=float)
    if not np.all(np.isfinite(values)) or specific_leaf_area < 0.0:
        raise ValueError("foliage carbon/SLA must be finite and SLA nonnegative")
    # dmTc=2.2, then source conversions /100 and /2 for one-sided area.
    return max(specific_leaf_area * foliage_carbon * 2.2 / 100.0 / 2.0, 0.0)


def visitc_irradiance_extinction(
    initial_extinction: float, solar_height_degrees: float
) -> float:
    """Transcribe ``irr_attn`` including its minimum sine factor of 0.3."""
    values = np.asarray((initial_extinction, solar_height_degrees), dtype=float)
    if not np.all(np.isfinite(values)) or initial_extinction < 0.0:
        raise ValueError("extinction inputs must be finite and nonnegative")
    sine_height = np.clip(math.sin(solar_height_degrees * 0.0174533), 0.3, 1.0)
    return float(initial_extinction / sine_height)


def visitc_canopy_conductance(
    parameters: VISITCCanopyConductanceParameters,
    *,
    fixed_350_ppm_co2: bool = False,
) -> VISITCCanopyConductance:
    """Transcribe ``f_canopy_cond`` for one plant canopy.

    The source-local ``gpp`` is returned as a diagnostic.  It has the same
    source-specific canopy-rate basis used by the Ball--Berry-like expression;
    it is not the ecosystem daily GPP flux reported by ``daily_scheme``.
    """
    p = parameters
    if p.photosynthetic_capacity > 0.0:
        cc1 = 2.0 * p.photosynthetic_capacity / p.radiation_extinction
        bb = (
            p.radiation_extinction
            * p.light_use_efficiency
            * p.canopy_top_ppfd
            / p.photosynthetic_capacity
        )
        cc2 = 1.0 + math.sqrt(1.0 + bb)
        cc3 = 1.0 + math.sqrt(
            1.0 + bb * math.exp(-p.radiation_extinction * p.leaf_area_index)
        )
        gpp = cc1 * math.log(cc2 / cc3)
    else:
        gpp = 0.0
    co2_value = 350.0 if fixed_350_ppm_co2 else p.atmospheric_co2
    co2_factor = 1.0 / (co2_value - p.co2_compensation_point)
    vpd_factor = 1.0 / (1.0 + p.vapor_pressure_deficit / p.vpd_scale)
    conductance = (
        p.minimum_stomatal_conductance * p.leaf_area_index
        + p.ball_berry_slope * co2_factor * vpd_factor * gpp
    )
    return VISITCCanopyConductance(
        gross_photosynthesis_proxy=gpp,
        conductance=conductance,
        co2_factor=co2_factor,
        vpd_factor=vpd_factor,
    )


def visitc_saturated_vapor_pressure(air_temperature: float) -> float:
    """Transcribe ``f_vap_pre_sat`` (Tetens equation), returning hPa."""
    if not np.isfinite(air_temperature):
        raise ValueError("air_temperature must be finite")
    if air_temperature > 0.0:
        value = 6.1078 * 10.0 ** (
            7.5 * air_temperature / (237.3 + air_temperature)
        )
    else:
        value = 6.1078 * 10.0 ** (
            9.5 * air_temperature / (265.3 + air_temperature)
        )
    return max(value, 0.0)


def visitc_saturation_vapor_pressure_slope(air_temperature: float) -> float:
    """Transcribe ``f_slope_vps``, returning hPa K-1."""
    if not np.isfinite(air_temperature):
        raise ValueError("air_temperature must be finite")
    absolute_temperature = 273.15 + air_temperature
    if absolute_temperature <= 0.0:
        raise ValueError("air_temperature must exceed absolute zero")
    if air_temperature > 0.0:
        numerator = 6.1078 * (2500.0 - 2.4 * air_temperature)
        exponential = 10.0 ** (
            7.5 * air_temperature / (237.3 + air_temperature)
        )
    else:
        numerator = 6.1078 * 2834.0
        exponential = 10.0 ** (
            9.5 * air_temperature / (265.3 + air_temperature)
        )
    denominator = 0.4615 * absolute_temperature**2
    return numerator / denominator * exponential


def visitc_air_density(
    air_temperature: float, air_pressure: float, vapor_pressure: float
) -> float:
    """Transcribe ``f_airdens``, returning kg m-3."""
    values = np.asarray((air_temperature, air_pressure, vapor_pressure))
    if not np.all(np.isfinite(values)) or air_pressure <= 0.0 or vapor_pressure < 0.0:
        raise ValueError("air-density inputs must be finite with positive pressure")
    absolute_temperature = 273.15 + air_temperature
    if absolute_temperature <= 0.0:
        raise ValueError("air_temperature must exceed absolute zero")
    return (
        1.293
        * 273.15
        / absolute_temperature
        * air_pressure
        / 1013.25
        * (1.0 - 0.378 * vapor_pressure / air_pressure)
    )


def visitc_aerodynamic_resistance(wind_speed: float) -> float:
    """Transcribe ``f_r_aero``, returning s m-1 including source clipping."""
    if not np.isfinite(wind_speed) or wind_speed < 0.0:
        raise ValueError("wind_speed must be finite and nonnegative")
    wind = max(wind_speed, 0.1)
    resistance = math.log(10.0) ** 2 / (0.41**2 * wind)
    return min(max(resistance, 0.1), 59.5)


def visitc_net_radiation(
    environment: VISITCPenmanMonteithEnvironment,
) -> VISITCNetRadiation:
    """Transcribe LAI cover and ``radiation.c::f_net_rad`` daily algebra."""
    params = environment.radiation
    lai = np.asarray(params.leaf_area_index)
    extinction0 = np.asarray(params.extinction_initial)
    extinction = np.asarray(params.extinction_radiation)
    albedo = np.asarray(params.albedo)

    tree_cover = 1.0 - math.exp(-extinction0[0] * lai[0])
    c3_cover = (
        (1.0 - tree_cover)
        * params.c3_understory_fraction
        * (1.0 - math.exp(-extinction0[1] * lai[1]))
    )
    c4_cover = (
        (1.0 - tree_cover)
        * params.c4_understory_fraction
        * (1.0 - math.exp(-extinction0[2] * lai[2]))
    )
    ground_cover = 1.0 - tree_cover - c3_cover - c4_cover
    covers = np.array((tree_cover, c3_cover, c4_cover, ground_cover))

    stefan_boltzmann = 5.6703e-8
    upward_longwave = (
        0.95 * (environment.surface_temperature + 273.15) ** 4 * stefan_boltzmann
    )
    atmospheric_emissivity = 0.53 + 0.06 * math.sqrt(environment.vapor_pressure)
    air_emission = (environment.air_temperature + 273.15) ** 4 * stefan_boltzmann
    downward_longwave = (
        (1.0 - environment.cloud_fraction) * atmospheric_emissivity * air_emission
        + environment.cloud_fraction * (air_emission - 9.0)
    )
    ecosystem_longwave = downward_longwave - upward_longwave
    partitioned_longwave = ecosystem_longwave * covers

    incoming_shortwave = max(environment.incoming_shortwave, 0.0)
    transmit = 0.1
    fsw_tree = 1.0 - math.exp(-extinction[0] * lai[0] * (1.0 - transmit))
    fsw_c3 = (
        (1.0 - fsw_tree)
        * params.c3_understory_fraction
        * (1.0 - math.exp(-extinction[1] * lai[1]) * (1.0 - transmit))
    )
    fsw_c4 = (
        (1.0 - fsw_tree)
        * params.c4_understory_fraction
        * (1.0 - math.exp(-extinction[2] * lai[2]) * (1.0 - transmit))
    )
    fsw_ground = 1.0 - fsw_tree - fsw_c3 - fsw_c4
    shortwave_fractions = np.array((fsw_tree, fsw_c3, fsw_c4, fsw_ground))
    absorbed_shortwave = incoming_shortwave * shortwave_fractions * (1.0 - albedo)

    # Preserve the native sign/order literally: rn_layer = short - rn_long_layer.
    net = absorbed_shortwave - partitioned_longwave
    ecosystem = float(absorbed_shortwave.sum() - ecosystem_longwave)
    return VISITCNetRadiation(
        tree=float(net[0]),
        c3=float(net[1]),
        c4=float(net[2]),
        ground=float(net[3]),
        ecosystem=ecosystem,
        cover_fractions=tuple(float(value) for value in covers),
        absorbed_shortwave=tuple(float(value) for value in absorbed_shortwave),
        partitioned_longwave=tuple(float(value) for value in partitioned_longwave),
    )


def _penman_monteith_flux(
    *,
    slope: float,
    net_radiation: float,
    heat_capacity: float,
    air_density: float,
    vapor_pressure_deficit: float,
    aerodynamic_resistance: float,
    surface_resistance: float,
    day_length: float,
) -> float:
    psychrometric_constant = 0.667
    latent_heat = 695.0
    numerator = slope * net_radiation + (
        heat_capacity
        * air_density
        * vapor_pressure_deficit
        / aerodynamic_resistance
    )
    denominator = slope + psychrometric_constant * (
        1.0 + surface_resistance / aerodynamic_resistance
    )
    return max(day_length * numerator / denominator / latent_heat, 0.0)


def visitc_penman_monteith_fluxes(
    environment: VISITCPenmanMonteithEnvironment,
    *,
    upper_soil_water: float,
    field_capacity_upper: float,
) -> VISITCPenmanMonteithFluxes:
    """Evaluate source-exact PM potentials from native direct dependencies."""
    if (
        not np.isfinite(upper_soil_water)
        or upper_soil_water < 0.0
        or not np.isfinite(field_capacity_upper)
        or field_capacity_upper <= 0.0
    ):
        raise ValueError("upper soil water/capacity must be finite and valid")
    saturated = visitc_saturated_vapor_pressure(environment.air_temperature)
    slope = visitc_saturation_vapor_pressure_slope(environment.air_temperature)
    density = visitc_air_density(
        environment.air_temperature,
        environment.air_pressure,
        environment.vapor_pressure,
    )
    aerodynamic = visitc_aerodynamic_resistance(environment.wind_speed)
    radiation = visitc_net_radiation(environment)
    common = dict(
        slope=slope,
        heat_capacity=0.2813,
        air_density=density,
        vapor_pressure_deficit=environment.vapor_pressure_deficit,
        aerodynamic_resistance=aerodynamic,
        day_length=environment.day_length,
    )
    interception = tuple(
        _penman_monteith_flux(
            **common, net_radiation=value, surface_resistance=0.0
        )
        for value in (radiation.tree, radiation.c3, radiation.c4)
    )

    conductance_conversion = 0.0224 / 1000.0
    ground_conductance = 1000.0 * upper_soil_water / field_capacity_upper + 100.0
    soil_resistance = 1.0 / (ground_conductance * conductance_conversion)
    soil_evaporation = _penman_monteith_flux(
        slope=slope,
        net_radiation=radiation.ground,
        # pm_evap assigns cp twice; the second source assignment is effective.
        heat_capacity=1014.0,
        air_density=density,
        vapor_pressure_deficit=environment.vapor_pressure_deficit,
        aerodynamic_resistance=aerodynamic,
        surface_resistance=soil_resistance,
        day_length=environment.day_length,
    )
    transpiration = []
    for conductance, net_radiation in zip(
        (
            environment.canopy_conductance_tree,
            environment.canopy_conductance_c3,
            environment.canopy_conductance_c4,
        ),
        (radiation.tree, radiation.c3, radiation.c4),
    ):
        if conductance > 0.0 and net_radiation > 0.0:
            canopy_resistance = 1.0 / (
                conductance * conductance_conversion
            )
            value = _penman_monteith_flux(
                **common,
                net_radiation=net_radiation,
                surface_resistance=canopy_resistance,
            )
        else:
            value = 0.0
        transpiration.append(value)
    return VISITCPenmanMonteithFluxes(
        saturated_vapor_pressure=saturated,
        saturation_vapor_pressure_slope=slope,
        air_density=density,
        aerodynamic_resistance=aerodynamic,
        soil_resistance=soil_resistance,
        net_radiation=radiation,
        potential_interception_tree=interception[0],
        potential_interception_c3=interception[1],
        potential_interception_c4=interception[2],
        potential_soil_evaporation=soil_evaporation,
        potential_transpiration_tree=transpiration[0],
        potential_transpiration_c3=transpiration[1],
        potential_transpiration_c4=transpiration[2],
    )


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


def visitc_hydrology_step_with_penman_monteith(
    state: VISITCHydrologyState,
    forcing: VISITCHydrologyForcing,
    parameters: VISITCHydrologyParameters,
    environment: VISITCPenmanMonteithEnvironment,
) -> tuple[VISITCHydrologyStep, VISITCPenmanMonteithFluxes]:
    """Run one source-order hydrology day with transcribed PM potentials.

    Native ``f_hydrology`` receives atmospheric and radiation diagnostics that
    were prepared earlier by ``f_loct_proc``.  This wrapper preserves that
    boundary, computes all seven potential fluxes, and then passes them to the
    already isolated store update.
    """
    if not np.isclose(forcing.air_temperature, environment.air_temperature):
        raise ValueError(
            "forcing and Penman-Monteith air temperatures must be identical"
        )
    potentials = visitc_penman_monteith_fluxes(
        environment,
        upper_soil_water=state.upper_soil_water,
        field_capacity_upper=parameters.field_capacity_upper,
    )
    prepared_forcing = replace(
        forcing,
        potential_interception_tree=potentials.potential_interception_tree,
        potential_interception_c3=potentials.potential_interception_c3,
        potential_interception_c4=potentials.potential_interception_c4,
        potential_soil_evaporation=potentials.potential_soil_evaporation,
        potential_transpiration_tree=potentials.potential_transpiration_tree,
        potential_transpiration_c3=potentials.potential_transpiration_c3,
        potential_transpiration_c4=potentials.potential_transpiration_c4,
    )
    return visitc_hydrology_step(state, prepared_forcing, parameters), potentials


__all__ = [
    "VISITCHydrologyStep",
    "VISITCCanopyRadiationParameters",
    "VISITCCanopyConductance",
    "VISITCCanopyConductanceParameters",
    "VISITCNetRadiation",
    "VISITCPenmanMonteithEnvironment",
    "VISITCPenmanMonteithFluxes",
    "VISITCHydrologyFluxes",
    "VISITCHydrologyForcing",
    "VISITCHydrologyParameters",
    "VISITCHydrologyState",
    "visitc_aerodynamic_resistance",
    "visitc_air_density",
    "visitc_canopy_conductance",
    "visitc_hydrology_step",
    "visitc_hydrology_step_with_penman_monteith",
    "visitc_irradiance_extinction",
    "visitc_leaf_area_index",
    "visitc_net_radiation",
    "visitc_penman_monteith_fluxes",
    "visitc_saturated_vapor_pressure",
    "visitc_saturation_vapor_pressure_slope",
]
