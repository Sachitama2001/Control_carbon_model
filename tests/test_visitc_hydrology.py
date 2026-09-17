import math

import numpy as np
import pytest

from control_carbon.visitc_hydrology import (
    VISITCCanopyRadiationParameters,
    VISITCCanopyConductanceParameters,
    VISITCPenmanMonteithEnvironment,
    VISITCHydrologyForcing,
    VISITCHydrologyParameters,
    VISITCHydrologyState,
    visitc_hydrology_step,
    visitc_aerodynamic_resistance,
    visitc_air_density,
    visitc_canopy_conductance,
    visitc_irradiance_extinction,
    visitc_leaf_area_index,
    visitc_penman_monteith_fluxes,
    visitc_saturated_vapor_pressure,
)
from control_carbon.visitc_source_map import VISITC_SOURCE_COMMIT


@pytest.fixture
def parameters():
    return VISITCHydrologyParameters(
        field_capacity_upper=100.0,
        field_capacity_whole=200.0,
        hydraulic_conductivity=0.0,
        tree_leaf_area_index=2.0,
        c3_leaf_area_index=1.0,
        c4_leaf_area_index=1.0,
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
    )


def test_zero_potential_fluxes_expose_baseflow_double_subtraction(parameters):
    result = visitc_hydrology_step(
        VISITCHydrologyState(0.0, 50.0, 100.0),
        VISITCHydrologyForcing(
            precipitation=0.0,
            air_temperature=10.0,
            deep_soil_temperature=5.0,
        ),
        parameters,
    )

    assert result.fluxes.baseflow == pytest.approx(0.1)
    assert result.fluxes.lower_bucket_runoff == pytest.approx(0.0)
    assert result.state.whole_soil_water == pytest.approx(99.8)
    assert result.water_budget_residual == pytest.approx(-0.1)
    assert result.water_budget_residual == pytest.approx(
        result.expected_source_residual
    )
    assert result.exposes_baseflow_double_subtraction
    assert result.source.commit == VISITC_SOURCE_COMMIT


def test_frozen_deep_soil_disables_baseflow_and_closes_budget(parameters):
    result = visitc_hydrology_step(
        VISITCHydrologyState(0.0, 50.0, 100.0),
        VISITCHydrologyForcing(
            precipitation=0.0,
            air_temperature=-5.0,
            deep_soil_temperature=0.0,
        ),
        parameters,
    )

    assert result.fluxes.baseflow == 0.0
    assert result.water_budget_residual == pytest.approx(0.0, abs=1e-12)
    assert not result.exposes_baseflow_double_subtraction


def test_interception_uses_source_quadratic_and_c4_availability(parameters):
    forcing = VISITCHydrologyForcing(
        precipitation=10.0,
        air_temperature=20.0,
        deep_soil_temperature=0.0,
        potential_interception_tree=0.3,
        potential_interception_c3=0.2,
        potential_interception_c4=0.1,
    )
    result = visitc_hydrology_step(
        VISITCHydrologyState(0.0, 40.0, 80.0), forcing, parameters
    )

    rain = result.fluxes.rainfall
    tree_supply = min(rain, 0.5)
    expected_tree = (
        tree_supply
        + 0.3
        - math.sqrt((tree_supply + 0.3) ** 2 - 4 * 0.85 * tree_supply * 0.3)
    ) / 1.7
    assert result.fluxes.interception_tree == pytest.approx(expected_tree)
    # Native C3 and C4 branches each subtract tree interception, not each other.
    assert result.fluxes.interception_c3 > 0.0
    assert result.fluxes.interception_c4 > 0.0
    assert result.fluxes.liquid_input_to_soil == pytest.approx(
        result.fluxes.rainfall - result.fluxes.interception + result.fluxes.thaw
    )


def test_qhb_uses_threefold_baseflow_rate(parameters):
    qhb = VISITCHydrologyParameters(**{**parameters.__dict__, "site_id": "QHB"})
    result = visitc_hydrology_step(
        VISITCHydrologyState(0.0, 50.0, 100.0),
        VISITCHydrologyForcing(0.0, 10.0, 5.0),
        qhb,
    )

    assert result.fluxes.baseflow == pytest.approx(0.3)
    assert result.water_budget_residual == pytest.approx(-0.3)


def test_hydrology_inputs_are_validated(parameters):
    with pytest.raises(ValueError, match="state"):
        VISITCHydrologyState(0.0, -1.0, 1.0)
    with pytest.raises(ValueError, match="forcing"):
        VISITCHydrologyForcing(-1.0, 10.0, 5.0)
    with pytest.raises(ValueError, match="field capacities"):
        VISITCHydrologyParameters(0.0, 100.0, 0.0)
    with pytest.raises(ValueError, match="fractions"):
        VISITCHydrologyParameters(10.0, 100.0, 0.0, c3_understory_fraction=1.1)
    with pytest.raises(ValueError, match="sum"):
        VISITCHydrologyParameters(
            10.0,
            100.0,
            0.0,
            c3_understory_fraction=0.7,
            c4_understory_fraction=0.6,
        )


def test_atmospheric_helpers_match_hand_computable_source_values():
    assert visitc_saturated_vapor_pressure(0.0) == pytest.approx(6.1078)
    assert visitc_air_density(0.0, 1013.25, 0.0) == pytest.approx(1.293)
    # The unclipped log-law value exceeds the native upper bound at 0.1 m/s.
    assert visitc_aerodynamic_resistance(0.0) == pytest.approx(59.5)


def test_penman_monteith_zero_conductance_disables_transpiration():
    radiation = VISITCCanopyRadiationParameters(
        leaf_area_index=(2.0, 1.0, 0.5),
        extinction_initial=(0.5, 0.6, 0.7),
        extinction_radiation=(0.45, 0.55, 0.65),
        albedo=(0.12, 0.18, 0.20, 0.15),
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
    )
    environment = VISITCPenmanMonteithEnvironment(
        air_temperature=15.0,
        surface_temperature=17.0,
        air_pressure=1000.0,
        vapor_pressure=12.0,
        vapor_pressure_deficit=8.0,
        wind_speed=2.5,
        day_length=12.0,
        incoming_shortwave=220.0,
        cloud_fraction=0.4,
        canopy_conductance_tree=0.0,
        canopy_conductance_c3=0.0,
        canopy_conductance_c4=0.0,
        radiation=radiation,
    )
    result = visitc_penman_monteith_fluxes(
        environment, upper_soil_water=50.0, field_capacity_upper=100.0
    )
    assert result.potential_transpiration_tree == 0.0
    assert result.potential_transpiration_c3 == 0.0
    assert result.potential_transpiration_c4 == 0.0
    assert result.soil_resistance == pytest.approx(1.0 / (600.0 * 0.0000224))
    assert sum(result.net_radiation.cover_fractions) == pytest.approx(1.0)


def test_canopy_conductance_preserves_lai_gpp_co2_and_vpd_dependencies():
    parameters = VISITCCanopyConductanceParameters(
        photosynthetic_capacity=15.0,
        radiation_extinction=0.5,
        light_use_efficiency=0.05,
        canopy_top_ppfd=1000.0,
        leaf_area_index=2.0,
        atmospheric_co2=410.0,
        co2_compensation_point=40.0,
        vapor_pressure_deficit=10.0,
        minimum_stomatal_conductance=10.0,
        ball_berry_slope=9.0,
        vpd_scale=15.0,
    )
    result = visitc_canopy_conductance(parameters)
    assert result.gross_photosynthesis_proxy > 0.0
    assert result.co2_factor == pytest.approx(1.0 / 370.0)
    assert result.vpd_factor == pytest.approx(0.6)
    assert result.conductance == pytest.approx(
        20.0 + 9.0 * result.gross_photosynthesis_proxy / 370.0 * 0.6
    )


def test_lai_and_irradiance_extinction_match_source_hand_calculations():
    assert visitc_leaf_area_index(2.0, 20.0) == pytest.approx(0.44)
    assert visitc_irradiance_extinction(0.5, 90.0) == pytest.approx(0.5)
    assert visitc_irradiance_extinction(0.5, 0.0) == pytest.approx(0.5 / 0.3)
