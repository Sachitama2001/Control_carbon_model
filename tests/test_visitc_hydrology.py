import math

import numpy as np
import pytest

from control_carbon.visitc_hydrology import (
    VISITCHydrologyForcing,
    VISITCHydrologyParameters,
    VISITCHydrologyState,
    visitc_hydrology_step,
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
