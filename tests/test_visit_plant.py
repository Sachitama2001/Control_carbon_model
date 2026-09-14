from pathlib import Path

import numpy as np
import pytest

from control_carbon.visit_native import (
    build_native_visit_plant_allocation_bridge,
    build_native_visit_plant_respiration_bridge,
    build_native_visit_plant_turnover_bridge,
    run_native_visit_plant_respiration,
    run_native_visit_plant_allocation,
    run_native_visit_plant_turnover,
)
from control_carbon.visit_plant import (
    VISIT_LITTERFALL_SOURCE,
    VISIT_PLANT_TURNOVER_SOURCE,
    VISITPlantStructuralState,
    VISITPlantAllocationParameters,
    VISITPlantRespirationParameters,
    VISITPlantTurnoverContext,
    VISITPlantTurnoverParameters,
    visit_plant_turnover_step,
    visit_plant_allocation,
    visit_plant_respiration,
)
from control_carbon.visit_source_map import VISIT_SOURCE_COMMIT


VISIT_SOURCE_ROOT = Path(__file__).resolve().parents[2] / "VISIT-matrix" / "visit_local"


@pytest.fixture(scope="module")
def native_plant_bridge(tmp_path_factory):
    if not VISIT_SOURCE_ROOT.is_dir():
        pytest.skip("authoritative VISIT source checkout is unavailable")
    return build_native_visit_plant_turnover_bridge(
        VISIT_SOURCE_ROOT,
        tmp_path_factory.mktemp("native-plant") / "visit_plant_turnover_bridge",
    )


@pytest.fixture(scope="module")
def native_respiration_bridge(tmp_path_factory):
    if not VISIT_SOURCE_ROOT.is_dir():
        pytest.skip("authoritative VISIT source checkout is unavailable")
    return build_native_visit_plant_respiration_bridge(
        VISIT_SOURCE_ROOT,
        tmp_path_factory.mktemp("native-respiration")
        / "visit_plant_respiration_bridge",
    )


@pytest.fixture(scope="module")
def native_allocation_bridge(tmp_path_factory):
    if not VISIT_SOURCE_ROOT.is_dir():
        pytest.skip("authoritative VISIT source checkout is unavailable")
    return build_native_visit_plant_allocation_bridge(
        VISIT_SOURCE_ROOT,
        tmp_path_factory.mktemp("native-allocation")
        / "visit_plant_allocation_bridge",
    )


@pytest.fixture
def state():
    return VISITPlantStructuralState(foliage=10.0, stem=20.0, root=30.0)


@pytest.fixture
def parameters():
    return VISITPlantTurnoverParameters(
        leaf_rate=0.0004,
        stem_rate=0.00009,
        root_rate=0.00035,
        deciduous_shedding_fraction=0.12,
    )


def context(**changes):
    values = {
        "season": 1,
        "crop_stage": 0,
        "lai": 1.0,
        "deep_soil_water_potential": -1.0,
        "air_temperature": 10.0,
        "photosynthesis_type": 3,
    }
    return VISITPlantTurnoverContext(**{**values, **changes})


def test_tky_baseline_turnover_matches_source_products(state, parameters):
    result = visit_plant_turnover_step(state, parameters, context())

    assert np.allclose(result.litterfall, [0.004, 0.0018, 0.0105])
    assert np.allclose(result.state.as_array(), state.as_array() - result.litterfall)
    assert result.regime == "baseline-turnover"


def test_deciduous_shedding_uses_season_and_lai_branches(state, parameters):
    shedding = visit_plant_turnover_step(
        state, parameters, context(season=3, lai=1.0)
    )
    final_shedding = visit_plant_turnover_step(
        state, parameters, context(season=3, lai=0.05)
    )

    assert shedding.litterfall[0] == pytest.approx(1.2)
    assert final_shedding.litterfall[0] == pytest.approx(10.0)
    assert shedding.regime == "seasonal-shedding"
    assert final_shedding.regime == "seasonal-final-shedding"


def test_crop_stage_five_overrides_all_three_turnover_rates(state, parameters):
    result = visit_plant_turnover_step(
        state, parameters, context(crop_stage=5)
    )

    assert np.allclose(result.litterfall, 0.7 * state.as_array())
    assert result.regime == "crop-harvest"


def test_drought_adds_source_leaf_shedding_fraction(state, parameters):
    result = visit_plant_turnover_step(
        state,
        parameters,
        context(deep_soil_water_potential=-149.0),
    )

    assert result.litterfall[0] == pytest.approx((0.0004 + 0.012) * state.foliage)
    assert result.regime == "drought-enhanced-turnover"


def test_turnover_provenance_is_pinned():
    assert VISIT_LITTERFALL_SOURCE.commit == VISIT_SOURCE_COMMIT
    assert VISIT_PLANT_TURNOVER_SOURCE.commit == VISIT_SOURCE_COMMIT


@pytest.mark.parametrize(
    "turnover_context",
    (
        context(),
        context(season=3),
        context(crop_stage=5),
        context(deep_soil_water_potential=-149.0),
    ),
)
def test_turnover_branches_match_native_plant_process(
    native_plant_bridge, state, parameters, turnover_context
):
    native = run_native_visit_plant_turnover(
        native_plant_bridge, state, parameters, turnover_context
    )
    python = visit_plant_turnover_step(state, parameters, turnover_context)

    assert np.allclose(native.state.as_array(), python.state.as_array())
    assert np.allclose(native.litterfall, python.litterfall)


def test_tky_plant_respiration_matches_source_equations(state):
    parameters = VISITPlantRespirationParameters(
        growth_foliage=0.4,
        growth_stem=0.3,
        growth_root=0.35,
        maintenance_foliage_15c=1.45,
        maintenance_stem_sapwood=0.15,
        maintenance_stem_heartwood=0.01,
        maintenance_root_fine=0.75,
        maintenance_root_coarse=0.25,
        q10_foliage_base=2.0,
        q10_stem_base=2.0,
        q10_root_base=2.0,
        size_stem=20.0,
        size_root=10.0,
    )
    allocations = np.asarray((0.2, 0.3, 0.4))

    result = visit_plant_respiration(
        state,
        parameters,
        surface_temperature=15.0,
        upper_soil_temperature=15.0,
        allocation_fluxes=allocations,
    )

    stem_power = 1.0 - 0.33334 * 20.0 / 40.0
    root_power = 1.0 - 0.33334 * 30.0 / 40.0
    sapwood = min(20.0**stem_power, 20.0)
    fine_root = min(30.0**root_power, 30.0)
    stem_specific = (0.15 * sapwood + 0.01 * (20.0 - sapwood)) / 20.00001
    root_specific = (0.75 * fine_root + 0.25 * (30.0 - fine_root)) / 30.00001
    expected_maintenance = np.asarray(
        (10.0 * 1.45 / 1000.0, 20.0 * stem_specific / 1000.0, 30.0 * root_specific / 1000.0)
    )

    assert np.allclose(result.q10, 2.0)
    assert np.allclose(result.maintenance, expected_maintenance)
    assert np.allclose(result.growth, [0.08, 0.09, 0.14])


def test_nonpositive_epp_has_no_growth_respiration(state):
    parameters = VISITPlantRespirationParameters(
        0.4, 0.3, 0.35, 1.45, 0.15, 0.01, 0.75, 0.25, 2.0, 2.0, 2.0, 20.0, 10.0
    )

    result = visit_plant_respiration(
        state,
        parameters,
        surface_temperature=10.0,
        upper_soil_temperature=8.0,
        allocation_fluxes=np.asarray((-0.1, -0.2, -0.3)),
        positive_epp=False,
    )

    assert np.allclose(result.growth, 0.0)
    assert np.all(result.maintenance > 0.0)


def test_tky_plant_respiration_matches_native_functions(
    native_respiration_bridge, state
):
    parameters = VISITPlantRespirationParameters(
        0.4, 0.3, 0.35, 1.45, 0.15, 0.01, 0.75, 0.25, 2.0, 2.0, 2.0, 20.0, 10.0
    )
    allocations = np.asarray((0.2, 0.3, 0.4))
    native = run_native_visit_plant_respiration(
        native_respiration_bridge,
        state,
        parameters,
        surface_temperature=10.0,
        upper_soil_temperature=8.0,
        allocation_fluxes=allocations,
    )
    python = visit_plant_respiration(
        state,
        parameters,
        surface_temperature=10.0,
        upper_soil_temperature=8.0,
        allocation_fluxes=allocations,
    )

    for field in (
        "q10",
        "specific_maintenance",
        "maintenance",
        "growth",
        "sapwood_fine_root",
        "heartwood_coarse_root",
    ):
        assert np.allclose(getattr(native, field), getattr(python, field))


def test_positive_epp_allocation_respects_lai_season_and_crop():
    parameters = VISITPlantAllocationParameters(0.075, 0.67, 110.0)

    above = visit_plant_allocation(
        parameters,
        lai=4.0,
        optimum_lai=3.0,
        season=1,
        crop_stage=0,
        epp=2.0,
        gpp=2.5,
        maintenance_respiration=np.zeros(3),
    )
    dormant = visit_plant_allocation(
        parameters,
        lai=2.9,
        optimum_lai=3.0,
        season=0,
        crop_stage=4,
        epp=2.0,
        gpp=2.5,
        maintenance_respiration=np.zeros(3),
    )

    assert np.allclose(above.fractions, [0.0, 0.67, 0.33, 0.0])
    assert above.total_translocation == pytest.approx(2.0)
    assert dormant.fractions[3] == pytest.approx(0.7)
    assert dormant.fractions[:3].sum() == pytest.approx(1.0)
    assert dormant.fluxes[:3].sum() == pytest.approx(0.3 * 2.0)
    assert dormant.total_translocation == pytest.approx(2.0)


def test_nonpositive_epp_allocation_uses_gpp_minus_maintenance():
    parameters = VISITPlantAllocationParameters(0.075, 0.67, 110.0)
    maintenance = np.asarray((0.1, 0.2, 0.3))

    result = visit_plant_allocation(
        parameters,
        lai=2.0,
        optimum_lai=3.0,
        season=1,
        crop_stage=0,
        epp=-0.2,
        gpp=0.4,
        maintenance_respiration=maintenance,
    )

    expected_fractions = np.asarray((0.075, 0.925 * 0.67, 0.925 * 0.33, 0.0))
    assert np.allclose(result.fractions, expected_fractions)
    assert np.allclose(result.fluxes[:3], expected_fractions[:3] * 0.4 - maintenance)
    assert result.total_translocation == pytest.approx(-0.2)
    assert result.regime == "nonpositive-epp"


@pytest.mark.parametrize(
    "case",
    (
        dict(lai=4.0, optimum_lai=3.0, season=1, crop_stage=0, epp=2.0, gpp=2.5),
        dict(lai=2.9, optimum_lai=3.0, season=0, crop_stage=4, epp=2.0, gpp=2.5),
        dict(lai=2.0, optimum_lai=3.0, season=1, crop_stage=0, epp=-0.2, gpp=0.4),
    ),
)
def test_allocation_branches_match_native_f_allocation(
    native_allocation_bridge, case
):
    parameters = VISITPlantAllocationParameters(0.075, 0.67, 110.0)
    maintenance = np.asarray((0.1, 0.2, 0.3))
    native = run_native_visit_plant_allocation(
        native_allocation_bridge,
        parameters,
        maintenance_respiration=maintenance,
        **case,
    )
    python = visit_plant_allocation(
        parameters,
        maintenance_respiration=maintenance,
        **case,
    )

    assert np.allclose(native.fractions, python.fractions)
    assert np.allclose(native.fluxes, python.fluxes)
    assert native.total_translocation == pytest.approx(python.total_translocation)