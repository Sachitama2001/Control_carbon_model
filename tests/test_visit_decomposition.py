import numpy as np
import pytest

from control_carbon.visit_decomposition import (
    DecompositionTemperatureMode,
    HumusBranchPolicy,
    VISIT_DECOMPOSITION_SOURCE,
    VISIT_STYPE_CALL_SOURCE,
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
    VISITSourceUndefinedError,
    visit_decomposition_scalars,
    visit_humus_decomposition_scalars,
    visit_litter_decomposition_scalar,
)
from control_carbon.visit_source_map import VISIT_SOURCE_COMMIT


@pytest.fixture
def environment():
    return VISITDecompositionEnvironment(
        soil_temperature_10cm=20.0,
        soil_temperature_deep=20.0,
        soil_water_upper=50.0,
        soil_water_lower=50.0,
        soil_aperture_upper=1.0,
        soil_aperture_whole=1.0,
        field_capacity_30cm=100.0,
        field_capacity=100.0,
    )


@pytest.fixture
def parameters():
    return VISITDecompositionParameters(kml=1.0, kmh=1.0, kmsl=1.0, kmsh=1.0)


def lloyd_taylor(temperature, reference_offset, temperature_offset):
    return 0.01 + np.exp(
        308.56 * (1.0 / reference_offset - 1.0 / (temperature + temperature_offset))
    )


def test_baseline_frl_and_frh_match_source_equations(environment, parameters):
    moisture_factor = 0.8 * 50.0 / (1.0 * 100.0 + 50.0) + 0.2
    expected = lloyd_taylor(20.0, 56.02, 46.02) * moisture_factor

    result = visit_decomposition_scalars(environment, parameters)

    assert result.litter == pytest.approx(expected)
    assert np.allclose(result.humus, expected)
    assert result.mode is DecompositionTemperatureMode.BASELINE


def test_pool_specific_mode_selects_litter_stype_one(environment, parameters):
    moisture_factor = 0.8 * 50.0 / (1.0 * 100.0 + 50.0) + 0.2
    expected = lloyd_taylor(20.0, 66.02, 56.02) * moisture_factor

    result = visit_litter_decomposition_scalar(
        environment, parameters, DecompositionTemperatureMode.POOL_SPECIFIC
    )

    assert result == pytest.approx(expected)


def test_source_policy_rejects_undefined_frh_stype_one_and_two(
    environment, parameters
):
    with pytest.raises(VISITSourceUndefinedError, match="uninitialized"):
        visit_humus_decomposition_scalars(
            environment,
            parameters,
            DecompositionTemperatureMode.POOL_SPECIFIC,
        )


def test_inferred_fix_maps_humus_stypes_to_three_source_formulas(
    environment, parameters
):
    moisture_factor = 0.8 * 50.0 / (1.0 * 100.0 + 50.0) + 0.2
    expected = np.asarray(
        [
            lloyd_taylor(20.0, 56.02, 46.02),
            lloyd_taylor(20.0, 50.02, 40.02),
            lloyd_taylor(20.0, 44.02, 34.02),
        ]
    ) * moisture_factor

    result = visit_humus_decomposition_scalars(
        environment,
        parameters,
        DecompositionTemperatureMode.POOL_SPECIFIC,
        branch_policy=HumusBranchPolicy.INFERRED_STYPE_FIX,
    )

    assert np.allclose(result, expected)


def test_temperature_threshold_uses_source_floor(environment, parameters):
    frozen = VISITDecompositionEnvironment(
        **{
            **environment.__dict__,
            "soil_temperature_10cm": -20.0,
            "soil_temperature_deep": -30.0,
        }
    )
    moisture_factor = 0.8 * 50.0 / (1.0 * 100.0 + 50.0) + 0.2

    result = visit_decomposition_scalars(frozen, parameters)

    assert result.litter == pytest.approx(0.01 * moisture_factor)
    assert np.allclose(result.humus, 0.01 * moisture_factor)


def test_stype_provenance_and_input_domain(environment):
    assert VISIT_STYPE_CALL_SOURCE.commit == VISIT_SOURCE_COMMIT
    assert VISIT_DECOMPOSITION_SOURCE.commit == VISIT_SOURCE_COMMIT
    with pytest.raises(ValueError, match="field capacities"):
        VISITDecompositionEnvironment(
            **{**environment.__dict__, "field_capacity": 0.0}
        )