import numpy as np
import pytest

from control_carbon.visit_soil import (
    VISIT_SOIL_PROCESS_SOURCE,
    VISITSoilParameters,
    soil_decomposition_fractions,
    visit_soil_daily_matrices,
    visit_soil_daily_step,
)
from control_carbon.visit_source_map import VISIT_SOURCE_COMMIT


@pytest.fixture
def soil_parameters():
    return VISITSoilParameters(
        sr_lf=1000.0,
        sr_lc=1000.0,
        sr_lr=1000.0,
        sr_ha=1000.0,
        sr_hi=1000.0,
        sr_hp=1000.0,
        f_co2_lf=0.5,
        f_co2_lc=0.6,
        f_co2_lr=0.7,
        f_hm_a=0.5,
        f_hm_i=0.3,
        f_hm_p=0.2,
    )


def test_visit_soil_one_day_hand_calculation(soil_parameters):
    state = np.arange(10.0, 100.0, 10.0)
    litter_inputs = np.ones(6)

    result = visit_soil_daily_step(
        state,
        litter_inputs,
        soil_parameters,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )

    assert np.allclose(result.fluxes.degradation, np.arange(1.0, 10.0))
    assert np.allclose(
        result.fluxes.microbial_respiration,
        [0.5, 1.2, 2.1, 2.0, 3.0, 4.2, 7.0, 8.0, 9.0],
    )
    assert np.allclose(result.fluxes.humus_formation.sum(axis=1), [4.0, 2.4, 1.6])
    assert result.fluxes.heterotrophic_respiration == pytest.approx(37.0)
    assert result.fluxes.unaccounted_carbon == pytest.approx(0.0)
    assert np.allclose(result.state, [10.0, 19.0, 28.0, 37.0, 46.0, 55.0, 67.0, 74.4, 82.6])
    assert result.state.sum() == pytest.approx(
        state.sum() + litter_inputs.sum() - result.fluxes.heterotrophic_respiration
    )
    assert not np.any(result.clipped_pools)


def test_visit_soil_matrix_matches_direct_algebra(soil_parameters):
    rng = np.random.default_rng(42)
    state = rng.uniform(1.0, 100.0, size=9)
    litter_inputs = rng.uniform(0.0, 1.0, size=6)
    humus_scalars = np.array([0.05, 0.1, 0.2])

    matrix_a, matrix_b = visit_soil_daily_matrices(
        soil_parameters, litter_scalar=0.1, humus_scalars=humus_scalars
    )
    result = visit_soil_daily_step(
        state,
        litter_inputs,
        soil_parameters,
        litter_scalar=0.1,
        humus_scalars=humus_scalars,
        clip_negative=False,
    )

    assert np.allclose(result.state, matrix_a @ state + matrix_b @ litter_inputs)


def test_source_humification_imbalance_is_reported_not_rejected(soil_parameters):
    imbalanced = VISITSoilParameters(
        **{
            **soil_parameters.__dict__,
            "f_hm_a": 0.4,
            "f_hm_i": 0.3,
            "f_hm_p": 0.2,
        }
    )
    result = visit_soil_daily_step(
        np.arange(10.0, 100.0, 10.0),
        np.ones(6),
        imbalanced,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )

    assert imbalanced.humification_total == pytest.approx(0.9)
    assert result.fluxes.unaccounted_carbon == pytest.approx(0.8)


def test_native_negative_pool_clipping_is_explicit(soil_parameters):
    state = np.ones(9)
    result = visit_soil_daily_step(
        state,
        np.zeros(6),
        soil_parameters,
        litter_scalar=2.0,
        humus_scalars=[2.0, 2.0, 2.0],
    )

    assert np.all(result.clipped_pools[:6])
    assert np.all(result.state >= 0.0)


def test_soil_inputs_and_provenance_are_validated(soil_parameters):
    assert np.allclose(
        soil_decomposition_fractions(soil_parameters, 0.1, [0.2, 0.3, 0.4]),
        [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4],
    )
    assert VISIT_SOIL_PROCESS_SOURCE.commit == VISIT_SOURCE_COMMIT
    with pytest.raises(ValueError, match="state"):
        visit_soil_daily_step(
            np.ones(8), np.ones(6), soil_parameters, 0.1, [0.1, 0.1, 0.1]
        )
    with pytest.raises(ValueError, match="partition"):
        VISITSoilParameters(
            **{**soil_parameters.__dict__, "f_co2_lf": 1.1}
        )