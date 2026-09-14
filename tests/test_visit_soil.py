import numpy as np
import pytest

from control_carbon.state_space import (
    discrete_modal_analysis,
    discrete_transfer_function_matrix,
)
from control_carbon.visit_soil import (
    VISIT_SOIL_PROCESS_SOURCE,
    VISITSoilParameters,
    compare_visit_soil_irf_to_daily_simulation,
    soil_decomposition_fractions,
    visit_soil_daily_matrices,
    visit_soil_daily_step,
    visit_soil_discrete_system,
    visit_soil_fixed_point,
    visit_soil_environmental_sensitivity,
    visit_soil_heterotrophic_respiration_irf,
)
from control_carbon.visit_decomposition import (
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
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


def test_soil_system_output_matches_direct_respiration(soil_parameters):
    state = np.arange(10.0, 100.0, 10.0)
    litter_inputs = np.ones(6)
    system = visit_soil_discrete_system(
        soil_parameters, litter_scalar=0.1, humus_scalars=[0.1, 0.1, 0.1]
    )
    states, outputs = system.forced_response(
        litter_inputs[np.newaxis, :], initial_state=state
    )
    direct = visit_soil_daily_step(
        state,
        litter_inputs,
        soil_parameters,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )

    assert np.allclose(states[1], direct.state)
    assert outputs[0, 0] == pytest.approx(
        direct.fluxes.heterotrophic_respiration
    )


def test_soil_fixed_point_balances_input_and_respiration(soil_parameters):
    litter_inputs = np.ones(6)
    result = visit_soil_fixed_point(
        soil_parameters,
        litter_inputs,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )
    daily = visit_soil_daily_step(
        result.state,
        litter_inputs,
        soil_parameters,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )

    assert np.allclose(result.state, [10.0] * 6 + [12.0, 7.2, 4.8])
    assert result.residual_norm < 1e-12
    assert result.stable
    assert np.allclose(daily.state, result.state)
    assert daily.fluxes.heterotrophic_respiration == pytest.approx(
        litter_inputs.sum()
    )


def test_respiration_irf_respects_native_daily_update_order(soil_parameters):
    response = visit_soil_heterotrophic_respiration_irf(
        soil_parameters,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
        n_steps=500,
    )

    assert np.allclose(response[0], 0.0)
    assert np.allclose(response[1], [0.05, 0.06, 0.07, 0.05, 0.06, 0.07])
    assert np.allclose(response.sum(axis=0), 1.0)


def test_irf_prediction_matches_direct_soil_trajectory(soil_parameters):
    perturbations = np.zeros((40, 6))
    perturbations[0, 0] = 1.0

    comparison = compare_visit_soil_irf_to_daily_simulation(
        soil_parameters,
        baseline_litter_inputs=np.ones(6),
        input_perturbations=perturbations,
        litter_scalar=0.1,
        humus_scalars=[0.1, 0.1, 0.1],
    )

    assert comparison.days.shape == (41,)
    assert np.array_equal(comparison.output_days, comparison.days[:-1])
    assert comparison.phase_axis_names == (
        "total litter carbon",
        "total humus carbon",
    )
    assert comparison.perturbed_phase.shape == (41, 2)
    assert np.allclose(comparison.state_error, 0.0)
    assert np.allclose(comparison.heterotrophic_respiration_error, 0.0)


def test_soil_qse_environmental_sensitivity_matches_resolved_fixed_points(
    soil_parameters,
):
    environment = VISITDecompositionEnvironment(
        soil_temperature_10cm=10.0,
        soil_temperature_deep=8.0,
        soil_water_upper=50.0,
        soil_water_lower=100.0,
        soil_aperture_upper=0.8,
        soil_aperture_whole=0.9,
        field_capacity_30cm=120.0,
        field_capacity=300.0,
    )
    decomposition_parameters = VISITDecompositionParameters(
        kml=0.22, kmh=0.1, kmsl=0.21, kmsh=0.07
    )

    sensitivity = visit_soil_environmental_sensitivity(
        soil_parameters,
        litter_inputs=np.ones(6),
        decomposition_parameters=decomposition_parameters,
        environment=environment,
    )

    assert np.max(np.abs(sensitivity.state_sensitivity_error)) < 1e-4
    assert np.max(
        np.abs(sensitivity.heterotrophic_respiration_sensitivity_error)
    ) < 1e-9
    assert np.max(
        np.abs(sensitivity.heterotrophic_respiration_sensitivity)
    ) < 1e-9
    upper_temperature = sensitivity.environment_names.index(
        "soil_temperature_10cm"
    )
    assert np.all(sensitivity.state_sensitivity[:6, upper_temperature] < 0.0)


def test_tky_soil_grouped_modes_reconstruct_rh_response():
    parameters = VISITSoilParameters(
        sr_lf=1.1,
        sr_lc=0.3,
        sr_lr=0.75,
        sr_ha=0.18,
        sr_hi=0.08,
        sr_hp=0.025,
        f_co2_lf=0.5,
        f_co2_lc=0.5,
        f_co2_lr=0.5,
        f_hm_a=0.6,
        f_hm_i=0.3,
        f_hm_p=0.1,
    )
    environment = VISITDecompositionEnvironment(
        soil_temperature_10cm=10.0,
        soil_temperature_deep=8.0,
        soil_water_upper=50.0,
        soil_water_lower=100.0,
        soil_aperture_upper=0.8,
        soil_aperture_whole=0.9,
        field_capacity_30cm=120.0,
        field_capacity=300.0,
    )
    decomposition_parameters = VISITDecompositionParameters(
        kml=0.22, kmh=0.1, kmsl=0.21, kmsh=0.07
    )
    from control_carbon.visit_decomposition import visit_decomposition_scalars

    scalars = visit_decomposition_scalars(environment, decomposition_parameters)
    system = visit_soil_discrete_system(
        parameters, scalars.litter, scalars.humus
    )

    modal = discrete_modal_analysis(system)

    assert modal.multiplicities.sum() == 9
    assert np.count_nonzero(modal.multiplicities == 2) == 3
    assert np.allclose(modal.state_participation.sum(axis=0), 1.0)
    assert np.allclose(modal.impulse_response(200), system.impulse_response(200))
    modal_dc_gain = modal.direct_term + np.sum(
        modal.residues / (1.0 - modal.poles[:, np.newaxis, np.newaxis]), axis=0
    )
    assert np.allclose(modal_dc_gain, discrete_transfer_function_matrix(system, 1.0))