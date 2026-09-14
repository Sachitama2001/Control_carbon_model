import json

import numpy as np
import pytest

from control_carbon.coupled_carbon_water import (
    COUPLED_STATE_NAMES,
    WATER_INCIDENCE_MATRIX,
    CoupledCarbonWaterForcing,
    CoupledCarbonWaterParameters,
    coupled_equilibrium,
    coupled_jacobian_blocks,
    coupled_model_provenance_manifest,
    decompose_capacity_change,
    evaluate_coupled_carbon_water,
    instantaneous_carbon_capacity,
    schur_complement_carbon_jacobian,
    simulate_coupled_carbon_water,
)


@pytest.fixture
def parameters():
    return CoupledCarbonWaterParameters()


@pytest.fixture
def forcing():
    return CoupledCarbonWaterForcing(precipitation=1.5)


@pytest.fixture
def equilibrium(parameters, forcing):
    result = coupled_equilibrium(
        forcing,
        [2.0, 35.0, 6.0, 120.0, 0.3, 7.0, 3.0, 155.0],
        parameters,
    )
    assert result.converged
    return result


def test_incidence_matrix_conserves_every_internal_water_flux():
    np.testing.assert_array_equal(
        np.ones(4) @ WATER_INCIDENCE_MATRIX, np.zeros(3)
    )


def test_rhs_has_explicit_carbon_and_water_budget_closure(parameters, forcing):
    state = np.array([2.0, 30.0, 7.0, 120.0, 0.4, 7.0, 3.0, 150.0])
    result = evaluate_coupled_carbon_water(state, forcing, parameters)

    assert result.derivative.shape == (8,)
    assert result.carbon_budget_residual == pytest.approx(0.0, abs=1e-14)
    assert result.water_budget_residual == pytest.approx(0.0, abs=1e-14)
    assert np.any(result.fluxes.water_internal != 0.0)


def test_carbon_matrix_uses_donor_columns_and_internal_turnover(parameters, forcing):
    state = np.array([2.0, 30.0, 7.0, 120.0, 0.4, 7.0, 3.0, 150.0])
    evaluation = evaluate_coupled_carbon_water(state, forcing, parameters)
    matrix = evaluation.fluxes.carbon_matrix

    assert np.all(np.diag(matrix) < 0.0)
    np.testing.assert_allclose(matrix[:, :3].sum(axis=0), 0.0)
    assert matrix[:, 3].sum() < 0.0
    np.testing.assert_allclose(
        evaluation.fluxes.plant_turnover, matrix[3, :3] * state[:3]
    )


def test_frozen_capacity_solves_carbon_balance(parameters, forcing):
    state = np.array([2.0, 30.0, 7.0, 120.0, 0.4, 7.0, 3.0, 150.0])
    capacity = instantaneous_carbon_capacity(state, forcing, parameters)

    np.testing.assert_allclose(
        capacity.carbon_matrix @ capacity.state + capacity.carbon_input,
        0.0,
        atol=1e-14,
    )
    assert capacity.total == pytest.approx(
        capacity.input_rate * capacity.ecosystem_transit_time
    )


def test_capacity_change_decomposition_is_exact(parameters, forcing):
    wet_state = np.array([2.0, 30.0, 7.0, 120.0, 0.7, 8.0, 3.2, 160.0])
    dry_state = wet_state.copy()
    dry_state[4:] *= 0.55
    baseline = instantaneous_carbon_capacity(wet_state, forcing, parameters)
    comparison = instantaneous_carbon_capacity(dry_state, forcing, parameters)
    result = decompose_capacity_change(baseline, comparison)

    assert result.reconstructed_change == pytest.approx(result.direct_change)
    assert result.productivity_effect != 0.0
    assert result.transit_time_effect != 0.0


def test_coupled_equilibrium_is_stable_and_matches_frozen_capacity(
    equilibrium, parameters, forcing
):
    assert equilibrium.residual_norm < 1e-9
    assert equilibrium.stability == "stable"
    assert np.all(equilibrium.eigenvalues.real < 0.0)
    capacity = instantaneous_carbon_capacity(
        equilibrium.state, forcing, parameters
    )
    np.testing.assert_allclose(equilibrium.state[:4], capacity.state, rtol=1e-8)


def test_jacobian_blocks_capture_both_coupling_directions(
    equilibrium, parameters, forcing
):
    blocks = coupled_jacobian_blocks(equilibrium.state, forcing, parameters)

    assert blocks.full.shape == (8, 8)
    np.testing.assert_allclose(blocks.full[:4, 4:], blocks.carbon_water)
    np.testing.assert_allclose(blocks.full[4:, :4], blocks.water_carbon)
    assert np.linalg.norm(blocks.carbon_water) > 0.0
    assert np.linalg.norm(blocks.water_carbon) > 0.0


def test_schur_complement_matches_quasistatic_elimination(
    equilibrium, parameters, forcing
):
    blocks = coupled_jacobian_blocks(equilibrium.state, forcing, parameters)
    effective = schur_complement_carbon_jacobian(blocks)
    carbon_anomaly = np.array([0.1, -0.2, 0.05, 0.3])
    water_anomaly = -np.linalg.solve(
        blocks.water_water, blocks.water_carbon @ carbon_anomaly
    )

    np.testing.assert_allclose(
        effective @ carbon_anomaly,
        blocks.carbon_carbon @ carbon_anomaly
        + blocks.carbon_water @ water_anomaly,
    )


def test_equilibrium_is_stationary_under_direct_integration(
    equilibrium, parameters, forcing
):
    trajectory = simulate_coupled_carbon_water(
        equilibrium.state, np.linspace(0.0, 20.0, 41), forcing, parameters
    )

    np.testing.assert_allclose(
        trajectory.states,
        np.broadcast_to(equilibrium.state, trajectory.states.shape),
        atol=2e-9,
    )


def test_drydown_and_rewetting_modify_water_then_carbon(
    equilibrium, parameters, forcing
):
    def protocol(time):
        precipitation = 0.15 if 20.0 <= time < 80.0 else forcing.precipitation
        return CoupledCarbonWaterForcing(precipitation=precipitation)

    times = np.linspace(0.0, 160.0, 321)
    trajectory = simulate_coupled_carbon_water(
        equilibrium.state, times, protocol, parameters, max_step=0.2
    )

    assert np.all(trajectory.states >= 0.0)
    assert trajectory.states[160, 7] < equilibrium.state[7]
    assert trajectory.states[160, 0] < equilibrium.state[0]
    assert trajectory.states[-1, 7] > trajectory.states[160, 7]


def test_time_step_refinement_is_converged(equilibrium, parameters, forcing):
    times = np.linspace(0.0, 20.0, 41)

    def protocol(time):
        return CoupledCarbonWaterForcing(
            precipitation=forcing.precipitation,
            potential_transpiration=2.2 if time >= 5.0 else 1.8,
        )

    coarse = simulate_coupled_carbon_water(
        equilibrium.state, times, protocol, parameters, max_step=0.2
    )
    fine = simulate_coupled_carbon_water(
        equilibrium.state, times, protocol, parameters, max_step=0.1
    )
    assert np.max(np.abs(coarse.states - fine.states)) < 2e-5


def test_provenance_is_serializable_and_labels_reduction():
    manifest = coupled_model_provenance_manifest()

    assert len(COUPLED_STATE_NAMES) == 8
    assert "reduced" in manifest["approximation_level"]
    assert "literature" in manifest
    assert json.loads(json.dumps(manifest)) == manifest


def test_invalid_inputs_are_rejected(parameters, forcing):
    with pytest.raises(ValueError, match="state"):
        evaluate_coupled_carbon_water(np.ones(7), forcing, parameters)
    with pytest.raises(ValueError, match="nonnegative"):
        evaluate_coupled_carbon_water(-np.ones(8), forcing, parameters)
    with pytest.raises(ValueError, match="sum to one"):
        CoupledCarbonWaterParameters(carbon_allocation=(0.2, 0.2, 0.2))
    with pytest.raises(ValueError, match="relative_step"):
        coupled_jacobian_blocks(np.ones(8), forcing, parameters, relative_step=0.0)
