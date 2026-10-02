import numpy as np
import pytest
from numpy.testing import assert_allclose

from control_carbon.cover_carbon_baseline import (
    CarbonParameters,
    carbon_budget,
    carbon_capacity,
    carbon_matrix,
    carbon_potential,
    carbon_rhs,
    interior_lifted_equilibria,
    jacobian,
    lifted_equilibrium,
    rhs,
)
from control_carbon.forest_grass_tipping import interior_roots


@pytest.fixture
def carbon_parameters():
    return CarbonParameters(
        production=[1.1, 0.8],
        allocation=[[0.5, 0.3, 0.2], [0.4, 0.35, 0.25]],
        turnover=[[0.25, 0.12, 0.08], [0.3, 0.15, 0.1]],
        transfer_emission=[[0.2, 0.1, 0.0], [0.0, 0.0, 0.0]],
        decomposition=[0.4, 0.5],
        humification=[0.3, 0.25],
        humus_loss=[0.02, 0.03],
    )


def test_cover_bistability_lifts_to_two_stable_carbon_equilibria(carbon_parameters):
    roots = interior_roots()
    lifted = interior_lifted_equilibria(carbon_parameters)

    assert_allclose([state[0] for state in lifted], roots, atol=0.0)
    assert len(lifted) == 3
    assert_allclose(rhs(lifted[0], carbon_parameters), 0.0, atol=2e-14)
    assert_allclose(rhs(lifted[1], carbon_parameters), 0.0, atol=2e-14)
    assert_allclose(rhs(lifted[2], carbon_parameters), 0.0, atol=2e-14)

    eigenvalues = [np.linalg.eigvals(jacobian(state, carbon_parameters))
                   for state in lifted]
    assert np.all(eigenvalues[0].real < 0)
    assert np.count_nonzero(eigenvalues[1].real > 0) == 1
    assert np.all(eigenvalues[2].real < 0)


def test_carbon_matrix_is_hurwitz_and_capacity_is_its_equilibrium(carbon_parameters):
    for g in [0.0, 0.2, 0.7, 1.0]:
        matrix = carbon_matrix(g, carbon_parameters)
        capacity = carbon_capacity(g, carbon_parameters)
        assert np.all(np.linalg.eigvals(matrix).real < 0)
        assert np.all(capacity >= 0)
        assert_allclose(carbon_rhs(g, capacity, carbon_parameters), 0.0, atol=2e-14)

    g = 0.37
    capacity = carbon_capacity(g, carbon_parameters)
    assert_allclose(carbon_potential(g, capacity, carbon_parameters), 0.0, atol=0.0)
    state = lifted_equilibrium(g, carbon_parameters)
    assert state.shape == (11,)


def test_carbon_budget_closes_and_distinguishes_external_losses(carbon_parameters):
    g = 0.42
    carbon = np.array([0.8, 0.5, 0.3, 1.2, 2.0, 0.4, 0.3, 0.2, 0.9, 1.5])
    budget = carbon_budget(g, carbon, carbon_parameters)

    assert budget["npp"] > 0
    assert budget["heterotrophic_respiration"] > 0
    assert budget["transition_emission"] > 0
    assert_allclose(
        budget["stock_change"],
        budget["npp"] - budget["heterotrophic_respiration"]
        - budget["transition_emission"],
        atol=2e-15,
    )
    assert abs(budget["residual"]) < 2e-15


def test_analytic_jacobian_matches_finite_difference(carbon_parameters):
    state = np.array([0.51, 0.7, 0.4, 0.2, 1.0, 1.6, 0.5, 0.3, 0.2, 0.8, 1.4])
    analytic = jacobian(state, carbon_parameters)
    step = 1e-6
    numeric = np.column_stack([
        (rhs(state + step * np.eye(11)[i], carbon_parameters)
         - rhs(state - step * np.eye(11)[i], carbon_parameters)) / (2 * step)
        for i in range(11)
    ])
    assert_allclose(analytic, numeric, rtol=2e-9, atol=2e-10)
    assert_allclose(analytic[0, 1:], 0.0, atol=0.0)


def test_nonnegative_boundary_signs_and_parameter_validation(carbon_parameters):
    zero_carbon = np.zeros(10)
    assert np.all(carbon_rhs(0.4, zero_carbon, carbon_parameters) >= 0)
    assert rhs(np.r_[0.0, zero_carbon], carbon_parameters)[0] > 0
    assert rhs(np.r_[1.0, zero_carbon], carbon_parameters)[0] == 0

    with pytest.raises(ValueError, match="lie in \\[0, 1\\]"):
        carbon_rhs(1.01, zero_carbon, carbon_parameters)
    with pytest.raises(ValueError, match="allocations"):
        CarbonParameters(
            production=[1.0, 1.0],
            allocation=[[0.5, 0.2, 0.2], [0.4, 0.3, 0.3]],
            turnover=np.ones((2, 3)),
            transfer_emission=np.zeros((2, 3)),
            decomposition=[0.2, 0.2],
            humification=[0.5, 0.5],
            humus_loss=[0.1, 0.1],
        )
