import numpy as np
import pytest

from control_carbon.state_space import (
    ContinuousLTI,
    DiscreteLTI,
    discrete_fixed_point,
    discrete_modal_analysis,
    discrete_transfer_function_matrix,
    impulse_response,
    transfer_function_matrix,
)


def test_first_order_transfer_function_and_irf():
    # dx/dt = -2 x + 3 u; y = 4 x
    system = ContinuousLTI(
        A=np.array([[-2.0]]),
        B=np.array([[3.0]]),
        C=np.array([[4.0]]),
        D=np.array([[0.0]]),
    )

    # G(s) = 12/(s+2)
    assert np.allclose(transfer_function_matrix(system, 1.0), [[4.0]])

    times = np.array([0.0, 0.5, 1.0])
    h = impulse_response(system, times)[:, 0, 0]
    assert np.allclose(h, 12.0 * np.exp(-2.0 * times))


def test_frequency_response_shape():
    system = ContinuousLTI(
        A=np.diag([-1.0, -0.1]),
        B=np.ones((2, 2)),
        C=np.ones((1, 2)),
        D=np.zeros((1, 2)),
    )
    response = system.frequency_response([0.01, 1.0, 10.0])
    assert response.shape == (3, 1, 2)


def test_discrete_first_order_transfer_and_impulse_response():
    # x[k+1] = 0.5 x[k] + 2 u[k]; y[k] = 3 x[k] + 4 u[k]
    system = DiscreteLTI(
        A=np.array([[0.5]]),
        B=np.array([[2.0]]),
        C=np.array([[3.0]]),
        D=np.array([[4.0]]),
        dt=1.0,
    )

    # G(z) = 6/(z-0.5) + 4
    assert np.allclose(discrete_transfer_function_matrix(system, 2.0), [[8.0]])
    assert np.allclose(system.impulse_response(4)[:, 0, 0], [4.0, 6.0, 3.0, 1.5])
    assert system.is_stable()
    assert np.allclose(system.poles(), [0.5])


def test_discrete_step_and_forced_response_agree():
    system = DiscreteLTI(
        A=np.array([[0.5]]),
        B=np.array([[2.0]]),
        C=np.array([[3.0]]),
        D=np.array([[4.0]]),
        dt=1.0,
    )

    states, outputs = system.forced_response(np.ones((4, 1)))
    assert np.allclose(states[:, 0], [0.0, 2.0, 3.0, 3.5, 3.75])
    assert np.allclose(outputs[:, 0], [4.0, 10.0, 13.0, 14.5])
    assert np.allclose(system.step_response(4)[:, 0, 0], outputs[:, 0])


def test_discrete_frequency_response_and_validation():
    system = DiscreteLTI(
        A=np.diag([0.5, 0.9]),
        B=np.ones((2, 2)),
        C=np.ones((1, 2)),
        D=np.zeros((1, 2)),
        dt=0.5,
    )

    response = system.frequency_response([0.01, 1.0, 10.0])
    assert response.shape == (3, 1, 2)
    with pytest.raises(ValueError, match="n_steps"):
        system.impulse_response(0)
    with pytest.raises(ValueError, match="inputs"):
        system.forced_response(np.ones((4, 1)))
    with pytest.raises(ValueError, match="dt"):
        DiscreteLTI(system.A, system.B, system.C, system.D, dt=0.0)


def test_discrete_fixed_point_reports_solution_and_stability():
    system = DiscreteLTI(
        A=np.array([[0.5]]),
        B=np.array([[2.0]]),
        C=np.array([[3.0]]),
        D=np.array([[4.0]]),
    )

    result = discrete_fixed_point(system, [1.0])

    assert np.allclose(result.state, [4.0])
    assert result.residual_norm == pytest.approx(0.0)
    assert result.condition_number == pytest.approx(1.0)
    assert result.stable
    assert result.converged


def test_discrete_fixed_point_rejects_nonunique_solution():
    system = DiscreteLTI(
        A=np.array([[1.0]]),
        B=np.array([[1.0]]),
        C=np.array([[1.0]]),
        D=np.array([[0.0]]),
    )

    with pytest.raises(ValueError, match="unique"):
        discrete_fixed_point(system, [1.0])


def test_grouped_discrete_modal_residues_reconstruct_impulse_response():
    system = DiscreteLTI(
        A=np.diag([0.5, 0.5, 0.9]),
        B=np.eye(3),
        C=np.eye(3),
        D=np.zeros((3, 3)),
    )

    modal = discrete_modal_analysis(system)

    assert np.allclose(modal.poles, [0.5, 0.9])
    assert np.array_equal(modal.multiplicities, [2, 1])
    assert np.allclose(modal.state_participation, [[1, 1, 0], [0, 0, 1]])
    assert np.allclose(modal.impulse_response(8), system.impulse_response(8))
    assert modal.decay_times[0] == pytest.approx(-1.0 / np.log(0.5))
    assert modal.reconstruction_error < 1e-12


def test_discrete_modal_analysis_rejects_jordan_block():
    system = DiscreteLTI(
        A=np.array([[0.5, 1.0], [0.0, 0.5]]),
        B=np.ones((2, 1)),
        C=np.ones((1, 2)),
        D=np.zeros((1, 1)),
    )

    with pytest.raises(ValueError, match="not diagonalizable"):
        discrete_modal_analysis(system)
