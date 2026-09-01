import numpy as np

from control_carbon.state_space import ContinuousLTI, impulse_response, transfer_function_matrix


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
