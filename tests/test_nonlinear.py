import numpy as np
import pytest

from control_carbon.nonlinear import (
    ContinuousNonlinear,
    compartment_rhs,
    continuous_equilibrium,
    simulate_continuous,
)


def linear_system():
    return ContinuousNonlinear(
        rhs=lambda state, driver: -2.0 * state + driver,
        jacobian=lambda state, driver: np.array([[-2.0]]),
        state_names=("carbon",), state_units=("Mg C ha^-1",),
        driver_name="input", driver_unit="Mg C ha^-1 day^-1",
        time_unit="day", approximation="model-agnostic derived analysis",
    )


def test_nonlinear_integrator_matches_analytic_forced_solution():
    times = np.linspace(0, 4, 41)
    result = simulate_continuous(linear_system(), [1.0], times, lambda time: time)
    expected = 1.25 * np.exp(-2 * times) + times / 2 - 0.25
    np.testing.assert_allclose(result.states[:, 0], expected, rtol=2e-9, atol=2e-10)
    np.testing.assert_array_equal(result.drivers, times)
    assert result.system.time_unit == "day"


def test_frozen_equilibrium_diagnostics():
    result = continuous_equilibrium(linear_system(), [0.0], 6.0)
    np.testing.assert_allclose(result.state, [3.0])
    np.testing.assert_allclose(result.eigenvalues, [-2.0])
    assert result.converged and result.stability == "stable"
    assert result.residual_norm < 1e-12


def test_compartment_sign_balance_and_boundary_positivity():
    transfer = [[-1.0, 0.0], [0.4, -1.0]]
    result = compartment_rhs([10, 20], transfer, [0.1, 0.2], [2, 1], [[1], [0]], [3])
    np.testing.assert_allclose(result, [1.0, -3.2])
    assert result.sum() == pytest.approx(3 - 0.6 * 2 - 4)
    boundary = compartment_rhs([0, 20], transfer, [0.1, 0.2], [2, 1], [[1], [0]], [3])
    assert boundary[0] >= 0
    with pytest.raises(ValueError, match="signed compartment"):
        compartment_rhs([1, 1], [[-1, 0], [1.2, -1]], [1, 1], [1, 1], [[1], [0]], [1])


@pytest.mark.parametrize("times", [[0], [1, 0], [0, 0], [0, np.nan]])
def test_invalid_times_rejected(times):
    with pytest.raises(ValueError, match="times"):
        simulate_continuous(linear_system(), [0], times, lambda time: 1.0)


def test_invalid_callback_and_solver_inputs_rejected():
    with pytest.raises(ValueError, match="driver"):
        simulate_continuous(linear_system(), [0], [0, 1], lambda time: np.nan)
    with pytest.raises(ValueError, match="max_step"):
        simulate_continuous(linear_system(), [0], [0, 1], lambda time: 0, max_step=0)
    with pytest.raises(ValueError, match="initial_state"):
        simulate_continuous(linear_system(), [0, 1], [0, 1], lambda time: 0)