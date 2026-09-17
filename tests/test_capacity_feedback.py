import numpy as np
import pytest

from control_carbon.capacity_feedback import (
    CubicCapacityFeedback,
    critical_ramp_duration,
    periodic_attractor,
    periodic_capacity,
    simulate_capacity_ramp,
)


def test_self_consistent_capacities_and_full_stability():
    model = CubicCapacityFeedback()
    for driver in np.linspace(0, 1, 11):
        low, threshold, high = model.equilibria(driver)
        for stock in (low, threshold, high):
            assert model.capacity(stock, driver) == pytest.approx(stock)
            assert model.rhs(stock, driver) == pytest.approx(0, abs=1e-14)
        assert model.jacobian(low, driver) < 0
        assert model.jacobian(threshold, driver) > 0
        assert model.jacobian(high, driver) < 0
        # The frozen matrix M=-relaxation is stable at every state; the
        # capacity feedback, not M alone, reverses the middle branch.
        assert -model.relaxation < 0


def test_capacity_is_positive_on_documented_domain():
    model = CubicCapacityFeedback()
    values = [model.capacity(x, driver)
              for x in np.linspace(0, 12, 301)
              for driver in np.linspace(0, 1, 101)]
    assert min(values) > 0


def test_storage_potential_is_not_always_contracting():
    model = CubicCapacityFeedback()
    # Near the unstable self-consistent capacity its squared magnitude grows.
    assert model.potential_energy_derivative(6.9, 1.0) > 0
    # Near a stable branch it declines.
    assert model.potential_energy_derivative(9.0, 1.0) < 0


def test_fast_and_slow_ramps_reach_different_nonzero_capacities():
    fast = simulate_capacity_ramp(20)
    slow = simulate_capacity_ramp(120)
    assert fast["outcome"] == "low"
    assert slow["outcome"] == "high"
    assert fast["stock"][-1] == pytest.approx(1.0, abs=1e-7)
    assert slow["stock"][-1] == pytest.approx(10.0, abs=1e-7)
    assert np.min(fast["stock"]) > 0


def test_critical_duration_and_numerical_refinement():
    bracket = critical_ramp_duration(iterations=10)
    refined = critical_ramp_duration(iterations=13)
    assert bracket[0] <= refined[0] < refined[1] <= bracket[1]
    coarse = simulate_capacity_ramp(120, max_step=1.0, rtol=1e-8)
    fine = simulate_capacity_ramp(120, max_step=.05, rtol=1e-11)
    assert coarse["outcome"] == fine["outcome"] == "high"
    np.testing.assert_allclose(coarse["stock"][-1], fine["stock"][-1], atol=1e-8)


def test_periodic_attractor_solves_forced_linear_equation():
    t = np.linspace(0, 2, 2001)
    relaxation = 1.3
    state = periodic_attractor(t, relaxation=relaxation)
    derivative = np.gradient(state, t)
    rhs = relaxation*(periodic_capacity(t)-state)
    np.testing.assert_allclose(derivative[2:-2], rhs[2:-2], atol=2e-4)
    assert np.max(np.abs(state-periodic_capacity(t))) > 1


@pytest.mark.parametrize("kwargs", [
    {"rate": 0}, {"relaxation": -1}, {"threshold_shift": 20},
])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        CubicCapacityFeedback(**kwargs)
