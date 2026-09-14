import numpy as np
import pytest

from control_carbon.nonlinear import continuous_equilibrium
from control_carbon.tipping import SaturatingFeedback, SmoothRamp, run_feedback_ramp


def test_feedback_analytic_roots_stability_and_jacobian():
    model = SaturatingFeedback()
    system = model.system("gain")
    roots = model.equilibria("gain", 3.0)
    np.testing.assert_allclose(roots, [0, (3 - np.sqrt(5)) / 2, (3 + np.sqrt(5)) / 2])
    for stock, stability in zip(roots, ["stable", "unstable", "stable"]):
        result = continuous_equilibrium(system, [stock], 3.0)
        assert result.converged and result.stability == stability
        assert result.residual_norm < 1e-12
        expected = -1.0 if stock == 0 else (1 - stock**2) / (1 + stock**2)
        assert result.eigenvalues[0] == pytest.approx(expected)
        step = 1e-5
        numerical = (system.evaluate([stock + step], 3) - system.evaluate([stock - step], 3)) / (2 * step)
        np.testing.assert_allclose(result.jacobian[:, 0], numerical, atol=1e-9)


def test_fold_and_affine_negative_control():
    model = SaturatingFeedback()
    fold = continuous_equilibrium(model.system("gain"), [1], 2)
    assert fold.converged and fold.stability == "nonhyperbolic"
    assert model.equilibria("gain", 1.99).size == 1
    assert model.equilibria("gain", 2.01).size == 3
    for gain in [0.5, 2, 3, 5]:
        result = continuous_equilibrium(model.system("gain", feedback=False), [0], gain)
        assert result.stability == "stable"
        np.testing.assert_allclose(result.state, model.equilibria("gain", gain, feedback=False))


def test_feedback_has_nonnegative_influx_and_exact_loss_balance():
    model = SaturatingFeedback(mortality=0.5, gain=4, half_saturation=2)
    system = model.system("gain")
    for stock in [0, 0.1, 2, 20]:
        influx = 4 * stock**2 / (4 + stock**2)
        assert influx >= 0
        assert system.evaluate([stock], 4)[0] == pytest.approx(influx - 0.5 * stock)
    assert system.evaluate([0], 4)[0] == 0


def test_smooth_ramp_holds_and_rate_definition():
    ramp = SmoothRamp(1, 10, 2, start_time=3)
    assert ramp(2) == 1
    assert ramp(3.25) == 5.5
    assert ramp(4) == 10
    assert ramp.end_time == 3.5
    assert ramp.max_driver_speed == 27


def test_bifurcation_loss_after_crossing_analytic_fold():
    result = run_feedback_ramp(SaturatingFeedback(), "gain", SmoothRamp(3, 1.5, 0.02))
    assert result.final_equilibria.size == 1
    assert result.final_basin == "low"
    assert result.distance_to_attractor < 1e-8
    assert result.final_residual < 1e-8


@pytest.mark.parametrize("rtol,atol,steps", [(1e-8, 1e-10, 50), (1e-10, 1e-12, 100)])
def test_rate_changes_final_basin_without_frozen_branch_loss(rtol, atol, steps):
    model = SaturatingFeedback()
    for half_saturation in np.linspace(1, 10, 31):
        roots = model.equilibria("half_saturation", half_saturation)
        assert roots.size == 3
        assert model.system("half_saturation").linearize([roots[-1]], half_saturation)[0, 0] < -0.7
    final_low_threshold = model.equilibria("half_saturation", 10)[1]
    assert model.equilibria("half_saturation", 1)[-1] < final_low_threshold
    for rate, basin in [(0.01, "high"), (10.0, "low")]:
        result = run_feedback_ramp(
            model, "half_saturation", SmoothRamp(1, 10, rate),
            rtol=rtol, atol=atol, steps_per_ramp=steps,
        )
        assert result.final_basin == basin
        assert result.distance_to_attractor < 1e-7
        assert result.final_residual < 1e-7
        assert result.hold.max_step == pytest.approx(25 / (steps * model.mortality))
        assert np.min(result.ramp.states) >= -atol
        assert np.min(result.hold.states) >= -atol


def test_longer_hold_preserves_rate_outcome():
    result = run_feedback_ramp(
        SaturatingFeedback(), "half_saturation", SmoothRamp(1, 10, 10), hold_time=80,
    )
    assert result.final_basin == "low"
    assert result.distance_to_attractor < 1e-10


@pytest.mark.parametrize("kwargs", [{"mortality": 0}, {"half_saturation": -1}, {"gain": np.nan}])
def test_invalid_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        SaturatingFeedback(**kwargs)


def test_invalid_protocol_and_boundary_classification():
    model = SaturatingFeedback()
    with pytest.raises(ValueError, match="driver"):
        model.system("rainfall")
    with pytest.raises(ValueError, match="domain"):
        model.equilibria("half_saturation", 0)
    with pytest.raises(ValueError, match="rate"):
        SmoothRamp(1, 2, 0)
    with pytest.raises(ValueError, match="initial forcing"):
        run_feedback_ramp(model, "gain", SmoothRamp(1, 3, 1))
    result = run_feedback_ramp(model, "gain", SmoothRamp(3, 2, 1))
    assert result.final_basin == "unresolved"
    assert np.isnan(result.distance_to_attractor)