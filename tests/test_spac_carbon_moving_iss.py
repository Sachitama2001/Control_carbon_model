from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_iss import equilibrium_tracking_diagnostics
from control_carbon.spac_carbon_moving_iss import (
    certify_moving_equilibrium, equilibrium_tangent_enclosure,
    moving_initial_energy_upper, moving_ramp_bounds, ramp_speed_upper,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_stability import _Interval
from control_carbon.spac_carbon_tracking import SmoothVPDRamp


@pytest.fixture
def model(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from prepare_manuscript_phase_portraits import model_family
    return RationalSPAC(model_family()["E4"])


@pytest.fixture
def certificate(model):
    return certify_moving_equilibrium(
        model, (.0024399, .0024401), reference_scales(model.base).state)


def test_tangent_interval_matches_implicit_full_jacobian_solve(model, certificate):
    c = certificate
    tangent_box = c.evidence["equilibrium_tangent_intervals_per_mpa"]
    for v in np.linspace(*c.vpd_band_mpa, 7):
        y = model.equilibria(v)[0]
        tangent = equilibrium_tracking_diagnostics(model, y, v)["tangent_per_mpa"]
        for value, bounds in zip(tangent, tangent_box):
            assert Decimal(bounds[0]) < Decimal.from_float(float(value)) < Decimal(bounds[1])
    _, slope = equilibrium_tangent_enclosure(model, _Interval(*c.evidence["leaf_band"]))
    assert slope.hi < 0


def test_common_metric_dissipation_across_band_and_tube(model, certificate):
    c = certificate
    lower, upper = np.array(c.evidence["state_box"], dtype=float).T
    for state, v in zip(np.random.default_rng(1020).uniform(lower, upper, (20, 6)),
                        np.linspace(*c.vpd_band_mpa, 20)):
        a = model.jacobian(state, v)*c.scales[None, :]/c.scales[:, None]
        assert np.linalg.eigvalsh(-a.T@c.metric-c.metric@a).min() > c.dissipation_lower_per_day
    assert c.max_vpd_rate_mpa_per_day > 3e-15


def test_nonlinear_moving_target_obeys_bound_and_capacity_output(model, certificate):
    c = certificate
    ramp = SmoothVPDRamp(.00244, .00244+5e-11, 100*YEAR_DAYS)
    y0 = model.equilibria(ramp.start_mpa)[0]
    initial = moving_initial_energy_upper(model, c, y0, ramp.start_mpa)
    times = np.linspace(0, 200*YEAR_DAYS, 121)
    bound = moving_ramp_bounds(c, ramp, times, initial)
    run = solve_ivp(lambda t, y: model.rhs(y, ramp.value_rate(t)[0]),
                    (0, times[-1]), y0, method="Radau", t_eval=times,
                    rtol=2e-12, atol=1e-14, max_step=YEAR_DAYS)
    assert run.success
    target = np.array([model.equilibria(ramp.value_rate(t)[0])[0] for t in times])
    error = (run.y.T-target)/c.scales
    metric_radius = np.sqrt(np.einsum("ni,ij,nj->n", error, c.metric, error))
    capacity_error = np.array([
        (y[2:]-model.capacity(y))/c.scales[2:] for y in run.y.T])
    assert np.all(metric_radius <= bound["metric_radius"]+1e-11)
    assert np.all(np.linalg.norm(capacity_error, axis=1)
                  <= bound["capacity_error_norm_scaled"]+1e-11)
    assert bound["metric_radius"][-1] < bound["metric_radius"][60]
    assert np.max(bound["metric_radius"]) < np.sqrt(c.energy_limit)


def test_rate_upper_handles_both_directions_and_zero_change():
    up, down = SmoothVPDRamp(.002, .003, 10), SmoothVPDRamp(.003, .002, 10)
    assert ramp_speed_upper(up) >= abs(up.value_rate(5)[1])
    assert ramp_speed_upper(up) == ramp_speed_upper(down)
    assert ramp_speed_upper(SmoothVPDRamp(.002, .002, 10)) < 1e-300


def test_fast_small_motion_is_certified_by_total_variation(certificate):
    ramp = SmoothVPDRamp(.00244, .00244+5e-11, 10.)
    assert ramp_speed_upper(ramp) > certificate.max_vpd_rate_mpa_per_day
    bound = moving_ramp_bounds(certificate, ramp, [0., 10., 1000.], 0.)
    assert bound["metric_radius"][0] == 0.
    assert bound["metric_radius"][2] < bound["metric_radius"][1]
    assert max(bound["metric_radius"]) < np.sqrt(certificate.energy_limit)/2


@pytest.mark.parametrize("reason", ["fast", "outside", "negative_time", "large_energy"])
def test_tracking_bound_refuses_uncertified_schedules(model, certificate, reason):
    ramp = SmoothVPDRamp(.00244, .00244+5e-11, 100*YEAR_DAYS)
    energy, times = 0., [0., 1.]
    if reason == "fast":
        ramp = SmoothVPDRamp(*certificate.vpd_band_mpa, 10)
    elif reason == "outside":
        ramp = SmoothVPDRamp(.002, .0021, 1e15)
    elif reason == "negative_time":
        times = [-1.]
    else:
        energy = 2*certificate.energy_limit
    with pytest.raises(ValueError):
        moving_ramp_bounds(certificate, ramp, times, energy)


@pytest.mark.parametrize("band", [(0., .001), (.003, .002), (.002439, .002441),
                                  (np.nan, .003)])
def test_invalid_or_unresolved_band_is_not_guaranteed(model, band):
    with pytest.raises(ValueError):
        certify_moving_equilibrium(model, band, reference_scales(model.base).state)
