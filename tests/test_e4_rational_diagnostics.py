import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_rational import RationalSPAC


@pytest.fixture
def scripts(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    return (importlib.import_module("analyze_spac_carbon_rational_e4_diagnostics"),
            importlib.import_module("animate_spac_carbon_rational_e4"))


def test_stability_plot_breaks_disjoint_intervals(scripts):
    diagnostic, _ = scripts
    fig, ax = plt.subplots()
    vpd = np.linspace(.001, .005, 5)
    leaves = np.array([.2, .1, .05, .01, .001])
    diagnostic.plot_equilibrium_segments(ax, vpd, leaves, np.array([-1, -1, 1, 1, -1]))
    stable, unstable = ax.lines
    np.testing.assert_allclose(stable.get_xdata(), vpd*1000)
    np.testing.assert_allclose(stable.get_ydata(), [.2, .1, np.nan, np.nan, .001],
                               equal_nan=True)
    np.testing.assert_allclose(unstable.get_ydata(), [np.nan, np.nan, .05, .01, np.nan],
                               equal_nan=True)
    plt.close(fig)


@pytest.mark.parametrize("ramp", [10., 10*YEAR_DAYS])
def test_pulses_are_continuous_with_zero_endpoint_slopes(scripts, ramp):
    diagnostic, _ = scripts
    pulse = diagnostic.pulse_vpd
    delta = diagnostic.VPD_PEAK-diagnostic.VPD_START
    hold = diagnostic.VPD_HOLD_DAYS
    for day, expected in [(0., diagnostic.VPD_START), (ramp, diagnostic.VPD_PEAK),
                          (ramp+hold, diagnostic.VPD_PEAK),
                          (2*ramp+hold, diagnostic.VPD_START)]:
        assert pulse(day, ramp) == pytest.approx(expected, abs=1e-15)
        epsilon = ramp*1e-7
        for side in [-1, 1]:
            slope = abs(pulse(day+side*epsilon, ramp)-pulse(day, ramp))/epsilon
            assert slope < 1e-5*delta/ramp
    assert pulse(ramp/2, ramp) == pytest.approx(
        (diagnostic.VPD_START+diagnostic.VPD_PEAK)/2)
    values = np.array([pulse(t, ramp) for t in np.linspace(0., 2*ramp+hold, 1001)])
    assert np.all(values >= diagnostic.VPD_START)
    assert np.all(values <= diagnostic.VPD_PEAK)


@pytest.mark.parametrize("day,ramp", [(np.nan, 10.), (0., 0.), (1., -1.)])
def test_invalid_forcing_arguments_raise(scripts, day, ramp):
    with pytest.raises(ValueError):
        scripts[0].pulse_vpd(day, ramp)


def test_movie_clock_resolves_both_fast_ramps_and_uses_common_horizon(scripts):
    diagnostic, movie = scripts
    days = movie.frame_days()
    assert days[0] == 0
    assert days[-1] == diagnostic.EARLY_YEARS*YEAR_DAYS
    assert np.all(np.diff(days) > 0)
    for start in [0., 10.+diagnostic.VPD_HOLD_DAYS]:
        ramp_days = days[(days >= start) & (days <= start+10.)]
        assert len(ramp_days) >= 40
        assert np.max(np.diff(ramp_days)) < .25
    assert len(days[(days >= 10.) & (days <= 15.)]) >= 60


def test_movie_interpolation_keeps_positive_carbon_and_rejects_extrapolation(scripts):
    interpolate = scripts[1].interpolate_state
    days = np.array([0., 2.])
    states = np.array([[-1., -2., 1., 2., 3., 4.], [-2., -4., 4., 8., 12., 16.]])
    np.testing.assert_allclose(interpolate(days, states, 1.), [-1.5, -3., 2., 4., 6., 8.])
    for day in [-1., 3., np.nan]:
        with pytest.raises(ValueError):
            interpolate(days, states, day)


def test_movie_enumerates_all_roots_and_colors_unstable_equilibria(scripts):
    diagnostic, movie = scripts
    models = diagnostic.model_family()
    e4 = RationalSPAC(models["E4"])
    records = movie.frozen_equilibria(e4, np.array([0., 10.]), 10.)
    assert [len(row["positive_equilibria"]) for row in records] == [1, 1]
    assert [row["positive_equilibria"][0]["stability"] for row in records] == [
        "stable", "unstable"]
    e0 = RationalSPAC(models["E0"])
    multiple = movie.frozen_equilibria(e0, np.array([0.]), 10.)
    assert len(multiple[0]["positive_equilibria"]) == 2


def test_fast_closure_event_and_conditional_decay_bound(scripts):
    diagnostic, _ = scripts
    model = RationalSPAC(diagnostic.model_family()["E4"])
    data = diagnostic.integrate_pulse(
        model, model.equilibria(diagnostic.VPD_START)[0], 10/YEAR_DAYS,
        np.asarray(reference_scales(model.base).state))
    assert data["status"] == "leaf_domain_exit"
    assert data["xylem_closure_day"] == pytest.approx(11.094697707, abs=2e-5)
    assert data["state_at_xylem_closure"][1] == pytest.approx(-3., abs=1e-12)
    later = data["day"] >= data["xylem_closure_day"]
    assert np.all(data["state"][later, 1] <= -3.+1e-9)
    bound = data["state_at_xylem_closure"][2]*np.exp(
        -.35*(data["day"][later]-data["xylem_closure_day"])/YEAR_DAYS)
    assert np.all(data["state"][later, 2] <= bound*(1+1e-8))
    assert data["exit_day"]/YEAR_DAYS == pytest.approx(13.94, abs=.01)


@pytest.mark.parametrize("psi", [-3.5, -2.5, -1.])
def test_zero_equilibrium_continuum_and_transverse_leaf_eigenvalue(scripts, psi):
    model = RationalSPAC(scripts[0].model_family()["E4"])
    state = np.r_[psi, psi, np.zeros(4)]
    np.testing.assert_allclose(model.rhs(state, .00244), 0., atol=1e-15)
    jac = model.jacobian(state, .00244)
    c = model.base.carbon
    leaf_eigenvalue = c.allocation[0]*c.assimilation_max*c.canopy_coefficient*max(
        0., 1+psi/4)-c.losses[0]
    assert jac[2, 2] == pytest.approx(leaf_eigenvalue)
    eigenvalues = np.linalg.eigvals(jac)
    assert np.min(abs(eigenvalues)) < 1e-12
    assert np.min(abs(eigenvalues-leaf_eigenvalue)) < 1e-12
    if psi > -20/9:
        assert leaf_eigenvalue > 0
