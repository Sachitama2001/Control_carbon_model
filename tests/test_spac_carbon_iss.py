from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from control_carbon.spac_carbon import YEAR_DAYS, carbon_matrix
from control_carbon.spac_carbon_iss import (
    assimilation_gradient, capacity_diagnostics, carbon_feedback_certificate,
    carbon_iss_certificate, componentwise_iss_bound, equilibrium_tracking_diagnostics,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_equilibria import stock_ratios


@pytest.fixture
def model(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from prepare_manuscript_phase_portraits import model_family
    return RationalSPAC(model_family()["E4"])


def test_linear_iss_constants_and_semigroup_bound(model):
    cert = carbon_iss_certificate(model.base.carbon)
    assert cert.l1_decay_per_day*YEAR_DAYS == pytest.approx(.005)
    assert cert.euclidean_prefactor == pytest.approx(7.1275926032)
    assert cert.euclidean_decay_per_day*YEAR_DAYS == pytest.approx(.01618487409)
    assert cert.lyapunov_residual < 1e-13
    for day in [0., 1., 10., 100., 1000., 10000.]:
        transition = expm(cert.matrix*day)
        assert np.linalg.norm(transition, 1) <= np.exp(-cert.l1_decay_per_day*day)+1e-13
        assert np.linalg.norm(transition, 2) <= (
            cert.euclidean_prefactor*np.exp(-cert.euclidean_decay_per_day*day)+1e-13)


def test_componentwise_sup_bound_is_exact_for_positive_constant_input(model):
    initial, input_rate = np.array([.01, .1, .02, .4]), np.array([.001, .002, .003, .004])
    matrix = carbon_matrix(model.base.carbon)
    solution = solve_ivp(lambda t, e: matrix@e+input_rate, (0, 1000), initial,
                         rtol=1e-11, atol=1e-12)
    bound = componentwise_iss_bound(model.base.carbon, 1000., initial, input_rate)
    np.testing.assert_allclose(bound, solution.y[:, -1], rtol=1e-10, atol=1e-11)
    np.testing.assert_array_equal(
        componentwise_iss_bound(model.base.carbon, 0., -initial, input_rate), initial)


@pytest.mark.parametrize("psi,leaf", [(-1.6, .07), (-3., .03), (-3.5, .02), (-4.2, .01)])
def test_analytic_capacity_rate_and_two_independent_error_identities(model, psi, leaf):
    y = np.r_[-.9, psi, leaf, 1.4, .27, 3.7]
    diag = capacity_diagnostics(model, y, .00244)
    assert diag["carbon_identity_residual"] < 1e-15
    assert diag["closed_error_identity_residual"] < 1e-14
    gradient = assimilation_gradient(model, y)
    direction = np.linalg.solve(-carbon_matrix(model.base.carbon),
                                np.r_[model.base.carbon.allocation, 0.])
    step = 1e-6
    for j in [1, 2]:
        offset = np.eye(6)[j]*step
        fd = (model.assimilation(y+offset)-model.assimilation(y-offset))/(2*step)
        assert gradient[j] == pytest.approx(fd, rel=2e-8, abs=1e-12)
    time_step = min(1e-5, 1e-5/max(1., np.max(np.abs(diag["rhs"]))))
    fd_capacity_rate = (model.capacity(y+time_step*diag["rhs"])
                        -model.capacity(y-time_step*diag["rhs"]))/(2*time_step)
    np.testing.assert_allclose(diag["capacity_rate"], fd_capacity_rate,
                               rtol=2e-8, atol=1e-9)
    np.testing.assert_allclose(diag["capacity_rate"], direction*(gradient@diag["rhs"]))


def test_stomatal_kink_uses_directional_path_rate_not_a_classical_gradient(model):
    state = np.r_[-1., -4., .03, 1., .2, 2.]
    with pytest.raises(ValueError, match="kink"):
        assimilation_gradient(model, state)
    diag = capacity_diagnostics(model, state, .003)
    assert diag["capacity_derivative_kind"] == "directional_at_stomatal_kink"
    np.testing.assert_array_equal(diag["capacity_rate"], np.zeros(4))
    np.testing.assert_array_equal(diag["rhs"][:2], np.zeros(2))
    assert diag["closed_error_identity_residual"] < 1e-15


def test_zero_carbon_boundary_has_finite_one_sided_derivative(model):
    y = np.r_[-1., -1., np.zeros(4)]
    diag = capacity_diagnostics(model, y, .00244)
    assert diag["leaf_feedback_ratio"] > 1
    assert np.all(np.isfinite(assimilation_gradient(model, y)))
    np.testing.assert_array_equal(diag["capacity_rate"], np.zeros(4))
    np.testing.assert_array_equal(diag["error"], np.zeros(4))


def test_carbon_balance_does_not_imply_full_equilibrium_or_invariant_zero_error(model):
    leaf = (model.base.carbon.allocation[0]*model.base.carbon.assimilation_max*.8
            /model.base.carbon.losses[0]-1/model.base.carbon.canopy_coefficient)
    state = np.r_[-.3, -.8, stock_ratios(model.base)*leaf]
    diag = capacity_diagnostics(model, state, .00244)
    np.testing.assert_allclose(diag["error"], 0., atol=1e-14)
    np.testing.assert_allclose(diag["rhs"][2:], 0., atol=1e-15)
    assert np.linalg.norm(diag["rhs"][:2]) > 1e-3
    assert np.linalg.norm(diag["error_rate"]) > 1e-3
    np.testing.assert_allclose(diag["error_rate"], -diag["capacity_rate"], atol=1e-14)


@pytest.mark.parametrize("rho", [0., .5625, .9, .95])
def test_uniform_carbon_feedback_comparison_certificate(model, rho):
    certificate = carbon_feedback_certificate(model, rho)
    weights = certificate["weights"]
    matrix = certificate["comparison_matrix"]
    alpha = certificate["decay_per_day"]
    assert np.all(weights > 0)
    assert alpha > 0
    assert certificate["water_rate_gain"] > 0
    assert certificate["leaf_decay_per_day"]*YEAR_DAYS == pytest.approx(.8*(1-rho))
    assert np.max((matrix@weights)/weights) == pytest.approx(-alpha)
    for leaf in [.01, .05, .2]:
        y = np.r_[-1., -3.5, leaf, 1., .2, 2.]
        diag = capacity_diagnostics(model, y, .00244)
        if diag["leaf_feedback_ratio"] <= rho:
            assert np.all(diag["carbon_feedback_matrix"] <= matrix+1e-15)


@pytest.mark.parametrize("bad", [-1., 1., np.inf, np.nan])
def test_invalid_feedback_certificate_fails(model, bad):
    with pytest.raises(ValueError):
        carbon_feedback_certificate(model, bad)


@pytest.mark.parametrize("vpd", [.00244, .003])
def test_full_equilibrium_coordinate_change_preserves_spectrum_and_tangent(model, vpd):
    y = model.equilibria(vpd)[0]
    d = equilibrium_tracking_diagnostics(model, y, vpd)
    ratio = d["leaf_feedback_ratio"]
    assert ratio == pytest.approx(1/(1+model.base.carbon.canopy_coefficient*y[2]))
    assert d["coordinate_determinant"] > 0
    np.testing.assert_allclose(d["tangent_per_mpa"], d["curve_tangent_per_mpa"],
                               rtol=1e-11, atol=1e-9)
    step = 1e-8
    fd = (model.equilibria(vpd+step)[0]-model.equilibria(vpd-step)[0])/(2*step)
    np.testing.assert_allclose(d["tangent_per_mpa"], fd, rtol=1e-7)
    full = np.sort_complex(np.linalg.eigvals(d["jacobian"]))
    transformed = np.sort_complex(np.linalg.eigvals(d["transformed_jacobian"]))
    np.testing.assert_allclose(full, transformed, rtol=1e-8, atol=1e-11)
    assert np.max(d["carbon_feedback_eigenvalues_per_day"].real) < 0
    assert (np.max(full.real) > 0) == (vpd == .003)


def test_error_coordinates_rejected_at_non_equilibrium(model):
    with pytest.raises(ValueError, match="equilibrium"):
        equilibrium_tracking_diagnostics(model, np.r_[-1., -2., .05, 1., .2, 3.], .00244)


@pytest.mark.parametrize("elapsed,initial,rate", [
    (-1., [0]*4, [0]*4), (0., [0]*3, [0]*4), (0., [0]*4, [-1]*4),
    (np.nan, [0]*4, [0]*4), (1., [np.inf]*4, [0]*4),
])
def test_invalid_sup_bounds_fail_explicitly(model, elapsed, initial, rate):
    with pytest.raises(ValueError):
        componentwise_iss_bound(model.base.carbon, elapsed, initial, rate)
