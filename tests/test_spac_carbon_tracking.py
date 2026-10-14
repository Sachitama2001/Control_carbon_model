from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_tracking import (
    FrozenBranch, SmoothVPDRamp, forcing_derivative, tracking_diagnostics, linearized_tracking,
    frozen_lyapunov_metric,
)
from control_carbon.spac_carbon_e4 import hopf_point, periodic_orbit
from control_carbon.spac_carbon_section import local_cycle_section, section_return

H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
REF = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
F = StructuralFeedback(REF[0], REF[1], .01, 1e-6, 1., 1., 1.)
M = PressureSPAC(H, C, F, 0., .101325)
S = reference_scales(M)


@pytest.mark.parametrize("strengths,bounds,v", [
    ((0, 0, 0), (-1.8, -.1), .0013),
    ((1, 0, 0), (-1.1, -.1), .0013),
    ((0, 1, 0), (-2.05, -.3), .0013),
    ((1, 1, 0), (-2.05, -.3), .0013),
    ((1, 1, 1), (-2.05, -.3), .0024),
])
def test_forcing_and_full_equilibrium_tangent_against_independent_differences(strengths, bounds, v):
    m = replace(M, feedback=replace(F, leaf_supply=strengths[0], leaf_demand=strengths[1],
                                  wood_supply=strengths[2]))
    branch = FrozenBranch(m, bounds)
    y, tangent, _ = branch.tangent(v, S)
    fv = forcing_derivative(y, v, m)
    step = 1e-8
    np.testing.assert_allclose(fv, (m.rhs(y, v+step)-m.rhs(y, v-step))/(2*step), rtol=1e-8, atol=1e-8)
    fd = (branch.state(v+step)-branch.state(v-step))/(2*step)
    np.testing.assert_allclose(tangent, fd, rtol=2e-7, atol=1e-5)
    np.testing.assert_allclose(m.jacobian(y, v)@tangent+fv, 0., atol=2e-10)
    np.testing.assert_allclose(m.rhs(y, v), 0., atol=1e-10)


def test_exact_error_equation_and_quadratic_remainder():
    branch = FrozenBranch(M, (-2.05, -.3))
    y = branch.state(.0024)
    direction = np.array([.2, -.3, .1, .05, -.04, .08])*np.array(S.state)
    rates = []
    for epsilon in [1e-3, 5e-4]:
        d = tracking_diagnostics(y+epsilon*direction, .0024, 1e-6, branch, S)
        rates.append(d["nonlinear_remainder_norm_scaled_per_day"])
        step = 1e-5
        fd = (M.rhs(y+epsilon*direction, .0024)/np.array(S.state)
              -(branch.state(.0024+step*1e-6)-branch.state(.0024-step*1e-6))/(2*step)/np.array(S.state))
        np.testing.assert_allclose(d["exact_error_rate_scaled"], fd, atol=1e-8)
    assert rates[0]/rates[1] == pytest.approx(4., rel=.01)


def test_linearized_constant_branch_matches_matrix_exponential_and_zero_error():
    branch = FrozenBranch(M, (-2.05, -.3))
    times = np.linspace(0, 100, 11)
    ramp = SmoothVPDRamp(.0024, .0024, 10.)
    e0 = np.array([1e-4, -1e-4, 2e-4, -1e-4, 1e-4, 2e-4])
    y = branch.state(.0024)
    scale = np.array(S.state)
    a = M.jacobian(y, .0024)*scale[None, :]/scale[:, None]
    result = linearized_tracking(times, ramp, branch, S, initial_error_scaled=e0)
    np.testing.assert_allclose(result["error_scaled"], [expm(a*t)@e0 for t in times], atol=1e-10)
    np.testing.assert_allclose(linearized_tracking(times, ramp, branch, S)["error_scaled"], 0., atol=1e-14)


def test_smooth_ramp_units_and_endpoint_derivatives():
    ramp = SmoothVPDRamp(.001, .002, 10.)
    assert ramp.value_rate(-1) == (.001, 0.)
    assert ramp.value_rate(11) == (.002, 0.)
    v, rate = ramp.value_rate(5)
    assert v == pytest.approx(.0015)
    assert rate == pytest.approx(.00015)
    step = 1e-4
    assert rate == pytest.approx((ramp.value_rate(5+step)[0]-ramp.value_rate(5-step)[0])/(2*step))
    with pytest.raises(ValueError):
        SmoothVPDRamp(-1., .002, 10.)
    with pytest.raises(ValueError):
        SmoothVPDRamp(.001, .002, 0.)
    with pytest.raises(ValueError):
        ramp.value_rate(np.nan)


def test_fold_and_out_of_branch_inputs_fail_without_selecting_another_root():
    e0 = replace(M, feedback=replace(F, leaf_supply=0., leaf_demand=0., wood_supply=0.))
    with pytest.raises(ValueError, match="fold"):
        FrozenBranch(e0, (-2.222, -.001))
    branch = FrozenBranch(M, (-2.05, -.3))
    with pytest.raises(ValueError, match="outside"):
        branch.state(.003)
    with pytest.raises(ValueError):
        tracking_diagnostics(branch.state(.0024), .0024, np.nan, branch, S)
    with pytest.raises(ValueError):
        linearized_tracking([1, 2], SmoothVPDRamp(.0024, .0024, 1.), branch, S)


def test_moving_lyapunov_metric_derivative_and_unstable_guard():
    branch = FrozenBranch(M, (-2.1, -.3))
    v = .0024
    d = frozen_lyapunov_metric(v, branch, S, step_mpa=5e-8)
    step = 1e-7
    fd = (frozen_lyapunov_metric(v+step, branch, S)["metric"]
          -frozen_lyapunov_metric(v-step, branch, S)["metric"])/(2*step)
    np.testing.assert_allclose(d["metric_derivative_per_mpa"], fd, rtol=1e-5, atol=1.)
    assert d["lyapunov_residual"] < 1e-6
    assert d["minimum_metric_eigenvalue_days"] > 0
    with pytest.raises(ValueError, match="stable"):
        frozen_lyapunov_metric(.00268, branch, S)


def test_small_ramp_nonlinear_linear_tracking_difference_is_quadratic():
    branch = FrozenBranch(M, (-2.05, -.3))
    times = np.r_[np.linspace(0., 10., 21), np.linspace(11., 100., 20)]
    initial = branch.state(.0024)
    differences = []
    for delta in [1e-6, 5e-7]:
        ramp = SmoothVPDRamp(.0024, .0024+delta, 10.)
        nonlinear = M.simulate(initial, times, lambda t: ramp.value_rate(t)[0], rtol=1e-11)
        linear = linearized_tracking(times, ramp, branch, S, rtol=1e-11)
        exact = np.array([(y-branch.state(ramp.value_rate(t)[0]))/np.array(S.state)
                          for t, y in zip(nonlinear["time_days"], nonlinear["state"])])
        differences.append(np.max(np.linalg.norm(exact-linear["error_scaled"], axis=1)))
    assert differences[0]/differences[1] == pytest.approx(4., rel=.05)


def test_outcome_labels_require_full_errors_and_a_resolved_hold_window(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from analyze_spac_carbon_tracking import sampled_outcome
    solution = {"status": "completed", "time_days": np.array([0., 300.])}
    rows = [{"day": day, "equilibrium_error_norm_scaled": 1e-5,
             "full_residual_norm_scaled_per_day": 1e-8} for day in np.linspace(270., 300., 31)]
    assert sampled_outcome(solution, rows, 300.) == "sampled_equilibrium_neighborhood_not_global_basin_proof"
    rows[0]["equilibrium_error_norm_scaled"] = .002
    assert sampled_outcome(solution, rows, 300.) == "finite_time_unresolved"
    with pytest.raises(RuntimeError, match="resolved"):
        sampled_outcome(solution, rows[-1:], 300.)
    solution["status"] = "vegetation_domain_exit"
    assert sampled_outcome(solution, [], 300.) == "vegetation_domain_exit_not_alternate_attractor"


@pytest.fixture(scope="module")
def cycle():
    return periodic_orbit(hopf_point((-2.01, -1.99), M), .002, M, tol=1e-9)


def test_poincare_phase_removal_and_nonlinear_return_derivative(cycle):
    section = local_cycle_section(cycle, M, S)
    basis, normal = section["section_basis_scaled"], section["normal_scaled"]
    np.testing.assert_allclose(normal@basis, 0., atol=1e-15)
    np.testing.assert_allclose(basis.T@basis, np.eye(3), atol=1e-15)
    expected = np.delete(cycle["core_floquet"], np.argmin(abs(cycle["core_floquet"]-1)))
    np.testing.assert_allclose(np.sort(section["multipliers"]), np.sort(expected.real), atol=2e-7)
    np.testing.assert_allclose(section["left_covectors"]@section["right_vectors"], np.eye(3), atol=1e-12)
    step = 2e-7
    direction = section["right_vectors"][:, 0]
    left = section_return(section, -step*direction, M, rtol=1e-12)
    right = section_return(section, step*direction, M, rtol=1e-12)
    fd = (right["coordinates"]-left["coordinates"])/(2*step)
    np.testing.assert_allclose(fd, section["return_derivative"]@direction, atol=1e-5)
    assert left["mode_coordinates"][0] < 0 < right["mode_coordinates"][0]
    assert right["section_residual"] < 1e-12
    np.testing.assert_allclose(section["stable_tangent_covector_scaled"]@basis@section["right_vectors"][:, 1:],
                               0., atol=1e-12)
    with pytest.raises(ValueError):
        section_return(section, [0., 0.], M)
    with pytest.raises(ValueError):
        local_cycle_section(cycle, M, S, rtol=0.)
