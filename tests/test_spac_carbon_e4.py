from dataclasses import replace

import numpy as np
import pytest

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, positive_carbon_equilibrium, hydraulic_factor, carbon_budget,
)
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_e4 import core_jacobian, core_tensors, hopf_point, periodic_orbit

H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
REF = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
F = StructuralFeedback(REF[0], REF[1], .01, 1e-6, 1., 1., 1.)
M = PressureSPAC(H, C, F, 0., .101325)


@pytest.mark.parametrize("exponent", [.5, 1., 2.])
def test_four_state_factorization_and_exact_derivatives(exponent):
    m = replace(M, feedback=replace(F, wood_exponent=exponent))
    y = np.array([-.8, -1.4, .18, 3.])
    full = m.jacobian(np.r_[y, .5, 4.], .002)
    a = core_jacobian(y, .002, m)
    expected = np.r_[np.linalg.eigvals(a), -np.array(C.losses[2:])]
    np.testing.assert_allclose(np.sort_complex(np.linalg.eigvals(full)), np.sort_complex(expected), atol=1e-12)
    b, c = core_tensors(y, .002, m)
    step = 1e-5
    for j, e in enumerate(np.eye(4)):
        fd2 = (core_jacobian(y+step*e, .002, m)-core_jacobian(y-step*e, .002, m))/(2*step)
        fd3 = (core_tensors(y+step*e, .002, m)[0]-core_tensors(y-step*e, .002, m)[0])/(2*step)
        np.testing.assert_allclose(b[:, :, j], fd2, rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(c[:, :, :, j], fd3, rtol=1e-8, atol=1e-8)


def test_e4_hopf_resolves_simple_pair_and_nonlinear_type():
    h = hopf_point((-2.05, -1.95), M)
    assert h["vpd_mpa"] == pytest.approx(.0026647529, abs=1e-10)
    assert M.air_diagnostic(h["vpd_mpa"], 25.)["admissible_at_diagnostic_temperature"]
    assert h["omega"] > .02
    assert np.count_nonzero(abs(h["eigenvalues"].real) < 1e-9) == 2
    assert abs(h["l1"]) > .01
    assert h["slope_refinement_per_mpa"][-1] > 400.


def test_invalid_smooth_domain():
    with pytest.raises(ValueError):
        core_tensors([-1., -3.1, .1, 2.], .002, M)
    with pytest.raises(ValueError):
        core_tensors([-1., -2., .1, 0.], .002, M)
    with pytest.raises(ValueError, match="complex-pair neighborhood"):
        hopf_point((-2.12, -1.9), M)


def test_wood_exponent_change_must_not_be_mislabelled_attracting_hopf():
    m = replace(M, feedback=replace(F, wood_exponent=1.25))
    with pytest.raises(ValueError, match="stable remaining eigenvalue"):
        hopf_point((-2.01, -1.99), m)
    folds = m.stationary_points((-2.222, -.001))
    assert len(folds) == 1
    state, vpd, _ = m.equilibrium_curve(-2.)
    assert max(np.linalg.eigvals(core_jacobian(state[:4], vpd, m)).real) > 1e-5


def test_e4_periodic_collocation_full_lift_and_budgets():
    h = hopf_point((-2.05, -1.95), M)
    cycle = periodic_orbit(h, .001, M, tol=1e-9)
    assert h["l1"] > 0
    assert cycle["vpd_mpa"] < h["vpd_mpa"]
    assert max(abs(cycle["core_floquet"])) > 1.0004
    assert min(abs(cycle["core_floquet"]-1)) < 1e-7
    assert cycle["independent_flow_max_state_error"] < 1e-7
    assert cycle["state"].shape == (401, 6)
    assert cycle["full_state_closure"] < 1e-8
    for y in cycle["state"][::40]:
        rhs = M.rhs(y, cycle["vpd_mpa"])
        erx, exl, e = M.fluxes(y, cycle["vpd_mpa"])
        assert H.cx*rhs[0]+H.cl*rhs[1] == pytest.approx(erx-e, abs=1e-15)
        assert sum(rhs[2:]) == pytest.approx(carbon_budget(y[2:], hydraulic_factor(y[1], H), C)["nep"], abs=1e-15)
    with pytest.raises(ValueError):
        periodic_orbit(h, 0., M)
