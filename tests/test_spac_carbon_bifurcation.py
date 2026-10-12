from dataclasses import replace
from decimal import Decimal, localcontext
from itertools import permutations

import numpy as np
import pytest

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium, carbon_budget,
)
from control_carbon.spac_carbon_two_way import (
    StructuralFeedback, two_way_rhs, two_way_jacobian, equilibrium_curve, structural_fluxes,
)
from control_carbon.spac_carbon_bifurcation import (
    core_rhs, core_jacobian, core_tensors, cubic_hurwitz, first_lyapunov, hopf_point, periodic_orbit,
)
from control_carbon.spac_carbon_stability import _Interval, _branch_intervals, certify_stable_path


H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
STOCKS = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
F = StructuralFeedback(STOCKS[0], STOCKS[1], .01, 1e-6, leaf_demand=1)


@pytest.mark.parametrize("supply", [0., .5, 1.])
def test_exact_core_full_spectra_and_higher_derivative_tensors(supply):
    f = replace(F, leaf_supply=supply)
    y = np.array([-.8, -1.4, .18])
    jac = core_jacobian(y, 0., .02, H, C, f)
    full = two_way_jacobian(np.r_[y, STOCKS[1:]], 0., .02, H, C, f)
    expected = np.r_[np.linalg.eigvals(jac), -np.array(C.losses[1:])]
    np.testing.assert_allclose(np.sort_complex(np.linalg.eigvals(full)), np.sort_complex(expected), atol=1e-12)
    b, c = core_tensors(y, 0., .02, H, C, f)
    step = 1e-5
    for j, e in enumerate(np.eye(3)):
        fd_b = (core_jacobian(y+step*e, 0., .02, H, C, f)
                -core_jacobian(y-step*e, 0., .02, H, C, f))/(2*step)
        fd_c = (core_tensors(y+step*e, 0., .02, H, C, f)[0]
                -core_tensors(y-step*e, 0., .02, H, C, f)[0])/(2*step)
        np.testing.assert_allclose(b[:, :, j], fd_b, rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(c[:, :, :, j], fd_c, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(b, b.swapaxes(1, 2), atol=0)
    for perm in permutations([1, 2, 3]):
        np.testing.assert_allclose(c, c.transpose(0, *perm), atol=0)
    np.testing.assert_allclose(cubic_hurwitz(jac)[:3], np.poly(jac)[1:], rtol=1e-9)


@pytest.mark.parametrize("a", [-1., 1.])
def test_lyapunov_formula_on_known_planar_normal_form(a):
    omega = 2.
    jac = np.array([[0., -omega, 0.], [omega, 0., 0.], [0., 0., -3.]])
    third = np.zeros((3, 3, 3, 3))
    third[0, 0, 0, 0], third[1, 1, 1, 1] = 6*a, 6*a
    for perm in set(permutations((0, 1, 1))):
        third[(0, *perm)] = 2*a
    for perm in set(permutations((0, 0, 1))):
        third[(1, *perm)] = 2*a
    normal = first_lyapunov(jac, np.zeros((3, 3, 3)), third)
    assert normal["l1"] == pytest.approx(2*a/omega)
    assert np.vdot(normal["p"], normal["q"]) == pytest.approx(1.)
    assert np.linalg.norm(normal["q"]) == pytest.approx(1.)


def test_interval_arithmetic_contains_high_precision_results():
    with localcontext() as ctx:
        ctx.prec = 80
        a, b = Decimal(1)/7, Decimal(3)/11
        ia, ib = _Interval(a), _Interval(b)
        for enclosure, value in [(ia+ib, a+b), (ia-ib, a-b), (ia*ib, a*b), (ia/ib, a/b),
                                 (ia.exp(), a.exp()), ((-ib).exp(), (-b).exp())]:
            assert enclosure.lo <= value <= enclosure.hi
        assert _Interval(-2, 3).exp().lo <= Decimal(-2).exp()
        assert _Interval(-2, 3).exp().hi >= Decimal(3).exp()
    with pytest.raises(ValueError):
        _Interval(1)/_Interval(-1, 1)
    with pytest.raises(ValueError):
        _Interval("NaN")


@pytest.mark.parametrize("supply,bounds", [(0., (.19, .30)), (1., (.145, .30))])
def test_continuous_path_certificate_covers_all_cells_and_float_hurwitz(supply, bounds):
    f = replace(F, leaf_supply=supply)
    report = certify_stable_path(bounds, (.01, .03), 0., H, C, f)
    assert report["cell_count"] > 1
    cells = report["cells"]
    assert Decimal(cells[0]["leaf_interval"][0]) == Decimal.from_float(bounds[0])
    assert Decimal(cells[-1]["leaf_interval"][1]) == Decimal.from_float(bounds[1])
    for left, right in zip(cells[:-1], cells[1:]):
        assert left["leaf_interval"][1] == right["leaf_interval"][0]
    assert all(Decimal(z) > 0 for z in report["strict_lower_bounds"].values())
    assert Decimal(report["Dprime_L_upper_bound"]) < 0
    # Midpoints also agree with the independently checked binary64 Jacobian.
    for cell in cells[::max(1, len(cells)//15)]:
        leaf = sum(map(float, cell["leaf_interval"]))/2
        exact = _branch_intervals(_Interval(leaf), 0., H, C, f)
        psi, d = float(exact["psi"].lo), float(exact["D"].lo)
        y = equilibrium_curve(psi, 0., H, C, f)[0]
        assert y[2] == pytest.approx(leaf, rel=1e-10)
        values = cubic_hurwitz(core_jacobian(y[:3], 0., d, H, C, f))
        for key, value in zip(["a1", "a2", "a3", "delta"], values):
            assert float(cell[key][0]) <= value <= float(cell[key][1])


@pytest.mark.parametrize("supply", [0., 1.])
def test_hopf_point_and_nonzero_transversality(supply):
    f = replace(F, leaf_supply=supply)
    hopf = hopf_point((-2.1, -1.9), 0., H, C, f)
    assert hopf["omega"] > .02
    assert abs(hopf["l1"]) > 1e-3
    assert abs(hopf["hurwitz"][-1]) < 1e-8
    assert hopf["eigenvalue_slope_D_refinement"][-1] > 0
    eigen = hopf["eigenvalues"]
    assert np.count_nonzero(abs(eigen.real) < 1e-8) == 2
    assert min(eigen.real) < -1


def test_certificate_does_not_pass_uncovered_or_unstable_paths():
    with pytest.raises(ValueError):
        certify_stable_path((.20, .30), (.01, .03), 0., H, C, F)
    with pytest.raises(RuntimeError, match="unresolved/nonstable"):
        certify_stable_path((.02, .30), (.01, .3), 0., H, C, F, max_depth=3)
    with pytest.raises(ValueError):
        certify_stable_path((.19, .30), (.01, .03), 0., H, C, replace(F, wood_supply=1.))
    with pytest.raises(ValueError):
        core_tensors([-1., -3.5, .2], 0., .02, H, C, F)
    with pytest.raises(ValueError):
        core_rhs([-1., -2., .2], 0., .02, H, C, replace(F, wood_supply=1.))


def test_subcritical_cycle_is_not_equilibrium_and_full_lift_preserves_budgets():
    f = replace(F, leaf_supply=1.)
    hopf = hopf_point((-2.1, -1.9), 0., H, C, f)
    cycle = periodic_orbit(hopf, .001, 0., H, C, f, leaf_bracket=(-2.1, -1.9), rtol=1e-10)
    assert hopf["l1"] > 0
    assert cycle["D"] < hopf["D"]
    assert abs(cycle["radial_floquet"]) > 1.0002
    assert cycle["full_state_closure"] < 1e-8
    np.testing.assert_allclose(cycle["full_floquet"][3:], np.exp(-np.array(C.losses[1:])*cycle["period_days"]))
    assert np.ptp(cycle["state"][:, 1]) > .001
    for y in cycle["state"][::40]:
        rhs = two_way_rhs(y, 0., cycle["D"], H, C, f)
        erx, exl, e = structural_fluxes(y, 0., cycle["D"], H, f)
        assert H.cx*rhs[0]+H.cl*rhs[1] == pytest.approx(erx-e, abs=1e-15)
        assert sum(rhs[2:]) == pytest.approx(carbon_budget(y[2:], hydraulic_factor(y[1], H), C)["nep"], abs=1e-15)
    with pytest.raises(ValueError):
        periodic_orbit(hopf, 0., 0., H, C, f, leaf_bracket=(-2.1, -1.9))


@pytest.mark.parametrize("supply", [0., 1.])
def test_formal_zero_leaf_limits_are_not_isolated_dead_attractors(supply):
    f = replace(F, leaf_supply=supply)
    eps = 1e-12
    wet = np.r_[0., 0., [eps]*4]
    dry = np.r_[-1., -3.5, [eps]*4]  # gx=0 but gs>0; FE tends to zero
    assert max(abs(two_way_rhs(wet, 0., .03, H, C, f))) < 1e-9
    assert max(abs(two_way_rhs(dry, 0., .03, H, C, f))) < 1e-9
    wet_eig = np.linalg.eigvals(two_way_jacobian(wet, 0., .03, H, C, f))
    dry_eig = np.linalg.eigvals(two_way_jacobian(dry, 0., .03, H, C, f))
    assert max(wet_eig.real) == pytest.approx(C.losses[0]*(2.25-1), rel=1e-6)
    assert np.count_nonzero(abs(dry_eig) < 1e-8) == 2
