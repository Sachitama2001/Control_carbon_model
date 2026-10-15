from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import expm

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, assimilation, carbon_matrix, hydraulic_factor, positive_carbon_equilibrium,
)
from control_carbon.spac_carbon_equilibria import stock_ratios
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_two_way import StructuralFeedback


H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3),
                     tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]), (.7, .9, .8), 10., .5, .25)
REF = positive_carbon_equilibrium(hydraulic_factor(transpiring_equilibria(0., .01, H)[0], H), C)


def variant(label):
    strengths = {"E0": (0, 0, 0), "E1": (1, 0, 0), "E2": (0, 1, 0),
                 "E3": (1, 1, 0), "E4": (1, 1, 1)}[label]
    return RationalSPAC(PressureSPAC(
        H, C, StructuralFeedback(REF[0], REF[1], .01, 1e-6, *strengths), 0., .101325))


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_rational_elimination_polynomial_and_capacity(label):
    m = variant(label)
    assert m.wet_leaf == pytest.approx(.25)
    assert m.opening(0.) == pytest.approx(4/9)
    assert m.equilibrium_polynomial(.002).degree() == (2 if label in ["E0", "E2"] else 3)
    for leaf in [1e-6, .005, .025, .1, .2, m.wet_leaf]:
        state, v, derivative = m.curve(leaf)
        np.testing.assert_allclose(m.rhs(state, v), 0., atol=1e-12)
        np.testing.assert_allclose(m.capacity(state), state[2:], atol=2e-14)
        polynomial = m.equilibrium_polynomial(v)
        assert abs(polynomial(leaf)) < 1e-14
        np.testing.assert_allclose(state[2:], stock_ratios(m.base)*leaf, atol=1e-14)
        if .005 <= leaf < .2:
            step = 1e-7
            fd = (m.curve(leaf+step)[1]-m.curve(leaf-step)[1])/(2*step)
            assert derivative == pytest.approx(fd, rel=2e-7, abs=1e-10)


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_positive_root_counts_and_discarded_zero_endpoint(label):
    m = variant(label)
    mark = m.landmarks()
    limit, fold = mark["zero_leaf_vpd_limit_mpa"], mark["fold_vpd_mpa"]
    if label in ["E0", "E1"]:
        cases = [(0, 1), (limit/2, 1), (limit, 1), ((limit+fold)/2, 2),
                 (fold, 1), (fold*1.1, 0)]
    elif label == "E4":
        cases = [(0, 1), (limit/2, 1), (limit, 0), (limit*1.1, 0)]
    else:
        cases = [(0, 1), (.001, 1), (.003, 1), (.02, 1), (.1, 1)]
    for v, count in cases:
        roots = m.equilibria(v)
        assert len(roots) == count
        for state in roots:
            assert state[2] > 0
            np.testing.assert_allclose(m.rhs(state, v), 0., atol=1e-10)
    if label == "E0":
        assert fold == pytest.approx(.0024720425531914894)
        assert mark["fold_leaf"] == pytest.approx(.025)


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_full_jacobian_matches_central_differences(label):
    m = variant(label)
    y = np.r_[-.7, -1.3, REF*np.array([.8, 1.1, .9, 1.2])]
    jac = m.jacobian(y, .002)
    for j in range(6):
        step = 1e-6*max(1., abs(y[j]))
        direction = np.eye(6)[j]*step
        fd = (m.rhs(y+direction, .002)-m.rhs(y-direction, .002))/(2*step)
        np.testing.assert_allclose(jac[:, j], fd, rtol=2e-7, atol=1e-8)
    np.testing.assert_allclose(jac[:2], m.base.jacobian(y, .002)[:2])
    assert not np.allclose(jac[2:, 2], m.base.jacobian(y, .002)[2:, 2])


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_shape_proof_and_numerical_derivative_signs(label):
    m = variant(label)
    bounds = m.shape_bounds()
    assert all(float(gap[0]) > 0 for gap in bounds["strict_gaps"].values())
    leaves = np.linspace(1e-6, m.wet_leaf, 101)
    signs = np.array([m.curve(leaf)[2] for leaf in leaves])
    if label in ["E2", "E3", "E4"]:
        assert np.all(signs < 0)
    else:
        assert np.sum(signs[:-1]*signs[1:] < 0) == 1


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_closure_invariance_budget_and_exact_carbon_tail(label):
    m = variant(label)
    y = np.r_[-2., -4.2, REF]
    np.testing.assert_array_equal(m.rhs(y, .002)[:2], [0., 0.])
    np.testing.assert_array_equal(m.capacity(y), np.zeros(4))
    assert m.assimilation(y) == 0
    run = m.simulate(y, np.linspace(0, 100, 11), .002, rtol=1e-10)
    exact = np.array([expm(carbon_matrix(C)*t)@REF for t in run["day"]])
    np.testing.assert_allclose(run["state"][:, 2:], exact, rtol=2e-8)
    wet = np.r_[-.7, -1.4, REF]
    gpp = m.assimilation(wet)
    respiration = np.dot((1-np.array(C.litter_fractions))*np.array(C.losses[:3]), REF[:3])
    assert sum(m.rhs(wet, .002)[2:]) == pytest.approx(gpp-respiration-C.losses[3]*REF[3])


def test_saturation_origin_slope_and_linearization_distinction():
    m = variant("E0")
    for leaf in [1e-8, .01, 10.]:
        y = np.r_[-.5, -1., leaf, 1., 1., 1.]
        h = hydraulic_factor(y[1], H)
        rational = m.assimilation(y)
        assert 0 < rational < C.assimilation_max*h
        assert rational <= assimilation(leaf, h, C)
        assert rational <= C.assimilation_max*h*C.canopy_coefficient*leaf
    y = np.r_[-.5, -1., 1e-10, 1., 1., 1.]
    assert m.assimilation(y)/y[2] == pytest.approx(
        C.assimilation_max*.75*C.canopy_coefficient, rel=1e-8)
    y[2] = 1e8
    assert m.assimilation(y) == pytest.approx(C.assimilation_max*.75, rel=1e-8)
    h = .7
    leaf_exact = C.allocation[0]*C.assimilation_max*h/C.losses[0]-1/C.canopy_coefficient
    y = np.r_[-.6, -4*(1-h), stock_ratios(m.base)*leaf_exact]
    assert m.rhs(y, .002)[2] == pytest.approx(0., abs=1e-16)
    assert C.allocation[0]*C.assimilation_max*h*C.canopy_coefficient-C.losses[0] > 0


def test_unsupported_parameters_and_invalid_inputs_fail_explicitly():
    m = variant("E0")
    with pytest.raises(ValueError):
        RationalSPAC(replace(m.base, root_potential_mpa=-.5))
    for bad in [0., -.1, .3, np.nan]:
        with pytest.raises(ValueError):
            m.curve(bad)
    with pytest.raises(ValueError):
        m.equilibria(-.001)
    with pytest.raises(ValueError):
        m.jacobian(np.r_[-1., -3., REF], .002)
    with pytest.raises(ValueError):
        m.simulate(np.r_[-1., -2., REF], np.array([0., 0.]), .002)


@pytest.mark.parametrize("label", ["E0", "E1", "E4"])
def test_adjacent_float_endpoint_and_fold_classification(label):
    m = variant(label)
    mark = m.landmarks()
    limit = mark["zero_leaf_vpd_limit_mpa"]
    below, above = np.nextafter(limit, 0.), np.nextafter(limit, np.inf)
    assert len(m.equilibria(below)) == 1
    assert len(m.equilibria(above)) == (0 if label == "E4" else 2)
    if label != "E4":
        fold = mark["fold_vpd_mpa"]
        assert len(m.equilibria(np.nextafter(fold, 0.))) == 2
        assert len(m.equilibria(np.nextafter(fold, np.inf))) == 0


@pytest.mark.parametrize("label,vpd", [("E0", .00244), ("E1", .0009),
                                      ("E2", .00244), ("E3", .0061), ("E4", .003)])
def test_polynomial_roots_independently_match_monotone_brackets(label, vpd):
    m = variant(label)
    algebraic = sorted(z.real for z in m.equilibrium_polynomial(vpd).roots()
                       if abs(z.imag) < 1e-12 and 0 < z.real <= m.wet_leaf)
    bracketed = sorted(y[2] for y in m.equilibria(vpd))
    np.testing.assert_allclose(algebraic, bracketed, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_xylem_only_closed_decay_bound(label):
    m = variant(label)
    state = np.r_[-1.8, -3.2, REF]
    derivative = m.rhs(state, .002)
    assert derivative[0] == 0
    assert derivative[1] <= 0
    sigma = C.losses[0]-C.allocation[0]*C.assimilation_max*C.canopy_coefficient/4
    assert sigma*YEAR_DAYS == pytest.approx(.35)
    assert derivative[2] <= -sigma*state[2]


@pytest.mark.parametrize("label", ["E0", "E1", "E2", "E3", "E4"])
def test_formal_zero_carbon_equilibria_are_not_counted_as_positive(label):
    m = variant(label)
    closed = np.r_[-1.8, -4.2, np.zeros(4)]
    np.testing.assert_array_equal(m.rhs(closed, .002), np.zeros(6))
    assert np.all(np.isfinite(m.jacobian(closed, .002)))
    if label in ["E0", "E1"]:
        opening = .7
        fl, _, _ = m.factors(np.zeros(6))
        psi = -4*(1-opening)
        state = np.r_[fl*psi/(1+fl), psi, np.zeros(4)]
        d = (4*H.gx_max/(3*H.gs_max)*(5-4*opening-1/opening)*2*fl/(1+fl))
        np.testing.assert_allclose(m.rhs(state, m.base.pressure(d)), 0., atol=1e-12)
    else:
        np.testing.assert_array_equal(m.rhs(np.zeros(6), .002), np.zeros(6))
        np.testing.assert_array_equal(m.rhs(np.r_[-1., -3.2, np.zeros(4)], .002), np.zeros(6))
    if label == "E4":
        np.testing.assert_array_equal(m.rhs(np.r_[-1., -1., np.zeros(4)], .002), np.zeros(6))
