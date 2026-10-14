from dataclasses import replace
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from control_carbon.manzoni_hydraulics import HydraulicParameters, conductances, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, carbon_matrix, hydraulic_factor, positive_carbon_equilibrium,
)
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_capacity import reference_scales, diagnose_capacity
from control_carbon.spac_carbon_equilibria import (
    stock_ratios, leaf_opening, positive_domain, leaf_curve, baseline_shape_proof,
    baseline_landmarks, positive_equilibria,
)

H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
PSI = transpiring_equilibria(0., .01, H)[0]
REF = positive_carbon_equilibrium(hydraulic_factor(PSI, H), C)
F = StructuralFeedback(REF[0], REF[1], .01, 1e-6)
M = PressureSPAC(H, C, F, 0., .101325)


def variant(label):
    strengths = {"E0": (0, 0, 0), "E1": (1, 0, 0), "E4": (1, 1, 1)}[label]
    return replace(M, feedback=replace(F, leaf_supply=strengths[0], leaf_demand=strengths[1],
                                      wood_supply=strengths[2]))


@pytest.mark.parametrize("label", ["E0", "E1", "E4"])
def test_elimination_matches_prior_curve_full_equations_and_derivatives(label):
    m = variant(label)
    top = positive_domain(m)["leaf_upper_closed_at_zero_vpd"]
    for leaf in [1e-10, 1e-6, .01, .08, .2, top]:
        y, v, derivative = leaf_curve(leaf, m)
        np.testing.assert_allclose(y[2:], stock_ratios(m)*leaf, atol=1e-14)
        np.testing.assert_allclose(m.rhs(y, v), 0., atol=1e-10)
        data = diagnose_capacity(y, v, m, reference_scales(m), equilibrium=y)
        assert data["potential_norm_scaled"] < 1e-12
        if leaf > 1e-6 and leaf < top:
            old, old_v, old_slope = m.equilibrium_curve(y[1])
            np.testing.assert_allclose(y, old, atol=1e-12)
            assert v == pytest.approx(old_v, abs=1e-13)
            assert derivative == pytest.approx(old_slope*(-2*H.psi50_s)*leaf_opening(leaf, m)[1], rel=1e-10)
            step = 1e-7
            fd = (leaf_curve(leaf+step, m)[1]-leaf_curve(leaf-step, m)[1])/(2*step)
            assert derivative == pytest.approx(fd, rel=1e-7, abs=1e-10)
    assert leaf_curve(top, m)[1] == 0.


@pytest.mark.parametrize("leaf", [0., 1e-12, 1e-8, 1e-5, .01])
def test_opening_derivative_avoids_zero_leaf_cancellation(leaf):
    opening, derivative = leaf_opening(leaf, M)
    with localcontext() as ctx:
        ctx.prec = 60
        r = Decimal.from_float(C.losses[0])/(Decimal.from_float(C.allocation[0])*Decimal.from_float(C.assimilation_max))
        q = Decimal.from_float(C.canopy_coefficient)
        z = q*Decimal.from_float(leaf)
        if leaf == 0:
            exact_h, exact_dh = r/q, r/2
        else:
            decay = (-z).exp()
            exact_h = r*Decimal.from_float(leaf)/(1-decay)
            exact_dh = r*(1-decay-z*decay)/(1-decay)**2
    assert opening == pytest.approx(float(exact_h), rel=1e-14)
    assert derivative == pytest.approx(float(exact_dh), rel=1e-12)


@pytest.mark.parametrize("label", ["E0", "E1", "E4"])
def test_global_root_counts_including_open_zero_and_closed_zero_vpd(label):
    m = variant(label)
    landmarks = baseline_landmarks(m)
    limit, fold = landmarks["zero_leaf_vpd_limit_mpa"], landmarks["fold_vpd_mpa"]
    if fold is None:
        cases = [(0., 1), (limit/2, 1), (limit, 0), (limit*1.1, 0)]
    else:
        cases = [(0., 1), (limit/2, 1), (limit, 1), ((limit+fold)/2, 2), (fold, 1), (fold*1.1, 0)]
    for v, count in cases:
        roots = positive_equilibria(v, m, reference_scales(m))
        assert len(roots["roots"]) == count
        for root in roots["roots"]:
            assert root["state"][2] > 0
            assert root["diagnostics"]["full_residual_norm_scaled_per_day"] < 1e-8
            np.testing.assert_allclose(root["diagnostics"]["capacity"], root["state"][2:], atol=1e-12)
        if fold is not None and v == fold:
            assert roots["roots"][0]["diagnostics"]["equilibrium_center_dimensions"] == 1
        if fold is not None and v == (limit+fold)/2:
            assert {root["diagnostics"]["equilibrium_type"] for root in roots["roots"]} == {
                "hyperbolic_sink", "hyperbolic_saddle",
            }
    assert positive_equilibria(limit, m, reference_scales(m))["zero_leaf_limit_hit"]


@pytest.mark.parametrize("label", ["E0", "E1", "E4"])
def test_positive_endpoint_neighbor_is_not_rounded_to_excluded_zero_leaf(label):
    m = variant(label)
    limit = baseline_landmarks(m)["zero_leaf_vpd_limit_mpa"]
    v = np.nextafter(limit, -np.inf if label == "E4" else np.inf)
    roots = positive_equilibria(v, m, reference_scales(m))["roots"]
    assert len(roots) == (1 if label == "E4" else 2)
    assert 0 < min(root["state"][2] for root in roots) < 1e-12
    assert all(root["diagnostics"]["full_residual_norm_scaled_per_day"] < 1e-8 for root in roots)


def test_e4_global_monotonic_margin_and_nonmonotone_hydraulic_base():
    m = variant("E4")
    proof = baseline_shape_proof(m)
    assert Decimal(proof["strict_gap_interval"][0]) > 0
    domain = positive_domain(m)
    leaves = np.geomspace(1e-10, domain["leaf_upper_closed_at_zero_vpd"], 101)
    assert all(leaf_curve(leaf, m)[2] < 0 for leaf in leaves)
    e0 = baseline_landmarks(M)
    assert e0["fold_vpd_mpa"] == pytest.approx(.0024720425531914894)
    e1 = baseline_landmarks(variant("E1"))
    assert e1["fold_vpd_mpa"] == pytest.approx(.001546784415086454, abs=1e-13)
    assert e1["fold_leaf"] > e0["fold_leaf"]


@pytest.mark.parametrize("label", ["E0", "E1", "E4"])
def test_fully_closed_water_invariance_capacity_zero_and_positive_carbon_decay(label):
    m = variant(label)
    y = np.r_[-1.7, -4.2, REF]
    data = diagnose_capacity(y, .00244, m, reference_scales(m))
    np.testing.assert_array_equal(data["capacity"], np.zeros(4))
    np.testing.assert_array_equal(data["rhs"][:2], np.zeros(2))
    np.testing.assert_allclose(data["rhs"][2:], carbon_matrix(C)@REF, atol=1e-15)
    future = expm(carbon_matrix(C)*750*YEAR_DAYS)@REF
    assert np.all(future > 0) and max(future) < 1e-10
    solution = m.simulate(y, np.linspace(0., 10*YEAR_DAYS, 41), lambda t: .00244, rtol=1e-10)
    np.testing.assert_allclose(solution["state"][:, :2], y[:2][None, :]+np.zeros((len(solution["time_days"]), 2)), atol=1e-13)
    exact = np.array([expm(carbon_matrix(C)*t)@REF for t in solution["time_days"]])
    np.testing.assert_allclose(solution["state"][:, 2:], exact, rtol=1e-8, atol=1e-10)


def test_xylem_closed_alone_is_not_zero_capacity_or_frozen_water():
    y = np.r_[-1.5, -3.5, REF]
    data = diagnose_capacity(y, .00244, M, reference_scales(M))
    gx, gs = conductances(y[1], H)
    assert gx == 0 and gs > 0
    assert data["rhs"][0] == 0 and data["rhs"][1] < 0
    assert np.all(data["capacity"] > 0)


def test_e4_zero_leaf_curve_derivative_has_a_finite_limit():
    m = variant("E4")
    values = [leaf_curve(leaf, m)[2] for leaf in [1e-12, 1e-10, 1e-8]]
    assert values[0] == pytest.approx(values[1], rel=1e-7)
    assert values[0] == pytest.approx(values[2], rel=1e-5)
    assert values[0] < 0


def test_e4_unique_root_can_recover_stability_in_formal_high_vpd_range():
    m = variant("E4")
    for leaf, status in [(.001, "unstable"), (2e-5, "stable")]:
        _, v, _ = leaf_curve(leaf, m)
        assert v > .003
        roots = positive_equilibria(v, m, reference_scales(m))["roots"]
        assert len(roots) == 1
        assert roots[0]["diagnostics"]["equilibrium_stability"] == status
        np.testing.assert_allclose(roots[0]["state"][2], leaf, atol=1e-13)


def test_xylem_closed_low_leaf_keeps_gs_positive_in_e4_before_domain_exit():
    m = variant("E4")
    stocks = stock_ratios(m)*(2*F.minimum_leaf)
    initial = np.r_[-1.7, -3.5, stocks]
    result = m.simulate(initial, np.linspace(0, 10*YEAR_DAYS, 101), lambda t: .00244, rtol=1e-10)
    assert result["status"] == "vegetation_domain_exit"
    assert np.all(result["state"][:, 1] <= -3.5)
    assert np.all(result["state"][:, 1] > -4.)
    assert conductances(result["state"][-1, 1], H)[1] > 0
    assert result["state"][-1, 2] == pytest.approx(F.minimum_leaf, rel=1e-8)


def test_flux_timeseries_handles_exact_kinks_without_fabricating_jacobians(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from analyze_spac_carbon_equilibria import flux_capacity_record
    eq = leaf_curve(.1, M)[0]
    y = np.r_[-1.5, -4., REF]
    row = flux_capacity_record("kink", 1., y, .00244, M, reference_scales(M), eq)
    assert row["jacobian_status"] == "undefined_at_conductance_kink"
    assert row["gx"] == row["gs"] == row["GPP"] == row["capacity_total"] == 0.
    assert row["identity_residual_max"] < 1e-15
    with pytest.raises(ValueError, match="kink"):
        M.jacobian(y, .00244)


def test_proof_hypotheses_and_domains_fail_explicitly():
    with pytest.raises(ValueError, match="beta1"):
        baseline_shape_proof(replace(variant("E4"), feedback=replace(variant("E4").feedback, wood_exponent=1.25)))
    with pytest.raises(ValueError, match="margin"):
        baseline_shape_proof(replace(variant("E4"), feedback=replace(variant("E4").feedback, leaf_floor=0.)))
    with pytest.raises(ValueError):
        positive_domain(replace(M, root_potential_mpa=-3.1))
    for leaf in [-.1, 0., 1., np.nan]:
        with pytest.raises(ValueError):
            leaf_curve(leaf, M)
    with pytest.raises(ValueError):
        positive_equilibria(-.001, M, reference_scales(M))
