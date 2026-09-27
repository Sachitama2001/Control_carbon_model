from dataclasses import replace
from pathlib import Path
import math
import numpy as np
import pytest

from control_carbon.temperature_nsc import Parameters
from control_carbon.temperature_nsc_visit import (
    FixedEnvironment, VisitTreeParameters, evaluate_variant, leaf_area,
    leaf_carbon_from_lai, living_tissue, load_tky_tree_parameters,
    optimum_lai_source, q10_acclimated, size_dependent_respiration,
    visit_canopy_gpp)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = Path("/mnt/d/ct/VISIT-matrix/visit_local/INPUT/Config_TKY.xlsx")
X = np.array([2., 80., 15., 100., .2, 5., 1., 1.])


def test_tky_loader_uses_tree_column():
    if not CONFIG.is_file():
        pytest.skip("pinned TKY workbook unavailable")
    p = load_tky_tree_parameters(CONFIG)
    assert p.sla == 150
    assert p.lue0 == .052
    assert p.kmci == 40
    assert p.cmpcd0 == 50
    assert p.alloc_ass == .1
    assert p.rmc_s == .033
    assert p.rmc_h == .004
    assert p.f_size_stem == 40


def test_lai_conversion_is_inverse_and_hand_computable():
    p = VisitTreeParameters(sla=150)
    assert leaf_area(2, p) == pytest.approx(3.3)
    assert leaf_carbon_from_lai(3.3, p) == pytest.approx(2)


def test_q10_source_algebra():
    assert q10_acclimated(2, 15) == 2
    assert q10_acclimated(2, 25) == pytest.approx(2*math.exp(-.09))


def test_optimum_lai_valid_and_invalid_branches():
    p, env = VisitTreeParameters(), FixedEnvironment()
    lai, status = optimum_lai_source(20, p, env)
    assert status == "valid" and lai >= 0
    frozen, status = optimum_lai_source(0, p, env)
    assert frozen == 0 and status == "psat_nonpositive"
    bad, status = optimum_lai_source(20, p, replace(env, day_length=1e-9))
    assert np.isnan(bad) and status == "source_denominator_nonpositive"


def test_visit_gpp_zero_leaf_and_saturating():
    p, env = VisitTreeParameters(), FixedEnvironment()
    assert visit_canopy_gpp(0, 20, p, env) == 0
    values = [visit_canopy_gpp(v, 20, p, env) for v in (1, 2, 100)]
    assert 0 < values[0] < values[1] < values[2]


def test_living_tissue_and_size_respiration():
    p = VisitTreeParameters()
    assert living_tissue(0, 40) == 0
    assert living_tissue(.1, 40) == pytest.approx(.1)
    assert living_tissue(100, 40) < 100
    r = size_dependent_respiration(X[:3], 20, p)
    assert np.all(r > 0)
    assert r[1]/X[1] < p.rmc_s/1000*q10_acclimated(p.qtc0,20)**.5


@pytest.mark.parametrize("variant", ["V0", "V1f", "V1t", "V2", "V3"])
def test_first_variants_have_no_nsc_mortality_and_close_budget(variant):
    base = replace(Parameters(), nsc_mortality=False, mu_max=0)
    dx, d = evaluate_variant(X, 20, base, VisitTreeParameters(), FixedEnvironment(), variant)
    assert d["nsc_mortality_enabled"] is False
    assert d["budget"] == pytest.approx(0, abs=1e-14)
    assert np.isfinite(dx).all()


def test_mortality_gate_rejects_accidental_reintroduction():
    with pytest.raises(ValueError, match="mortality disabled"):
        evaluate_variant(X, 20, Parameters(), VisitTreeParameters(), FixedEnvironment(), "V1t")
