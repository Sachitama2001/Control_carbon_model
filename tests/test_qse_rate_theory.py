"""Checks for algebraic proof certificates; no parameter-search experiments."""
import numpy as np
import pytest
from control_carbon.qse_rate_theory import (
    COMOVING_FOLD_RATE,crossing_displacement,critical_rate,cubic_field,
    bounded_input_coefficients,bounded_production,bounded_field,smooth_certificates,
)


def test_finite_ramp_critical_rate_is_not_comoving_fold():
    r=critical_rate()
    assert r>COMOVING_FOLD_RATE
    assert crossing_displacement(r)==pytest.approx(1.5,abs=1e-10)
    assert crossing_displacement(r*.9)>1.5
    assert crossing_displacement(r*1.1)<1.5
    with pytest.raises(ValueError):
        critical_rate(1)


def test_bounded_carbon_balance_factorization_and_positive_domain():
    for lam in np.linspace(0,1.5,9):
        p,a,h2=bounded_input_coefficients(lam)
        assert min(p,a,h2)>0
        for x in np.r_[0,np.geomspace(1e-6,1e3,50),lam+np.array([1,2,3])]:
            production=bounded_production(x,lam)
            assert 0<=production<=p+a
            assert bounded_field(x,lam)==pytest.approx(production-x,abs=1e-10)
        assert bounded_field(0,lam)==0
        assert cubic_field(1,lam)>=0 and cubic_field(4.5,lam)<=0


def test_simple_root_stabilities_are_two_negative_one_positive():
    for lam in (0,.75,1.5):
        for j,sign in ((1,-1),(2,1),(3,-1)):
            root=lam+j
            h=1e-5
            slope=(bounded_field(root+h,lam)-bounded_field(root-h,lam))/(2*h)
            assert sign*slope>0


def test_smooth_ramp_invariant_tube_and_fast_bound():
    for bounded in (False,True):
        c=smooth_certificates(bounded=bounded)
        assert c["net_stock_change_fast"]==-.5
        assert c["net_stock_change_slow"]==1.5
        assert c["fast_duration_strict_upper"]<c["slow_duration_lower"]
        for lam in np.linspace(0,1.5,9):
            for y in np.linspace(.5,1,9):
                x=lam+2+y
                q=(x/((1+x)*(bounded_input_coefficients(lam)[2]+x*x)) if bounded else 1)
                assert c["q_lower"]<=q<=c["q_upper"]


def test_total_potential_can_be_positive_while_carbon_flux_is_negative():
    p=np.array([1.,-.2]);loss=np.array([1.,10.])
    assert p.sum()>0 and loss@p<0


def test_upper_qse_gap_and_pointwise_capacity_are_distinct():
    x,lam=1.5,0
    assert (lam+3)-x>0
    assert bounded_production(x,lam)-x<0
