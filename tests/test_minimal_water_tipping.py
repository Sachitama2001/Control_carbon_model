import numpy as np
import pytest
from scipy.integrate import solve_ivp

from control_carbon.minimal_water_tipping import (
    MinimalTipping, advance_relative, crossing_log_shift, dynamic_jacobian,
    dynamic_rhs, finite_critical_rate, moving_rhs, scalar_rhs,
    simulate_drying, simulate_ramp,
)


def test_roots_and_frozen_stability():
    model = MinimalTipping()
    lo, hi = model.roots
    assert (lo, hi) == pytest.approx((.8, 1.25))
    for scale in [1.,1.5,2.]:
        for x in [0,lo,hi]:
            assert scalar_rhs(scale*x,scale,model) == pytest.approx(0,abs=1e-14)
        # f'(lambda*x*) = m(1-x*^2)/(1+x*^2), independent of lambda.
        assert model.mortality*(1-hi*hi)/(1+hi*hi) < 0
        assert model.mortality*(1-lo*lo)/(1+lo*lo) > 0


def test_rate_criterion_and_small_shift_no_tip():
    m = MinimalTipping()
    rc = finite_critical_rate(m)
    assert rc == pytest.approx(.00507617318728073,rel=1e-8)
    assert rc > m.comoving_fold_rate
    assert crossing_log_shift(rc,m) == pytest.approx(np.log(2))
    assert np.isinf(finite_critical_rate(MinimalTipping(final_scale=1.5)))
    assert crossing_log_shift(rc*.9,m) > np.log(2)
    assert crossing_log_shift(rc*1.1,m) < np.log(2)


def test_ramp_basin_threshold_not_only_endpoint_transient():
    m = MinimalTipping()
    rc=finite_critical_rate(m)
    below=simulate_ramp(.95*rc,m,hold_years=1000)
    above=simulate_ramp(1.05*rc,m,hold_years=1000)
    assert below['margin']>0 and above['margin']<0
    assert below['carbon'][-1] == pytest.approx(2.5,rel=1e-7)
    assert above['carbon'][-1] < 1e-8
    at=simulate_ramp(rc,m,hold_years=0)
    assert abs(at['margin']) < 1e-8


def test_month_year_flow_semigroup_and_direct_ode():
    m=MinimalTipping()
    x=m.roots[1]
    year=advance_relative(x,1.,.01,m)
    for _ in range(12): x=advance_relative(x,1/12,.01,m)
    assert x == pytest.approx(year,abs=1e-11)
    direct=solve_ivp(lambda t,y:[scalar_rhs(y[0],np.exp(.01*t),m)],
                     [0,1],[m.roots[1]],rtol=1e-12,atol=1e-14,first_step=.01)
    assert direct.success
    assert year*np.exp(.01) == pytest.approx(direct.y[0,-1],abs=1e-10)
    monthly=simulate_ramp(.01,m,update_years=1/12,hold_years=20)
    annual=simulate_ramp(.01,m,update_years=1.,hold_years=20)
    assert monthly['margin'] == pytest.approx(annual['margin'],abs=1e-9)
    assert monthly['carbon'][-1] == pytest.approx(annual['carbon'][-1],abs=1e-9)


def test_water_carbon_budgets_and_coupling():
    # Physical reconstruction C0=100 MgC/ha, eta=.02 MgC/ha/mm, ell=12/y.
    C0, eta, ell, m, a, lam, c=100.,.02,12.,.1,2.05,1.4,.9
    P0=a*m*C0/eta
    W=P0*lam/(ell*(lam*lam+c*c))
    uptake=ell*c*c*W
    loss=ell*lam*lam*W
    assert P0*lam-uptake-loss == pytest.approx(0,abs=1e-10)
    npp=eta*uptake; death=m*C0*c; dead=50.; kd=.03
    live_rhs=C0*scalar_rhs(c,lam)
    dead_rhs=death-kd*dead
    assert live_rhs+dead_rhs == pytest.approx(npp-kd*dead,abs=1e-12)
    J=dynamic_jacobian([1.2,.7])
    assert J[0,1]>0 and J[1,0]<0
    for j in range(2):
        h=np.eye(2)[j]*1e-6
        fd=(dynamic_rhs(np.array([1.2,.7])+h,2.1)-dynamic_rhs(np.array([1.2,.7])-h,2.1))/(2e-6)
        np.testing.assert_allclose(fd,J[:,j],rtol=1e-8)
    assert dynamic_rhs([0.,.5],2.1)[0] == 0
    assert dynamic_rhs([1.,0.],2.1)[1] > 0


def test_dynamic_rainfall_only_outcomes_and_frozen_branch():
    for a in np.linspace(2.002,10,50):
        c=(a+np.sqrt(a*a-4))/2
        np.testing.assert_allclose(dynamic_rhs([c,1/c],a),0,atol=1e-12)
        assert np.max(np.linalg.eigvals(dynamic_jacobian([c,1/c])).real)<0
    fast=simulate_drying(100.)
    fine=simulate_drying(100.,rtol=2e-11)
    slow=simulate_drying(1000.)
    assert fast['certified_bare'] and fine['certified_bare']
    assert not slow['certified_bare']
    target=(2.002+np.sqrt(2.002**2-4))/2
    assert slow['state'][0] == pytest.approx(target,abs=1e-8)
    np.testing.assert_allclose(fast['state'],fine['state'],atol=1e-8)


def test_smooth_rainfall_protocol_matters_but_tipping_persists():
    assert simulate_drying(10.,profile="smoothstep")['certified_bare']
    smooth=simulate_drying(100.,profile="smoothstep")
    assert not smooth['certified_bare']
    assert smooth['state'][0] == pytest.approx(1.04573253849,abs=1e-8)


@pytest.mark.parametrize('kwargs',[{'productivity':2.},{'mortality':0.},{'final_scale':1.},{'mortality':np.nan}])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError): MinimalTipping(**kwargs)
