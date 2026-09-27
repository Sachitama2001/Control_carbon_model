from dataclasses import replace
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest

from control_carbon.temperature_nsc import Parameters
from control_carbon.temperature_nsc_visit import (VisitTreeParameters, FixedEnvironment,
    visit_physiology,optimum_lai_source,visit_canopy_gpp,leaf_area)
from control_carbon.nsc_allocation import Mobilization, source_allocation, evaluate_allocation
from control_carbon.nsc_allocation_analysis import (donor_balance, reconstruct, frozen_audit,
    exclusion_bounds,rescue_audit)
from control_carbon.nsc_allocation_followups import native_positive_storage_step,rescue_increment

P = replace(Parameters(), nsc_mortality=False, mu_max=0)
V = VisitTreeParameters()
E = FixedEnvironment()
M = Mobilization((.03,.001,.003))


def test_piecewise_fractions_and_equality():
    # SLA chosen so LAI equals carbon mass exactly.
    v=replace(V,sla=200/2.2)
    for lai,expected,branch in ((3.1,0,"above"),(3,.1,"near"),
                                (2.95,.1,"near"),(2,.05,"far")):
        c,b=source_allocation(lai,3,1,v)
        assert c[0]==pytest.approx(expected)
        assert c.sum()==pytest.approx(1)
        assert b==branch


@pytest.mark.parametrize("variant",["V4","V5","V7_V4"])
def test_conservation_and_positive_orthant(variant):
    rng=np.random.default_rng(81)
    for _ in range(15):
        x=10**rng.uniform(-7,2,8)
        dx,d=evaluate_allocation(x,20,P,V,E,M,variant)
        assert d["budget"]==pytest.approx(0,abs=1e-13)
        assert d["investment"].sum()==pytest.approx(d["donor"].sum())
        for i in range(8):
            z=x.copy(); z[i]=0
            f,_=evaluate_allocation(z,20,P,V,E,M,variant)
            assert f[i]>=0


def test_quadratic_balance_and_scalar_elimination():
    assert donor_balance(0,2,3,4)==0
    n=donor_balance(2+3/5,2,3,4)
    assert n==pytest.approx(1)
    for variant in ("V4","V5"):
        x=reconstruct(.2,.05,20,P,V,M,variant)
        dx,_=evaluate_allocation(x,20,P,V,E,M,variant,.05)
        assert np.max(np.abs(dx[[0,1,2,3,5,6,7]]))<1e-12


def test_admissible_roots_satisfy_unforced_equations():
    report=frozen_audit(20,P,V,E,M,"V4",points=100)
    for row in report["roots"]:
        if row["admissible"]:
            dx,_=evaluate_allocation(row["state"],20,P,V,E,M,"V4")
            assert np.max(np.abs(dx))<1e-9


def test_original_c_allocation(tmp_path):
    source=Path("/mnt/d/ct/VISIT-matrix/visit_local")
    if not source.is_dir() or not shutil.which("gcc"):
        pytest.skip("native source/compiler unavailable")
    harness=tmp_path/"audit.c"
    harness.write_text('''#include <stdio.h>
#include "structure.h"
#include "prototype.h"
int main(void) {
 struct Pchar p={0}; struct Pmas m={0}; struct Pflx f={0};
 p.sla=150; p.alloc_ass=.1; p.alloc_abg=.67; p.season=1; p.opt_lai=3;
 double cases[]={0,2.0,2.9,3.0,3.1};
 for(int i=0;i<5;i++) {m.lai=cases[i]; f.epp=1;
 f_allocation(&p,&m,&f);
 printf("%.17g %.17g %.17g\\n",f.tpf,f.tpc,f.tpr);}
 return 0;
}
''')
    binary=tmp_path/"audit"
    subprocess.run(["gcc","-ffunction-sections","-Wl,--gc-sections","-I",str(source),str(harness),str(source/"allocation.c"),
                    "-lm","-o",str(binary)],check=True,capture_output=True)
    native=np.loadtxt(subprocess.check_output([str(binary)],text=True).splitlines())
    py=np.array([source_allocation(lai,3,1,V)[0] for lai in (0,2,2.9,3,3.1)])
    np.testing.assert_allclose(native,py,rtol=1e-14,atol=1e-15)


def test_sequential_storage_overshoot_and_budget():
    x=np.array([1.,2.,3.,.9]); a=np.array([.1,.2,.3]);r=np.array([.04,.036,.096])
    y,route=native_positive_storage_step(x,a,r,1)
    assert route==["storage","structure"]
    assert y[3]==pytest.approx(1.1)
    assert y.sum()-x.sum()==pytest.approx(a.sum()-r.sum())
    # A small donor can become negative in the source sequence: report it.
    y,_=native_positive_storage_step([0,.001,.001,0],a,r,10)
    assert y[1]<0 and y[2]<0


def test_rescue_is_structural_conservative_and_donor_bounded():
    for stem in (0,1e-12,.1,100):
        for root in (0,1e-12,.1,100):
            d=rescue_increment(0,stem,root,V)
            assert d.sum()==pytest.approx(0,abs=1e-15)
            assert np.min(np.array([0,stem,root])+d)>=0
    assert not rescue_increment(1,100,100,V).any()


def test_native_target_and_gpp_given_physiology(tmp_path):
    source=Path("/mnt/d/ct/VISIT-matrix/visit_local")
    if not source.is_dir() or not shutil.which("gcc"):
        pytest.skip("native source/compiler unavailable")
    harness=tmp_path/"audit_target.c"
    harness.write_text('''#include <stdio.h>
#include "structure.h"
#include "prototype.h"
static struct Grid g; static struct Loct l; static struct Pchar p; static struct Pmas m;
int main(void) {
 p.sla=150; p.eK=.48; p.ppfd_t=1000; p.rmf=1.3; p.lf=8.2e-5; p.rgf=.4;
 p.qTf0=p.qTc0=p.qTr0=2; l.doy=1; l.daylen[1]=12; m.lai=3.3;
 while(scanf("%lf %lf %lf",&l.tmp_sfc,&p.psat,&p.lue)==3) {
 f_q10_ar(&l,&p); f_opt_lai(&g,&l,&p);
 printf("%.17g %.17g\\n",p.opt_lai,f_gpp(&g,&l,&p,&m));
 }
 return 0;
}
''')
    binary=tmp_path/"audit_target"
    subprocess.run(["gcc","-ffunction-sections","-Wl,--gc-sections","-I",str(source),
                    str(harness),str(source/"ecophysiology.c"),str(source/"photosynthesis.c"),
                    "-lm","-o",str(binary)],check=True,capture_output=True)
    inputs=[];expected=[]
    for t in (0,10,20,30,34):
        v=visit_physiology(t,V,E)
        inputs.append(f'{t} {v["psat"]:.17g} {v["lue"]:.17g}')
        expected.append([optimum_lai_source(t,V,E)[0],visit_canopy_gpp(2,t,V,E)])
    output=subprocess.check_output([str(binary)],input="\n".join(inputs)+"\n",text=True)
    np.testing.assert_allclose(np.loadtxt(output.splitlines()),expected,rtol=1e-13,atol=1e-14)


def test_exclusion_bound_matches_concave_gpp():
    bounds=exclusion_bounds(20,P,V,E)
    assert bounds["positive_equilibria_excluded"]
    assert bounds["gpp_origin_slope"]==pytest.approx(visit_canopy_gpp(1e-5,20,V,E)/1e-5,rel=1e-5)
    assert bounds["gpp_upper_bound"]==pytest.approx(visit_canopy_gpp(1000,20,V,E))
    for leaf in (1e-12,1e-30,1e-100):
        assert visit_canopy_gpp(leaf,20,V,E)/leaf==pytest.approx(bounds["gpp_origin_slope"],rel=1e-10)


def test_rescue_boundary_is_not_a_classical_root():
    report=rescue_audit(20,P,V,E,M,points=100)
    assert len(report["boundaries"])==1
    q=report["boundaries"][0]
    assert q["attracting"] and q["far_branch_valid"]
    assert q["residual"]<1e-10
    assert q["normal_below"]>0>q["normal_above"]
    # Setting exactly critical leaf picks rescue OFF in the native rule.
    x=np.array(q["state"]); x[0]=report["critical_leaf"]
    dx,_=evaluate_allocation(x,20,P,V,E,M,"V7_V4")
    assert dx[0]<0
