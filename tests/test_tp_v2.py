import ctypes
from dataclasses import replace
from pathlib import Path
import subprocess

import numpy as np
import pytest
from numpy.testing import assert_allclose

from control_carbon import temperature_precipitation as v1
from control_carbon.tp_processes import (COMMIT, ProcessOptions, process_provenance,
    visit_gpp, visit_temperature, visit_quadratic_supply, visit_c3_soil_scalar)
from control_carbon.tp_experiments_v2 import Model, Solver, SpinupSettings, integrate, spinup


def test_scalar_hand_values_and_limits():
    assert visit_temperature(27, 0, 27, 45) == 1
    assert visit_temperature(0, 0, 27, 45) == 0
    assert visit_temperature(70, 0, 27, 45) == 0
    a = (30-45)*30
    assert_allclose(visit_temperature(30, 0, 27, 45), a/(a-9))
    assert visit_gpp(1, 12, .5, .04, 1000, 0) == 0
    assert visit_gpp(0, 12, .5, .04, 1000, 3) == 0
    assert visit_quadratic_supply(0, 4) == 0
    assert_allclose(visit_quadratic_supply(3, 3, 1), 3)
    assert_allclose(visit_quadratic_supply(2, 4, .85), (6-np.sqrt(36-4*.85*8))/(2*.85))
    assert visit_c3_soil_scalar(0, 240, .3) == .05
    assert_allclose(visit_c3_soil_scalar(240, 240, 1), .525)


def test_native_c_photosynthesis_and_water_limiter(tmp_path):
    source = Path(__file__).resolve().parents[2]/"VISIT-matrix/visit_local"
    if not (source/"photosynthesis.c").exists():
        pytest.skip("pinned old VISIT source unavailable")
    assert subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip() == COMMIT
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "photosynthesis.c", "hydro_balance.c", "definition.h", "structure.h"], check=True)
    # Actual source lines from the soil-evaporation block, not a new C formula.
    native_lines = "\n".join((source/"hydro_balance.c").read_text().splitlines()[128:133])
    cfile = tmp_path/"bridge.c"
    cfile.write_text('''#include <stdio.h>
#include <math.h>
#include "structure.h"
#include "prototype.h"
double probe_gpp(double psat,double hours,double ek,double lue,double ppfd,double lai){
 struct Grid grid={0}; struct Loct loct={0}; struct Pchar p={0}; struct Pmas m={0};
 p.psat=psat;p.eK=ek;p.lue=lue;p.ppfd_t=ppfd;loct.doy=0;loct.daylen[0]=hours;m.lai=lai;
 return f_gpp(&grid,&loct,&p,&m);
}
double probe_tem(double t){
 struct Grid g={0};struct Loct l={0};struct Pchar p={0};
 g.fieldcap=1;l.soilwtr_h=1;l.tmp_sfc=t;l.tmp10_soil=10;
 p.phototype=3;p.ci=280;p.topt0=24.2;p.tmin=0;p.tmax=45;p.pmax=1;
 f_pc_sat(&g,&l,&p);return p.psat;
}
double probe_supply(double s,double e){
 struct Mass m={0};struct Loct l={0};struct Mass *mass=&m;struct Loct *loct=&l;
 double aa,bb,cc;m.sw30=s;l.pm_evpr=e;
''' + native_lines + '\nreturn l.evpr;\n}\n')
    lib = tmp_path/"bridge.so"
    subprocess.run(["gcc", "-shared", "-fPIC", "-O2", "-I", str(source), str(cfile),
                    str(source/"photosynthesis.c"), "-lm", "-o", str(lib)], check=True, capture_output=True)
    dll = ctypes.CDLL(str(lib))
    for name, argc in (("probe_gpp", 6), ("probe_tem", 1), ("probe_supply", 2)):
        f = getattr(dll, name)
        f.argtypes = [ctypes.c_double]*argc
        f.restype = ctypes.c_double
    for t in (0, 10, 27, 31, 44, 45):
        assert_allclose(visit_temperature(t, 0, 27, 45), dll.probe_tem(t), atol=1e-14)
    for psat in (0, 1, 10):
        for lai in (0, .2, 3, 10):
            args = (psat, 12, .5, .04, 1000, lai)
            assert_allclose(visit_gpp(*args), dll.probe_gpp(*args), atol=1e-14)
    for s, e in ((0, 4), (1, 3), (100, .1), (3, 3)):
        assert_allclose(visit_quadratic_supply(s, e), dll.probe_supply(s, e), atol=2e-14)


def variants():
    return [ProcessOptions(), ProcessOptions(photo_temperature="visit_tem"),
            ProcessOptions(photo_canopy="visit_monsi"), ProcessOptions(photo_water="visit_soil_with_plant_gate"),
            ProcessOptions(evap_supply="visit_quadratic"),
            ProcessOptions(drainage="visit_baseflow_plus_excess"),
            ProcessOptions(photo_temperature="visit_tem", photo_canopy="visit_monsi", evap_supply="visit_quadratic")]


@pytest.mark.parametrize("options", variants())
def test_each_closure_budget_nonnegative_and_capacity_boundaries(options):
    model = Model(control=v1.Controls(mortality="additive"), processes=options)
    rng = np.random.default_rng(5)
    for _ in range(4):
        c = rng.uniform(1, 200, 4)
        x = np.r_[c, v1.capacity(c, model.parameters)*rng.uniform(.1, 1, 3), rng.uniform(0, 360)]
        climate = v1.Climate(31, 2)
        assert_allclose(model.budgets(x, climate), 0, atol=1e-10)
        for k in range(8):
            y = x.copy(); y[k] = 0
            assert model.rhs(0, y, climate)[k] >= -1e-12
        for k in range(3):
            y = x.copy(); y[k+4] = v1.capacity(c, model.parameters)[k]
            dx = model.rhs(0, y, climate)
            assert model.parameters.omega[k]*dx[k]-dx[k+4] >= -1e-12
        y = x.copy(); y[7] = model.parameters.soil_sat
        assert model.rhs(0, y, climate)[7] <= 1e-12
    for ref in process_provenance(options)["sources"].values():
        assert ref["commit"] == COMMIT and ref["units"] and ref["native_order"] and ref["adaptation"]


def test_default_matches_v1_without_changing_equations():
    x = v1.equilibrium().state
    for climate in (v1.Climate(), v1.Climate(31, 2)):
        assert_allclose(Model().rhs(0, x, climate), v1.rhs(0, x, climate), atol=1e-14)


def test_guard_rejects_negative_seed_and_audits_internal_steps():
    model = Model()
    x = v1.equilibrium().state
    bad = x.copy(); bad[0] = -1e-14
    with pytest.raises(ValueError, match="nonnegative"):
        integrate(model, bad, [0, 100], v1.Climate())
    result = integrate(model, x, [0, 100], v1.Climate(40, 0))
    assert result.audit["inspected_points"] > result.audit["accepted_steps"] > 2
    assert min(result.audit["minima"].values()) >= 0
    assert result.audit["mass_added_by_clipping"] == 0


def test_guard_retries_and_stops_if_rhs_drives_a_pool_negative():
    class BrokenModel:
        parameters = v1.Parameters()

        def rhs(self, time, state, climate):
            return np.array([-1., 0., 0., 0., 0., 0., 0., 0.])

    seed = np.array([1., 40., 8., 80., 0., 0., 0., 120.])
    with pytest.raises(RuntimeError, match="negative pool"):
        integrate(BrokenModel(), seed, [0, 3], v1.Climate(), Solver(retries=1))


def test_spinup_converges_from_nonequilibrium_and_records_failure():
    p = replace(v1.Parameters(), gpp_max=.8, turnover=(.02, .0008, .006),
                respiration=(.006, .0003, .0015), soil_respiration=.002)
    model = Model(parameters=p)
    seed = np.array([3., 40., 8., 80., .4, 2., .8, 120.])
    fail = spinup(model, seed, v1.Climate(), SpinupSettings(max_years=1, chunk_years=1))
    assert not fail["converged"]
    result = spinup(model, seed, v1.Climate(), SpinupSettings(max_years=100, chunk_years=10))
    assert result["converged"]
    assert_allclose(result["states"][0], seed)
    assert np.linalg.norm(result["states"][-1]-seed) > 1
    eq = model.equilibrium(v1.Climate(), result["states"][-1])
    assert np.max(abs(result["states"][-1]-eq.state)/np.maximum(1, abs(eq.state))) < 1e-5
    assert all(c["passes"] for c in result["checks"][-3:])


def test_unknown_choice_and_unphysical_settings_fail():
    with pytest.raises(ValueError):
        ProcessOptions(photo_canopy="misspelled")
    with pytest.raises(ValueError):
        ProcessOptions(supply_time_days=0)
    with pytest.raises(ValueError):
        Solver(max_step=30)
