from datetime import date

import pytest

from control_carbon.visit_workbook import TKYModel, SOURCE_DEFAULT
from control_carbon.visit_workbook_data import (
    DATA_DEFAULT, corrected_temperature, load_baseline, read_native_snapshot,
    read_rows, summer_forcing,
)

NATIVE = SOURCE_DEFAULT.parent.parent/'Control_carbon_model/artifacts/visit_native_20000715/run_pinned2013/TKY_native_snapshot_20000715.tsv'


def test_temperature_and_hand_computable_forcing():
    assert corrected_temperature(273.15) == pytest.approx(-8.3206)
    row = [2000, 196, 273.15, 0, 0, 0, 1/86400, 0, 50,
           273.15, 273.15, 273.15, 0, 0, 0, 3, 4, 0, 100000]
    v = summer_forcing(row, date(2000, 7, 15), 0, 0, 3.141592653)
    assert v['z_doy'] == 196
    assert v['z_rain'] == pytest.approx(.9626)
    assert v['z_wind'] == 5
    assert v['z_cloud'] == .5
    assert v['z_vp'] == 0
    assert v['z_vpd'] == pytest.approx(6.1078)


@pytest.mark.skipif(not DATA_DEFAULT.exists() or not SOURCE_DEFAULT.exists(), reason='local VISIT data required')
def test_real_snapshot_does_not_invent_stocks_or_change_parameters():
    model = TKYModel()
    parameters = {k: n.value for k, n in model.graph.nodes.items() if k.startswith('p_')}
    baseline = load_baseline(model)
    assert (baseline['before'], baseline['after']) == (70, 46)
    assert all(model.graph.nodes[s.name].value is None for s in model.states)
    assert parameters == {k: n.value for k, n in model.graph.nodes.items() if k.startswith('p_')}
    assert model.graph.nodes['z_rain'].value == pytest.approx(7.73468352)
    assert model.graph.nodes['z_depo_no3'].value == pytest.approx(13.076)
    assert model.graph.nodes['h_t_flush'].value == 54
    for name in ('z_tl','z_th','h_t_gc','h_g_shed','z_ngas2_aa','z_curry_undefined'):
        assert model.graph.nodes[name].value is None
    assert len(baseline['restart']) == 21
    assert all(':' in r[3] for r in baseline['imports'])
    assert all(not isinstance(v, str) for k, v in model.graph.evaluate().items() if k.startswith('q_'))


def test_numeric_parser_accepts_cr_and_rejects_duplicate(tmp_path):
    # Test fixtures may use write_bytes; production importer is read-only.
    path = tmp_path/'data.txt'
    path.write_bytes(b'2000 0 1\r2000 1 2\r')
    assert read_rows(path, 3)[2000, 1] == (2, (2000., 1., 2.))
    path.write_bytes(b'2000 0 1\n2000 0 2\n')
    with pytest.raises(ValueError, match='duplicate'):
        read_rows(path, 3)


@pytest.mark.skipif(not NATIVE.exists() or not DATA_DEFAULT.exists(), reason='native snapshot required')
def test_native_snapshot_populates_all_defined_inputs_and_matches_day_map():
    model=TKYModel()
    baseline=load_baseline(model,native_snapshot=NATIVE)
    assert (baseline['before'],baseline['after'])==(70,2)
    assert {n.name for n in model.graph.nodes.values()
            if n.expression is None and n.value is None}=={'z_ngas2_aa','z_curry_undefined'}
    assert model.graph.nodes['z_vp'].value==pytest.approx(.35322809525173571)
    assert model.graph.nodes['x_t_nsc'].value==pytest.approx(5.0636156963329864)
    assert model.graph.nodes['h_t_gc'].value==pytest.approx(389.09685500300714)
    assert baseline['native_max_scaled_error']<1e-12
    assert len(baseline['native_validation'])==37
    rows=read_native_snapshot(NATIVE,'2000-07-15')
    assert rows['after_loct']['x_w_sw']==pytest.approx(104.4413838408682)
