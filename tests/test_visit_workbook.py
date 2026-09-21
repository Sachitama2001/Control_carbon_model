"""Source-order, matrix, missing-value and OOXML regression tests."""
import ast
from pathlib import Path
import shutil

import numpy as np
import pytest

openpyxl=pytest.importorskip('openpyxl')
from control_carbon.visit_workbook import Graph, TKYModel, SOURCE_DEFAULT
from control_carbon.visit_workbook_export import Exporter, prepare_matrix
from control_carbon.visit_workbook_validation import artificial_inputs, validate


@pytest.fixture
def model():
    if not (SOURCE_DEFAULT/'INPUT/Config_TKY.xlsx').is_file():
        pytest.skip('User-provided TKY source/config unavailable')
    return TKYModel()


def test_missing_is_not_zero_and_if_is_lazy():
    g=Graph();g.add('input_x','','','',value=None)
    g.add('twice','','','',expression='2*input_x')
    g.add('inactive','','','',expression='IF(1==1,7,input_x/0)')
    assert g.evaluate()=={'input_x':None,'twice':None,'inactive':7}
    assert 'ISNUMBER(input_x)' in g.excel('twice')
    assert g.evaluate({'input_x':3})['twice']==6
    with pytest.raises(ValueError):g.evaluate({'twice':5})
    with pytest.raises(ValueError):g.add('future','','','',expression='missing_name+1')


def test_config_defaults_and_unknown_states(model):
    assert len(model.states)==37
    assert len(model.fluxes)==132
    assert model.graph.nodes['p_fc30'].value==pytest.approx(104.81027281148704)
    assert model.graph.nodes['p_fc'].value==pytest.approx(704.059727188513)
    assert model.graph.nodes['p_t_lc0'].value==.0000078
    values=model.graph.evaluate()
    assert values['z_rain'] is None
    assert values['x_t_fol'] is None
    assert values['m_rhs_t_fol'] is None
    assert values['p_t_pmax']==13
    assert all(values[s.name] is None for s in model.states)
    for name,node in model.graph.nodes.items():
        if node.source.startswith('Config_TKY.xlsx!'):
            _,sheet,addr=node.source.split('!')
            assert node.value==model.config_book[sheet][addr].value


@pytest.mark.parametrize('case',['summer','flush','autumn','winter','dry','zero_upper','zero_leaf'])
def test_factorization_budgets_and_stage_endpoints(model,case):
    prepare_matrix(model)
    v=model.graph.evaluate(artificial_inputs(model,case))
    x,b,mu,a,xi,k=model.matrix_values(v)
    rhs=np.array([v['m_rhs_'+s.name[2:]] for s in model.states])
    np.testing.assert_allclose(b@mu+a@xi@k@x,rhs,rtol=1e-12,atol=1e-10)
    for s in model.states:
        assert abs(v['mx_error_'+s.name[2:]])<1e-9
        if s.name.startswith('x_s_'):endpoint='s_end_'+s.name[4:]
        elif s.name.startswith('x_w_'):endpoint='w_'+s.name[4:]+'_end'
        elif s.group=='C':endpoint='c_'+s.name[2:3]+'_end_'+s.name[4:]
        elif s.group=='N':endpoint='n_end_'+s.name[2:]
        else:endpoint='d_'+s.name[2:]+'_end'
        assert v['m_end_'+s.name[2:]]==pytest.approx(v[endpoint],abs=1e-9)
    for group in ('c','n','w'):assert abs(v[f'b_{group}_error'])<1e-9
    assert v['b_w_vs_reported']==pytest.approx(-v['w_base']+v['w_swclip']+v['w_dwclip'],abs=1e-10)
    if case=='flush':assert v['b_c_vs_nep']==pytest.approx(-.1*v['c_t_emerge'],abs=1e-10)
    if case=='zero_upper':
        assert v['m_regular_w_sw']==0
        assert v['m_xi_w_sw']==0
        assert v['m_res_w_sw']!=0


def test_zero_baseline_k_routes_to_affine_residual(model):
    prepare_matrix(model)
    inputs=artificial_inputs(model)
    inputs['p_t_lf0']=0
    v=model.graph.evaluate(inputs)
    assert v['m_regular_t_fol']==0
    assert isinstance(v['mx_error_t_fol'],(int,float))
    assert abs(v['mx_error_t_fol'])<1e-10


def test_native_c_comparison(model):
    if shutil.which('gcc') is None:pytest.skip('gcc unavailable')
    result=validate(model)
    assert result['行列恒等式_最大絶対誤差']<1e-9
    assert result['C原典比較数']>=1300
    assert result['C原典比較_最大規格化誤差']<1e-10
    assert result['C原典未比較_定義域等']==[]


def test_xlsx_cached_formulas_no_external_links_and_missing_inputs(model,tmp_path):
    exporter=Exporter(model)
    dest=tmp_path/'matrix.xlsx'
    result=exporter.build(dest)
    assert result['formula_cells']>19000
    formulas=openpyxl.load_workbook(dest,data_only=False)
    values=openpyxl.load_workbook(dest,data_only=True)
    assert not formulas.calculation.iterate
    assert not formulas._external_links
    assert len(formulas['11_B'][5])==80
    for name in ('x_t_fol','z_rain','z_ngas2_aa'):
        sheet,addr=exporter.locations[name]
        assert values[sheet][addr].value is None
    sheet,addr=exporter.locations['p_t_lf0']
    assert values[sheet][addr].value==.000082
    sheet,addr=exporter.locations['m_rhs_t_fol']
    assert formulas[sheet][addr].data_type=='f'
    assert values[sheet][addr].value=='未設定/定義域'
    assert formulas.sheetnames[-1]=='99_不明と要確認'
    assert not any(row[2].value=='補助／分岐監査対象' for row in formulas['91_過程網羅'].iter_rows(min_row=5))


def test_formula_lengths_nesting_and_sources(model):
    prepare_matrix(model)
    def depth(t):
        children=list(ast.iter_child_nodes(t))
        return int(isinstance(t,ast.Call))+(max(map(depth,children)) if children else 0)
    for name,node in model.graph.nodes.items():
        if node.expression is None:continue
        assert len(model.graph.excel(name))<=8192
        # IFERROR wrapper and the guarded leaf references add at most 3 levels.
        assert depth(node.tree)+3<=64,(name,depth(node.tree))
        assert node.source
