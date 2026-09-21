"""Artificial verification cases; never default forest/meteorology values.

The C comparison executes native radiation, SLA, atmospheric relations,
hydrology, physiology, phenology and the daily C/N path. Pressure and cover
fractions are supplied as slice inputs. It does not validate file ingestion,
the full driver, unresolved diagnostics or all branches.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from .visit_workbook import TKYModel
from .visit_workbook_export import prepare_matrix


def artificial_inputs(model,case='summer'):
    x={'z_doy':180,'z_ta':20,'z_ts':22,'z_tl':17,'z_th':12,'z_rain':4,
       'z_cloud':.4,'z_wind':2,'z_vp':15,'z_vpd':8,'z_co2':400,'z_ch4':1.8,
       'z_depo_no3':20,'z_depo_nh4':25,'z_tsoil_mean':10,
       'z_rain_annual':1800,
       # Undefined diagnostic aa/Curry intentionally remain unset.
       'x_w_snow':0,'x_w_sw':45,'x_w_dw':500,'x_doc':5,'x_casa':.1}
    for p,scale in [('t',1),('g',.05),('v',.01)]:
        for organ,value in [('fol',2),('stm',60),('rot',20),('nsc',3)]:x[f'x_{p}_{organ}']=value*scale
        x[f'x_{p}_ncan']=30000*scale;x[f'x_{p}_nstr']=100000*scale
        for key,value in [('gc',120),('gdd',900),('cdd',0),('flush',60),('shed',0)]:x[f'h_{p}_{key}']=value
    for pool in ('tf','tc','tr','gf','gc','gr'):x['x_s_'+pool]=2
    for pool in ('ha','hi','hp'):x['x_s_'+pool]=40
    for pool,value in [('mic',300),('lit',50000),('hum',100000),('no3',2000),('nh4',3000)]:x['x_n_'+pool]=value
    if case=='flush':x.update(z_doy=120,z_ta=14,z_ts=15,h_t_gdd=280,h_t_flush=4,x_t_fol=.08)
    elif case=='autumn':x.update(z_doy=280,z_ta=7,z_ts=8,z_tl=10,z_th=9,h_t_cdd=-80,h_t_flush=145,h_t_shed=15)
    elif case=='winter':x.update(z_doy=20,z_ta=-8,z_ts=-5,z_tl=-2,z_th=-1,z_rain=3,x_w_snow=80,h_t_gdd=0,h_t_flush=0,h_t_shed=80,x_t_fol=.01)
    elif case=='dry':x.update(z_rain=0,x_w_sw=8,x_w_dw=120,z_vpd=20)
    elif case=='zero_upper':x.update(x_w_sw=0,z_rain=10)
    elif case=='zero_leaf':x.update(x_t_fol=0)
    elif case!='summer':raise ValueError(case)
    return x


def native_bridge(model,directory):
    """Generate a named-input executable without editing authoritative sources."""
    assigns=[];keys=[]
    def setv(lhs,key):
        keys.append(key);assigns.append(f'  if(scanf("%lf", &v)!=1) return 2; {lhs}=v;')
    static=['#include <stdio.h>','#include <math.h>','#include <string.h>',
            '#include "structure.h"','#include "prototype.h"',
            'short SA_PARA=0; short SA_PARA_EN=0; double SA_PARA_VAR=0;',
            'static struct Grid g; static struct Loct l; static struct Echar e;',
            'static struct Mass m; static struct Flux f;',
            'int main(void){double v; strcpy(g.site_id,"TKY"); g.veg_type=4; l.hour=24; l.time=200;',
            'l.adyear=2000; l.CO2y=1750; l.funder_c3=1; l.funder_c4=0;']
    grid={'lat':'lat','lon':'lon','alt':'topo','texture':'stexture','fc30':'fieldcap30','fc':'fieldcap','rootdepth':'soildpth','hydcond':'hyd_cond','density':'bulkdens','ph':'soil_ph','sand':'sand_frac','clay':'clay_frac'}
    for key,field in grid.items():setv('g.'+field,'p_'+key)
    loct={'z_doy':'doy','z_ta':'tmp_2m','z_ts':'tmp_sfc','z_tl':'tmp10_soil','z_th':'tmp200_soil','z_rain':'prate_sfc','z_cloud':'tcdc_clm','z_wind':'wnd_10m','z_co2':'aCO2','z_ch4':'atm_ch4_a1[0]','z_vp':'vp','z_vpd':'vpd','e_slope':'slope_vps','e_density':'air_dns','e_ra':'r_aero','e_dl':'daylen[l.doy]',
          'e_rn_t':'rn_tree','e_rn_g':'rn_c3','e_rn_v':'rn_c4','e_rn_s':'rn_ground','e_ppfd':'ppfd_h[l.hour]','e_ppfdc':'ppfd_c_h[l.hour]',
          'x_w_sw':'soilwtr_l','x_w_dw':'soilwtr_h','x_casa':'m_casa_pre','z_depo_no3':'depo_no3','z_depo_nh4':'depo_nh4'}
    for key,field in loct.items():setv('l.'+field,key)
    setv('g.stmp10cm_av','z_tsoil_mean')
    setv('l.prate_ann','z_rain_annual')
    setv('l.air_prsr','e_pressure')
    setv('l.comp_over1','c_comp1');setv('l.comp_over2','c_comp2')
    for p,field in [('t','tree'),('g','c3'),('v','c4'),('s','ground')]:setv('l.fcover_'+field,'e_cov_'+p)
    for key,field in [('snow','snwa'),('sw','sw30'),('dw','sww')]:setv('m.'+field,'x_w_'+key)
    output={}
    for p,field in [('t','tree'),('g','c3'),('v','c4')]:
        for key,node in model.graph.nodes.items():
            if key.startswith('p_'+p+'_'):setv('e.'+field+'.'+key[len('p_'+p+'_'):],key)
        for key,target in [('gdd','gdd'),('cdd','cdd'),('flush','day_flush'),('shed','day_shed')]:setv(f'e.{field}.{target}',f'h_{p}_{key}')
        setv(f'e.{field}.gc',f'h_{p}_gc')
        for key,target in [('fol','fol'),('stm','stm'),('rot','rot'),('nsc','nsch_storage'),('ncan','n_canopy'),('nstr','n_storage')]:setv(f'm.{field}.{target}',f'x_{p}_{key}')
        for key,target in [('sla','sla'),('ek','eK')]:setv(f'e.{field}.{target}',f'e_{p}_{key}')
        setv(f'm.{field}.lai',f'e_{p}_lai');setv(f'e.{field}.ppfd_t',f'pht_{p}_ppfd')
        for key,target in [('fol','fol'),('stm','stm'),('rot','rot'),('nsc','nsch_storage')]:output[f'c_{p}_end_{key}']=f'm.{field}.{target}'
        for key,target in [('ncan','n_canopy'),('nstr','n_storage')]:output[f'n_end_{p}_{key}']=f'm.{field}.{target}'
        for key,target in [('gpp','gpp'),('npp','npp'),('lf','lf'),('lc','lc'),('lr','lr'),('rmf','rfm'),('rmc','rcm'),('rmr','rrm'),('rgf','rfg'),('rgc','rcg'),('rgr','rrg'),('tpf','tpf'),('tpc','tpc'),('tpr','tpr')]:output[f'c_{p}_{key}']=f'f.{field}.{target}'
        for key,target in [('psat6','psat'),('gc','gc'),('optlai','opt_lai'),('lue6','lue'),('ci6','ci')]:output[f'pht_{p}_{key}']=f'e.{field}.{target}'
        for key,target in [('season','season'),('gdd','gdd'),('cdd','cdd'),('flush','day_flush'),('shed','day_shed')]:output[f'q_{p}_{key}']=f'e.{field}.{target}'
        for key in ('mass','photo'):output[f'd_{p}_ch4{key}']=f'f.{field}.emit_ch4_kirschbaum_{key}'
    for key,node in model.graph.nodes.items():
        if key.startswith('p_s_'):setv('e.soil.'+key[4:],key)
    for pool in ('tf','tc','tr','gf','gc','gr','ha','hi','hp'):
        field='msl_'+pool[-1] if pool[0]=='h' else 'ltr_'+pool
        setv('m.soil.'+field,'x_s_'+pool);output['s_end_'+pool]='m.soil.'+field
    for key,field in [('mic','n_mcrb'),('lit','n_lttr'),('hum','n_hums'),('no3','n_no3'),('nh4','n_nh4')]:
        setv('m.soil.'+field,'x_n_'+key);output['n_end_n_'+key]='m.soil.'+field
    setv('m.soil.doc','x_doc')
    lines=static+assigns+[
        'l.atm_ch4_a1[0]*=1000; f_soil_saxton(&g);',
        'f_sla_change(&g,&l,&e,&m);',
        'm.tree.lai=lai_mass(&m.tree,&e.tree); m.c3.lai=lai_mass(&m.c3,&e.c3); m.c4.lai=lai_mass(&m.c4,&e.c4);',
        'l.soldec[l.doy]=f_soldec(&g,&l); l.daylen[l.doy]=f_daylen(&g,&l);',
        'for(l.hour=0;l.hour<48;l.hour++){ l.hangle=-180+7.5*l.hour; l.solhgt_h[l.doy][l.hour]=f_solhgt(&g,&l); l.toprad_h[l.doy][l.hour]=f_toprad(&g,&l); f_ppfd(&g,&l); } l.hour=24;',
        'e.tree.eK=irr_attn(&g,&l,&e.tree); e.c3.eK=irr_attn(&g,&l,&e.c3); e.c4.eK=irr_attn(&g,&l,&e.c4);',
        'l.snow_acc=m.snwa; e.soil.albedo=albedo_soil(&l,&e.soil); f_net_rad(&g,&l,&e,&m);',
        'l.air_dns=f_airdens(&g,&l); l.slope_vps=f_slope_vps(&g,&l); l.r_aero=f_r_aero(&g,&l);',
        'f_hydrology(&g,&l,&e,&m);',
        'l.soilwtr_l=m.sw30; l.soilwtr_h=m.sww;',
        'l.soilappr_l=fmax(0,fmin(1,(g.fieldcap30-m.sw30)/g.fieldcap30));',
        'l.soilappr_w=fmax(0,fmin(1,(g.fieldcap-m.sww)/g.fieldcap));',
        'l.wfps=fmax(0.05,m.sw30/300/(1-g.bulkdens/2.65));',
        'l.pot_total_l=-0.05-0.478*pow(fmax(0.2,m.sw30)/g.fieldcap30,-5.39);',
        'l.pot_total_h=-1-0.478*pow(fmax(0.2,m.sww)/g.fieldcap,-5.39);',
        'f_ecophysiology(&g,&l,&e.tree,&m.tree);',
        'f_ecophysiology(&g,&l,&e.c3,&m.c3);',
        'f_ecophysiology(&g,&l,&e.c4,&m.c4);',
        'daily_scheme(&g,&l,&e,&m,&f); f_erosion_rusle(&g,&l,&m,&f);']
    output.update({'w_snow_end':'m.snwa','w_sw_end':'m.sw30','w_dw_end':'m.sww','w_aet':'l.aet','w_ro1':'l.ro1','w_ro2':'l.ro2','w_evap':'l.evpr','w_pet':'l.pet','d_doc_end':'m.soil.doc','d_casa_end':'l.m_casa','s_rh':'f.soil.hr'})
    output.update({'e_dl':'l.daylen[l.doy]','e_dec':'l.soldec[l.doy]','e_t_sla':'e.tree.sla','e_rn_t':'l.rn_tree','e_rn_g':'l.rn_c3','e_rn_s':'l.rn_ground','e_ppfd':'l.ppfd_h[24]','e_ppfdc':'l.ppfd_c_h[24]','e_vps':'f_vap_pre_sat(&g,&l)','e_density':'l.air_dns','e_slope':'l.slope_vps','e_ra':'l.r_aero','d_ero_c':'f.soil.erosion_carbon'})
    for name in ('isopr','monotrp','methanl','acetone','actaldhd','frmardhd','formacd','acetacd','co'):output['d_voc_'+name]='f.voc_'+name+'_g97'
    for key,field in [('d_ridg_ch4','ch4oxy_ridg'),('d_dg_ch4','ch4oxy_delgrosso'),('d_casa_ch4','ch4oxy_casa'),('d_casa_no','d_no_casa'),('d_casa_n2o','d_n2o_casa'),('d_casa_n2','d_n2_casa')]:output[key]='f.soil.'+field
    for h in range(48):output[f'e_sw{h:02d}']=f'l.sfcrad_h[{h}]'
    for key,field in output.items():lines.append(f'printf("RESULT {key} %.17g\\n",(double)({field}));')
    lines.append('return 0;}')
    path=Path(directory)/'visit_workbook_bridge.c';path.write_text('\n'.join(lines))
    exe=path.with_suffix('')
    sources=['radiation','leaf_aging','erosion','hydro_balance','hydro_flows','ecophysiology','photosynthesis','daily_scheme','phenology','agriculture','plant_proc','litterfall','respiration','allocation','soil_proc','decomposition','doc','soil_physics','n_flux','n_budget','n2o_emit','vocemit','ch4_emit','ch4_oxy']
    cmd=['gcc','-O2','-fcommon','-ffunction-sections','-fdata-sections','-I',str(model.source),str(path)]+[str(model.source/(s+'.c')) for s in sources]+['-Wl,--gc-sections','-lm','-o',str(exe)]
    p=subprocess.run(cmd,text=True,capture_output=True)
    if p.returncode:raise RuntimeError(p.stderr)
    return exe,keys,output


def validate(model=None):
    m=model or TKYModel();prepare_matrix(m)
    matrix_errors=[];native_errors=[];comparisons=0;unresolved=[]
    with tempfile.TemporaryDirectory(prefix='visit-xlsx-native-') as directory:
        exe,keys,out=native_bridge(m,directory)
        for case in ('summer','flush','autumn','winter','dry','zero_upper','zero_leaf'):
            v=m.graph.evaluate(artificial_inputs(m,case))
            x,b,u,a,xi,k=m.matrix_values(v)
            rhs=np.array([v['m_rhs_'+s.name[2:]] for s in m.states],dtype=float)
            matrix_errors.append(float(np.max(np.abs(b@u+a@xi@k@x-rhs))))
            for s in m.states:
                assert isinstance(v['mx_error_'+s.name[2:]],(int,float)),(case,s.name,v['mx_error_'+s.name[2:]])
            text=' '.join(str(v[key]) for key in keys)
            result=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True)
            differences=[]
            for line in result.stdout.splitlines():
                if not line.startswith('RESULT '):continue
                _,key,raw=line.split();actual=float(raw);expected=v[key]
                if not isinstance(expected,(int,float)) or not math.isfinite(actual):
                    unresolved.append((case,key,expected,raw));continue
                comparisons+=1
                err=abs(actual-expected)/(1+abs(actual))
                native_errors.append(err)
                if err>1e-9:differences.append((key,expected,actual,err))
            if differences:print(case,json.dumps(differences,ensure_ascii=False))
    return {'人工ケース数':7,'行列恒等式_最大絶対誤差':max(matrix_errors),
            'C原典比較数':comparisons,'C原典比較_最大規格化誤差':max(native_errors),
            'C原典未比較_定義域等':unresolved}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=validate();print(json.dumps(result,ensure_ascii=False,indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
