"""Read-only, provenance-preserving import of a TKY summer forcing snapshot.

The workbook's pinned equations/parameters remain authoritative. Historical
output stocks are references, never silently substituted for summer stocks.
"""
from datetime import date, timedelta
import csv
import hashlib
import math
from pathlib import Path

from .visit_workbook import Graph

DATA_DEFAULT = Path('/mnt/d/VISIT/honban/point/visitb_v1.1.0')
RAIN_FACTORS = (2.3636, 1.9079, 1.3716, 1.3368, 1.1556, 1.1687,
                .9626, .9846, 2.0962, 1.6507, 1.7851, 1.8242)


def read_rows(path, width, keys=2):
    """Accept LF/CRLF/CR records; reject malformed/duplicate records."""
    result = {}
    for line, text in enumerate(path.read_text().splitlines(), 1):
        if not text.strip():
            continue
        row = tuple(map(float, text.split()))
        if len(row) != width or not all(math.isfinite(v) for v in row):
            raise ValueError(f'{path}:{line}: invalid numeric record')
        key = tuple(int(v) for v in row[:keys])
        if tuple(row[:keys]) != key or key in result:
            raise ValueError(f'{path}:{line}: invalid/duplicate key')
        result[key] = (line, row)
    return result


def corrected_temperature(kelvin):
    return (kelvin - 273.15) * 1.1079 - 8.3206


def summer_forcing(row, reference, latitude, altitude, pi):
    """Native TKY preprocessing; soil temperatures require unknown snow."""
    ta = corrected_temperature(row[2])
    q = row[5]
    raw_t = row[2]-273.15
    # Pinned reader tests Celsius; supplied visitb tests Kelvin instead.
    exponent = 7.5*raw_t/(237.3+raw_t) if raw_t >= 0 else 9.5*raw_t/(265.3+raw_t)
    year = reference.year
    doy = (reference-date(year, 1, 1)).days
    co2 = (1904299 - 3322.4242*year + 1.6541596*year**2
           + 1.3362655*year**3/10000 - 3.0828809*year**4/10000000
           + 6.2121261*year**5/100000000000 + 1.6*latitude/85
           + math.exp(.04*latitude)/2*math.sin(doy/365*2*pi))
    return dict(z_doy=float(doy), z_ta=ta, z_ts=corrected_temperature(row[9]),
                z_rain=row[6]*86400*RAIN_FACTORS[reference.month-1],
                z_cloud=row[8]/100, z_wind=math.hypot(row[15], row[16]),
                # location_proc.c has a missing "==0" after one strcmp;
                # for TKY this makes the CEAMIP branch execute.
                z_vp=(6.1078*10**(7.5*ta/(237.3+ta) if ta >= 0 else 9.5*ta/(265.3+ta)))*q,
                z_vpd=max(0, 6.1078*10**exponent-row[18]*.01*q/(.622+.378*q)),
                z_co2=co2)


def reconstruct_history(model, climate, reference):
    """Start unknown; retain only histories resolved by actual reset branches."""
    graph = Graph()
    for n in model.graph.nodes.values():
        if n.expression is None or n.name.startswith('q_'):
            expression = n.expression
            if n.name == 'q_t_gdd':
                # Equivalent nested IF resolves the cold reset even when CDD
                # is still unknown; do not require an irrelevant AND operand.
                expression = 'IF(q_t_gdd1<580,q_t_gdd1,IF(q_t_cdd1<(-190),0,q_t_gdd1))'
            graph.add(n.name, n.label, n.unit, n.sheet, expression, n.value)
    history = {n: None for n in graph.nodes if n.startswith('h_') and not n.endswith('_gc')}
    start = date(reference.year-1, 1, 1)
    for offset in range((reference-start).days):
        current = start+timedelta(days=offset)
        day = (current-date(current.year, 1, 1)).days
        _, row = climate[current.year, day]
        values = graph.evaluate(dict(history, z_doy=day,
                                     z_ta=corrected_temperature(row[2]),
                                     z_ts=corrected_temperature(row[9])))
        history = {n: values['q_'+n[2:]] for n in history}
    return {n: v for n, v in history.items() if isinstance(v, (int, float))}


def read_native_snapshot(path, reference):
    """Read the three source-order checkpoints emitted by diagnostic ansis.c."""
    reference = date.fromisoformat(str(reference))
    doy = (reference-date(reference.year, 1, 1)).days
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        rows={}
        for record in reader:
            phase=record.pop('phase')
            year=int(record.pop('year')); row_doy=int(record.pop('doy'))
            if year!=reference.year or row_doy!=doy or phase in rows:
                raise ValueError(f'Unexpected native snapshot row: {year}/{row_doy}/{phase}')
            values={key:float(value) for key,value in record.items()}
            if not all(math.isfinite(value) for value in values.values()):
                raise ValueError(f'Non-finite native snapshot value: {phase}')
            rows[phase]=values
    if set(rows)!={'start','after_loct','end'}:
        raise ValueError(f'Native snapshot phases are incomplete: {set(rows)}')
    return rows


def apply_native_snapshot(model, metadata, path, reference):
    """Populate X/h and post-location forcing, then compare the native day map."""
    path=Path(path)
    rows=read_native_snapshot(path,reference)
    nodes=model.graph.nodes
    adopted={}
    for name,node in nodes.items():
        phase=None
        if name.startswith(('x_','h_')) and name in rows['start']:
            phase='start'
        elif name.startswith('z_') and name in rows['after_loct'] and name!='z_rain_annual':
            phase='after_loct'
        if phase is None:
            continue
        node.value=rows[phase][name]
        node.source=f'{path}:{phase}'
        node.note=(f'{reference} native C診断実行の{phase}。startは当日f_loct_proc前、'
                   'after_loctは水文・当日環境設定後。更新順序を混同しない。')
        adopted[name]=(name,node.value,node.unit,f'{path}:{phase}',node.note)
    metadata['imports']=[row for row in metadata['imports'] if row[0] not in adopted]
    metadata['imports'].extend(adopted.values())
    order={name:i for i,name in enumerate(nodes)}
    metadata['imports'].sort(key=lambda row:order[row[0]])
    metadata['references']=[row for row in metadata['references'] if '積雪補正前' not in row[0]]
    metadata['native_rows']=[]
    for name in rows['start']:
        unit=nodes[name].unit if name in nodes else '原典内部単位'
        adopted_at=('X/hはstart、気象zはafter_loctを採用'
                    if name in adopted else '相間比較用')
        metadata['native_rows'].append((name,rows['start'][name],rows['after_loct'][name],
                                        rows['end'][name],unit,adopted_at))
    values=model.graph.evaluate()
    comparison=[]
    for state in model.states:
        predicted=values['m_end_'+state.name[2:]]
        actual=rows['end'][state.name]
        error=predicted-actual
        scaled=abs(error)/max(1,abs(actual))
        comparison.append((state.name,rows['start'][state.name],actual,predicted,error,scaled,state.unit))
    metadata['native_validation']=comparison
    metadata['native_max_scaled_error']=max(row[5] for row in comparison)
    metadata['native_snapshot']=str(path)
    metadata['manifest'].append((str(path),hashlib.sha256(path.read_bytes()).hexdigest(),
                                 '固定ソースをFLUX_SCHEME=0で実行した3時点診断出力'))
    for filename,label in [('diagnostic_changes.patch','nativeソースへの全変更'),
                           ('native_manifest.sha256','実行物・入力・出力のSHA256一覧')]:
        audit=path.parent/filename
        if audit.is_file():
            metadata['manifest'].append((str(audit),hashlib.sha256(audit.read_bytes()).hexdigest(),label))
    metadata['after']=sum(n.expression is None and n.value is None for n in nodes.values())
    model.issues=[issue for issue in model.issues if not issue[0].startswith('DATA')]
    model.issue('NATIVE01','native実行とソース補修',
                '2013固定ソースを使用。FLUX_SCHEME=0と診断出力のみ変更。16バイオーム表はparameter_VISITc_16.xlsxの出力用sheetから復元し、コードが要求するparameter_S1b.txt名で使用。',
                '08・09・98とdiagnostic_changes.patchを参照')
    model.issue('NATIVE02','TKY水蒸気圧分岐の比較漏れ',
                'location_proc.cのstrcmp(grid->site_id,"CEAMIP_TMK")に==0がなく、TKYでもvp=vps×specific_humidityが実行される。native値0.3532 hPaを採用。',
                '意図式に修正する場合は別のモデル版として比較')
    model.issue('NATIVE03','実行時保水容量',
                'Configのfc30=64.21 mmはSaxton式で104.810272811 mmへ上書き。全層808.87 mmから差し引いた704.059727189 mmが深層バケット。ブックを実行値へ修正。',
                'init_site.c::f_init_site, soil_physics.c::f_soil_saxton')
    return rows


def load_baseline(model, root=DATA_DEFAULT, reference='2000-07-15', native_snapshot=None):
    root = Path(root)
    reference = date.fromisoformat(reference)
    if not 1980 <= reference.year <= 2003:
        raise ValueError('Importer requires direct daily N deposition (1980–2003)')
    inp, out = root/'INPUT', root/'OUTPUT'
    climate = read_rows(inp/'ext_28-73.dat', 19)
    deposition = read_rows(inp/'ndepo_1980-2003_tky.txt', 6)
    ghg = read_rows(inp/'AtmGHG_timeseries.txt', 13, 1)
    doy = (reference-date(reference.year, 1, 1)).days
    line, row = climate[reference.year, doy]
    nline, nr = deposition[reference.year, doy]
    gline, gr = ghg[reference.year,]
    nodes = model.graph.nodes
    before = sum(n.expression is None and n.value is None for n in nodes.values())
    values = summer_forcing(row, reference, nodes['p_lat'].value,
                            nodes['p_alt'].value, nodes['c_pi'].value)
    values.update(z_depo_nh4=(nr[2]+nr[4])*1000,
                  z_depo_no3=(nr[3]+nr[5])*1000, z_ch4=gr[5]/1000)
    # Pinned visit_local uses 1980–2009, unlike supplied visitb (1980–1999).
    values['z_tsoil_mean'] = sum(climate[y, d][1][10]-273.15
                                for y in range(1980, 2010) for d in range(365))/(30*365)
    days = (date(reference.year+1, 1, 1)-date(reference.year, 1, 1)).days
    values['z_rain_annual'] = sum(climate[reference.year, d][1][6]*86400
                                *RAIN_FACTORS[(date(reference.year, 1, 1)+timedelta(days=d)).month-1]
                                for d in range(days))
    histories = reconstruct_history(model, climate, reference)
    values.update(histories)
    metadata = {'date': str(reference), 'before': before, 'imports': [], 'references': [], 'manifest': []}
    rules = {
        'z_doy': '暦日→0始まり年内日', 'z_ta': '(列3−273.15)×1.1079−8.3206',
        'z_ts': '(列10−273.15)×1.1079−8.3206',
        'z_rain': '列7×86400×月補正係数（7月0.9626）',
        'z_cloud': '列9/100', 'z_wind': 'sqrt(列16²+列17²)',
        'z_vp': '飽和水蒸気圧×列6。location_proc.cのCEAMIP_TMKに対するstrcmp比較漏れによりTKYでもこの分岐。意図式ではなくnative実行式',
        'z_vpd': 'max(0,補正前気温の飽和水蒸気圧−列19×0.01×列6/(0.622+0.378×列6))',
        'z_co2': 'atm_co2.c::atmco2_trend 歴史期間の多項式＋緯度・季節項。GHG表CO2ではない',
        'z_ch4': '列6 ppbv /1000→ppmv; initialize.c:142–154',
        'z_depo_nh4': '(列3+列5)×1000; initialize.c:199–208',
        'z_depo_no3': '(列4+列6)×1000; initialize.c:199–208',
        'z_tsoil_mean': '1980–2009各年のdoy0–364、列11−273.15の平均。旧visit_local/location_init.c:76–126を採用',
        'z_rain_annual': '当年全日列7×86400×各月補正係数の和。旧visit_local/location_proc.c:276–306はサイト補正後にprate_annへ加算（提供visitbは補正前）。年次診断専用',
    }
    for name, value in values.items():
        source = f'{inp}/ext_28-73.dat:{line}'
        if name.startswith('z_depo_'):
            source = f'{inp}/ndepo_1980-2003_tky.txt:{nline}'
        elif name == 'z_ch4':
            source = f'{inp}/AtmGHG_timeseries.txt:{gline}'
        elif name in histories:
            source = f'{inp}/ext_28-73.dat:{climate[reference.year-1, 0][0]}–{climate[reference.year, doy-1][0]}'
        elif name == 'z_tsoil_mean':
            source = f'{inp}/ext_28-73.dat:1980–2009、各年先頭365日'
        elif name == 'z_rain_annual':
            source = f'{inp}/ext_28-73.dat:{climate[reference.year, 0][0]}–{climate[reference.year, days-1][0]}'
        rule = rules.get(name, '前年1月1日の履歴は未知として旧ブックのq_*を前日まで逐次適用。リセットにより一意に定まった値のみ採用。実出力の再現ではない')
        node = nodes[name]
        assert node.expression is None and node.value is None
        node.value = float(value)
        node.source += ' | データ: '+source
        node.note = f'{reference}用の固定スナップショット。{rule}。日付・Config変更時は再生成。'
        metadata['imports'].append((name, node.value, node.unit, source, rule))
    for name, index in [('z_tl', 10), ('z_th', 11)]:
        nodes[name].note = '積雪水当量が未提供のため補正後地温を確定できない。06の補正前地温を参照。'
        metadata['references'].append((name+'（積雪補正前）', corrected_temperature(row[index]), 'degC',
                                       f'ext_28-73.dat:{line} 列{index+1}', '積雪S>0.2ならa=S/(3+S)、Tl=(1−a)Tl、Th=2a+(1−a)Th。S不明につき未採用'))
    restart = list(map(float, (out/'TKY_restart.txt').read_text().split()))
    mapping = ['x_w_snow','x_w_sw','x_w_dw'] + [f'x_{p}_{o}' for p in ('t','g','v') for o in ('fol','stm','rot')] + ['x_s_'+p for p in ('tf','tc','tr','gf','gc','gr','ha','hi','hp')]
    if len(restart) != len(mapping):
        raise ValueError('Unexpected restart layout')
    metadata['restart'] = [(name, value, nodes[name].unit, f'{out}/TKY_restart.txt:1 列{i}',
                            'spinup終了時。夏の状態には未採用。INPUT/spinup.c:234–256')
                           for i, (name, value) in enumerate(zip(mapping, restart), 1)]
    daily = read_rows(out/'TKY_110127_daily.txt', 5)
    dline, dr = daily[reference.year, doy]
    for index, label in enumerate(('GPP', '生態系呼吸', 'NEP'), 2):
        metadata['references'].append((label, dr[index], 'gC m-2 day-1',
                                       f'{out}/TKY_110127_daily.txt:{dline} 列{index+1}',
                                       'INPUT/ansis.c::output_ansis_daily。別版の出力参考値。現在ブックの計算結果とは比較未検証'))
    for path in sorted(inp.iterdir()) + sorted(out.iterdir()):
        if path.is_file() and (path.suffix in ('.c', '.h') or path.name in
                ('ext_28-73.dat','ndepo_1980-2003_tky.txt','AtmGHG_timeseries.txt',
                 'setting.txt','site_TKY.txt','parameter_TKY.txt') or path.parent == out):
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            old = model.source/path.name
            comparison = ('同一' if old.is_file() and hashlib.sha256(old.read_bytes()).hexdigest() == digest
                          else '別ファイル・式やパラメータの置換には使用しない')
            metadata['manifest'].append((str(path), digest, comparison))
    metadata['after'] = before-len(values)
    model.baseline = metadata
    model.issues = [issue for issue in model.issues if issue[0] not in ('INPUT01','INPUT03','INPUT04','INPUT05')]
    model.issue('DATA01', '夏の状態・履歴の不足', '日次出力はGPP/呼吸/NEPのみ。restart21量はspinup終了時で夏の状態ではない。未採用。', '基準日開始時の全37状態と前日gc等。99の個別空欄参照')
    model.issue('DATA02', '異なるソース版', '式・パラメータは旧visit_local/Configを維持。visitbデータのみ利用。長期地温平均は旧版1980–2009を採用（visitbは1980–1999）。', '06–08参照。visitb出力と本ブックの全モデル同値は主張しない')
    model.issue('DATA03', '気象補正と版の違い', '旧版は年降水量をTKY補正後に積算、提供visitbは補正前。飽差の水/氷分岐は旧版が℃、visitbがK。旧版規約を採用。vpと飽差は別補正経路のため整合性に注意。', '旧location_proc.c:276–306, init_site.c:353–362。vpが補正後の飽和水蒸気圧を超えるが独自補正は加えない。')
    if native_snapshot:
        apply_native_snapshot(model,metadata,native_snapshot,reference)
    return metadata
