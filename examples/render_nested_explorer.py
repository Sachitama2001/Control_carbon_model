"""Server-rendered equations with native HTML disclosure; no JavaScript."""
import ast
from html import escape


def render_nested(data):
    nodes, equations = data['nodes'], data['equations']
    used = set()
    formulas = {
        'x': 'x = (x炭素, x窒素, x水, x補助)ᵀ\nx翌日 = x + Δt fΔ(x, h, u)\nΔt = 1 day',
        'carbon': 'x炭素 = (x樹木, xC3林床, xC4仮想群, x土壌炭素)ᵀ',
        'plant': 'x各植生群 = (C葉, C幹, C根, C_NSC)ᵀ',
        'soil': 'x土壌炭素 = (Ctf, Ctc, Ctr, Cgf, Cgc, Cgr, Cha, Chi, Chp)ᵀ',
        'nitrogen': 'x窒素 = (x植物窒素, x土壌窒素)ᵀ',
        'plant_n': 'x植物窒素 = (N樹木葉群, N樹木貯蔵, NC3葉群, NC3貯蔵, NC4葉群, NC4貯蔵)ᵀ',
        'soil_n': 'x土壌窒素 = (N微生物, Nリター, N腐植, N硝酸態, Nアンモニア態)ᵀ',
        'water': 'x水 = (W積雪, W上層土壌, W深層土壌)ᵀ',
        'aux': 'x補助 = (DOC濃度, CASA水分記憶)ᵀ',
        'B': 'fΔ = B(x,h,u) μ(x,h,u) + A(x,h,u) ξ(x,h,u) K x',
        'allocation': 'EPP = GPP − R維持\n(tpf, tpc, tpr)ᵀ = (b葉, b幹, b根)ᵀ EPP\n（非作物・EPP > 0。成長呼吸とNSCへの振分けは後続処理）',
        'dynamic': 'b = b(L, L*, EPP, season)\nL = γ C葉\nL* = f_opt_lai(光・温度・生理係数)',
        'mu': 'μ = (41本の系外／有効入力, 37状態の残差)ᵀ\nBμ = 有効入力と残差の状態別の和',
        'A': 'Aⱼⱼ = −1\nAᵢⱼ = Tᵢⱼ / Lⱼ  （i ≠ j、通常列）\nLⱼ：状態jからの総控除、Tᵢⱼ：jからiへの移行',
        'xi': 'ξⱼⱼ = Lⱼ / (Kⱼⱼ xⱼ)  （通常列）\nxⱼ = 0 または Lⱼ = 0 または Kⱼⱼ = 0 の列は、有効残差としてBμ側へ移す',
        'K': 'K = diag(k₁, …, k₃₇)\n各状態の基準率は下の「状態別の係数」で確認',
    }
    children = {k: list(n['children']) for k,n in nodes.items()}
    children['tree'].append('nsc')
    children['nsc'] = ['nsc_flush', 'nsc_storage']
    children['soil_n'].append('microbes')
    children['microbes'] = []
    children['allocation'] = ['fixed', 'dynamic', 'nsc_storage']
    children['x'] += ['history', 'aggregation']

    def box(text):
        return '<pre class="equation">'+escape(text)+'</pre>'

    def refs(expr):
        if not expr:
            return []
        return sorted({n.id for n in ast.walk(ast.parse(expr, mode='eval'))
                       if isinstance(n, ast.Name) and n.id in equations})

    def collect(key):
        if key in used:
            return
        used.add(key)
        for ref in refs(equations[key]['expression']):
            collect(ref)

    def equation_detail(key, depth=1):
        collect(key)
        n=equations[key]
        expr=n['expression']
        definition=expr if expr is not None else str(n['value']) if n['value'] is not None else '外部入力／現在値（未設定）'
        body=box(key+' = '+definition)
        body+='<p class="unit">'+escape(n['unit'])+'</p>'
        if n['note']:
            body+='<p>'+escape(n['note'])+'</p>'
        body+='<p class="source">出典：'+escape(n['source'])+'</p>'
        if depth and expr:
            body+='<div class="nest">'+''.join(equation_detail(r,depth-1) for r in refs(expr))+'</div>'
        if not depth and expr:
            body+='<p class="source">さらに深い定義：'+ ' · '.join('<a href="#eq-'+escape(r)+'">'+escape(r)+'</a>' for r in refs(expr))+'</p>'
        return '<details class="process"><summary>'+escape(n['label'])+'<small>'+escape(key)+'</small></summary><div class="inside">'+body+'</div></details>'

    def source_detail(key):
        s=data['sources'][key]
        return '<details class="source-detail"><summary>原典：'+escape(s['symbol'])+'</summary><div class="inside"><p>'+escape(s['repository']+' / '+s['commit'])+'</p><p>'+escape(s['path']+' · '+s['level'])+'</p><p>更新順：'+escape(s['order'])+'</p><pre class="source-code">'+escape(s['excerpt'])+'</pre></div></details>'

    def tree(key, path=()):
        if key in path:
            return ''
        n=nodes[key]
        body=('<span id="open-'+escape(key)+'"></span>' if not path else '')+'<p class="description">'+escape(n['summary'])+'</p>'
        formula=n['formula'] or formulas.get(key,'')
        if key in ('tree','floor','c4'):
            formula='x = (C葉, C幹, C根, C_NSC)ᵀ'
        if formula:
            body+=box(formula)
        if n['unit']:
            body+='<p class="unit">'+escape(n['unit'])+'</p>'
        if 'state' in n:
            suffix=key[2:]
            body+=box(key+'翌日 = '+key+' + Δt × f_'+suffix+'\nf_'+suffix+' = '+equations['m_rhs_'+suffix]['expression'])
            body+='<div class="nest">'+''.join(equation_detail(r,1) for r in refs(equations['m_rhs_'+suffix]['expression']))+'</div>'
        for p in n['paragraphs']:
            body+='<p>'+escape(p)+'</p>'
        body+='<div class="nest">'+''.join(tree(c,path+(key,)) for c in children[key])+'</div>'
        if key in ('K','xi'):
            prefix='m_k_' if key=='K' else 'm_xi_'
            body+='<details><summary>状態別の係数</summary><div class="inside nest">'+''.join(equation_detail(prefix+s[2:],1) for s,v in nodes.items() if 'state' in v)+'</div></details>'
        body+=''.join(source_detail(s) for s in n['sources'])
        return '<details class="topic" id="node-'+'-'.join(path+(key,))+'"><summary><span>'+escape(n['title'])+'</span><small>'+escape(n['badge'])+'</small></summary><div class="inside">'+body+'</div></details>'

    out=''.join(tree(key) for key in ('x','B','mu','A','xi','K','time'))
    glossary=[]
    for key in sorted(used):
        n=equations[key]
        value=n['expression'] if n['expression'] is not None else str(n['value']) if n['value'] is not None else '未設定の入力'
        # Definitions are visible once this single appendix is opened; fragment
        # navigation can reveal ancestors in browsers supporting hidden-until-found.
        glossary.append('<section class="definition" id="eq-'+escape(key)+'"><h3>'+escape(n['label'])+'</h3>'+box(key+' = '+value)+'<p class="unit">'+escape(n['unit'])+'</p><p class="source">'+escape(n['source'])+'</p></section>')
    out+='<details class="topic glossary"><summary>式中の記号と定義<small>参照された全依存式</small></summary><div class="inside"><p>各式の識別子を見出しの式に併記しています。深い定義へのリンクが移動しないビューアでは、この付録を開いて参照してください。</p>'+''.join(glossary)+'</div></details>'
    return out
