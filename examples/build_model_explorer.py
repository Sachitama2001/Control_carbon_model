"""Build a self-contained, read-only explanation of the audited TKY pathway.

Run from the repository root. No new ecosystem equations are evaluated here.
State inventory is extracted from the existing workbook adapter; explanations
and source excerpts are embedded so the exported HTML also works via file://.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from control_carbon.visit_workbook import SOURCE_COMMIT, TKYModel


def build_data(source: Path) -> dict:
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if revision != SOURCE_COMMIT:
        raise ValueError("The explorer requires the pinned old VISIT revision")
    changes = subprocess.check_output(
        ["git", "-C", str(source), "status", "--porcelain", "--", "."], text=True)
    if changes.strip():
        raise ValueError("Source checkout contains changes; audit them before exporting")
    model = TKYModel(source=source)
    nodes = {}

    def node(key, title, summary, parent="home", children=(), paragraphs=(),
             formula="", sources=(), related=(), badge="仕組み", unit="", aliases=""):
        nodes[key] = dict(id=key, title=title, summary=summary, parent=parent,
                          children=list(children), paragraphs=list(paragraphs),
                          formula=formula, sources=list(sources), related=list(related),
                          badge=badge, unit=unit, aliases=aliases)

    node("home", "一つの式から、生態系の仕組みへ。",
         "貯留量が変わる理由を、入力・配分・移行からたどります。気になる記号や箱を選んでください。",
         parent=None, children=["x", "B"], badge="概要",
         paragraphs=["この資料は旧VISITのTKY・日次経路を案内します。最初に見える式は共通の概念表現です。VISITの具体的な説明では、1日ごとの増分 Δx/Δt に切り替わります。",
                     "数値入力は不要です。箱を開く操作は説明の展開であり、計算モデルの状態数や仮定を変更しません。"],
         related=["time", "guide"])
    node("x", "x ｜ 何が蓄えられているか", "炭素、窒素、水。モデルが翌日へ持ち越す量の一覧です。",
         children=["carbon", "nitrogen", "water", "aux"], badge="内訳", unit="物質ごとに異なる",
         paragraphs=["選択した監査ブックは37状態を持ちます。物質ごとに単位が違うため、すべてを足した総量は作りません。",
                     "前日の生理量や展葉の履歴も更新に必要です。37状態だけで閉じた自律ODEではありません。"], related=["history", "aggregation"])
    node("carbon", "炭素の箱", "植物の体と貯蔵炭素、地面に戻ったリター・腐植に分けます。", "x",
         ["plant", "soil"], unit="Mg C ha⁻¹", badge="内訳",
         paragraphs=["葉・幹・根は構造炭素、NSCは非構造性の貯蔵炭素です。どちらも炭素ですが、担う役割と移行の仕方が異なります。"])
    node("plant", "植物炭素", "樹木・C3林床・C4仮想群それぞれに、葉・幹・根・NSCがあります。", "carbon",
         ["tree", "floor", "c4"], badge="内訳", unit="Mg C ha⁻¹",
         paragraphs=["TKYのC4面積は0ですが、原典の窒素処理に残る群を監査ブックは保持しています。群ごとの炭素をそのまま足して実サイトの総炭素としないでください。"], related=["nsc"])
    for key, title, prefix in [("tree", "樹木", "t"), ("floor", "C3林床", "g"), ("c4", "C4仮想群", "v")]:
        node(key, title, "器官と貯蔵の違いを開いて見ます。", "plant",
             [f"x_{prefix}_{part}" for part in ("fol", "stm", "rot", "nsc")], badge="内訳", unit="Mg C ha⁻¹")
    node("soil", "土壌炭素", "6つのリタープールと3つの腐植プールを区別します。", "carbon",
         [f"x_s_{p}" for p in ("tf", "tc", "tr", "gf", "gc", "gr", "ha", "hi", "hp")],
         badge="内訳", unit="Mg C ha⁻¹", related=["microbes"],
         paragraphs=["この監査ブックの9炭素プールに、独立した微生物炭素の箱はありません。活性腐植を微生物炭素と読み替えることはしません。"])
    node("nitrogen", "窒素の箱", "植物の葉群・貯蔵窒素と、土壌の有機・無機窒素を分けます。", "x",
         ["plant_n", "soil_n"], badge="内訳", unit="g N ha⁻¹（原典収支の規約）")
    node("plant_n", "植物窒素", "この表現では、炭素の葉・幹・根と同じ分割ではありません。", "nitrogen",
         [f"x_{p}_n{o}" for p in ("t", "g", "v") for o in ("can", "str")], badge="内訳")
    node("soil_n", "土壌窒素", "微生物、リター、腐植、硝酸態、アンモニア態の5つです。", "nitrogen",
         [f"x_n_{p}" for p in ("mic", "lit", "hum", "no3", "nh4")], badge="内訳", related=["microbes"])
    node("water", "水の貯留", "積雪水当量、上層土壌水、深層土壌水の3状態です。", "x",
         ["x_w_snow", "x_w_sw", "x_w_dw"], badge="内訳", unit="mm（水深換算）",
         paragraphs=["葉・幹・根の独立した植物水貯留は、この選択経路の監査ブックにはありません。プロジェクト内の4炭素＋4水SPACモデルは別モデルです。",
                     "水収支には蒸発散・流出・融雪などが入ります。水の移動をすべて炭素の配分率と同じ解釈にはしません。"], sources=["hydrology"], related=["models"])
    node("aux", "補助状態", "濃度や水分記憶は、炭素量・水量そのものと分けて扱います。", "x",
         ["x_doc", "x_casa"], badge="内訳", related=["history"])
    node("history", "前日からの履歴", "気孔コンダクタンス、積算温度、展葉・落葉開始からの日数も必要です。", "x",
         paragraphs=["日次経路では、水収支が当日の生理特性更新より前に評価されます。前日の気孔コンダクタンスを現在の状態から勝手に再計算して置き換えません。",
                     "表示上の37状態とは別に、樹木・C3・C4それぞれの5履歴を監査ブックが入力として保持します。"], sources=["daily", "location"], related=["time"])
    node("nsc", "NSC ｜ 使うために蓄える炭素", "参照している旧VISITには、nsch_storage が実装されています。", "plant",
         ["x_t_nsc", "nsc_flush", "nsc_storage"], aliases="非構造性炭水化物 non structural carbohydrate",
         paragraphs=["各植生群が一つのNSC貯蔵を持ちます。樹木の葉・幹・根に別々のNSC状態を置く、器官別モデルではありません。",
                     "展葉時には貯蔵から葉へ炭素が移ります。幹・根へ向かう正の配分も、貯蔵量が容量に達していない場合はNSCへ回されます。"], sources=["flush", "storage"], related=["models", "allocation"])
    node("nsc_flush", "NSCから展葉へ", "TKYの展葉期には、その時点のNSCの一部を新しい葉に使います。", "nsc",
         formula="emerge = 0.03 × NSC\n葉への増加 = emerge − 0.10 × emerge\nNSCからの控除 = emerge",
         unit="Mg C ha⁻¹／選択した1日ステップ",
         paragraphs=["veg_type=4、season=2の分岐です。これは日次の移行量であり、任意の時間刻みへそのまま使う連続時間速度ではありません。",
                     "0.10の構築コストに関連する後段の呼吸計上も原典に存在します。この抜粋だけで、原典全体の炭素収支が正しいと認定するものではありません。"], sources=["flush"], related=["time"])
    node("nsc_storage", "成長より先に貯蔵する分岐", "幹向け、次いで根向けの配分を、現在のNSC量と容量で振り分けます。", "nsc",
         paragraphs=["幹向けの正の配分を処理した後に、更新されたNSCで根向けの条件を判定します。両者を同時に判定する式への置換は同じ処理ではありません。",
                     "容量までをminで充填する規則ではないため、1ステップで容量を超えることもあります。成長呼吸は後続の原典処理も含めて読む必要があります。"], sources=["storage"], related=["allocation"])
    node("microbes", "微生物はどこにいるか", "この監査範囲には、微生物窒素の独立した状態があります。", "soil_n",
         children=["x_n_mic"], aliases="microbe microbial 微生物炭素",
         paragraphs=["微生物窒素 n_mcrb は、有機化による流入と微生物由来の損失で更新されます。",
                     "一方、監査ブックの9つの土壌炭素状態には独立した微生物炭素を含めていません。原典の構造体に名称があることと、選択経路の状態として使うことは区別します。"], sources=["nitrogen"], related=["soil"])
    node("B", "B ｜ どの箱へ配るか", "入力や有効流量と、その受け取り先を結び付けます。",
         children=["allocation", "B_actual"], badge="仕組み",
         paragraphs=["最初は「炭素を葉・幹・根へどう配るか」を例に考えます。固定配分と動的配分は、異なる仮定として比較できます。",
                     "既存ブックのBは、この説明例より広い役割を持ちます。分解の実際も別に確認できます。"])
    node("allocation", "植物への炭素配分", "葉に投資するか、幹・根へ回すか。選択したVISITではLAIに応じて割合が変わります。", "B",
         ["dynamic", "compare"], aliases="allocation アロケーション",
         paragraphs=["正のEPPは、GPPから維持呼吸を引いた炭素です。配分後には成長呼吸・NSC貯蔵などの処理もあるため、配分量がそのまま最終的な器官の増加量になるとは限りません。"], related=["nsc_storage"], sources=["allocation"])
    node("compare", "固定配分と動的配分を比べる", "表示の細かさではなく、配分を決める仮定の違いです。", "allocation",
         ["fixed", "dynamic"], badge="仮定を比べる",
         paragraphs=["固定割合の例は理解のための概念モデルです。VISITの実際の動的配分を、説明を閉じるだけで固定割合へ変更することはありません。"])
    node("fixed", "固定割合の配分", "入力を、あらかじめ決めた割合で葉・幹・根へ分ける説明用のモデルです。", "compare",
         formula="b = (b葉, b幹, b根)ᵀ\nbᵢ ≥ 0,  Σbᵢ = 1\n各器官への入力 = bᵢ μ",
         badge="説明用の仮定", unit="b：無次元",
         paragraphs=["μが増減しても配分率bは変わりません。ここではTKYの採用値を設定せず、実サイトの配分率とも主張しません。"], related=["dynamic"])
    node("dynamic", "LAIに応じた動的配分", "実際の葉面積と目標の葉面積を比べ、葉・幹・根への割合を決めます。", "allocation",
         ["lai", "target", "branches"], sources=["allocation"],
         paragraphs=["実LAIが目標を超えると、正のEPPの葉への配分はゼロになります。目標以下では不足分と利用可能なEPPを比較します。",
                     "この規則には切替境界があります。単純に「不足するほど葉へ多く」という滑らかな関数ではありません。"], related=["compare"])
    node("lai", "実LAI ｜ 現在の葉面積", "葉の構造炭素量を、SLAを使って地面積当たりの葉面積へ換算します。", "dynamic",
         formula="L = γ C葉\nγ = 2.2 × SLA / 200", unit="L：m² m⁻²、C葉：Mg C ha⁻¹、SLA：cm² gDM⁻¹",
         paragraphs=["これは選択ソースの数値換算です。NSCを葉の構造炭素に加えてからLAIへ換算することはしません。"], sources=["lai"], related=["x_t_fol"])
    node("target", "目標LAI ｜ どこまで葉を持つか", "光獲得と葉の維持・更新コストから、原典の式が目標値を計算します。", "dynamic",
         ["target_formula"], sources=["target"],
         paragraphs=["光、光合成能力、温度、比葉面積、更新・呼吸の係数が関係します。目標値を計算する関数と、実際の配分を決める関数は別です。",
                     "「最適」はこの式の仮定の下での意味です。森林一般の最適性や、常にその葉量を実現できることを意味しません。"], related=["branches", "lai"])
    node("target_formula", "目標LAIの具体式", "f_opt_lai の代数を、原典と対応する記号で示します。", "target",
         formula="aₘ = (rmf / 1000) qTc^((T−15)/10) × 2.2 × 10000 / SLA\na𝗀 = lf × 2.2 × 10000 / SLA × (1+rgf)\na = aₘ + a𝗀\nc = P × [Ph / (Ph−24a) − 1]\nL* = log(max(1, kqI / c)) / k   （P > 0）\nL* = 0   （P ≤ 0）",
         unit="L*：m² m⁻²（LAI）",
         paragraphs=["P=psat、h=日長、k=eK、q=lue、I=ppfd_t。Pはµmol CO₂ m⁻² s⁻¹、hはhour、Iはµmol photon m⁻² s⁻¹です。",
                     "葉のコストですがコードはqTcを使用します。qTfへ読み替えません。rmf、lfなどは、この関数より前の生理計算で設定されます。",
                     "Ph=24a、c=0、k=0、SLA=0等では特異性があります。aはソースのコスト中間量で、Pと時間の基準に対する次元整合性は未解決です。この画面は代数の説明であり数値予測は行いません。"], sources=["target", "physiology"])
    node("branches", "配分の条件分岐", "非作物・活動期・EPP>0の分岐を先に読みます。", "dynamic",
         formula="E = EPP,  ΔC = (L*−L) / γ\na_f = alloc_ass,  a_s = alloc_abg\n\nL > L*                       → b葉 = 0\nL ≤ L*, ΔC ≤ a_f E           → b葉 = a_f\nL ≤ L*, ΔC > a_f E           → b葉 = min(ΔC/E, 0.05)\n\nb幹 = (1−b葉) a_s\nb根 = (1−b葉) (1−a_s)",
         unit="b：無次元、E・ΔC：Mg C ha⁻¹（Eは1日分）",
         paragraphs=["この整理はa_f>0の範囲です。a_f=0では元コードの除算を別途確認する必要があります。L=L*でもb葉=a_fが残り、直上でゼロになります。0.05は全分岐共通の上限ではありません。",
                     "休眠期season=0では目標以下の枝で葉配分に0.7を掛けます。EPP≤0ではGPPの固定配分から器官別維持呼吸を引く別式です。作物の穀物配分は今回のTKY樹木の説明範囲外です。",
                     "原典のコメントにmonthlyとある箇所もありますが、選択した呼び出し経路は日次です。"], sources=["allocation"], related=["nsc_storage", "time"])
    node("B_actual", "監査ブックでのBμ", "37状態に対し、41本の系外・有効入力と37本の残差を対応付けます。", "B",
         formula="B：37 × 78\nμ：78成分（41有効入力＋37残差）",
         paragraphs=["植物のμは正の器官純入力を含み、GPPそのものではありません。NSC動員・内部再配分や負の配分は流量台帳に別途記録されます。",
                     "供給状態がゼロなどでAξKの除算ができない列は、有効残差としてBμ側で表現します。これを新たな大気からの炭素入力とは解釈しません。",
                     "分解は一意ではありません。これは既存ブックの規約です。"], related=["mu", "A"], sources=["matrix"])
    for key, title, summary, body in [
        ("mu", "μ ｜ 入力の大きさ", "どれだけの入力・有効流量があるかを表します。", "降水、窒素沈着、植物への炭素入力では単位が異なります。既存ブックは内部の有効残差も含むため、すべてが系外入力ではありません。"),
        ("A", "A ｜ 移行と損失", "ある箱から出る量が、どの箱へ移るかを表します。", "炭素の典型的な規約は供給側が列、受け取り側が行です。原典の負の流量・残差がある場合まで、正の移行率行列であることを保証しません。"),
        ("xi", "ξ ｜ 環境・状態の影響", "温度や水分などによって、基準の速度が変わります。", "ξ(t)という短い表記でも、実際には状態・気象・履歴に依存します。原典を再構成した有効係数と、単純な温度倍率を同一視しません。"),
        ("K", "K ｜ 基準の速度", "turnoverや分解の基準を置く行列です。", "既存ブックでは窒素・水・補助状態のK=1/dayは代数上の規格化です。生物学的な寿命が1日という意味ではありません。"),
    ]:
        node(key, title, summary, paragraphs=[body, "この試作では概要までを収録しています。"], related=["B_actual"], badge="概要", sources=["matrix"])
    node("time", "概念の微分式と、VISITの日次更新", "このモデルの具体的な計算は、順序を持つ1日ごとの更新です。",
         formula="fΔ = [FΔ(x,h,u) − x] / Δt\nx翌日 = x + Δt fΔ     （Δt = 1 day）\nfΔ = Bμ + AξKx",
         paragraphs=["日次写像の増分を、この形で書き直しています。この右辺を任意の刻みでODE積分しても、原典の日次結果と一致するとは限りません。",
                     "元の式の短い表記では省略した状態x・履歴h・気象uへの依存が、係数にも含まれます。AξKは一般に全体系のヤコビアンではありません。"], sources=["daily", "matrix"])
    node("aggregation", "箱をまとめる意味", "説明の集約は、元のモデルを計算し直しません。", "x",
         formula="y = P x\ndy/dt = P f(x,u)",
         paragraphs=["同じ物質・面積基準で重複しない状態だけを合計できます。yだけで閉じた方程式が得られるかは別問題です。",
                     "この試作では未設定の状態に数値を入れず、合計値も捏造しません。C・N・水やDOC濃度をまとめて加算しません。"])
    node("models", "別モデルでの表現", "同じプロジェクトでも、持っている状態と仮定は異なります。",
         badge="仮定を比べる", paragraphs=["旧VISIT・TKY：各植生群に単一NSC、選択した水貯留は雪と2層土壌水。今回の主な説明対象です。",
         "温度–NSC理論モデル：4構造炭素＋4可動・易分解性炭素。葉・幹・根の器官別NSCがあります。土壌の可動炭素は植物NSCではありません。",
         "SPACモデル：4炭素＋4水。葉・幹・根・土壌の水貯留とポテンシャル輸送を持つ、独立した連続モデルです。",
         "この試作は比較の説明までです。クリックによって別モデルの計算へ切り替わる機能はありません。"], related=["nsc", "water"])
    node("guide", "この資料の読み方", "知りたい記号を選び、一つずつ説明を開いてください。",
         paragraphs=["「内訳」は同じモデルの箱を開く操作、「仕組み」は同じ過程の詳しい説明、「仮定を比べる」は異なる定式化の比較です。",
                     "各ページの原典を開くと、ネット接続なしでソースの抜粋を読めます。外部のGitHubリンクを開く場合だけネット接続が必要です。",
                     "リンクを保存すると同じ項目へ戻れます。印刷は現在の説明と原典を出力します。全ページの一括印刷ではありません。",
                     "既存の収支・native照合の結果を参照していますが、この閲覧資料はシミュレーターではありません。"])

    soil_labels = dict(tf="樹木・葉リター", tc="樹木・幹リター", tr="樹木・根リター",
                       gf="林床・葉リター", gc="林床・幹リター", gr="林床・根リター",
                       ha="活性腐植", hi="中間腐植", hp="安定腐植")
    for s in model.states:
        if s.name.startswith(("x_t_", "x_g_", "x_v_")):
            prefix, part = s.name.split("_")[1:]
            parent = "plant_n" if s.group == "N" else dict(t="tree", g="floor", v="c4")[prefix]
        elif s.name.startswith("x_s_"):
            parent = "soil"
        elif s.group == "N":
            parent = "soil_n"
        elif s.group == "W":
            parent = "water"
        else:
            parent = "aux"
        title = soil_labels.get(s.name[4:], s.label) if s.name.startswith("x_s_") else s.label
        summary = "監査ブックに実在する状態です。現在値はこの資料では設定していません。"
        related = []
        if s.name.endswith("_nsc"):
            summary = "この植生群が持つ単一の非構造性炭素貯蔵です。器官別に分割されたNSCではありません。"
            related = ["nsc", "nsc_flush", "nsc_storage"]
        elif s.name.endswith("_fol"):
            summary = "葉の構造を構成する炭素です。LAIの計算に使われます。"
            related = ["lai", "allocation"]
        elif s.name == "x_n_mic":
            summary = "微生物に保持された窒素。炭素の活性腐植とは別の状態です。"
            related = ["microbes"]
        node(s.name, title, summary, parent=parent, unit=s.unit, badge="状態",
             paragraphs=[f"対応ID：{s.name}", f"既存ブックの状態出典：{s.source}",
                         "単位は監査ブックの規約を保持しています。初期化用の小さな種値を、現在の森林の値として表示していません。"],
             sources=["state_inventory"], related=related, aliases=s.name)
        nodes[s.name]["state"] = asdict(s)

    # Source text is embedded verbatim, with numbered lines and content hashes.
    specs = {
        "allocation": (source / "allocation.c", 16, 110, "f_allocation", "配分率→配分量。呼吸・NSC処理より前"),
        "target": (source / "ecophysiology.c", 270, 300, "f_opt_lai", "生理係数更新後、配分より前"),
        "physiology": (source / "ecophysiology.c", 60, 88, "f_ecophysiology", "呼吸・更新率→目標LAI→貯蔵容量"),
        "lai": (source / "ecophysiology.c", 108, 127, "lai_mass", "構造炭素からLAIを診断"),
        "flush": (source / "plant_proc.c", 16, 66, "plant_process", "展葉→後続の日次植物過程"),
        "storage": (source / "plant_proc.c", 200, 242, "plant_process", "葉更新→幹／貯蔵→根／貯蔵→救済再配分"),
        "nitrogen": (source / "n_budget.c", 85, 116, "n_budget", "日次窒素流量の計算後に状態更新"),
        "hydrology": (source / "hydro_balance.c", 16, 62, "f_hydrology", "日次の水収支。抜粋は関数の冒頭のみ"),
        "daily": (source / "daily_scheme.c", 15, 85, "daily_scheme", "選択した日次経路。抜粋の呼び出し順を参照"),
        "location": (source / "location_proc.c", 1, 55, "f_loct_proc", "水収支は当日の生理更新より前。抜粋は冒頭"),
        "matrix": (ROOT / "docs/visit_matrix_workbook.md", 30, 70, "行列の具体的な構成", "source-order daily effective-rate representation"),
        "state_inventory": (ROOT / "src/control_carbon/visit_workbook.py", 328, 350, "TKYModel._inputs_states", "監査ブックの状態登録。native更新関数ではない"),
    }
    refs = {}
    for key, (path, start, end, symbol, order) in specs.items():
        text = path.read_text(errors="replace")
        native = path.parent == source
        rel = "visit_local/" + path.name if native else str(path.relative_to(ROOT))
        refs[key] = dict(path=rel, symbol=symbol, start=start, end=end, order=order,
                         repository="Sachitama2001/VISIT-matrix" if native else "Control_carbon_model（ローカル監査資料）",
                         commit=SOURCE_COMMIT if native else "生成manifestのファイルhashで固定",
                         level="native nonlinear discrete VISIT map（抜粋）" if native else "source-order transcription / 監査資料",
                         sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         excerpt="\n".join(f"{i}: {line}" for i, line in enumerate(text.splitlines(), 1) if start <= i <= end),
                         url=f"https://github.com/Sachitama2001/VISIT-matrix/blob/{SOURCE_COMMIT}/{rel}#L{start}" if native else "")
    for item in nodes.values():
        for target in item["children"] + item["related"]:
            if target not in nodes:
                raise ValueError(f"Broken reference {item['id']} → {target}")
        for target in item["sources"]:
            if target not in refs:
                raise ValueError(f"Missing source {target}")
    return dict(model_id="visit-tky-daily-3285bd8", revision=revision, nodes=nodes,
                sources=refs, state_count=len(model.states), config_sha256=model.config_hash,
                equations={k: dict(label=n.label, expression=n.expression, unit=n.unit, source=n.source, note=n.note, value=n.value) for k,n in model.graph.nodes.items()},
                scope="旧VISIT / TKY / FLUX_SCHEME=0 / N_CYCLE=1", version="0.2")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT.parent / "VISIT-matrix/visit_local")
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts/model_explorer/index.html")
    args = parser.parse_args()
    data = build_data(args.source.resolve())
    template = ROOT / "docs/model_explorer/template.html"
    from render_nested_explorer import render_nested
    result = template.read_text().replace("__NESTED_CONTENT__", render_nested(data))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(result)
    manifest = dict(model_id=data["model_id"], source_commit=data["revision"],
                    state_count=data["state_count"], content_count=len(data["nodes"]),
                    config_sha256=data["config_sha256"],
                    build_inputs={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [Path(__file__), ROOT / "examples/render_nested_explorer.py", template, ROOT / "src/control_carbon/visit_workbook.py"]},
                    html_sha256=hashlib.sha256(result.encode()).hexdigest(),
                    source_files={r["path"]: r["sha256"] for r in data["sources"].values()})
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f"Built {args.output} ({len(data['nodes'])} explanations, {data['state_count']} states)")


if __name__ == "__main__":
    main()
