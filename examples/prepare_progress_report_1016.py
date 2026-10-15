"""Render presentation assets from unchanged, independently audited SPAC results."""

import argparse
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from scipy.linalg import expm

from analyze_spac_carbon_capacity import encode, require, write_csv
from control_carbon.manzoni_hydraulics import HydraulicParameters
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, carbon_matrix
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_e4 import hopf_point

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "b8275212c11ef1e05ddbf57d052735aa9f20f171"
ARCHIVES = {
    "scalar": ROOT/"artifacts/spac_carbon/scalar_equilibria_final_v5_20261014/results.json",
    "capacity": ROOT/"artifacts/spac_carbon/stage_c_capacity_final_v2_20261013/results.json",
    "E4": ROOT/"artifacts/spac_carbon/stage_b3_e4_final_v2_20261013/results.json",
}
REPORTS = [
    "spac_carbon_scalar_equilibrium_results.md", "spac_carbon_capacity_definition_ledger.md",
    "spac_carbon_stage_c_results.md", "spac_carbon_stage_b3_e4_results.md",
    "spac_carbon_stage_c23_results.md", "manzoni2014_reproduction.md",
]
TEST_FILES = ["test_spac_carbon_equilibria.py", "test_spac_carbon_capacity.py",
              "test_spac_carbon_e4.py", "test_spac_carbon_mpa.py"]
COLORS = {"E0": "#0072B2", "E1": "#D55E00", "E4": "#009E73"}
CASES = ["wet_initial", "dry_initial", "both_closed_initial", "xylem_closed_low_leaf"]
CASE_COLORS = ["#009E73", "#D55E00", "#CC79A7", "#0072B2"]
CASE_LABELS = ["wet: leaf=-1.5 MPa", "dry: leaf=-2.5 MPa",
               "both closed: leaf=-4.2 MPa", "xylem closed: leaf=-3.5 MPa, low L"]

# Japanese prose is intentional: this is the requested speaker-facing evidence note.
CLAIMS = [
    ("A1", "解析的証明", "capacityは写像、収束先とは別",
     "Capacity is a pointwise diagnostic, not automatically an attracting equilibrium.",
     "炭素収支をdC/dt=u+MCと書き、Mが可逆ならK=-M^{-1}u、P=K-Cを現在の状態で計算できる。固定したuとHurwitz行列Mの線形系ではKが一意な安定平衡になる。一方、状態に依存する入力や係数を持つ場合も定義式は有効だが、C=Kの自己整合性と水理収支の停止を別に確かめる必要がある。P=0は炭素収支が止まることと同値でも、水理状態の停止は意味しない。Luoらが非線形性を除外したという主張はせず、瞬間の診断値を長期貯留先と読むための条件を問う。"),
    ("A2", "解析的証明", "6状態の平衡条件を葉量一変数に帰着",
     "All positive conducting equilibria reduce to one scalar equation in leaf carbon.",
     "正の炭素平衡では葉の収支からGPP=k_L L/b_Lを得て、木部・根・土壌の平衡式へ代入するとW=a_W L、R=a_R L、S=a_S Lとなる。葉収支を光合成式と照合してh(L)を求め、気孔式を反転して葉水ポテンシャルを得る。木部水収支で木部ポテンシャルを消し、最後の葉水収支がD=D_curve(L)になる。通水conductanceと葉量が正であることが条件で、ゼロ炭素や閉鎖集合は別に扱う。この消去は平衡条件の整理であって、時間発展を一次元に近似したものではない。"),
    ("A3", "解析的証明", "基準E0/E1は最大二根、E4は最大一根",
     "Baseline E0 and E1 have at most two positive roots; E4 has at most one.",
     "基準unsafe係数ではh(L)は厳密に増え、正の通水平衡の葉量域は0<L<=L_wetとなる。E0のVPD曲線はh=1/2で一つだけ最大点を持つ。E1では供給係数を掛けた曲線の対数導関数を比較して、一意な最大点を示せる。E4では供給・需要・木部の係数を整理し、増加項より減少項が大きいことを解析的不等式と外向き丸めの定数評価で示した。その結果、最大二根と最大一根の分類は単なる走査結果ではない。L=0の開端点は正根に数えず、別の形質・木部指数への一般化もしていない。"),
    ("A4", "解析的証明", "不安定な平衡でもK=C*",
     "Equilibrium capacity equals stored carbon even when the equilibrium is unstable.",
     "全結合平衡ではMC*+u*=0なので、K(y*,v)=-M^{-1}u*=C*が成立する。これは安定性を使わない恒等式であり、不安定な平衡でも成立する。同じ固定VPDに正根が二つある場合には、水理停止とC=Kを同時に満たす自己整合的capacity候補が二つある。瞬間capacity写像が二値になるという意味ではない。図では各根を全6状態Jacobianで分類し、安定枝と不安定枝を別の線種で示す。候補の数、安定な貯留先の数、実際の軌道が到達する先は別の問題として説明する。"),
    ("N1", "数値検証", "一意なE4正根も不安定化する",
     "A unique positive equilibrium can lose stability through a numerically verified subcritical Hopf.",
     "E4では葉だけでなく木部炭素も水理へ戻るため、水理・葉・木部の4状態の応答時間を残して全6状態Jacobianを調べる必要がある。主VPD範囲の共役固有値対の横断と正の第一Lyapunov係数は、既存局所解析を現在のコードで再評価して照合した。近傍の既知周期軌道は不安定で、別の吸引貯留先ではない。一意な正根やK=C*だけでは安定性を保証できないことが結論である。Hopfと係数は精度を確認した局所数値結果で、区間認証や大域的な回復境界の証明にはしていない。"),
    ("A5", "解析的証明", "気孔閉鎖ではKが先にゼロになる",
     "Stomatal shutdown eliminates capacity immediately; stored carbon decays only under persistent input loss.",
     "気孔conductanceがゼロなら開度hもゼロとなり、GPP、炭素入力u、capacity Kはその時点でゼロになる。現在の炭素ストックCまで直ちにゼロになるわけではない。両conductanceが閉じる領域では水理フラックスも消え、水ポテンシャルが固定されるため、炭素収支はC(t)=exp(Mt)C(0)として持続する。MのHurwitz性から長期にはストックも消失するが、土壌への内部移動により各プールの単調減少は一般には言えない。図の退出後の破線は、この条件が成り立つ両閉鎖例の炭素部分だけに用いる。"),
    ("A6", "解析的証明", "木部閉鎖だけでも炭素衰退は起こり得る",
     "Xylem closure can drive carbon loss while stomatal conductance remains positive.",
     "木部だけが閉じた領域では気孔がまだ開き、瞬間GPPとKは正になり得る。ただし非負VPDでは葉水ポテンシャルが増えず、この領域は前向き不変となる。今回のunsafe係数では最大開度が同化維持に必要な開度より低く、葉量は正の減衰率で指数的に上から評価できる。入力も減るため木部・根・土壌とKも漸近的に消える。E4では蒸散面積が葉量と共に減り、正の気孔開度を極限に残せる。したがって完全気孔閉鎖は炭素消失の必要条件ではなく、この結論は基準係数と形式的ODE極限に限定する。"),
    ("N2", "数値検証", "同じ固定VPDでも初期状態で帰結が分かれる",
     "At the same fixed VPD, initial states approach recovery or leave the positive-leaf domain.",
     "固定VPDでE0とE4を比較し、各モデルの検証済み正平衡を基にした炭素初期値と、異なる初期葉水ポテンシャルを用いた。湿潤例は100年後に正平衡近傍へ接近し、乾燥例と両閉鎖例は約14年で葉量下限に達する。小葉量・木部閉鎖例はより早く退出し、E4ではその時点でも気孔conductanceが正に残る。現在HEADで既存runnerを再実行し、全初期値とsolver精緻化を照合した。この少数例は機構を見せるためのもので、退出を死亡や別吸引子と呼んだり、大域basinを確定したりはしない。"),
    ("D1", "計算領域に限定", "正枝と形式的極限を区別する",
     "Positive equilibria, numerical domain exits, and formal zero-carbon limits must remain distinct.",
     "数学的な正通水枝は正葉量の開区間から湿潤端点までであり、数値積分のL_minとは異なる。積分はL_minで停止し、それより先の全水理–炭素軌道は表示しない。ゼロ炭素拡張には非孤立の平衡集合が現れるので、孤立した死亡sinkと解釈できない。高VPDのE4端点と第2固有値横断は25°Cの飽和水蒸気圧を超える形式的結果として分ける。極小葉量で許容値内の固有値は真のゼロ固有値の証明ではない。25°CはVPDの気象整合性診断であり、温度依存生理や死亡閾値には用いていない。"),
    ("U1", "未解決", "局所安定性から大域予測への接続",
     "Connecting local stability to global recovery boundaries and changing-climate tracking remains unresolved.",
     "今回整理できたのは、capacityの計算、自己整合的平衡の存在、全結合系の局所安定性、閉鎖条件下の炭素減衰を区別する枠組みである。局所の不安定周期軌道や安定方向が大域的な回復・退出境界をどこまで構成するかは未解決である。環境が時間変化するときの追従も、固定外力の安定性だけからは保証されない。今後は同一条件で局所構造と大域帰結をつなぐ検証、係数の実測根拠、必要なら退出後の生理を検討する。既知モデルの結合や標準解析手法を使ったこと自体の一般的新規性は主張せず、文献比較を別途行う。"),
]

SLIDES = [
    (1, "Title", "We distinguish capacity, self-consistent equilibria, stability, and realized carbon storage.", "01", "A1,U1",
     "題名のmultistabilityは複数の安定吸引先を連想させる。正根二つの一方は不安定なので、equilibriaとcapacityを軸にした題名を推奨。"),
    (2, "Table of contents", "The argument proceeds from definitions to equilibria, stability, and shutdown.", "", "A1,A2,A4,A5",
     "何を知りたいか、使う式、計算結果、まだ言えないことの順に案内する。"),
    (3, "Previous research", "Matrix notation separates carbon inputs, transfers, and losses.", "01", "A1",
     "既存研究の導入。先行研究の内容と今回の合成モデルの成果を分ける。"),
    (4, "Land carbon cycle", "Carbon storage reflects the balance between inputs, transfers, and losses.", "", "A1",
     "陸域炭素収支の背景。観測・地球規模の炭素動向を今回の低次元モデルの検証結果とはしない。"),
    (5, "Carbon equilibrium", "An equilibrium requires every modeled state to stop changing.", "01", "A1",
     "全状態の停止を定義する。既存スライドの『Earth system is closed』というCox論文の要約は未監査なので原文確認まで断定を避ける。"),
    (6, "Matrix approach", "Capacity is defined pointwise whenever the carbon balance matrix is invertible.", "01", "A1",
     "u、M、K、Pを順に示す。Luoが非線形モデルを禁止したという研究ギャップを作らない。"),
    (7, "Carbon storage capacity", "Capacity is an attracting equilibrium only under additional conditions.", "01", "A1,A4",
     "固定u・Hurwitz Mの線形系では一意な安定平衡。非線形系では写像、固定点、吸引先の区別が必要。"),
    (8, "Research gap", "Computing capacity does not establish a self-consistent hydraulic-carbon equilibrium.", "01", "A1,A2",
     "非線形性があるというだけではなく、Kの計算から全平衡への論理の不足を問いにする。"),
    (9, "Research gap", "Local stability and actual convergence are separate questions.", "04", "A4,N1,U1",
     "『非線形条件はfrontier』という一般論を避ける。標準手法で解ける部分と大域接続の未解決部分を分ける。"),
    (10, "Research gap", "Zero capacity does not imply instantaneously zero carbon stocks.", "05", "A5",
     "入力の停止と既存ストックの減衰の時間順序を示す。"),
    (11, "Research question", "When can capacity be interpreted as a long-term carbon storage target?", "01", "A1,A4,U1",
     "中心となる問いを一文に限定。存在、安定性、閉鎖下の減衰を回答の三層として予告する。"),
    (12, "Methods", "We couple four carbon pools to two hydraulic states under fixed VPD.", "02", "A2,D1",
     "root ψ=0、合成係数、VPDのみを動かす設定を明示。今の発表の代表例は固定VPDであり、新しいramp実験ではない。"),
    (13, "SPAC states", "Hydraulic storage follows supply and transpiration flux balances.", "02", "A2",
     "木部と葉の二つの水ポテンシャルを示す。土壌水収支はなく、根境界を固定している。"),
    (14, "SPAC conductance", "Xylem and stomatal closure occur at different leaf water potentials.", "05", "A5,A6",
     "unsafeでは木部−3 MPa、気孔−4 MPa。閉木部から直ちにGPP=0とは言えない。"),
    (15, "SPAC feedback", "Leaf, demand, and wood feedbacks change the equilibrium structure.", "03", "A2,A3",
     "E0、E1、E4の係数表を用いる。新しい生理式を導入したのではなく、既存feedbackの比較。"),
    (16, "Nonlinear analysis", "Equilibrium elimination classifies roots; full Jacobians classify stability.", "02", "A2,A3,N1",
     "catastrophe theory一語でまとめず、変数消去、導関数の符号、全Jacobianの順に具体化する。"),
    (17, "Model equations", "Six equilibrium conditions reduce exactly to VPD as a function of leaf carbon.", "02", "A2",
     "消去の順序を見せる。正通水・L>0の仮定と、動力学を一次元化しないことを口頭でも添える。"),
    (18, "Equilibria", "E0 and E1 permit two positive roots; baseline E4 permits only one.", "03", "A3,A4",
     "foldと2.44 kPaを指す。E1にはこの代表VPDで正根がない。E4の高VPD端点は気象域外。"),
    (19, "Stability", "A unique positive equilibrium need not be an attracting storage target.", "04", "A4,N1",
     "このページはK=C*と全Jacobianの安定性だけに限定。追加の19Bページで図05とA5の一文結論を使い、閉鎖時の入力停止とストック減衰を説明する。"),
    (20, "Contribution", "We separate capacity calculation, equilibrium existence, coupled stability, and carbon loss.", "01,04", "A1,A3,A4,A5,U1",
     "既存水・炭素モデルを結合しただけの新規性ではなく、長期貯留解釈の条件を三層で整理したことを述べる。"),
    (21, "Future plans", "Global boundaries, climate tracking, and physiological validation require further evidence.", "", "U1,D1",
     "大域境界、非自律追従、実測係数、退出後生理は今後の課題。本作業では追加しない。"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    with path.open(newline="") as stream:
        return [{**row, "source_csv_line": i} for i, row in enumerate(csv.DictReader(stream), 2)]


def column(rows, key):
    return np.array([float(row[key]) for row in rows])


def model_from_evidence(raw, label="E4"):
    cfg = raw["baseline_config"]
    c = cfg["carbon"]
    carbon = CarbonParameters(c["assimilation_max_per_year"]/YEAR_DAYS, tuple(c["allocation"]),
                             tuple(v/YEAR_DAYS for v in c["losses_per_year"]), tuple(c["litter_fractions"]),
                             c["specific_leaf_area"], c["carbon_dry_fraction"], c["extinction"])
    refs = raw["state_scales"]["state"]
    feedback = StructuralFeedback(refs[2], refs[3], cfg["leaf_floor"], cfg["minimum_leaf"],
                                  *cfg["variants"][label], wood_exponent=cfg["wood_exponent"])
    return PressureSPAC(HydraulicParameters(**cfg["hydraulic"]), carbon, feedback,
                        cfg["root_potential"], raw["settings"]["atmospheric_pressure_mpa"])


def audit(raw):
    old, checks, inputs = {}, {}, {}
    for label, path in ARCHIVES.items():
        old[label] = json.loads(path.read_text())
        mismatches = []
        for name, expected in old[label]["source_sha256"].items():
            source = Path(name)
            require(source.is_file(), f"Archived source is missing: {source}")
            actual = digest(source)
            if actual != expected:
                mismatches.append(name)
            inputs[str(source)] = actual
        require(not mismatches, f"Archived evidence differs from current sources: {mismatches}")
        inputs[str(path)] = digest(path)
        checks[label] = {"archive": str(path), "source_hashes_checked": len(old[label]["source_sha256"]),
                         "source_mismatches": mismatches}
    require(raw["run_git_HEAD"] == subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "Fresh scalar evidence does not come from current HEAD")
    for key in ["settings", "baseline_config", "state_scales"]:
        require(raw[key] == old["scalar"][key], f"Scalar configuration changed: {key}")
    for label in COLORS:
        for key in ["zero_leaf_vpd_limit_mpa", "fold_leaf", "fold_vpd_mpa"]:
            a, b = raw["landmarks_and_shape_proofs"][label][key], old["scalar"]["landmarks_and_shape_proofs"][label][key]
            require(a == b, f"Scalar landmark differs from verified archive: {label}/{key}")
    require(len(raw["outcomes"]) == 8, "Expected the eight verified initial states")
    for new, prior in zip(raw["outcomes"], old["scalar"]["outcomes"]):
        require((new["model"], new["case"], new["status"]) == (prior["model"], prior["case"], prior["status"]),
                "Fixed-VPD outcome classification differs")
        require(new["initial"] == prior["initial"], "Initial states differ from verified JSON")
        require(abs(new["last_day"]-prior["last_day"]) < .1, "Exit-time comparison failed")
    model = model_from_evidence(raw)
    hopf = hopf_point((-2.01, -1.99), model)
    require(abs(hopf["vpd_mpa"]-old["E4"]["hopf"]["vpd_mpa"]) < 1e-11, "E4 Hopf VPD audit failed")
    require(abs(hopf["l1"]-old["E4"]["hopf"]["l1"]) < 1e-7 and hopf["l1"] > 0, "E4 l1 audit failed")
    counter_path = ARCHIVES["capacity"].parent/"counterexamples.csv"
    counters = read_rows(counter_path)
    examples = [r for r in counters if float(r["potential_norm_scaled"]) < 1e-10
                and float(r["water_residual_norm_scaled_per_day"]) > 1e-4]
    require(examples, "Archived P=0 / moving-water evidence was not found")
    inputs[str(counter_path)] = digest(counter_path)
    checks["capacity"]["carbon_stationary_water_moving_csv_keys"] = [
        {"case": row["case"], "source_csv_line": row["source_csv_line"]} for row in examples]
    checks["E4"]["current_hopf_recalculation"] = {
        key: hopf[key] for key in ["vpd_mpa", "omega", "l1", "state", "slope_refinement_per_mpa"]}
    checks["E4"]["periodic_results"] = "source-hash audit of archived collocation/Floquet evidence; not rerun"
    return old, checks, inputs, model


def save_figure(fig, output, name, title, footer):
    fig.suptitle(title, fontsize=15, y=.985)
    fig.text(.5, .015, footer, ha="center", va="bottom", fontsize=9)
    multiple = len(fig.axes) > 1
    fig.subplots_adjust(left=.10 if multiple else .125, right=.97 if multiple else .9,
                        bottom=.18 if multiple else .16, top=.80 if multiple else .86,
                        hspace=.55, wspace=.48 if len(fig.axes) == 10 else .32)
    for suffix in ["png", "pdf"]:
        fig.savefig(output/f"{name}.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def box(ax, xy, width, height, text, *, color="#F1F4F6", fontsize=13):
    ax.add_patch(FancyBboxPatch(xy, width, height, boxstyle="round,pad=.012",
                              facecolor=color, edgecolor="#54626F", linewidth=1.))
    ax.text(xy[0]+width/2, xy[1]+height/2, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="->", mutation_scale=15,
                                linewidth=1.2, color="#54626F"))


def conceptual_figures(output, raw):
    fig, ax = plt.subplots(figsize=(12.8, 6.4))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    box(ax, (.035, .56), .4, .29,
        "Fixed linear coefficients\n"
        r"$\dot C=u+MC,\quad K=-M^{-1}u$" "\n"
        "Fixed u; Hurwitz M\nUnique attracting equilibrium C*=K")
    box(ax, (.565, .56), .4, .29,
        "State-dependent carbon balance\n"
        r"$K(y,v)=-M(y,v)^{-1}u(y,v)$" "\n"
        "Pointwise diagnostic\nNot automatically an attracting target")
    arrow(ax, (.435, .7), (.565, .7))
    box(ax, (.05, .2), .9, .21,
        r"$P=K-C=-M^{-1}\dot C$" "\n"
        r"$P=0\ \Longleftrightarrow\ \dot C=0$" "\n"
        "Carbon balance stops; hydraulic balance need not stop.")
    arrow(ax, (.765, .56), (.765, .41))
    save_figure(fig, output, "01_concept_capacity", "Schematic: computing capacity is not proving convergence",
                "Conceptual diagram, not experimental data. Nonlinear capacity requires invertible M; coupled stability is a separate question.")
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    steps = [
        ("Four carbon balances", r"$W=a_WL,\ R=a_RL,\ S=a_SL$" "\n"
         r"$a_W=21.3333,\quad a_R=4,\quad a_S=54$"),
        ("Leaf balance\n+ canopy assimilation", r"$h(L)=rL/(1-e^{-qL}),\quad r=k_L/(b_Lp)$"),
        ("Invert stomatal conductance", r"$\psi(L)=-4[1-h(L)]$"),
        ("Eliminate xylem potential", r"$x=G_l\psi/(G_r+G_l),\quad G_{\rm eff}=G_rG_l/(G_r+G_l)$"),
        ("Remaining leaf-water balance", r"$D=D_{\rm curve}(L)=G_{\rm eff}(-\psi)/(g_sF_E)$"),
    ]
    for i, (label, equation) in enumerate(steps):
        y = .83-i*.17
        box(ax, (.025, y), .26, .115, label, fontsize=11)
        box(ax, (.34, y), .625, .115, equation, fontsize=13)
        arrow(ax, (.285, y+.0575), (.34, y+.0575))
        if i < 4:
            arrow(ax, (.65, y-.003), (.65, y-.05))
    top = raw["landmarks_and_shape_proofs"]["E0"]["domain"]["leaf_upper_closed_at_zero_vpd"]
    save_figure(fig, output, "02_scalar_elimination", "Exact equilibrium elimination, not one-dimensional dynamics",
                f"Synthetic unsafe baseline; root=0 MPa; 0<L<={top:.4f} kg C/m2; gx,gs>0.\n"
                "L=0 is excluded. v(kPa)=1000 p_atm D; p_atm=0.101325 MPa. Formal closed/zero-carbon sets are treated separately.")


def equilibrium_figures(output, rows, raw):
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 5.1), sharey=True)
    for ax, (label, color) in zip(axes, COLORS.items()):
        selected = [r for r in rows if r["case"] == label]
        ax.plot(column(selected, "vpd_mpa")*1000, column(selected, "L"), color=color, lw=2)
        landmark = raw["landmarks_and_shape_proofs"][label]
        ax.plot(1000*landmark["zero_leaf_vpd_limit_mpa"], 0., "o", mfc="white", mec=color, ms=7)
        if landmark["fold_leaf"] is not None:
            v = 1000*landmark["fold_vpd_mpa"]
            ax.plot(v, landmark["fold_leaf"], "o", color=color)
            ax.annotate(f"fold: {v:.4g} kPa", (v, landmark["fold_leaf"]),
                        xytext=(v+.18, landmark["fold_leaf"]+.06), fontsize=10,
                        arrowprops={"arrowstyle": "->", "color": color})
        ax.axvline(2.44, color="#555555", ls="--", lw=1)
        ax.axvspan(3.167489, 5.05, color="#EAEAEA", zorder=-1)
        ax.set(xlabel="VPD v (kPa)", xlim=(0., 5.05), ylim=(-.012, .41),
               title=f"{label}: at most {'one' if label == 'E4' else 'two'} positive roots")
        ax.grid(alpha=.22)
        if label == "E1":
            ax.text(2.52, .28, "2.44 kPa:\nno positive root", fontsize=10)
        elif label == "E0":
            zoom = ax.inset_axes([.57, .55, .40, .33])
            zoom.plot(column(selected, "vpd_mpa")*1000, column(selected, "L"), color=color, lw=1.4)
            zoom.axvline(2.44, color="#555555", ls="--", lw=.9)
            zoom.plot(1000*landmark["zero_leaf_vpd_limit_mpa"], 0., "o", mfc="white", mec=color, ms=4)
            zoom.set(xlim=(2.395, 2.48), ylim=(-.003, .105), title="E0 two-root window")
            zoom.set_xticks([2.40, 2.44, 2.47])
            zoom.tick_params(labelsize=7)
            zoom.title.set_fontsize(9)
        elif label == "E4":
            ax.annotate("open zero-leaf limit\n4.807 kPa (formal)", (4.806749409, 0.),
                        xytext=(3.25, .13), fontsize=10, arrowprops={"arrowstyle": "->"})
    axes[0].set_ylabel("Equilibrium leaf carbon L (kg C/m2)")
    fig.legend(handles=[Line2D([0], [0], color="#555555", ls="--", label="representative v=2.44 kPa"),
                        Line2D([0], [0], marker="o", mfc="white", color="#555555", ls="none", label="L=0 excluded")],
               loc="upper center", bbox_to_anchor=(.5, .91), ncol=2, fontsize=10)
    save_figure(fig, output, "03_equilibrium_curves", "Scalar equilibrium curves: different feedbacks change the number of positive roots",
                "Data evaluated by the unchanged model. Synthetic unsafe baseline; root=0 MPa; positive L and gx,gs.\n"
                "Shading: v>3.167 kPa exceeds saturation at diagnostic 25 C (not a temperature response or mortality threshold).\n"
                "E0 inset uses the same units: v in kPa, L in kg C/m2.")
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 8.))
    keys = ["K_L", "K_W", "K_R", "K_S", "capacity_total", "capacity_norm"]
    labels = ["Leaf K_L", "Wood K_W", "Root K_R", "Soil K_S", "Total capacity sum(K_i)", "Capacity Euclidean norm"]
    for ax, key, title in zip(axes.ravel(), keys, labels):
        for label, color in COLORS.items():
            selected = [r for r in rows if r["case"] == label]
            values = column(selected, key)
            for status, style in [("stable", "-"), ("unstable", "--"), ("nonhyperbolic_within_tolerance", ":")]:
                mask = np.array([r["equilibrium_stability"] == status for r in selected])
                ax.plot(column(selected, "vpd_mpa")*1000, np.where(mask, values, np.nan),
                        color="grey" if status == "nonhyperbolic_within_tolerance" else color,
                        ls=style, lw=1.8)
                if status == "nonhyperbolic_within_tolerance":
                    ax.scatter(column(selected, "vpd_mpa")[mask]*1000, values[mask], color="grey", s=4)
        ax.axvspan(3.167489, 5.05, color="#EAEAEA", zorder=-1)
        ax.axvline(1000*raw["hopf_E4_vpd_mpa"], color="#777777", ls=":", lw=.9)
        ax.set(title=title, xlabel="VPD v (kPa)", ylabel="kg C/m2", xlim=(0., 5.05))
        ax.grid(alpha=.22)
    handles = [Line2D([0], [0], color=color, label=label) for label, color in COLORS.items()]
    handles += [Line2D([0], [0], color="black", ls="-", label="stable"),
                Line2D([0], [0], color="black", ls="--", label="unstable"),
                Line2D([0], [0], color="grey", ls=":", label="within spectral tolerance")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .935), ncol=6, fontsize=10)
    save_figure(fig, output, "04_capacity_branches", "K=C* at every equilibrium, including unstable equilibria",
                "Classification: full 6-state Jacobian, tolerance 1e-10/day; E4 stability is sampled, not interval-certified.\n"
                "Vertical dotted line: E4 local Hopf near 2.665 kPa. Shaded region is formal high VPD; extremely small L may be unresolved.\n"
                "Positive branches avoid conductance kinks; grey marks include folds and near-zero small-L spectral modes.")


def time_data(raw, trajectories, closed, model):
    full, tails = [], []
    minimum_leaf = raw["baseline_config"]["minimum_leaf"]
    matrix = carbon_matrix(model.carbon)
    for outcome_index, outcome in enumerate(raw["outcomes"]):
        case = f"{outcome['model']}:{outcome['case']}"
        selected = [r for r in trajectories if r["case"] == case]
        require(selected, f"Missing trajectory: {case}")
        require(max(column(selected, "day")) <= outcome["last_day"], "Numerical trajectory extended beyond termination")
        require(min(column(selected, "L")) >= minimum_leaf*(1-1e-8), "Numerical trajectory crossed L_min")
        require(float(selected[-1]["day"]) == outcome["last_day"], "Terminal point is missing")
        for row in selected:
            stocks = [float(row[c]) for c in ["L", "W", "R", "S"]]
            k = [float(row[f"K_{c}"]) for c in ["L", "W", "R", "S"]]
            full.append({**row, "model": outcome["model"], "initial_case": outcome["case"],
                         "year": float(row["day"])/YEAR_DAYS, "C_norm": float(np.linalg.norm(stocks)),
                         "K_norm": float(np.linalg.norm(k)), "phase": "full_model_until_termination",
                         "source_file": "raw_scalar/trajectories.csv"})
        if outcome["case"] != "both_closed_initial":
            continue
        exact_event = expm(matrix*outcome["last_day"])@np.array(outcome["initial"][2:])
        require(np.max(abs(exact_event-outcome["final_state"][2:])) < 1e-8,
                "Exact closed extension disagrees at the terminal point")
        tail = [{"model": outcome["model"], "day": outcome["last_day"],
                 **dict(zip(["L", "W", "R", "S"], exact_event)),
                 "stock_total": float(sum(exact_event)), "capacity_total": 0.,
                 "in_integrator_leaf_domain": bool(exact_event[0] > minimum_leaf),
                 "interpretation": "analytic_carbon_subsystem_limit_not_post_exit_full_model_simulation",
                 "source_file": f"raw_scalar/results.json#/outcomes/{outcome_index}/initial + expm(M*t_exit)",
                 "source_csv_line": ""}]
        tail.extend({**r, "source_file": "raw_scalar/analytic_closed_carbon.csv"} for r in closed
                    if r["model"] == outcome["model"] and float(r["day"]) > outcome["last_day"])
        for row in tail:
            stocks = [float(row[c]) for c in ["L", "W", "R", "S"]]
            tails.append({**row, "initial_case": "both_closed_initial", "year": float(row["day"])/YEAR_DAYS,
                          "C_norm": float(np.linalg.norm(stocks)), "K_norm": 0., "GPP": 0.,
                          "phase": "exact_carbon_subsystem_only_not_post_exit_full_model"})
    require(len({r["model"] for r in tails}) == 2, "Both E0 and E4 need a closed-carbon extension")
    return full, tails


def shutdown_figure(output, full, tails, raw):
    fig, axes = plt.subplots(2, 5, figsize=(17.5, 8.2))
    titles = ["Xylem", "Stomata", "GPP", "K and C", "Leaf carbon"]
    units = ["gx (m water/day/MPa)", "gs (m water/day)", "GPP (kg C/m2/day)", "Euclidean norm (kg C/m2)", "L (kg C/m2)"]
    for i, label in enumerate(["E0", "E4"]):
        for case, color in zip(CASES, CASE_COLORS):
            rows = [r for r in full if r["model"] == label and r["initial_case"] == case]
            times = column(rows, "year")
            for j, key in enumerate(["gx", "gs", "GPP", "K_norm", "L"]):
                axes[i, j].plot(times, column(rows, key), color=color, lw=1.65)
            axes[i, 3].plot(times, column(rows, "C_norm"), color=color, ls=":", lw=1.7)
            outcome = next(o for o in raw["outcomes"] if o["model"] == label and o["case"] == case)
            if outcome["status"] == "vegetation_domain_exit":
                for ax, key in [(axes[i, 3], "C_norm"), (axes[i, 4], "L")]:
                    ax.plot(times[-1], float(rows[-1][key]), "o", color=color, ms=4)
        extension = [r for r in tails if r["model"] == label and float(r["year"]) <= 100]
        axes[i, 3].plot(column(extension, "year"), column(extension, "C_norm"),
                        color=CASE_COLORS[2], ls="--", lw=1.8)
        axes[i, 3].plot(column(extension, "year"), column(extension, "K_norm"),
                        color=CASE_COLORS[2], ls="--", lw=1.2)
        axes[i, 4].plot(column(extension, "year"), column(extension, "L"),
                        color=CASE_COLORS[2], ls="--", lw=1.8)
        for j, ax in enumerate(axes[i]):
            ax.set_xscale("symlog", linthresh=.01/YEAR_DAYS)
            ax.set(xlim=(0., 100), xlabel="Time (year)", ylabel=units[j],
                   title=f"{label}: {titles[j]}")
            ax.set_xticks([0., 1e-4, 1e-2, 1., 100.])
            ax.set_xticklabels(["0", r"$10^{-4}$", r"$10^{-2}$", "1", "100"], fontsize=8)
            ax.tick_params(axis="y", labelsize=9)
            ax.xaxis.label.set_fontsize(9)
            ax.yaxis.label.set_fontsize(9)
            ax.title.set_fontsize(11)
            ax.grid(alpha=.2)
        axes[i, 2].set_yscale("symlog", linthresh=1e-10)
        axes[i, 2].set_yticks([0., 1e-9, 1e-6, 1e-3])
        axes[i, 3].set_yscale("symlog", linthresh=1e-6)
        axes[i, 3].set_ylim(0., 30)
        axes[i, 3].set_yticks([0., 1e-4, 1e-2, 1., 10.])
        for ax in [axes[i, 2], axes[i, 3]]:
            ax.spines["bottom"].set_zorder(0)
        axes[i, 4].set_yscale("log")
        axes[i, 4].set_ylim(1e-8, 1.)
        axes[i, 4].axhline(raw["baseline_config"]["minimum_leaf"], color="#555555", ls="-.", lw=1)
    handles = [Line2D([0], [0], color=c, label=n) for c, n in zip(CASE_COLORS, CASE_LABELS)]
    handles += [Line2D([0], [0], color="black", ls="-", label="K / full-model quantities"),
                Line2D([0], [0], color="black", ls=":", label="full-model C"),
                Line2D([0], [0], color=CASE_COLORS[2], ls="--", label="exact carbon-only tail; both closed only"),
                Line2D([0], [0], color="black", marker="o", ls="none", label="L_min exit; not mortality")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .94), ncol=4, fontsize=9)
    save_figure(fig, output, "05_shutdown_vs_stock", "At fixed 2.44 kPa, capacity can disappear before carbon stocks",
                "Synthetic unsafe baseline; root=0 MPa; unchanged log-carbon Radau. All numerical trajectories stop at L_min=1e-6 kg C/m2.\n"
                "Dashed tails: exp(Mt)C0 for both-closed cases only, NOT post-exit full-system simulation. Other cases are not extrapolated.\n"
                "Time, GPP and norm axes are symlog: exact zero is retained. Norms are unscaled 4-pool Euclidean norms.\n"
                "Leaf is on a log axis; carbon-only leaf tails below 1e-8 are outside the plot.")


def audited_values(raw, hopf):
    rows = []
    def add(key, value, unit, source, evidence):
        rows.append({"key": key, "value": value, "unit": unit, "source": source, "evidence": evidence})
    for label, landmark in raw["landmarks_and_shape_proofs"].items():
        add(f"{label}_zero_leaf_limit", 1000*landmark["zero_leaf_vpd_limit_mpa"], "kPa",
            f"raw_scalar/results.json#/landmarks_and_shape_proofs/{label}/zero_leaf_vpd_limit_mpa *1000", "analytic shape; numerical constant")
        if landmark["fold_leaf"] is not None:
            add(f"{label}_fold_vpd", 1000*landmark["fold_vpd_mpa"], "kPa",
                f"raw_scalar/results.json#/landmarks_and_shape_proofs/{label}/fold_vpd_mpa *1000", "analytic unique fold; numerical location")
            add(f"{label}_fold_L", landmark["fold_leaf"], "kg C/m2",
                f"raw_scalar/results.json#/landmarks_and_shape_proofs/{label}/fold_leaf", "numerical location")
    for key in ["vpd_mpa", "l1", "omega"]:
        name = "vpd_kpa" if key == "vpd_mpa" else key
        conversion = " *1000" if key == "vpd_mpa" else ""
        add(f"E4_hopf_{name}", hopf[key]*1000 if key == "vpd_mpa" else hopf[key],
            {"vpd_mpa": "kPa", "l1": "coordinate-dependent", "omega": "1/day"}[key],
            f"07_provenance.json#/archive_audit/E4/current_hopf_recalculation/{key}{conversion}", "local numerical verification")
    add("E4_hopf_L", hopf["state"][2], "kg C/m2",
        "07_provenance.json#/archive_audit/E4/current_hopf_recalculation/state/2", "local numerical verification")
    add("E4_second_spectral_crossing", 1000*raw["E4_detected_stability_crossings"][0]["vpd_mpa"], "kPa",
        "raw_scalar/results.json#/E4_detected_stability_crossings/0/vpd_mpa *1000", "formal high VPD; nonlinear type unresolved")
    for name, value in zip(["L", "W", "R", "S"], raw["stock_ratios"]):
        add(f"stock_ratio_{name}", value, "dimensionless", f"raw_scalar/results.json#/stock_ratios/{['L','W','R','S'].index(name)}", "baseline algebra")
    add("L_wet", raw["landmarks_and_shape_proofs"]["E0"]["domain"]["leaf_upper_closed_at_zero_vpd"], "kg C/m2",
        "raw_scalar/results.json#/landmarks_and_shape_proofs/E0/domain/leaf_upper_closed_at_zero_vpd", "positive domain; wet endpoint included")
    add("L_min", raw["baseline_config"]["minimum_leaf"], "kg C/m2", "raw_scalar/results.json#/baseline_config/minimum_leaf", "numerical exit, not zero carbon")
    add("xylem_closure_psi", 2*raw["baseline_config"]["hydraulic"]["psi50_x"], "MPa",
        "raw_scalar/results.json#/baseline_config/hydraulic/psi50_x *2", "baseline algebra")
    add("stomatal_closure_psi", 2*raw["baseline_config"]["hydraulic"]["psi50_s"], "MPa",
        "raw_scalar/results.json#/baseline_config/hydraulic/psi50_s *2", "baseline algebra")
    model = model_from_evidence(raw)
    sigma = model.carbon.losses[0]-model.carbon.allocation[0]*model.carbon.assimilation_max*model.carbon.canopy_coefficient/4
    add("closed_xylem_leaf_decay_bound", sigma*YEAR_DAYS, "1/year",
        "baseline carbon coefficients: (kL-bL*p*q/4)*365.25", "baseline analytical inequality")
    add("representative_vpd", 1000*raw["representative_vpd_mpa"], "kPa", "raw_scalar/results.json#/representative_vpd_mpa *1000", "fixed forcing")
    for i, outcome in enumerate(raw["outcomes"]):
        prefix = f"{outcome['model']}:{outcome['case']}"
        add(prefix+":last_year", outcome["last_day"]/YEAR_DAYS, "year",
            f"raw_scalar/results.json#/outcomes/{i}/last_day /365.25", outcome["status"])
        add(prefix+":final_gs", outcome["final_gs"], "m water/day",
            f"raw_scalar/results.json#/outcomes/{i}/final_gs", "numerical; zero is not proof of finite-time exact closure")
    return rows


def write_notes(output, deck, raw, facts):
    with zipfile.ZipFile(deck) as archive:
        paths = sorted([p for p in archive.namelist() if re.fullmatch(r"ppt/slides/slide\d+\.xml", p)],
                       key=lambda p: int(re.search(r"slide(\d+)", p).group(1)))
        require(len(paths) == len(SLIDES), "Slide count changed; update the explicit slide mapping first")
        original = {}
        for i, name in enumerate(paths, 1):
            elements = ET.fromstring(archive.read(name)).findall(".//{http://schemas.openxmlformats.org/drawingml/2006/main}t")
            original[i] = " | ".join(e.text for e in elements if e.text)
    for key, _, _, english, body in CLAIMS:
        require(200 <= len(body) <= 400, f"Evidence paragraph {key} has {len(body)} characters, expected 200--400")
        require(len(english.split()) <= 20, f"English claim {key} exceeds 20 words")
    slide_rows = []
    for number, title, english, figures, claims, note in SLIDES:
        require(len(english.split()) <= 20, f"Slide {number} conclusion exceeds 20 words")
        slide_rows.append({"original_slide": number, "section": title, "english_conclusion": english,
                           "english_word_count": len(english.split()), "suggested_figure_numbers": figures,
                           "claim_keys": claims, "speaker_note_ja": note, "original_text": original[number]})
    write_csv(output/"data/slide_notes.csv", slide_rows)
    write_csv(output/"data/audited_values.csv", facts)
    initial_rows = [{"key": f"{r['model']}:{r['case']}", "model": r["model"], "case": r["case"],
                     **dict(zip(["x_mpa", "leaf_mpa", "L", "W", "R", "S"], r["initial"])),
                     "source": f"raw_scalar/results.json#/outcomes/{i}/initial"}
                    for i, r in enumerate(raw["outcomes"])]
    write_csv(output/"data/initial_states.csv", initial_rows)
    text = [
        "# 10/16発表用：主張と根拠の整理", "",
        "capacityを計算できることと、自己整合的な平衡が存在することを分けて説明できる。",
        "平衡でK=C*でも、全結合系が安定とは限らないことを結果図で示せる。",
        "閉鎖時の即時capacityゼロと遅いストック減衰を、死亡の主張なしに区別できる。", "",
        "## 前提と監査", "",
        f"- 基準HEAD：`{BASELINE}`。PPT原本も改訂版も編集していない。",
        "- 係数は合成。root ψ=0固定、VPDのみforcing。温度は25°Cの湿度診断のみ。",
        "- 発表図のVPDはkPa、内部D=v(MPa)/0.101325。炭素kg C/m²、元solver時間day。",
        "- 正通水域は0<L≤L_wet≈.3840 kg C/m²。L=0除外、L_min=1e−6は積分退出境界。",
        "- 旧scalar/capacity/E4の入力15/18/10ハッシュは現在ソースと一致。E4 Hopf/l1は再計算。",
        "- 周期collocation、Floquet、非自律ramp全体は今回再実行していない。旧結果のハッシュ監査と関連テストを用いた。",
        "- raw JSONのleaf_xtol=1e−14は表記が広すぎる。実際の正根探索はfloatの最小正規化値、fold探索は1e−14。",
        "  元runnerは変更せず、この差を[再現記録](07_provenance.json)に残した。", "",
        "## 主張：各段落200〜400字", "",
    ]
    evidence_sources = {
        "A1": ("spac_carbon_capacity_definition_ledger.md", "capacity counterexamples.csv: case + source_csv_line in provenance"),
        "A2": ("spac_carbon_scalar_equilibrium_results.md", "stock_ratio_*, L_wet"),
        "A3": ("spac_carbon_scalar_equilibrium_results.md", "E0_fold_vpd, E1_fold_vpd, E4_zero_leaf_limit"),
        "A4": ("spac_carbon_scalar_equilibrium_results.md", "04_capacity_branches.csv: model + vpd_kpa + L + source_csv_line"),
        "N1": ("spac_carbon_stage_b3_e4_results.md", "E4_hopf_*; archived E4 results.json#/cycles"),
        "A5": ("spac_carbon_scalar_equilibrium_results.md", "05_carbon_only_extension.csv: model + day + phase"),
        "A6": ("spac_carbon_scalar_equilibrium_results.md", "closed_xylem_leaf_decay_bound, E4:xylem_closed_low_leaf:final_gs"),
        "N2": ("spac_carbon_scalar_equilibrium_results.md", "*:last_year; initial_states.csv key"),
        "D1": ("spac_carbon_scalar_equilibrium_results.md", "L_min, E4_second_spectral_crossing, E4_zero_leaf_limit"),
        "U1": ("spac_carbon_stage_c23_results.md", "local/current scope only; no new quantitative global claim"),
    }
    for key, category, title, english, body in CLAIMS:
        report, sources = evidence_sources[key]
        text += [f"### {key} — {category}：{title}", "", body, "",
                 f"**English conclusion ({len(english.split())} words):** {english}", "",
                 f"根拠：[研究ノート](../../../docs/{report})。数値・CSV行キー：`{sources}`。",
                 f"段落文字数：{len(body)}。", ""]
    text += ["## 監査済み数値と行キー", "",
             "発表では3〜4桁を目安にする。厳密値とJSON pointerは[data/audited_values.csv](data/audited_values.csv)。", "",
             "| key | value | unit | 根拠の範囲 |", "|---|---:|---|---|"]
    text += [f"| `{r['key']}` | {r['value']:.12g} | {r['unit']} | {r['evidence']} |" for r in facts]
    text += ["", "## PPT原本21枚への対応", "",
             "対象：[ProgressReport1016.pptx](../../../docs/ProgressReport1016.pptx)。改訂版の番号対応は別途確認。",
             "原文・話す内容・一文結論の機械可読版は[data/slide_notes.csv](data/slide_notes.csv)。", ""]
    for row in slide_rows:
        text += [f"### Slide {row['original_slide']}: {row['section']}",
                 f"**{row['english_conclusion']}**",
                 f"図：{row['suggested_figure_numbers'] or '文章のみ'}。主張：`{row['claim_keys']}`。",
                 row["speaker_note_ja"], ""]
    text += ["## 一枚一結論にするための追加案", "",
             "- **19B: Shutdown versus stock** — A5の英語結論をそのまま使い、図05でKゼロとC減衰の時間差を示す。",
             "- **任意の19C: Xylem-only closure** — A6の英語結論を使い、図05の青い小葉量例を指して気孔開度が正でも葉量が減ることを説明する。",
             "- 元21枚の番号対応は維持したメモであり、PPTへページを追加したわけではない。", ""]
    text += ["## 再解釈：短い言い換え", "",
             "- capacityの定義が消えるのではない。Kの値がゼロになる場合と、正の自己整合的候補が消える場合を分ける。",
             "- 複数capacityは複数の自己整合的全平衡であって、瞬間写像が多価という意味ではない。",
             "- 一意性だけでは、安定性、basin内の初期状態、変動外力への追従を保証しない。",
             "- hydraulic failureは今回のモデル内の低通水・衰退の数学的挙動。死亡・再生生理の実証ではない。",
             "- 一般的新規性には文献比較が必要。今回の価値は存在、全結合安定性、閉鎖の不変性を分けたこと。", "",
             "## 図の注意", "",
             "- 01は模式図。02は既存式からの厳密な平衡消去。03〜05は再生成データからの描画。",
             "- 04の灰色は固有値が許容値内。真の非双曲性の証明ではない。E4の全枝安定性は標本評価。",
             "- 05はK/Cの4プールEuclideanノルム。総量sumとは別。GPP・ノルムのsymlog軸は厳密ゼロと小葉量例の正入力を保持する。",
             "- 05の破線延長は両閉鎖だけ。木部のみ閉鎖例にはGPPが残るので線形延長を適用しない。",
             "- 元スライドのCox論文のclosed-system解釈と、広い新規性は原文比較を要する。",
             "- ignoredな成果物はローカルで存在確認済み。GitHubから閲覧できるとは約束しない。", ""]
    (output/"06_summary_evidence.md").write_text("\n".join(text), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="New unchanged-scalar runner output")
    parser.add_argument("--output", type=Path, required=True, help="Bundle root; numbered assets must not exist")
    parser.add_argument("--test-log", type=Path, required=True, help="Log from the declared four-file pytest command")
    args = parser.parse_args()
    output, source = args.output.resolve(), args.source.resolve()
    require(source == output/"raw_scalar", "Keep raw evidence inside the bundle as raw_scalar")
    require(not any(output.glob("0*")), "Numbered presentation assets already exist; use a new bundle")
    require(not (output/"data").exists(), "Presentation data already exists; do not overwrite")
    for key, _, _, english, body in CLAIMS:
        require(200 <= len(body) <= 400 and len(english.split()) <= 20, f"Evidence text constraints failed: {key}, {len(body)} characters")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    baseline_diff = subprocess.check_output(["git", "--no-pager", "diff", "--name-only", BASELINE, "HEAD"],
                                            cwd=ROOT, text=True).splitlines()
    worktree_status = subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True)
    raw = json.loads((source/"results.json").read_text())
    for name, expected in raw["source_sha256"].items():
        require(digest(Path(name)) == expected, f"Fresh scalar input changed: {name}")
    old, checks, inputs, model = audit(raw)
    deck_paths = [ROOT/"docs/ProgressReport1016.pptx", ROOT/"docs/ProgressReport1016_revised.pptx"]
    deck_hashes = {str(p): digest(p) for p in deck_paths if p.exists()}
    test_text = args.test_log.read_text()
    match = re.search(r"(\d+) passed in ([\d.]+)s", test_text)
    require(match is not None and "failed" not in test_text.lower() and "error" not in test_text.lower(),
            "Current pytest log does not show a clean successful run")
    curves = read_rows(source/"leaf_curves.csv")
    trajectories = read_rows(source/"trajectories.csv")
    closed = read_rows(source/"analytic_closed_carbon.csv")
    for row in curves:
        k = np.array([float(row[f"K_{c}"]) for c in ["L", "W", "R", "S"]])
        require(float(row["full_residual_norm_scaled_per_day"]) < 1e-8, "Curve row fails full equilibrium residual")
        require(np.max(abs(k-np.array([float(row[c]) for c in ["L", "W", "R", "S"]]))) < 1e-12,
                "Capacity does not equal equilibrium stocks")
        row.update({"vpd_kpa": float(row["vpd_mpa"])*1000, "capacity_total": float(sum(k)),
                    "capacity_norm": float(np.linalg.norm(k)), "source_file": "raw_scalar/leaf_curves.csv"})
    full, tails = time_data(raw, trajectories, closed, model)
    facts = audited_values(raw, checks["E4"]["current_hopf_recalculation"])
    fact_map = {row["key"]: row["value"] for row in facts}
    for key, expected in [("E0_fold_vpd", 2.472042553), ("E1_fold_vpd", 1.546784415),
                          ("E4_zero_leaf_limit", 4.806749409), ("E4_hopf_vpd_kpa", 2.664752905)]:
        require(abs(fact_map[key]-expected) < 1e-8, f"Slide handoff numerical audit failed: {key}")
    (output/"data").mkdir()
    write_csv(output/"data/03_equilibrium_curves.csv", curves)
    write_csv(output/"data/04_capacity_branches.csv", curves)
    write_csv(output/"data/05_shutdown_vs_stock.csv", full)
    write_csv(output/"data/05_carbon_only_extension.csv", tails)
    shutil.copyfile(ARCHIVES["capacity"].parent/"counterexamples.csv", output/"data/capacity_counterexamples.csv")
    if args.test_log.resolve() != output/"validation.log":
        shutil.copyfile(args.test_log, output/"validation.log")
    inputs.update({str(p): digest(p) for p in source.iterdir() if p.is_file()})
    inputs.update({str(ROOT/"docs"/name): digest(ROOT/"docs"/name) for name in REPORTS})
    inputs.update(deck_hashes)
    for path in [Path(__file__).resolve(), ROOT/"docs/slide_agent_handoff_1016.md",
                 *[ROOT/"tests"/name for name in TEST_FILES], ROOT/"src/control_carbon/spac_carbon_tracking.py"]:
        inputs[str(path)] = digest(path)
    inputs[str(output/"validation.log")] = digest(output/"validation.log")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False})
    conceptual_figures(output, raw)
    equilibrium_figures(output, curves, raw)
    shutdown_figure(output, full, tails, raw)
    write_notes(output, deck_paths[0], raw, facts)
    relative_output = output.relative_to(ROOT).as_posix()
    tests = " ".join("tests/"+name for name in TEST_FILES)
    renderer = Path(__file__).relative_to(ROOT).as_posix()
    commands = [
        f"python examples/analyze_spac_carbon_equilibria.py --output {relative_output}/raw_scalar",
        f"python -m pytest -q {tests} > {relative_output}/validation.log",
        f"python {renderer} --source {relative_output}/raw_scalar --output {relative_output} --test-log {relative_output}/validation.log",
    ]
    readme = [
        "# 10/16 presentation assets", "",
        "PPT本体は変更せず、既存式の再計算から作った貼り込み用素材。各図は300 dpi PNGとPDF。",
        f"HEAD: `{head}`。基準HEADとの変更ファイルはprovenanceに記録。合成係数・root ψ=0・固定VPD。", "",
        "## Files", "",
        "- 01_concept_capacity: 模式図。固定線形capacityと状態依存写像。",
        "- 02_scalar_elimination: 正通水平衡の6状態消去。",
        "- 03_equilibrium_curves: E0/E1/E4のモデル評価曲線。",
        "- 04_capacity_branches: 全6状態スペクトルで分けたK=C*の各成分・総量・ノルム。",
        "- 05_shutdown_vs_stock: 8例のgx/gs/GPP/||K||/||C||/Lと両閉鎖例のみの解析的延長。",
        "- [06_summary_evidence.md](06_summary_evidence.md): 200〜400字の根拠別説明、原本21枚の一文結論・注意。",
        "- [07_provenance.json](07_provenance.json): hash、単位、solver、監査・限定、出力一覧。",
        "- data/: 図データ、8初期値、監査済み数値と元CSV行キー、スライド対応。raw_scalar/: 元runnerの12成果物。",
        "- validation.log: 今回実行した関連4ファイルのpytestログ。", "",
        "## Reproduce", "",
        "リポジトリrootで実行。下記の出力先を新規パスに置き換え、過去成果物は上書きしない。",
        "rendererだけに対して全モデルの再計算を期待しない。必ず元runner、試験、rendererの順。", "",
        "```bash", *commands, "```", "",
        f"今回は `{match.group(1)} passed`。旧報告の236を今回の実行結果として転記していない。",
        "logをこのbundle自身のvalidation.logに置いた再実行では、ファイルコピーは不要。",
        "VPD=kPa（図）、D=mol/mol（内部）、ψ=MPa、C/K=kg C/m²、time=day（solver）/year（時系列図）。",
        "0<L≤.3840 kg C/m²、conductance正の平衡枝。L=0除外。退出はL_min、死亡ではない。",
        "元runner JSONのleaf_xtolは古い包括表記。正根xtolとfold xtolはprovenanceで区別。",
        "ignoredな画像類のGitHub閲覧を保証しない。PPTへコピーする際はこのローカルbundleを使用する。", "",
    ]
    (output/"README.md").write_text("\n".join(readme), encoding="utf-8")
    for name, expected in deck_hashes.items():
        require(digest(Path(name)) == expected, f"User PowerPoint changed during asset preparation: {name}")
    for name, expected in inputs.items():
        require(digest(Path(name)) == expected, f"Evidence input changed during rendering: {name}")
    files = sorted(p for p in output.rglob("*") if p.is_file())
    listing = [{"path": p.relative_to(output).as_posix(), "bytes": p.stat().st_size, "sha256": digest(p)} for p in files]
    listing.append({"path": "07_provenance.json", "sha256": None, "policy": "no self-referential checksum"})
    ignored = subprocess.run(["git", "check-ignore", str(output/"03_equilibrium_curves.png")], cwd=ROOT,
                             capture_output=True, text=True)
    require(ignored.returncode in [0, 1], f"git check-ignore failed: {ignored.stderr}")
    provenance = {
        "generated_at": datetime.now().astimezone().isoformat(), "HEAD": head, "requested_baseline": BASELINE,
        "baseline_diff": baseline_diff, "worktree_status_at_start": worktree_status,
        "source_sha256": inputs, "archive_audit": checks, "commands": commands,
        "script": str(Path(__file__).resolve()), "settings": raw["settings"], "baseline_config": raw["baseline_config"],
        "solver": {**raw["solver_settings"], "positive_leaf_root_xtol": float(np.finfo(float).tiny),
                   "fold_leaf_xtol": 1e-14, "E4_spectral_crossing_leaf_xtol": 1e-14,
                   "E4_hopf_leaf_potential_xtol_mpa": 1e-13,
                   "scalar_brentq_rtol": float(4*np.finfo(float).eps),
                   "metadata_warning": "raw leaf_xtol=1e-14 does not describe positive_equilibria; original runner preserved"},
        "tests": {"command": f"python -m pytest -q {tests}", "files": TEST_FILES, "passed": int(match.group(1)),
                  "elapsed_seconds": float(match.group(2)), "log": "validation.log", "scope": "four targeted existing files, not full repository suite"},
        "axes": {"vpd": "kPa", "internal_D": "mol/mol=v(MPa)/0.101325", "water_potential": "MPa",
                 "conductance_xylem": "m water/day/MPa", "conductance_stomatal": "m water/day",
                 "GPP": "kg C/m2/day", "carbon": "kg C/m2", "time_plot": "year=day/365.25",
                 "norm": "unscaled Euclidean 2-norm of four carbon pools, not total carbon"},
        "data_ranges": {"vpd_kpa": [float(min(column(curves, "vpd_mpa"))*1000), float(max(column(curves, "vpd_mpa"))*1000)],
                        "curve_L_kgC_m2": [float(min(column(curves, "L"))), float(max(column(curves, "L")))],
                        "full_model_year": [0., 100.], "carbon_only_csv_year": [float(min(column(tails, "year"))), float(max(column(tails, "year")))],
                        "time_figure_year": [0., 100.], "leaf_axis_kgC_m2": [1e-8, 1.],
                        "norm_axis_kgC_m2": [0., 30.], "norm_scale": {"scale": "symlog", "linthresh_kgC_m2": 1e-6},
                        "GPP_scale": {"scale": "symlog", "linthresh_kgC_m2_day": 1e-10},
                        "zero_display": "GPP/norm symlog axes retain zero; leaf log axis omits tails below its positive lower bound"},
        "definitions": {"positive_leaf_domain": "0<L<=L_wet; L=0 excluded; both conductances positive",
                        "minimum_leaf": raw["baseline_config"]["minimum_leaf"],
                        "extension": "both-closed initial states only; carbon-subsystem exp(Mt)C0, not full-model continuation"},
        "initial_states": raw["outcomes"], "numeric_claim_sources": "data/audited_values.csv",
        "current_checks": {"curve_scaled_full_residual_max_per_day": float(max(column(curves, "full_residual_norm_scaled_per_day"))),
                           "curves": len(curves), "full_trajectory_rows": len(full), "carbon_only_tail_rows": len(tails),
                           "figure_pairs": 5, "evidence_character_counts": {key: len(body) for key, _, _, _, body in CLAIMS}},
        "unverified": ["global basin/recovery boundary", "autonomous stable-cycle attractor", "true R-tipping",
                       "new physiological mortality regime", "E4 continuous full-branch stability certificate",
                       "second high-VPD crossing nonlinear type", "coefficient calibration", "broad literature novelty",
                       "Cox closed-Earth interpretation in original slide 5"],
        "preserved_powerpoint_sha256": deck_hashes, "artifacts_gitignored": ignored.returncode == 0,
        "output_files": listing,
    }
    (output/"07_provenance.json").write_text(json.dumps(provenance, default=encode, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(f"Saved five 300-dpi PNG/PDF pairs and evidence bundle to {output}")
    print(f"Archive hashes 15/18/10 matched; current E4 l1={checks['E4']['current_hopf_recalculation']['l1']:.9g}; {match.group(1)} current tests passed.")


if __name__ == "__main__":
    main()
