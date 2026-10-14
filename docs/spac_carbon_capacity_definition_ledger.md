# Stage C：capacity定義・仮定・実装の台帳

2026-10-13。[研究計画](SPAC_nonlinear_carbon_storage_capacity_research_plan.md)の
WP-C0/C1。基準コミットは `2ef9438c080d8183edc1efeff6e3a3b9fad62eb1`。
既存Stage A/Bの方程式・設定・報告・成果物は変更しない。

## 1. 最初に固定する符号と単位

炭素状態は `C=(L,W,R,S)`、水理状態は `z=(psi_x,psi_leaf)`、
全状態の実装順は `y=(z,C)`。光合成GPPをここでは **A** と書き、
potentialの **P** と混同しない。

\[
\dot C=u+MC,\quad K=-M^{-1}u,\quad
P=K-C=-M^{-1}\dot C,\quad \dot C=M(C-K).
\]

Luo型のcapacityはK、potentialはP。この符号は
[既存の理論報告](capacity_qse_attractor_exploration.tex)と一致する。
本研究はLuoが状態依存を禁止したと主張するものではない。

| 量 | 単位 |
|---|---|
| C、K、P、条件付き炭素根 | kg C/m²地表面積 |
| u、GPP A、炭素RHS | kg C/m²/day |
| M、Jacobian固有値 | day⁻¹（混合単位のJacobian成分は別） |
| z | MPa、非正 |
| 外向けVPD v | MPa、非負、内部D=v/p_atm |
| 時間、p_atm | day、0.101325 MPa |
| 水フラックス、容量 | m水/day、m水/MPa |

25℃はRH整合性診断のみ。新規実験の主範囲は0～0.003 MPa。
旧速度経路の終点0.00303975 MPaは主範囲を超える履歴実験として保全する。

## 2. 五種類の対象を混ぜない

| 種類 | 定義 | 必要な確認 | 現在の根拠 |
|---|---|---|---|
| I：瞬間capacity写像 K(y,v) | 現在の入力・係数を凍結した値 | M可逆、uの定義と符号 | 恒等式は【解析的証明】 |
| II：水理条件付き炭素根 C*(z,v) | C=K、zを固定 | 全条件付き根、炭素方向の安定性 | 凹canopyのR0分類は【解析的証明】 |
| III：全結合凍結平衡 y*(v) | C=KかつH=0 | 全Jacobian、枝・領域の区別 | 根とスペクトルは数値；限定E2/E3枝は区間保証 |
| IV：凍結系の吸引集合 A(v) | 安定平衡・安定周期解など | 吸引性、吸引域 | E4の既知周期解は不安定、IVの実例ではない |
| V：非自律実軌道 y(t) | v(t)で駆動される解 | 枝の移動速度、誤差、吸引域 | I～IIIとの一致は保証されない |

条件付きの **ゼロ炭素根** は炭素方程式の数学的根である。
E0～E4のPressureSPACは正葉量を要求するため、それを全系に代入しない。
L=0での葉水容量・木部通水・死亡・再生の意味は未定義である。
「正根なし」「正根がL_min以下」「今回の枝区間に根なし」も別に記録する。

## 3. 実装位置と依存関係

| 量／処理 | 既存実装 | Stage Cでの利用 |
|---|---|---|
| u=bA、A=p h(1-exp(-qL)) | [assimilation / carbon_rhs](../src/control_carbon/spac_carbon.py) | 配分b=(bL,bW,bR,0)、GPP。FEを二重乗算しない |
| M | [carbon_matrix](../src/control_carbon/spac_carbon.py) | 固定・下三角Metzler。正のkにより可逆Hurwitz |
| K | [frozen_coefficient_capacity](../src/control_carbon/spac_carbon.py) | 既存solveを再利用し、inverseを直接作らない |
| 条件付き正根・R0 | [positive_carbon_equilibrium / reproduction_number](../src/control_carbon/spac_carbon.py) | ゼロ根も併記し、条件付き吸引性を分類 |
| H、全RHS、Jacobian | [two_way_rhs / two_way_jacobian](../src/control_carbon/spac_carbon_two_way.py) | E0～E4で同じ炭素式、feedbackは水フラックスのみ |
| v→D、平衡枝、積分 | [PressureSPAC](../src/control_carbon/spac_carbon_mpa.py) | 唯一のMPa変換入口。旧APIへMPaを直接渡さない |
| E4周期解 | [periodic_orbit](../src/control_carbon/spac_carbon_e4.py) | 能動4状態＋受動R/Sの周期応答。capacity診断を後付け |
| 正炭素の多根・potential増加反例 | [CubicCapacityFeedback](../src/control_carbon/capacity_feedback.py) | 数学的反例として再利用、森林モデルとは区別 |
| K、P、尺度付き誤差・残差 | [新診断器](../src/control_carbon/spac_carbon_capacity.py) | 元のRHSを変更しない純関数 |

既存E0～E4ではMとuにVPDの**直接依存はない**。
vは水理を介して間接的にGPPを変える。構造炭素はHへ戻る。
この診断器を、温度依存M・sigmoid・任意の状態依存Mへ無条件に転用しない。

## 4. 命題と証明の範囲

### 恒等式と固定点【解析的証明】

M可逆なら `M(C-K)=MC+u`。したがって `P=0` と `dot C=0` は同値だが、
`H=0` は含意しない。全系残差F=(H,dot C)を別に測る。
`1^T P=0` と `P=0` も同値でない。プール間の符号相殺がある。

### 平衡でのJacobian【解析的証明】

一般の滑らかな状態依存Mについて、方向deltaのcapacity微分は
\[
DK[\delta]=-M^{-1}\{DM[\delta]K+Du[\delta]\}.
\]
既存TeXのcapacity微分式の括弧内は、この**行列ベクトル積**として読む。
条件付き根ではC-K=0なので
\[
J_{CC,*}=M_*[I-D_CK|_*].
\]
現在の固定Mモデルでは、非平衡でも
`J_Cy=M(E_C-D_y K)`。E_Cは全状態から炭素を取り出す4×6行列。
実装では既存の解析Jacobianを再利用し、Kの独立中央差分で検査する。

MのHurwitz性、J_CCのHurwitz性、全JのHurwitz性、非自律追従性は別問題。
特にE4の正条件付き炭素根ではJ_CCが安定でも、全系はHopfの後で不安定になる。

### 条件付き根の完全分類【解析的証明】

`R0=bL p h q/kL`。
凹canopyによりゼロ根は常に存在し、R0>1で正根がちょうど一つ。
R0>1のゼロ根は不安定、正根は正葉量領域で吸引的。
R0<1ではゼロ根が指数安定、R0=1では非双曲的だが非負側で漸近安定。
これらは**水理を固定した**命題で、全系の吸引性ではない。

### potentialの限界【解析的証明】

`Pdot=DK F-dot C`（固定v・固定M）。
一般に `d(||P||²)/dt` の符号は定まらない。
スカラー反例では `d(P²/2)/dt=k(K'-1)P²`。
K'>1の点で正になる。単に総炭素量がcapacityと一致しても収束の証拠ではない。

### 追従と境界【解析的導出／仮説・未実施】

枝が滑らかでJ可逆なら `y*'(v)=-J^-1 F_v`、
`edot=F(y*+e,v)-F(y*,v)-y*'(v) vdot`。
一様な安定性、局所領域内の保持、吸引域からの余裕は別途必要。
foldやHopf近傍の長い時定数を、短期の見かけの回復で置き換えない。
E4の不安定周期解の安定多様体は**境界候補**であり、
回復／L_min退出を分ける大域境界という結論はまだない。

## 5. 再現・履歴監査

既存再現コマンドは以下（出力先は毎回新規パス）：

```bash
python examples/analyze_spac_carbon.py --output NEW_STAGE_A_PATH
python examples/analyze_spac_carbon_two_way.py --output NEW_B1_PATH
python examples/analyze_spac_carbon_bifurcation.py --output NEW_B2_PATH
python examples/review_spac_carbon_units.py --output NEW_PRESSURE_PATH
python examples/analyze_spac_carbon_e4.py --output NEW_B3_PATH
```

ローカル保存を確認したのは `stage_a_final_20261009`、
`stage_b1_final_20261009`、`stage_b2_verified_20261012`、
`pressure_review_final_v2_20261013`、`stage_b3_e4_final_v2_20261013`。
いずれも `artifacts/spac_carbon/` 内、Git無視下の成果物。
GitHubに掲載されているとは扱わない。

Stage Aのmanifest内7個、B1の7個、B3の10個の保存ソースhashは現行ソースと一致する。
B2とpressure reviewの保存hashは、B3で一般化された
`spac_carbon_bifurcation.py` **一つだけ**が不一致。
過去結果を改ざんせず、現行での回帰試験と区別する。
今回、旧campaign全部を再実行したとは報告しない。
Stage CのJSONには現行SHA・入力hash・旧成果物の監査結果・solver設定を保存する。
新規試験数、数値受入条件、成果物は[Stage C結果](spac_carbon_stage_c_results.md)へ記録する。
