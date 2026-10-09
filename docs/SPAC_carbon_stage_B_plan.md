# SPAC–炭素双方向結合 Stage B：研究計画・実装指針

2026-10-09。対象リポジトリ：`Sachitama2001/Control_carbon_model`、基準コミット `32acaacc3f115172694abd661af6b5763e958420`。本書は既存Stage Aの独立拡張であり、VISIT再現ではない。Leuning/FvCBの厳密な導入は後回しにする。

## 1. 研究目的と現状

目的：葉・木部・根の炭素貯留が水理流束を変え、水理状態が光合成入力を変える双方向結合のもとで、平衡点、安定性、hydraulic fold、carbon storage capacity、B-/R-tippingを調べる。

Stage Aは `y=(psi_x,psi_l,L,W,R,S)` の6状態。水理→炭素のみの片方向結合で、Jacobianの `J_HC=0`。したがって炭素量は水理foldを動かさない。Stage AのR0閾値・capacity写像・温度パルスのno-R結果は保全する。

## 2. Manzoniの2g_xの意味と拡張の注意

既存式：

```
E_rx = 2*g_x(psi_l)*(psi_r-psi_x)
E_xl = 2*g_x(psi_l)*(psi_x-psi_l)
E_atm = g_s(psi_l)*D
```

2つの等しい抵抗 `1/(2g_x)` の直列合成が `1/g_x` となる対称分割。等しい長さの経路という解釈は自然だが、**原著が葉量による経路長変化を導出したわけではない**。葉量が増えると葉脈の並列経路数・総通水断面積・蒸散面積が変わる可能性があり、経路長とは区別する。以下の倍率は本研究の仮説。

一般化：

```
E_rx = G_rx(psi_l,L,W,R)*(psi_r-psi_x)
E_xl = G_xl(psi_l,L,W,R)*(psi_x-psi_l)
E_atm = G_atm(psi_l,L)*D
C_x*psi_xdot = E_rx-E_xl
C_l*psi_ldot = E_xl-E_atm
```

### A：葉量→木部–葉conductance（ユーザー提案）

`G_rx=2g_x`、`G_xl=2g_x*F_L(L)`、`F_L(L_ref)=1`。最初は `F_L=(L+L_floor)/(L_ref+L_floor)`。葉量に伴う並列通水経路増加の仮説。`L_floor` は正則化であって生存閾値ではない。

定常時の直列合成は

```
G_eff=G_rx*G_xl/(G_rx+G_xl)=2*g_x*F_L/(1+F_L)
E_supply=G_eff*(psi_r-psi_l)
```

**`G_eff=g_x*F_L` ではない。** `F_L=1` なら厳密に `g_x` に戻る。最初の解析・単体テストとする。

### B：葉量→群落蒸散需要（必須対照）

`G_atm=g_s(psi_l)*F_E(L)`、第一候補 `F_E=LAI/LAI_ref=L/L_ref`。`LAI=SLA*L/f_C`。元の `g_s` は地表面積基準の基準群落conductanceと解釈する。葉面積基準conductanceとLAIを二重計上しない。葉量が増すと供給だけでなく蒸散需要も増えるため、**Aのみ・Bのみ・A+B**を比較する。飽和型 `F_E` は感度解析に回す。

### C：木部炭素→供給能力

`G_rx=2*g_x*F_W(W)`、例えば `F_W=(W/W_ref)^beta_W`。木質炭素量と辺材断面積は同一ではないため仮定を明記。A+B検証後に実装。

### D：根炭素→土壌–根conductance

根ポテンシャルを境界値から代数節点に変える：`E_sr=G_sr(R,psi_soil)*(psi_soil-psi_r)`、`E_sr=E_rx`。根に容量を置かなければDAEまたは代数消去。降水・土壌水分の動的状態は後続段階。

### E：炭素量→水理容量

`C_l(L)`、`C_x(W)` は候補だが、`C(L)*psi_dot=flux` だけでは組織増減時の水収支を破り得る。水貯留 `H_i(psi_i,C_i)` を導入して `dH_i/dt=H_psi*psi_dot+H_C*C_dot=net_flux` とする。流束feedbackを検証した後に実装。

### F：損傷・修復とNSC

embolism recovery、新規辺材形成、水理疲労は異なる過程。現段階では採用せず、根拠ある状態変数が用意できた段階で検討。

## 3. 最初に実装するStage B1

状態・炭素収支はStage Aを維持し、まず容量固定で流束だけを変える。

```
E_rx  = 2*g_x(psi_l)*F_W(W)*(psi_r-psi_x)
E_xl  = 2*g_x(psi_l)*F_L(L)*(psi_x-psi_l)
E_atm = g_s(psi_l)*F_E(L)*D
C_x*psi_xdot = E_rx-E_xl
C_l*psi_ldot = E_xl-E_atm
C_dot = M(T)*C + b*p(T)*h(psi_l)*(1-exp(-q*L))
h(psi_l)=g_s(psi_l)/g_s,max
```

最初は `F_W=1`。光合成の葉面積効果はすでに `1-exp(-qL)` に入るため、`F_E` をGPPへ再度掛けない。`L_ref,W_ref` はStage A湿潤正平衡で一度決め、シナリオ間で固定する。

モデル比較：E0=全Fが1、E1=F_Lのみ、E2=F_Eのみ、E3=両方、E4=さらにF_W。E1とE2でfoldが逆方向に動く可能性を含め、**効果の符号を先取りしない**。

## 4. 平衡点と全Jacobian

`G_rx,G_xl>0` の定常状態では

```
psi_x*=(G_rx*psi_r+G_xl*psi_l*)/(G_rx+G_xl)
G_eff(psi_l,L,W)*(psi_r-psi_l)=G_atm(psi_l,L)*D
```

固定した `psi_l,T` に対する炭素平衡のゼロ根と正根を調べ、`W*=b_W*P*/k_W` 等を代入して自己整合的な平衡を探す。複数初期値・分枝追跡が必須。

```
J = [[J_HH,J_HC],
     [J_CH,J_CC]]
```

Stage Bでは `J_HC != 0`。`J_CC` 可逆ならfold候補のゼロ固有値は `det(J_HH-J_HC*J_CC^-1*J_CH)=0` で調べられる。ただし単純ゼロ固有値・横断性・二次非退化性を確認しなければfoldの証明ではない。Stage AのR0は水理を凍結した炭素系の閾値であり、双方向全系の安定性を単独で決めない。

## 5. 作業パッケージ

**WP0：保全。** main SHAを固定。Stage Aのコード・テスト・結果を再実行。単位はMPa、day、kg C/m²、Dはmol/mol。全パラメータの出典と合成値を区別。

**WP1：流束feedback。** `spac_carbon_two_way.py` を独立作成。E0–E4を実装。F=1でStage AのRHS・Jacobian・平衡に一致すること。水収支 `C_x*psi_xdot+C_l*psi_ldot=E_rx-E_atm`、炭素収支 `d(sum C)/dt=GPP-Ra-Rh`、直列抵抗恒等式、J_HC非ゼロ、有限差分Jacobianをテスト。

**WP2：分枝・吸引域。** D・根ポテンシャル・feedback強度を走査。多初期値根探索、pseudo-arclength continuation、6×6固有値、foldとR0=1の区別、吸引域、capacity写像と実貯留の時間差を可視化。計算誤差と刻み依存を確認。

**WP3：ゼロ葉量の物理的定義。** `L=0` でも葉水ポテンシャルと固定容量が残る問題を明示。第一段階は `L>L_min` の植生存続領域に限定し、ゼロ平衡を実在する死亡状態と呼ばない。必要ならレジーム切替・保存則付き容量・残留通水を独立variantとして導入。線形Manzoniの完全非通水域は非孤立平衡集合である。

**WP4：tipping判定。** B-tippingは凍結系のfold横断。R-tippingは**全経路で安定凍結平衡が存続**することを確認した上で、同じ始終点・同じ経路で速度のみ変更し、最終外力で長く保持した吸引先を比較。起こらなければno-R条件を解析。温度とLeuningはその後。

**WP5：生理学的高度化。** 辺材・根のallometry、葉/木部容量、土壌水分、Leuning/FvCB/Tuzetの閉包、実測による較正を順に検討。

## 6. 最低限の実験・図・受入テスト

| ID | 比較 | 問い |
|---|---|---|
| E0 | Stage A | 厳密な基準線 |
| E1 | F_Lのみ | 葉量が供給を変えるか |
| E2 | F_Eのみ | 葉量が蒸散需要を変えるか |
| E3 | F_L+F_E | 供給・需要の競合でfoldはどう動くか |
| E4 | F_L+F_E+F_W | 木部炭素feedbackの追加効果 |
| E5 | D continuation | 分枝・固有値・capacityの変化 |
| E6 | 同一D・異なる初期値 | 吸引域と履歴依存 |
| E7 | 同一経路slow/fast | 真のR-tipping探索とB対照 |
| E8 | L_floor,基準値,指数の感度 | 仮定への頑健性 |

必須図：F関数とLAI、D–psi_l分枝、D–総炭素平衡、全Jacobian固有値、吸引域、capacityと実ストック、slow/fast最終保持。全図に単位・モデルvariant・合成パラメータを明記。

受入テスト：F=1極限、直列抵抗、水・炭素収支、非負性、単位変換、Jacobian差分、根残差・固有値、continuation収束、kink・ゼロ葉量の例外、R-tipping判定の凍結枝存続。数値クリップで人工的な炭素消失を作らない。

## 7. コーディングAIへの依頼文

> `docs/spac_carbon_stage_a.md`、`src/control_carbon/manzoni_hydraulics.py`、`src/control_carbon/spac_carbon.py` と関連テストを読み、Stage Aを変更せずStage B1を独立実装してください。まずF_Lのみ、F_Eのみ、両方を比較し、全6×6 Jacobian、平衡分枝、fold、吸引域、capacity写像と実炭素の違いを検証してください。最優先は直列抵抗の正しい扱い、単位、水・炭素収支、F=1極限の回帰テストです。Leuningを実装せず、森林死亡やR-tippingを先取りして主張しないでください。日本語の短い結果報告、再現コマンド、合成パラメータと出典、未解決点を残してください。

推奨ファイル：`src/control_carbon/spac_carbon_two_way.py`、`tests/test_spac_carbon_two_way.py`、`examples/analyze_spac_carbon_two_way.py`、`configs/spac_carbon_stage_b_synthetic.json`、`docs/spac_carbon_stage_b_results.md`。既存成果物を上書きしない。

## 8. 出典と限界

- Manzoni, Katul & Porporato (2014), DOI https://doi.org/10.1002/2013WR015236：2貯留水理・fold。F_L/F_E/F_Wは本研究独自の追加仮説。
- GitHub `docs/manzoni2014_reproduction.md`：2g_x、水収支、非通水退化、単位。
- GitHub `docs/spac_carbon_stage_a.md`：片方向結合、R0、capacity、温度対照。
- Leuning (1995)、Tuzet et al. (2003)：後続の気孔–同化–水理閉包。Stage B1では未導入。

**最重要の科学的注意：葉量依存E_xlだけでは供給側に偏る。蒸散面積のLAI依存を必ず対照に入れ、双方向feedbackの符号を実験で決める。**
